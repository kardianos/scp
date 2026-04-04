package main

import (
	"context"
	"fmt"
	"log"
	"os"
	"os/signal"
	"strings"
	"syscall"
	"time"
)

func main() {
	log.SetFlags(log.Ltime | log.Lshortfile)
	log.SetOutput(os.Stderr)

	ctx, cancel := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer cancel()

	// Load persisted state from disk.
	state, err := LoadState()
	if err != nil {
		log.Printf("scp-runner: warning: load state: %v (starting fresh)", err)
		state = newEmptyState()
	}

	// Create MCP server with no executor initially (sim_setup will configure them).
	server := NewMCPServer(nil)
	server.state = state

	// Start background monitor and link it to the server.
	monitor := NewMonitor(server, "/tmp/scp-runner")
	monitor.state = state
	server.monitor = monitor
	go monitor.Run(ctx)

	// Attempt to recover from persisted state.
	recoverFromState(ctx, server, state)

	log.Println("scp-runner: MCP server starting on stdio")
	if err := server.Serve(ctx); err != nil {
		if ctx.Err() != nil {
			log.Println("scp-runner: shutdown")
			return
		}
		log.Fatalf("scp-runner: %v", err)
	}
}

// recoverFromState attempts to reconnect to previously running remote instances
// and resume monitoring of any runs that were in progress.
func recoverFromState(ctx context.Context, server *MCPServer, state *RunnerState) {
	state.mu.Lock()
	instanceNames := make([]string, 0, len(state.Instances))
	for name := range state.Instances {
		instanceNames = append(instanceNames, name)
	}
	state.mu.Unlock()

	for _, name := range instanceNames {
		recoverInstance(ctx, server, state, name)
	}
}

// recoverInstance attempts to recover a single named instance.
func recoverInstance(ctx context.Context, server *MCPServer, state *RunnerState, name string) {
	pi := state.GetInstance(name)
	if pi == nil || pi.Connection == nil {
		return
	}

	inst := pi.Connection
	log.Printf("scp-runner: recovering instance %q (%d at %s:%d, %s)", name, inst.ID, inst.Host, inst.Port, inst.GPUName)

	remote := NewRemoteExecutor(RemoteConfig{WorkDir: "/tmp/scp-runner/" + name})
	if server.monitor != nil {
		remote.OnRunDone = server.monitor.Wakeup
	}

	instanceName := name // capture for closures
	remote.OnRunComplete = func(id string, info RunInfo) {
		server.handleRunDone(instanceName, id, info)
	}
	remote.OnDownloadDone = func(id, remotePath, localPath string, err error) {
		if err != nil {
			server.sendChannelEvent("download_failed", map[string]string{
				"instance":    instanceName,
				"download_id": id,
				"remote_path": remotePath,
				"local_path":  localPath,
				"error":       err.Error(),
			}, fmt.Sprintf("[%s] Download FAILED: %s -> %s: %v", instanceName, remotePath, localPath, err))
		} else {
			server.sendChannelEvent("download_complete", map[string]string{
				"instance":    instanceName,
				"download_id": id,
				"remote_path": remotePath,
				"local_path":  localPath,
			}, fmt.Sprintf("[%s] Download complete: %s -> %s", instanceName, remotePath, localPath))
		}
	}

	// Initialize the vast client (needed for teardown and instance checks).
	if err := remote.Setup(ctx); err != nil {
		log.Printf("scp-runner: recovery [%s]: vast setup failed: %v (clearing state)", name, err)
		state.ClearInstance(name)
		return
	}

	// Check if the instance is still alive via the Vast API before trying SSH.
	if remote.vast != nil {
		instances, err := remote.vast.ShowInstances(ctx)
		if err != nil {
			log.Printf("scp-runner: recovery [%s]: cannot list instances: %v (clearing state)", name, err)
			state.ClearInstance(name)
			return
		}
		found := false
		for _, vi := range instances {
			if vi.ID == inst.ID {
				if vi.Status == "running" {
					found = true
					// Update host/port in case they changed
					if vi.SSHHost != "" && vi.SSHPort > 0 {
						inst.Host = vi.SSHHost
						inst.Port = vi.SSHPort
					}
				} else {
					log.Printf("scp-runner: recovery [%s]: instance %d status=%s (not running, clearing state)", name, inst.ID, vi.Status)
					state.ClearInstance(name)
					return
				}
				break
			}
		}
		if !found {
			log.Printf("scp-runner: recovery [%s]: instance %d not found in Vast API (destroyed, clearing state)", name, inst.ID)
			state.ClearInstance(name)
			return
		}
	}

	// Try SSH connect with a short timeout.
	connectCtx, connectCancel := context.WithTimeout(ctx, 30*time.Second)
	err := remote.Connect(connectCtx, inst.Host, inst.Port)
	connectCancel()
	if err != nil {
		log.Printf("scp-runner: recovery [%s]: SSH connect to %s:%d failed: %v (clearing state)", name, inst.Host, inst.Port, err)
		state.ClearInstance(name)
		return
	}

	// Restore instance metadata.
	remote.instanceID = inst.ID
	remote.gpuName = inst.GPUName

	// Restore the last-built binary path.
	if pi.Binary != "" {
		remote.lastBinary = pi.Binary
	}

	// Set as active executor.
	server.setExecutor(name, remote)

	log.Printf("scp-runner: recovery [%s]: connected to instance %d", name, inst.ID)

	// Check and resume any runs that were in progress.
	recoveredAny := false
	for id, run := range pi.Runs {
		if run.Status != "running" && run.Status != "starting" {
			continue
		}
		recoverRun(ctx, remote, state, name, id, run)
		recoveredAny = true
	}

	// Trigger an immediate notification check so recovered runs get reported.
	if recoveredAny && server.monitor != nil {
		server.monitor.Wakeup()
	}
}

// recoverRun checks if a previously running simulation is still alive on the remote
// and either resumes monitoring or marks it as complete/failed.
func recoverRun(ctx context.Context, remote *RemoteExecutor, state *RunnerState, instanceName string, id string, pRun *PersistRun) {
	log.Printf("scp-runner: recovery [%s]: checking run %s", instanceName, id)

	// Check if the process is still running by looking for the exit file.
	exitFile := fmt.Sprintf("/root/run_%s.exit", id)
	logFile := fmt.Sprintf("/root/run_%s.log", id)
	exitOut, _ := remote.sshRun(ctx, fmt.Sprintf("cat %s 2>/dev/null", exitFile))
	exitCode := strings.TrimSpace(exitOut)

	if exitCode != "" {
		// Process already exited. Determine status.
		finalStatus := "complete"
		if exitCode != "0" {
			finalStatus = "failed"
		}
		log.Printf("scp-runner: recovery [%s]: run %s already finished (exit=%s)", instanceName, id, exitCode)

		// Do a final download if auto_download was configured.
		if pRun.AutoDownload != "" && len(pRun.OutputFiles) > 0 {
			log.Printf("scp-runner: recovery [%s]: doing final download for run %s to %s", instanceName, id, pRun.AutoDownload)
			for _, f := range pRun.OutputFiles {
				remotePath := pRun.RemoteDir + "/" + f
				if err := remote.rsyncFile(ctx, remotePath, pRun.AutoDownload); err != nil {
					log.Printf("scp-runner: recovery [%s]: download %s: %v", instanceName, f, err)
				}
			}
		}

		state.UpdateRunStatus(instanceName, id, finalStatus)
		// Also store in the executor's run map so RunStatus works.
		rr := &remoteRun{
			info: RunInfo{
				ID:     id,
				Status: RunState(finalStatus),
			},
		}
		// Grab last diag line.
		diagOut, _ := remote.sshRun(ctx, fmt.Sprintf("tail -1 %s 2>/dev/null", logFile))
		rr.info.LastDiag = strings.TrimSpace(diagOut)
		if finalStatus == "failed" {
			tailOut, _ := remote.sshRun(ctx, fmt.Sprintf("tail -20 %s 2>/dev/null", logFile))
			rr.info.Error = fmt.Sprintf("exit code %s\n%s", exitCode, strings.TrimSpace(tailOut))
		}
		remote.runs.Store(id, rr)
		return
	}

	// No exit file yet. Check if the process is still alive by looking at the log file.
	sizeOut, err := remote.sshRun(ctx, fmt.Sprintf("stat -c%%s %s 2>/dev/null", logFile))
	if err != nil {
		// No log file at all means the process never ran or was cleaned up.
		log.Printf("scp-runner: recovery [%s]: run %s has no log file, marking failed", instanceName, id)
		state.UpdateRunStatus(instanceName, id, "failed")
		rr := &remoteRun{
			info: RunInfo{
				ID:     id,
				Status: RunFailed,
				Error:  "no log file found during recovery",
			},
		}
		remote.runs.Store(id, rr)
		return
	}

	// Process appears to still be running. Find its PID and resume monitoring.
	log.Printf("scp-runner: recovery [%s]: run %s still alive (log size=%s), resuming monitor", instanceName, id, strings.TrimSpace(sizeOut))

	runCtx, cancel := context.WithCancel(ctx)
	adCtx, adCancel := context.WithCancel(context.Background())
	rr := &remoteRun{
		info: RunInfo{
			ID:     id,
			Status: RunRunning,
		},
		cancel:         cancel,
		notifyInterval: 60 * time.Second, // default for recovered runs
		autoDownload:   pRun.AutoDownload,
		remoteDir:      pRun.RemoteDir,
	}

	// Reconstruct output file paths for auto-download.
	if pRun.AutoDownload != "" && len(pRun.OutputFiles) > 0 {
		for _, f := range pRun.OutputFiles {
			rr.outputFiles = append(rr.outputFiles, pRun.RemoteDir+"/"+f)
		}
	}

	remote.runs.Store(id, rr)

	// Start auto-download if configured.
	remote.StartAutoDownload(adCtx, rr)

	// Resume the monitoring goroutine.
	go func() {
		defer adCancel()
		defer cancel()
		defer func() {
			// Wait for auto-download final sync.
			adCancel()
			if rr.adDone != nil {
				<-rr.adDone
			}
			if remote.OnRunDone != nil {
				remote.OnRunDone()
			}
		}()

		startWall := time.Now()
		ticker := time.NewTicker(5 * time.Second)
		defer ticker.Stop()

		for {
			select {
			case <-runCtx.Done():
				rr.mu.Lock()
				rr.info.Status = RunCancelled
				rr.info.WallSecs = time.Since(startWall).Seconds()
				rr.mu.Unlock()
				state.UpdateRunStatus(instanceName, id, "cancelled")
				return
			case <-ticker.C:
			}

			// Check exit file.
			exitOut, _ := remote.sshRun(runCtx, fmt.Sprintf("cat %s 2>/dev/null", exitFile))
			exitCode := strings.TrimSpace(exitOut)

			// Read recent log lines.
			diagOut, _ := remote.sshRun(runCtx, fmt.Sprintf("tail -20 %s 2>/dev/null", logFile))
			lines := strings.Split(strings.TrimSpace(diagOut), "\n")

			rr.mu.Lock()
			rr.info.WallSecs = time.Since(startWall).Seconds()
			if len(lines) > 0 && lines[len(lines)-1] != "" {
				rr.info.LastDiag = lines[len(lines)-1]
				for i := len(lines) - 1; i >= 0; i-- {
					if strings.HasPrefix(lines[i], "t=") {
						rest := strings.TrimPrefix(lines[i], "t=")
						if idx := strings.Index(rest, " E"); idx > 0 {
							tStr := strings.TrimSpace(rest[:idx])
							if tVal, parseErr := parseFloat(tStr); parseErr == nil {
								rr.info.SimTime = tVal
							}
						}
						break
					}
				}
			}

			if exitCode != "" {
				// Process exited.
				if exitCode == "0" {
					rr.info.Status = RunComplete
				} else {
					rr.info.Status = RunFailed
					tailOut, _ := remote.sshRun(runCtx, fmt.Sprintf("tail -20 %s 2>/dev/null", logFile))
					rr.info.Error = fmt.Sprintf("exit code %s\n%s", exitCode, strings.TrimSpace(tailOut))
				}
				rr.mu.Unlock()
				state.UpdateRunStatus(instanceName, id, string(rr.info.Status))
				return
			}
			rr.mu.Unlock()
		}
	}()
}

// parseFloat is a small helper to avoid importing strconv everywhere.
func parseFloat(s string) (float64, error) {
	var f float64
	_, err := fmt.Sscanf(s, "%f", &f)
	return f, err
}
