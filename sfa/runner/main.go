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
		state = &RunnerState{Runs: make(map[string]*PersistRun)}
	}

	// Create MCP server with no executor initially (sim_setup will configure one).
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

// recoverFromState attempts to reconnect to a previously running remote instance
// and resume monitoring of any runs that were in progress.
func recoverFromState(ctx context.Context, server *MCPServer, state *RunnerState) {
	if state.Instance == nil {
		return
	}

	inst := state.Instance
	log.Printf("scp-runner: recovering instance %d at %s:%d (%s)", inst.ID, inst.Host, inst.Port, inst.GPUName)

	remote := NewRemoteExecutor(RemoteConfig{WorkDir: "/tmp/scp-runner"})
	if server.monitor != nil {
		remote.OnRunDone = server.monitor.Wakeup
	}

	// Initialize the vast client (needed for teardown).
	if err := remote.Setup(ctx); err != nil {
		log.Printf("scp-runner: recovery: vast setup failed: %v (clearing state)", err)
		state.ClearInstance()
		return
	}

	// Try SSH connect with a short timeout.
	connectCtx, connectCancel := context.WithTimeout(ctx, 30*time.Second)
	err := remote.Connect(connectCtx, inst.Host, inst.Port)
	connectCancel()
	if err != nil {
		log.Printf("scp-runner: recovery: SSH connect to %s:%d failed: %v (clearing state)", inst.Host, inst.Port, err)
		state.ClearInstance()
		return
	}

	// Restore instance metadata.
	remote.instanceID = inst.ID
	remote.gpuName = inst.GPUName

	// Restore the last-built binary path.
	if state.Binary != "" {
		remote.lastBinary = state.Binary
	}

	// Set as active executor.
	server.mu.Lock()
	server.executor = remote
	server.mu.Unlock()

	log.Printf("scp-runner: recovery: connected to instance %d", inst.ID)

	// Check and resume any runs that were in progress.
	recoveredAny := false
	for id, run := range state.Runs {
		if run.Status != "running" && run.Status != "starting" {
			continue
		}
		recoverRun(ctx, remote, state, id, run)
		recoveredAny = true
	}

	// Trigger an immediate notification check so recovered runs get reported.
	if recoveredAny && server.monitor != nil {
		server.monitor.Wakeup()
	}
}

// recoverRun checks if a previously running simulation is still alive on the remote
// and either resumes monitoring or marks it as complete/failed.
func recoverRun(ctx context.Context, remote *RemoteExecutor, state *RunnerState, id string, pRun *PersistRun) {
	log.Printf("scp-runner: recovery: checking run %s", id)

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
		log.Printf("scp-runner: recovery: run %s already finished (exit=%s)", id, exitCode)

		// Do a final download if auto_download was configured.
		if pRun.AutoDownload != "" && len(pRun.OutputFiles) > 0 {
			log.Printf("scp-runner: recovery: doing final download for run %s to %s", id, pRun.AutoDownload)
			for _, f := range pRun.OutputFiles {
				remotePath := pRun.RemoteDir + "/" + f
				if err := remote.rsyncFile(ctx, remotePath, pRun.AutoDownload); err != nil {
					log.Printf("scp-runner: recovery: download %s: %v", f, err)
				}
			}
		}

		state.UpdateRunStatus(id, finalStatus)
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
	// We check if the log file exists and is growing by comparing sizes.
	sizeOut, err := remote.sshRun(ctx, fmt.Sprintf("stat -c%%s %s 2>/dev/null", logFile))
	if err != nil {
		// No log file at all means the process never ran or was cleaned up.
		log.Printf("scp-runner: recovery: run %s has no log file, marking failed", id)
		state.UpdateRunStatus(id, "failed")
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
	log.Printf("scp-runner: recovery: run %s still alive (log size=%s), resuming monitor", id, strings.TrimSpace(sizeOut))

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
				state.UpdateRunStatus(id, "cancelled")
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
				state.UpdateRunStatus(id, string(rr.info.Status))
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
