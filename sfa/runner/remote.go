package main

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"net"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"

	"golang.org/x/crypto/ssh"
)

// RemoteExecutor runs simulations on a remote GPU via SSH.
type RemoteExecutor struct {
	vast       *VastClient
	instanceID int
	sshHost    string
	sshPort    int
	gpuName    string
	sshClient  *ssh.Client
	startTime  time.Time
	workDir    string // local work dir for downloads
	lastBinary     string                            // set after successful Build
	OnRunDone      func()                            // called when any run reaches terminal state
	OnRunComplete  func(id string, info RunInfo)     // called with run details on completion
	OnDownloadDone func(id, remote, local string, err error) // called when a download completes or fails

	runs      sync.Map // map[string]*remoteRun
	downloads sync.Map // map[string]*DownloadInfo
	mu        sync.Mutex
}

type remoteRun struct {
	info           RunInfo
	cancel         context.CancelFunc
	notifyInterval time.Duration
	autoDownload   string   // local directory for auto-download
	outputFiles    []string // remote output files to sync
	remoteDir      string   // remote directory containing outputs (e.g. /root)
	adDone         chan struct{} // closed when auto-download goroutine exits
	mu             sync.Mutex
}

// RemoteConfig holds parameters for creating a remote executor.
type RemoteConfig struct {
	GPUFilter string // e.g. "gpu_name=Tesla_V100 num_gpus=1 rentable=true disk_space>=20"
	Image     string // e.g. "nvidia/cuda:12.2.0-devel-ubuntu22.04"
	DiskGB    int
	Onstart   string // startup script
	WorkDir   string // local directory for downloaded results
}

const (
	defaultImage = "nvidia/cuda:12.2.0-devel-ubuntu22.04"
	defaultDiskGB    = 100
	defaultOnstart   = "apt-get update && apt-get install -y libzstd-dev rsync && ldconfig"
)

func NewRemoteExecutor(cfg RemoteConfig) *RemoteExecutor {
	return &RemoteExecutor{
		workDir: cfg.WorkDir,
	}
}

func (r *RemoteExecutor) Type() ExecType { return ExecRemote }

func (r *RemoteExecutor) Setup(ctx context.Context) error {
	var err error
	r.vast, err = NewVastClient()
	if err != nil {
		return fmt.Errorf("vast client: %w", err)
	}

	r.startTime = time.Now()
	if err := os.MkdirAll(r.workDir, 0755); err != nil {
		return fmt.Errorf("create workdir: %w", err)
	}
	return nil
}

// Provision searches for a GPU offer, creates an instance, waits for it to be ready, and connects.
func (r *RemoteExecutor) Provision(ctx context.Context, gpuFilter map[string]string, diskGB int) error {
	if r.vast == nil {
		return fmt.Errorf("vast client not initialized — call Setup first")
	}
	if len(gpuFilter) == 0 {
		gpuFilter = map[string]string{"gpu_name": "Tesla_V100"}
	}
	if diskGB <= 0 {
		diskGB = defaultDiskGB
	}

	// Convert map to legacy filter string for SearchOffers.
	filterStr := gpuFilterMapToString(gpuFilter)

	// Check for existing running instances first.
	instances, err := r.vast.ShowInstances(ctx)
	if err == nil && len(instances) > 0 {
		for _, inst := range instances {
			if inst.Status == "running" && inst.SSHHost != "" && inst.SSHPort > 0 {
				fmt.Fprintf(os.Stderr, "scp-runner: reusing existing instance %d (%s:%d)\n", inst.ID, inst.SSHHost, inst.SSHPort)
				r.instanceID = inst.ID
				r.gpuName = inst.GPUName
				return r.Connect(ctx, inst.SSHHost, inst.SSHPort)
			}
		}
	}

	// Search for offers.
	offers, err := r.vast.SearchOffers(ctx, filterStr)
	if err != nil {
		return fmt.Errorf("search offers: %w", err)
	}
	if len(offers) == 0 {
		return fmt.Errorf("no GPU offers found for filter: %v", gpuFilter)
	}

	// Try offers in order (cheapest first). Offers are ephemeral — the first
	// one may be snatched between search and create, so retry with the next.
	maxTries := len(offers)
	if maxTries > 5 {
		maxTries = 5
	}
	var instanceID int
	var lastErr error
	for i := 0; i < maxTries; i++ {
		offer := offers[i]
		fmt.Fprintf(os.Stderr, "scp-runner: trying offer %d/%d: %s (%d GB VRAM, %s, $%.3f/hr) [%d offers]\n",
			i+1, maxTries, offer.GPUName, offer.GPUMemMB/1024, offer.Geolocation, offer.DPHTot, len(offers))

		instanceID, lastErr = r.vast.CreateInstance(ctx, offer.ID, defaultImage, diskGB, defaultOnstart)
		if lastErr == nil {
			r.instanceID = instanceID
			r.gpuName = offer.GPUName
			break
		}
		fmt.Fprintf(os.Stderr, "scp-runner: offer %d failed (%v), trying next...\n", offer.ID, lastErr)
	}
	if lastErr != nil && r.instanceID == 0 {
		return fmt.Errorf("create instance: all %d offers failed, last error: %w", maxTries, lastErr)
	}

	// Wait for instance to become ready with SSH info.
	fmt.Fprintf(os.Stderr, "scp-runner: waiting for instance %d to start...\n", instanceID)
	deadline := time.Now().Add(5 * time.Minute)
	for {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}
		if time.Now().After(deadline) {
			return fmt.Errorf("instance %d did not become ready within 5 minutes", instanceID)
		}

		instances, err := r.vast.ShowInstances(ctx)
		if err != nil {
			time.Sleep(5 * time.Second)
			continue
		}
		for _, i := range instances {
			if i.ID == instanceID && i.Status == "running" && i.SSHHost != "" && i.SSHPort > 0 {
				fmt.Fprintf(os.Stderr, "scp-runner: instance ready at %s:%d\n", i.SSHHost, i.SSHPort)
				if err := r.Connect(ctx, i.SSHHost, i.SSHPort); err != nil {
					return err
				}
				return r.waitOnstart(ctx)
			}
		}
		time.Sleep(5 * time.Second)
	}
}

// waitOnstart polls until the onstart script has finished by checking for libzstd and nvcc.
func (r *RemoteExecutor) waitOnstart(_ context.Context) error {
	fmt.Fprintf(os.Stderr, "scp-runner: waiting for onstart script to finish...\n")
	deadline := time.Now().Add(3 * time.Minute)
	for {
		// Check if libzstd and nvcc are available (installed by onstart).
		out, err := r.sshRun(context.Background(), "ldconfig -p | grep -q libzstd && which nvcc && echo READY")
		if err == nil && strings.Contains(out, "READY") {
			fmt.Fprintf(os.Stderr, "scp-runner: onstart complete, build tools ready\n")
			return nil
		}
		if time.Now().After(deadline) {
			return fmt.Errorf("onstart did not complete within 3 minutes (last: %s)", strings.TrimSpace(out))
		}
		time.Sleep(5 * time.Second)
	}
}

// Connect establishes SSH to an existing instance.
// Uses its own 3-minute timeout independent of the caller's context,
// since MCP tool calls may have short deadlines but SSH needs time.
func (r *RemoteExecutor) Connect(_ context.Context, host string, port int) error {
	r.sshHost = host
	r.sshPort = port

	signer, err := loadSSHKey()
	if err != nil {
		return fmt.Errorf("load ssh key: %w", err)
	}

	config := &ssh.ClientConfig{
		User:            "root",
		Auth:            []ssh.AuthMethod{ssh.PublicKeys(signer)},
		HostKeyCallback: ssh.InsecureIgnoreHostKey(),
		Timeout:         15 * time.Second,
	}

	addr := net.JoinHostPort(host, strconv.Itoa(port))

	// Retry SSH connection with a fixed 3-minute deadline.
	deadline := time.Now().Add(3 * time.Minute)
	var client *ssh.Client
	for attempt := 0; ; attempt++ {
		client, err = ssh.Dial("tcp", addr, config)
		if err == nil {
			break
		}
		if time.Now().After(deadline) {
			return fmt.Errorf("ssh dial %s after %d attempts: %w", addr, attempt+1, err)
		}
		time.Sleep(5 * time.Second)
	}

	r.sshClient = client
	return nil
}

func (r *RemoteExecutor) Teardown(ctx context.Context) error {
	if r.sshClient != nil {
		r.sshClient.Close()
	}
	if r.instanceID > 0 && r.vast != nil {
		return r.vast.DestroyInstance(ctx, r.instanceID)
	}
	return nil
}

func (r *RemoteExecutor) sshRun(ctx context.Context, cmd string) (string, error) {
	if r.sshClient == nil {
		return "", fmt.Errorf("ssh not connected")
	}
	session, err := r.sshClient.NewSession()
	if err != nil {
		return "", fmt.Errorf("new session: %w", err)
	}
	defer session.Close()

	// Use a channel to respect context cancellation.
	type result struct {
		out string
		err error
	}
	ch := make(chan result, 1)
	go func() {
		out, err := session.CombinedOutput(cmd)
		ch <- result{string(out), err}
	}()

	select {
	case <-ctx.Done():
		session.Signal(ssh.SIGTERM)
		return "", ctx.Err()
	case res := <-ch:
		return res.out, res.err
	}
}

func (r *RemoteExecutor) Status(ctx context.Context) (*ExecutorStatus, error) {
	s := &ExecutorStatus{
		Type:     ExecRemote,
		Location: fmt.Sprintf("%s:%d", r.sshHost, r.sshPort),
		GPUUtil:  -1,
		Uptime:   time.Since(r.startTime).Round(time.Second).String(),
	}

	// Count active runs.
	r.runs.Range(func(_, v any) bool {
		rr := v.(*remoteRun)
		rr.mu.Lock()
		switch rr.info.Status {
		case RunStarting, RunRunning:
			s.ActiveRuns++
		}
		rr.mu.Unlock()
		return true
	})

	// GPU info via nvidia-smi.
	out, err := r.sshRun(ctx, "nvidia-smi --query-gpu=name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits")
	if err == nil {
		parts := strings.Split(strings.TrimSpace(out), ", ")
		if len(parts) >= 4 {
			s.GPUName = parts[0]
			s.GPUUtil, _ = strconv.Atoi(strings.TrimSpace(parts[1]))
			s.GPUMemUsedMB, _ = strconv.ParseInt(strings.TrimSpace(parts[2]), 10, 64)
			s.GPUMemTotalMB, _ = strconv.ParseInt(strings.TrimSpace(parts[3]), 10, 64)
		}
	}

	// Disk info.
	out, err = r.sshRun(ctx, "df -BG --output=used,avail /root | tail -1")
	if err == nil {
		fields := strings.Fields(strings.TrimSpace(out))
		if len(fields) >= 2 {
			used, _ := strconv.ParseFloat(strings.TrimSuffix(fields[0], "G"), 64)
			avail, _ := strconv.ParseFloat(strings.TrimSuffix(fields[1], "G"), 64)
			s.DiskUsedGB = used
			s.DiskFreeGB = avail
		}
	}

	return s, nil
}

func (r *RemoteExecutor) Build(ctx context.Context, sources []string, cmd string) (*BuildResult, error) {
	hash, err := hashBuild(sources, cmd)
	if err != nil {
		return nil, fmt.Errorf("hash sources: %w", err)
	}

	binName := fmt.Sprintf("scp_sim_%s", hash)
	remoteBin := "/root/" + binName

	// Check remote cache.
	_, err = r.sshRun(ctx, fmt.Sprintf("test -f %s", remoteBin))
	if err == nil {
		r.mu.Lock()
		r.lastBinary = remoteBin
		r.mu.Unlock()
		return &BuildResult{Status: "cached", Binary: remoteBin, Cached: true}, nil
	}

	// Upload sources preserving directory structure under /root/src/.
	// This allows relative includes (e.g. "../format/sfa.h") to work.
	var cuFiles []string
	for _, src := range sources {
		remotePath := "/root/src/" + src
		remoteDir := "/root/src/" + filepath.Dir(src)
		r.sshRun(ctx, fmt.Sprintf("mkdir -p %s", remoteDir))
		if err := r.Upload(ctx, src, remotePath); err != nil {
			return nil, fmt.Errorf("upload %s: %w", src, err)
		}
		ext := filepath.Ext(src)
		if ext == ".cu" || ext == ".c" {
			cuFiles = append(cuFiles, remotePath)
		}
	}

	// Build.
	if cmd == "" {
		// Auto-detect: use nvcc for .cu, gcc for .c
		if len(cuFiles) > 0 && filepath.Ext(cuFiles[0]) == ".cu" {
			cmd = fmt.Sprintf("nvcc -O3 -arch=sm_70 -o %s %s -lzstd -lm -lpthread",
				remoteBin, strings.Join(cuFiles, " "))
		} else {
			cmd = fmt.Sprintf("gcc -O3 -march=native -fopenmp -o %s %s -lzstd -lm",
				remoteBin, strings.Join(cuFiles, " "))
		}
	} else {
		cmd = strings.ReplaceAll(cmd, "${OUTPUT}", remoteBin)
	}

	out, err := r.sshRun(ctx, cmd)
	if err != nil {
		return &BuildResult{Status: "failed", Error: fmt.Sprintf("%v: %s", err, out)}, nil
	}

	r.mu.Lock()
	r.lastBinary = remoteBin
	r.mu.Unlock()
	return &BuildResult{Status: "built", Binary: remoteBin}, nil
}

func (r *RemoteExecutor) Run(ctx context.Context, config string, id string, notifyInterval time.Duration) error {
	if _, loaded := r.runs.Load(id); loaded {
		return fmt.Errorf("run %s already exists", id)
	}

	runCtx, cancel := context.WithCancel(ctx)
	rr := &remoteRun{
		info: RunInfo{
			ID:     id,
			Status: RunStarting,
		},
		cancel:         cancel,
		notifyInterval: notifyInterval,
	}
	r.runs.Store(id, rr)

	go r.executeRemoteRun(runCtx, rr, config)
	return nil
}

// RunWithAutoDownload starts a run and configures auto-download if autoDownload is set.
func (r *RemoteExecutor) RunWithAutoDownload(ctx context.Context, config string, id string, notifyInterval time.Duration, autoDownload string) error {
	if _, loaded := r.runs.Load(id); loaded {
		return fmt.Errorf("run %s already exists", id)
	}

	runCtx, cancel := context.WithCancel(ctx)
	rr := &remoteRun{
		info: RunInfo{
			ID:     id,
			Status: RunStarting,
		},
		cancel:         cancel,
		notifyInterval: notifyInterval,
		autoDownload:   autoDownload,
		remoteDir:      "/root",
	}

	// Parse output file names from config.
	if autoDownload != "" && isConfigContent(config) {
		var files []string
		if out := parseConfigValue(config, "output"); out != "" {
			files = append(files, "/root/"+out)
		}
		if diag := parseConfigValue(config, "diag_file"); diag != "" {
			files = append(files, "/root/"+diag)
		}
		rr.outputFiles = files
	}

	r.runs.Store(id, rr)

	go r.executeRemoteRun(runCtx, rr, config)
	return nil
}

func (r *RemoteExecutor) executeRemoteRun(ctx context.Context, rr *remoteRun, config string) {
	// We use a separate context for auto-download so we can signal it to do a
	// final sync when the run finishes, while still allowing it to complete.
	adCtx, adCancel := context.WithCancel(context.Background())
	defer adCancel()
	defer rr.cancel()
	defer func() {
		// Signal auto-download to do final sync, then wait for it.
		adCancel()
		if rr.adDone != nil {
			<-rr.adDone
		}
		if r.OnRunComplete != nil {
			rr.mu.Lock()
			info := rr.info
			rr.mu.Unlock()
			r.OnRunComplete(info.ID, info)
		}
		if r.OnRunDone != nil {
			r.OnRunDone()
		}
	}()
	startWall := time.Now()

	// If config looks like file content (key=value), write it remotely and use lastBinary.
	runCmd := config
	if isConfigContent(config) {
		r.mu.Lock()
		binPath := r.lastBinary
		r.mu.Unlock()
		if binPath == "" {
			rr.mu.Lock()
			rr.info.Status = RunFailed
			rr.info.Error = "no binary built yet — call sim_build first"
			rr.mu.Unlock()
			return
		}

		cfgPath := fmt.Sprintf("/root/run_%s.cfg", rr.info.ID)
		// Write config file to remote.
		escaped := strings.ReplaceAll(config, "'", "'\\''")
		_, err := r.sshRun(ctx, fmt.Sprintf("cat > %s << 'CFGEOF'\n%s\nCFGEOF", cfgPath, escaped))
		if err != nil {
			rr.mu.Lock()
			rr.info.Status = RunFailed
			rr.info.Error = fmt.Sprintf("write config: %v", err)
			rr.mu.Unlock()
			return
		}

		// Parse TotalTime from config.
		for _, line := range strings.Split(config, "\n") {
			line = strings.TrimSpace(line)
			if strings.HasPrefix(line, "T") {
				parts := strings.SplitN(line, "=", 2)
				if len(parts) == 2 && strings.TrimSpace(parts[0]) == "T" {
					if t, err := strconv.ParseFloat(strings.TrimSpace(parts[1]), 64); err == nil {
						rr.mu.Lock()
						rr.info.TotalTime = t
						rr.mu.Unlock()
					}
				}
			}
		}

		runCmd = fmt.Sprintf("%s %s", binPath, cfgPath)
	}

	// Start simulation with nohup so it survives SSH drops.
	// Write exit code to a .exit file so we can check it from another session.
	logFile := fmt.Sprintf("/root/run_%s.log", rr.info.ID)
	exitFile := fmt.Sprintf("/root/run_%s.exit", rr.info.ID)
	cmd := fmt.Sprintf("nohup bash -c '%s; echo $? > %s' > %s 2>&1 & echo $!", runCmd, exitFile, logFile)

	out, err := r.sshRun(ctx, cmd)
	if err != nil {
		rr.mu.Lock()
		rr.info.Status = RunFailed
		rr.info.Error = fmt.Sprintf("start: %v: %s", err, out)
		rr.mu.Unlock()
		return
	}

	pid := strings.TrimSpace(out)
	rr.mu.Lock()
	rr.info.Status = RunRunning
	rr.mu.Unlock()

	// Start auto-download goroutine if configured.
	r.StartAutoDownload(adCtx, rr)

	// Monitor loop: check if process is alive and tail diag.
	ticker := time.NewTicker(5 * time.Second)
	defer ticker.Stop()

	for {
		select {
		case <-ctx.Done():
			// Kill remote process.
			r.sshRun(context.Background(), fmt.Sprintf("kill %s 2>/dev/null", pid))
			rr.mu.Lock()
			rr.info.Status = RunCancelled
			rr.info.WallSecs = time.Since(startWall).Seconds()
			rr.mu.Unlock()
			return
		case <-ticker.C:
		}

		// Check if process is still running.
		_, err := r.sshRun(ctx, fmt.Sprintf("kill -0 %s 2>/dev/null", pid))
		alive := err == nil

		// Read recent log lines: last line for LastDiag, last "t=" line for SimTime.
		diagOut, _ := r.sshRun(ctx, fmt.Sprintf("tail -20 %s 2>/dev/null", logFile))
		lines := strings.Split(strings.TrimSpace(diagOut), "\n")

		rr.mu.Lock()
		rr.info.WallSecs = time.Since(startWall).Seconds()
		if len(lines) > 0 && lines[len(lines)-1] != "" {
			rr.info.LastDiag = lines[len(lines)-1]
			// Find the last "t=" progress line for SimTime.
			for i := len(lines) - 1; i >= 0; i-- {
				if strings.HasPrefix(lines[i], "t=") {
					rest := strings.TrimPrefix(lines[i], "t=")
					if idx := strings.Index(rest, " E"); idx > 0 {
						if t, err := strconv.ParseFloat(strings.TrimSpace(rest[:idx]), 64); err == nil {
							rr.info.SimTime = t
						}
					}
					break
				}
			}
		}

		if !alive {
			// Process is gone. Check exit code from a wrapper file written at start.
			// Since `wait` doesn't work across SSH sessions, we check if the log
			// contains error indicators or if the nohup process wrote an exit marker.
			exitOut, _ := r.sshRun(ctx, fmt.Sprintf("cat /root/run_%s.exit 2>/dev/null", rr.info.ID))
			exitCode := strings.TrimSpace(exitOut)
			if exitCode == "" {
				// No exit file — assume success if log has content, failure otherwise.
				logCheck, _ := r.sshRun(ctx, fmt.Sprintf("wc -l < %s 2>/dev/null", logFile))
				lines := strings.TrimSpace(logCheck)
				if lines != "" && lines != "0" {
					rr.info.Status = RunComplete
				} else {
					rr.info.Status = RunFailed
					rr.info.Error = "process exited with no output"
				}
			} else if exitCode == "0" {
				rr.info.Status = RunComplete
			} else {
				// Non-zero exit — grab last 20 lines of log for error context.
				tailOut, _ := r.sshRun(ctx, fmt.Sprintf("tail -20 %s 2>/dev/null", logFile))
				tail := strings.TrimSpace(tailOut)
				if tail != "" {
					rr.info.Error = fmt.Sprintf("exit code %s\n%s", exitCode, tail)
				} else {
					rr.info.Error = fmt.Sprintf("exit code %s", exitCode)
				}
				rr.info.Status = RunFailed
			}
			rr.mu.Unlock()
			return
		}
		rr.mu.Unlock()
	}
}

// rsyncFile does a single rsync --append-verify of a remote file to a local directory.
// Returns nil on success, error on failure.
func (r *RemoteExecutor) rsyncFile(ctx context.Context, remotePath, localDir string) error {
	if err := os.MkdirAll(localDir, 0755); err != nil {
		return fmt.Errorf("mkdir %s: %w", localDir, err)
	}
	cmd := fmt.Sprintf("rsync -az --partial --append-verify -e 'ssh -o StrictHostKeyChecking=no -p %d' root@%s:%s %s/",
		r.sshPort, r.sshHost, remotePath, localDir)
	out, err := execOutput(ctx, "bash", "-c", cmd)
	if err != nil {
		return fmt.Errorf("rsync %s: %v: %s", remotePath, err, out)
	}
	return nil
}

// autoDownloadLoop periodically rsyncs output files to a local directory.
// It also scans for .sfp part files, downloads them, verifies, and deletes
// the remote copy to free disk space.
// Runs until the run's context is cancelled (run done or cancelled), then
// does one final sync. Closes rr.adDone when it exits.
func (r *RemoteExecutor) autoDownloadLoop(ctx context.Context, rr *remoteRun) {
	defer close(rr.adDone)

	ticker := time.NewTicker(60 * time.Second)
	defer ticker.Stop()

	syncAll := func() {
		rr.mu.Lock()
		files := rr.outputFiles
		localDir := rr.autoDownload
		remoteDir := rr.remoteDir
		rr.mu.Unlock()

		// Sync the main output files (sfa, tsv)
		for _, f := range files {
			dlCtx, cancel := context.WithTimeout(context.Background(), 5*time.Minute)
			if err := r.rsyncFile(dlCtx, f, localDir); err != nil {
				fmt.Fprintf(os.Stderr, "scp-runner: auto-download %s: %v\n", f, err)
			}
			cancel()
		}

		// Scan for completed .sfp part files (not .tmp) and download+delete them
		r.syncAndCleanSFP(context.Background(), remoteDir, localDir)
	}

	for {
		select {
		case <-ctx.Done():
			// Run is done. Do one final sync.
			fmt.Fprintf(os.Stderr, "scp-runner: auto-download final sync for run %s\n", rr.info.ID)
			syncAll()
			return
		case <-ticker.C:
			syncAll()
		}
	}
}

// syncAndCleanSFP lists .sfp files on the remote (excluding .tmp), downloads
// each one, verifies the local size matches, and deletes the remote copy.
func (r *RemoteExecutor) syncAndCleanSFP(ctx context.Context, remoteDir, localDir string) {
	// List .sfp files (completed only — .sfp.tmp excluded by the glob)
	listCmd := fmt.Sprintf("ls -1 %s/*.sfp 2>/dev/null || true", remoteDir)
	listCtx, cancel := context.WithTimeout(ctx, 30*time.Second)
	defer cancel()
	out, err := r.sshRun(listCtx, listCmd)
	if err != nil || strings.TrimSpace(out) == "" {
		return
	}

	for _, remotePath := range strings.Split(strings.TrimSpace(out), "\n") {
		remotePath = strings.TrimSpace(remotePath)
		if remotePath == "" || strings.HasSuffix(remotePath, ".tmp") {
			continue
		}

		// Get remote file size
		sizeCtx, sizeCancel := context.WithTimeout(ctx, 10*time.Second)
		sizeOut, err := r.sshRun(sizeCtx, fmt.Sprintf("stat -c%%s %s", remotePath))
		sizeCancel()
		if err != nil {
			continue
		}
		remoteSize, _ := strconv.ParseInt(strings.TrimSpace(sizeOut), 10, 64)
		if remoteSize == 0 {
			continue
		}

		// Download
		dlCtx, dlCancel := context.WithTimeout(ctx, 10*time.Minute)
		err = r.rsyncFile(dlCtx, remotePath, localDir)
		dlCancel()
		if err != nil {
			fmt.Fprintf(os.Stderr, "scp-runner: sfp download %s: %v\n", remotePath, err)
			continue
		}

		// Verify local size matches
		baseName := remotePath[strings.LastIndex(remotePath, "/")+1:]
		localPath := filepath.Join(localDir, baseName)
		info, err := os.Stat(localPath)
		if err != nil || info.Size() != remoteSize {
			fmt.Fprintf(os.Stderr, "scp-runner: sfp verify failed %s (remote=%d local=%d)\n",
				baseName, remoteSize, func() int64 { if info != nil { return info.Size() }; return -1 }())
			continue
		}

		// Delete remote copy to free disk
		rmCtx, rmCancel := context.WithTimeout(ctx, 10*time.Second)
		_, err = r.sshRun(rmCtx, fmt.Sprintf("rm -f %s", remotePath))
		rmCancel()
		if err != nil {
			fmt.Fprintf(os.Stderr, "scp-runner: sfp delete %s: %v\n", remotePath, err)
		} else {
			fmt.Fprintf(os.Stderr, "scp-runner: sfp downloaded+deleted %s (%d bytes)\n", baseName, remoteSize)
		}
	}
}

// StartAutoDownload launches the auto-download goroutine for an existing run.
// Used both for new runs and for recovery of previously running runs.
func (r *RemoteExecutor) StartAutoDownload(runCtx context.Context, rr *remoteRun) {
	rr.mu.Lock()
	if rr.autoDownload == "" || len(rr.outputFiles) == 0 {
		rr.mu.Unlock()
		return
	}
	rr.adDone = make(chan struct{})
	rr.mu.Unlock()
	go r.autoDownloadLoop(runCtx, rr)
}

func (r *RemoteExecutor) RunStatus(id string) *RunInfo {
	v, ok := r.runs.Load(id)
	if !ok {
		return nil
	}
	rr := v.(*remoteRun)
	rr.mu.Lock()
	defer rr.mu.Unlock()
	info := rr.info
	return &info
}

func (r *RemoteExecutor) RunCancel(id string) error {
	v, ok := r.runs.Load(id)
	if !ok {
		return fmt.Errorf("run %s not found", id)
	}
	rr := v.(*remoteRun)
	rr.cancel()
	return nil
}

func (r *RemoteExecutor) Upload(ctx context.Context, localPath, remotePath string) error {
	addr := fmt.Sprintf("%s:%d", r.sshHost, r.sshPort)
	cmd := fmt.Sprintf("scp -o StrictHostKeyChecking=no -P %d %s root@%s:%s",
		r.sshPort, localPath, r.sshHost, remotePath)
	_ = addr
	out, err := execOutput(ctx, "bash", "-c", cmd)
	if err != nil {
		return fmt.Errorf("scp upload: %v: %s", err, out)
	}
	return nil
}

func (r *RemoteExecutor) Download(ctx context.Context, remotePath, localPath string) (string, error) {
	id := fmt.Sprintf("dl-%d", time.Now().UnixNano())
	di := &DownloadInfo{
		ID:         id,
		Status:     DLRunning,
		RemotePath: remotePath,
		LocalPath:  localPath,
	}
	r.downloads.Store(id, di)

	go r.doDownload(ctx, di)
	return id, nil
}

func (r *RemoteExecutor) doDownload(ctx context.Context, di *DownloadInfo) {
	// Get remote file size.
	out, err := r.sshRun(ctx, fmt.Sprintf("stat -c%%s %s 2>/dev/null", di.RemotePath))
	if err == nil {
		di.BytesTotal, _ = strconv.ParseInt(strings.TrimSpace(out), 10, 64)
	}

	if err := os.MkdirAll(filepath.Dir(di.LocalPath), 0755); err != nil {
		di.Status = DLFailed
		di.Error = fmt.Sprintf("mkdir: %v", err)
		if r.OnDownloadDone != nil {
			r.OnDownloadDone(di.ID, di.RemotePath, di.LocalPath, err)
		}
		return
	}

	// Use rsync with --append-verify.
	cmd := fmt.Sprintf("rsync -avz --partial --append-verify -e 'ssh -o StrictHostKeyChecking=no -p %d' root@%s:%s %s",
		r.sshPort, r.sshHost, di.RemotePath, di.LocalPath)

	rsync := newTrackedExec(ctx, "bash", "-c", cmd)

	// Monitor progress in background.
	go func() {
		scanner := bufio.NewScanner(rsync.stderr)
		for scanner.Scan() {
			select {
			case <-ctx.Done():
				return
			default:
			}
			// rsync --progress output parsing could go here.
		}
	}()

	if err := rsync.cmd.Wait(); err != nil {
		di.Status = DLFailed
		di.Error = err.Error()
		if r.OnDownloadDone != nil {
			r.OnDownloadDone(di.ID, di.RemotePath, di.LocalPath, err)
		}
		return
	}

	// Verify download.
	stat, err := os.Stat(di.LocalPath)
	if err != nil {
		di.Status = DLFailed
		di.Error = fmt.Sprintf("stat local: %v", err)
		if r.OnDownloadDone != nil {
			r.OnDownloadDone(di.ID, di.RemotePath, di.LocalPath, err)
		}
		return
	}

	di.BytesDone = stat.Size()
	if di.BytesTotal > 0 && di.BytesDone != di.BytesTotal {
		di.Status = DLFailed
		di.Error = fmt.Sprintf("size mismatch: got %d, want %d", di.BytesDone, di.BytesTotal)
		if r.OnDownloadDone != nil {
			r.OnDownloadDone(di.ID, di.RemotePath, di.LocalPath, fmt.Errorf("%s", di.Error))
		}
		return
	}
	di.Status = DLComplete
	if r.OnDownloadDone != nil {
		r.OnDownloadDone(di.ID, di.RemotePath, di.LocalPath, nil)
	}
}

func (r *RemoteExecutor) DownloadStatus(id string) *DownloadInfo {
	v, ok := r.downloads.Load(id)
	if !ok {
		return nil
	}
	return v.(*DownloadInfo)
}

func (r *RemoteExecutor) ListFiles(ctx context.Context, pattern string) ([]FileInfo, error) {
	out, err := r.sshRun(ctx, fmt.Sprintf("stat -c '%%n %%s %%Y' %s 2>/dev/null", pattern))
	if err != nil {
		return nil, fmt.Errorf("list files: %w", err)
	}

	var files []FileInfo
	for _, line := range strings.Split(strings.TrimSpace(out), "\n") {
		fields := strings.Fields(line)
		if len(fields) < 3 {
			continue
		}
		size, _ := strconv.ParseInt(fields[1], 10, 64)
		epoch, _ := strconv.ParseInt(fields[2], 10, 64)
		files = append(files, FileInfo{
			Path:    fields[0],
			Size:    size,
			ModTime: time.Unix(epoch, 0).Format(time.RFC3339),
		})
	}
	return files, nil
}

func (r *RemoteExecutor) Exec(ctx context.Context, cmd string, timeout time.Duration) (string, error) {
	if timeout > 0 {
		var cancel context.CancelFunc
		ctx, cancel = context.WithTimeout(ctx, timeout)
		defer cancel()
	}
	return r.sshRun(ctx, cmd)
}

// loadSSHKey reads the user's SSH private key.
func loadSSHKey() (ssh.Signer, error) {
	home, err := os.UserHomeDir()
	if err != nil {
		return nil, fmt.Errorf("home dir: %w", err)
	}

	keyFiles := []string{
		filepath.Join(home, ".ssh", "id_ed25519"),
		filepath.Join(home, ".ssh", "id_rsa"),
	}

	for _, kf := range keyFiles {
		data, err := os.ReadFile(kf)
		if err != nil {
			continue
		}
		signer, err := ssh.ParsePrivateKey(data)
		if err != nil {
			continue
		}
		return signer, nil
	}
	return nil, fmt.Errorf("no SSH key found in %v", keyFiles)
}

// trackedExec wraps an exec.Cmd with accessible stderr.
type trackedExec struct {
	cmd    *exec.Cmd
	stderr io.ReadCloser
}

func newTrackedExec(ctx context.Context, name string, args ...string) *trackedExec {
	cmd := exec.CommandContext(ctx, name, args...)
	stderr, _ := cmd.StderrPipe()
	cmd.Start()
	return &trackedExec{cmd: cmd, stderr: stderr}
}
