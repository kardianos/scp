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

	runs      sync.Map // map[string]*remoteRun
	downloads sync.Map // map[string]*DownloadInfo
	mu        sync.Mutex
}

type remoteRun struct {
	info   RunInfo
	cancel context.CancelFunc
	mu     sync.Mutex
}

// RemoteConfig holds parameters for creating a remote executor.
type RemoteConfig struct {
	GPUFilter string // e.g. "gpu_name=Tesla_V100 num_gpus=1 rentable=true disk_space>=20"
	Image     string // e.g. "nvidia/cuda:12.2.0-devel-ubuntu22.04"
	DiskGB    int
	Onstart   string // startup script
	WorkDir   string // local directory for downloaded results
}

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

// Connect establishes SSH to an existing instance.
func (r *RemoteExecutor) Connect(ctx context.Context, host string, port int) error {
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

	// Retry connection with context awareness.
	var client *ssh.Client
	for attempt := 0; attempt < 30; attempt++ {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
		}
		client, err = ssh.Dial("tcp", addr, config)
		if err == nil {
			break
		}
		time.Sleep(5 * time.Second)
	}
	if err != nil {
		return fmt.Errorf("ssh dial %s: %w", addr, err)
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
	hash, err := hashSources(sources)
	if err != nil {
		return nil, fmt.Errorf("hash sources: %w", err)
	}

	binName := fmt.Sprintf("scp_sim_%s", hash)
	remoteBin := "/root/" + binName

	// Check remote cache.
	_, err = r.sshRun(ctx, fmt.Sprintf("test -f %s", remoteBin))
	if err == nil {
		return &BuildResult{Status: "cached", Binary: remoteBin, Cached: true}, nil
	}

	// Upload sources.
	for _, src := range sources {
		remotePath := "/root/" + filepath.Base(src)
		if err := r.Upload(ctx, src, remotePath); err != nil {
			return nil, fmt.Errorf("upload %s: %w", src, err)
		}
	}

	// Build.
	if cmd == "" {
		remoteNames := make([]string, len(sources))
		for i, s := range sources {
			remoteNames[i] = "/root/" + filepath.Base(s)
		}
		cmd = fmt.Sprintf("nvcc -O3 -arch=sm_70 -o %s %s -lzstd -lm",
			remoteBin, strings.Join(remoteNames, " "))
	} else {
		cmd = strings.ReplaceAll(cmd, "${OUTPUT}", remoteBin)
	}

	out, err := r.sshRun(ctx, cmd)
	if err != nil {
		return &BuildResult{Status: "failed", Error: fmt.Sprintf("%v: %s", err, out)}, nil
	}

	return &BuildResult{Status: "built", Binary: remoteBin}, nil
}

func (r *RemoteExecutor) Run(ctx context.Context, config string, id string) error {
	if _, loaded := r.runs.Load(id); loaded {
		return fmt.Errorf("run %s already exists", id)
	}

	runCtx, cancel := context.WithCancel(ctx)
	rr := &remoteRun{
		info: RunInfo{
			ID:     id,
			Status: RunStarting,
		},
		cancel: cancel,
	}
	r.runs.Store(id, rr)

	go r.executeRemoteRun(runCtx, rr, config)
	return nil
}

func (r *RemoteExecutor) executeRemoteRun(ctx context.Context, rr *remoteRun, config string) {
	defer rr.cancel()
	startWall := time.Now()

	// Start simulation with nohup so it survives SSH drops.
	logFile := fmt.Sprintf("/root/run_%s.log", rr.info.ID)
	cmd := fmt.Sprintf("nohup bash -c '%s' > %s 2>&1 & echo $!", config, logFile)

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

		// Read last diag line.
		diagOut, _ := r.sshRun(ctx, fmt.Sprintf("tail -1 %s 2>/dev/null", logFile))
		diagLine := strings.TrimSpace(diagOut)

		rr.mu.Lock()
		rr.info.WallSecs = time.Since(startWall).Seconds()
		if diagLine != "" {
			rr.info.LastDiag = diagLine
			if fields := strings.Fields(diagLine); len(fields) > 0 {
				if t, err := strconv.ParseFloat(fields[0], 64); err == nil {
					rr.info.SimTime = t
				}
			}
		}

		if !alive {
			// Check exit status.
			exitOut, _ := r.sshRun(ctx, fmt.Sprintf("wait %s 2>/dev/null; echo $?", pid))
			exitCode := strings.TrimSpace(exitOut)
			switch exitCode {
			case "0", "":
				rr.info.Status = RunComplete
			default:
				rr.info.Status = RunFailed
				rr.info.Error = fmt.Sprintf("exit code %s", exitCode)
			}
			rr.mu.Unlock()
			return
		}
		rr.mu.Unlock()
	}
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
		return
	}

	// Verify download.
	stat, err := os.Stat(di.LocalPath)
	if err != nil {
		di.Status = DLFailed
		di.Error = fmt.Sprintf("stat local: %v", err)
		return
	}

	di.BytesDone = stat.Size()
	if di.BytesTotal > 0 && di.BytesDone != di.BytesTotal {
		di.Status = DLFailed
		di.Error = fmt.Sprintf("size mismatch: got %d, want %d", di.BytesDone, di.BytesTotal)
		return
	}
	di.Status = DLComplete
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
