package main

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

// LocalExecutor runs simulations on the local machine.
type LocalExecutor struct {
	workDir   string
	hasNVCC   bool
	hasGCC    bool
	startTime time.Time

	runs      sync.Map // map[string]*localRun
	downloads sync.Map // map[string]*DownloadInfo
	mu        sync.Mutex
}

type localRun struct {
	info   RunInfo
	cancel context.CancelFunc
	mu     sync.Mutex
}

func NewLocalExecutor(workDir string) *LocalExecutor {
	return &LocalExecutor{
		workDir:   workDir,
		startTime: time.Now(),
	}
}

func (l *LocalExecutor) Type() ExecType { return ExecLocal }

func (l *LocalExecutor) Setup(ctx context.Context) error {
	if err := os.MkdirAll(l.workDir, 0755); err != nil {
		return fmt.Errorf("create workdir: %w", err)
	}
	_, err := exec.LookPath("nvcc")
	l.hasNVCC = err == nil
	_, err = exec.LookPath("gcc")
	l.hasGCC = err == nil
	if !l.hasGCC {
		return fmt.Errorf("gcc not found in PATH")
	}
	return nil
}

func (l *LocalExecutor) Teardown(_ context.Context) error { return nil }

func (l *LocalExecutor) Status(ctx context.Context) (*ExecutorStatus, error) {
	s := &ExecutorStatus{
		Type:     ExecLocal,
		Location: "localhost",
		GPUUtil:  -1,
		Uptime:   time.Since(l.startTime).Round(time.Second).String(),
	}

	// Count active runs.
	l.runs.Range(func(_, v any) bool {
		r := v.(*localRun)
		r.mu.Lock()
		switch r.info.Status {
		case RunStarting, RunRunning:
			s.ActiveRuns++
		}
		r.mu.Unlock()
		return true
	})

	// CPU info.
	s.CPUUtil = -1 // simplified; full /proc/stat parsing omitted for brevity
	s.RAMTotalMB = int64(runtime.NumCPU()) * 1024 // rough estimate

	// GPU info if available.
	if l.hasNVCC {
		out, err := execOutput(ctx, "nvidia-smi", "--query-gpu=name,utilization.gpu,memory.used,memory.total",
			"--format=csv,noheader,nounits")
		if err == nil {
			parts := strings.Split(strings.TrimSpace(out), ", ")
			if len(parts) >= 4 {
				s.GPUName = parts[0]
				s.GPUUtil, _ = strconv.Atoi(strings.TrimSpace(parts[1]))
				s.GPUMemUsedMB, _ = strconv.ParseInt(strings.TrimSpace(parts[2]), 10, 64)
				s.GPUMemTotalMB, _ = strconv.ParseInt(strings.TrimSpace(parts[3]), 10, 64)
			}
		}
	}

	// Disk info.
	out, err := execOutput(ctx, "df", "-BG", "--output=used,avail", l.workDir)
	if err == nil {
		lines := strings.Split(strings.TrimSpace(out), "\n")
		if len(lines) >= 2 {
			fields := strings.Fields(lines[len(lines)-1])
			if len(fields) >= 2 {
				used, _ := strconv.ParseFloat(strings.TrimSuffix(fields[0], "G"), 64)
				avail, _ := strconv.ParseFloat(strings.TrimSuffix(fields[1], "G"), 64)
				s.DiskUsedGB = used
				s.DiskFreeGB = avail
			}
		}
	}

	return s, nil
}

func (l *LocalExecutor) Build(ctx context.Context, sources []string, cmd string) (*BuildResult, error) {
	hash, err := hashSources(sources)
	if err != nil {
		return nil, fmt.Errorf("hash sources: %w", err)
	}

	// Determine output binary name from cmd or derive one.
	binName := fmt.Sprintf("scp_sim_%s", hash)
	binPath := filepath.Join(l.workDir, binName)

	// Check cache.
	if _, err := os.Stat(binPath); err == nil {
		return &BuildResult{Status: "cached", Binary: binPath, Cached: true}, nil
	}

	// Build.
	if cmd == "" {
		cmd = fmt.Sprintf("gcc -O3 -march=native -fopenmp -o %s %s -lzstd -lm",
			binPath, strings.Join(sources, " "))
	} else {
		// Replace placeholder output path if present.
		cmd = strings.ReplaceAll(cmd, "${OUTPUT}", binPath)
	}

	out, err := execOutput(ctx, "bash", "-c", cmd)
	if err != nil {
		return &BuildResult{Status: "failed", Error: fmt.Sprintf("%v: %s", err, out)}, nil
	}

	return &BuildResult{Status: "built", Binary: binPath}, nil
}

func (l *LocalExecutor) Run(ctx context.Context, config string, id string) error {
	if _, loaded := l.runs.Load(id); loaded {
		return fmt.Errorf("run %s already exists", id)
	}

	runCtx, cancel := context.WithCancel(ctx)
	lr := &localRun{
		info: RunInfo{
			ID:     id,
			Status: RunStarting,
		},
		cancel: cancel,
	}
	l.runs.Store(id, lr)

	go l.executeRun(runCtx, lr, config)
	return nil
}

func (l *LocalExecutor) executeRun(ctx context.Context, lr *localRun, config string) {
	defer lr.cancel()
	startWall := time.Now()

	// Parse config to find binary and args.
	// Expected format: binary_path [args...] or binary_path -c config_file
	parts := strings.Fields(config)
	if len(parts) == 0 {
		lr.mu.Lock()
		lr.info.Status = RunFailed
		lr.info.Error = "empty config"
		lr.mu.Unlock()
		return
	}

	cmd := exec.CommandContext(ctx, parts[0], parts[1:]...)
	cmd.Dir = l.workDir

	stderr, err := cmd.StderrPipe()
	if err != nil {
		lr.mu.Lock()
		lr.info.Status = RunFailed
		lr.info.Error = fmt.Sprintf("stderr pipe: %v", err)
		lr.mu.Unlock()
		return
	}

	if err := cmd.Start(); err != nil {
		lr.mu.Lock()
		lr.info.Status = RunFailed
		lr.info.Error = fmt.Sprintf("start: %v", err)
		lr.mu.Unlock()
		return
	}

	lr.mu.Lock()
	lr.info.Status = RunRunning
	lr.mu.Unlock()

	// Tail stderr for diagnostic output.
	go l.tailDiag(ctx, lr, stderr)

	// Wait for process.
	err = cmd.Wait()
	wall := time.Since(startWall).Seconds()

	lr.mu.Lock()
	defer lr.mu.Unlock()
	lr.info.WallSecs = wall

	switch {
	case ctx.Err() != nil:
		lr.info.Status = RunCancelled
	case err != nil:
		lr.info.Status = RunFailed
		lr.info.Error = err.Error()
	default:
		lr.info.Status = RunComplete
	}

	// Gather output files.
	matches, _ := filepath.Glob(filepath.Join(l.workDir, "*.sfa"))
	lr.info.OutputFiles = matches
}

func (l *LocalExecutor) tailDiag(ctx context.Context, lr *localRun, r io.Reader) {
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		select {
		case <-ctx.Done():
			return
		default:
		}
		line := scanner.Text()
		lr.mu.Lock()
		lr.info.LastDiag = line
		// Parse sim_time from first tab-separated field.
		if fields := strings.Fields(line); len(fields) > 0 {
			if t, err := strconv.ParseFloat(fields[0], 64); err == nil {
				lr.info.SimTime = t
			}
		}
		lr.mu.Unlock()
	}
}

func (l *LocalExecutor) RunStatus(id string) *RunInfo {
	v, ok := l.runs.Load(id)
	if !ok {
		return nil
	}
	lr := v.(*localRun)
	lr.mu.Lock()
	defer lr.mu.Unlock()
	info := lr.info
	return &info
}

func (l *LocalExecutor) RunCancel(id string) error {
	v, ok := l.runs.Load(id)
	if !ok {
		return fmt.Errorf("run %s not found", id)
	}
	lr := v.(*localRun)
	lr.cancel()
	return nil
}

func (l *LocalExecutor) Upload(_ context.Context, localPath, remotePath string) error {
	// Local executor: copy file to workdir.
	dst := filepath.Join(l.workDir, filepath.Base(remotePath))
	return copyFile(localPath, dst)
}

func (l *LocalExecutor) Download(_ context.Context, remotePath, localPath string) (string, error) {
	id := fmt.Sprintf("dl-%d", time.Now().UnixNano())
	di := &DownloadInfo{
		ID:         id,
		Status:     DLRunning,
		RemotePath: remotePath,
		LocalPath:  localPath,
	}

	stat, err := os.Stat(remotePath)
	if err != nil {
		di.Status = DLFailed
		di.Error = err.Error()
		l.downloads.Store(id, di)
		return id, fmt.Errorf("stat %s: %w", remotePath, err)
	}
	di.BytesTotal = stat.Size()

	if err := copyFile(remotePath, localPath); err != nil {
		di.Status = DLFailed
		di.Error = err.Error()
		l.downloads.Store(id, di)
		return id, fmt.Errorf("copy: %w", err)
	}

	di.BytesDone = di.BytesTotal
	di.Status = DLComplete
	l.downloads.Store(id, di)
	return id, nil
}

func (l *LocalExecutor) DownloadStatus(id string) *DownloadInfo {
	v, ok := l.downloads.Load(id)
	if !ok {
		return nil
	}
	di := v.(*DownloadInfo)
	return di
}

func (l *LocalExecutor) ListFiles(ctx context.Context, pattern string) ([]FileInfo, error) {
	if !filepath.IsAbs(pattern) {
		pattern = filepath.Join(l.workDir, pattern)
	}
	matches, err := filepath.Glob(pattern)
	if err != nil {
		return nil, fmt.Errorf("glob: %w", err)
	}
	var files []FileInfo
	for _, m := range matches {
		stat, err := os.Stat(m)
		if err != nil {
			continue
		}
		if stat.IsDir() {
			continue
		}
		files = append(files, FileInfo{
			Path:    m,
			Size:    stat.Size(),
			ModTime: stat.ModTime().Format(time.RFC3339),
		})
	}
	return files, nil
}

func (l *LocalExecutor) Exec(ctx context.Context, cmd string, timeout time.Duration) (string, error) {
	if timeout > 0 {
		var cancel context.CancelFunc
		ctx, cancel = context.WithTimeout(ctx, timeout)
		defer cancel()
	}
	return execOutput(ctx, "bash", "-c", cmd)
}

// execOutput runs a command and returns combined output.
func execOutput(ctx context.Context, name string, args ...string) (string, error) {
	cmd := exec.CommandContext(ctx, name, args...)
	out, err := cmd.CombinedOutput()
	return string(out), err
}

// copyFile copies src to dst.
func copyFile(src, dst string) error {
	if err := os.MkdirAll(filepath.Dir(dst), 0755); err != nil {
		return fmt.Errorf("mkdir: %w", err)
	}
	in, err := os.Open(src)
	if err != nil {
		return fmt.Errorf("open src: %w", err)
	}
	defer in.Close()

	out, err := os.Create(dst)
	if err != nil {
		return fmt.Errorf("create dst: %w", err)
	}
	defer out.Close()

	if _, err := io.Copy(out, in); err != nil {
		return fmt.Errorf("copy: %w", err)
	}
	return out.Close()
}
