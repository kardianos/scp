package main

import (
	"context"
	"time"
)

// ExecType identifies local vs remote execution.
type ExecType string

const (
	ExecLocal  ExecType = "local"
	ExecRemote ExecType = "remote"
)

// RunState tracks simulation lifecycle.
type RunState string

const (
	RunStarting  RunState = "starting"
	RunRunning   RunState = "running"
	RunComplete  RunState = "complete"
	RunFailed    RunState = "failed"
	RunCancelled RunState = "cancelled"
)

// DownloadState tracks file transfer lifecycle.
type DownloadState string

const (
	DLRunning  DownloadState = "running"
	DLComplete DownloadState = "complete"
	DLFailed   DownloadState = "failed"
)

// ExecutorStatus holds system resource information.
type ExecutorStatus struct {
	Type          ExecType `json:"type"`
	Location      string   `json:"location"`
	GPUName       string   `json:"gpu_name,omitempty"`
	GPUUtil       int      `json:"gpu_util"`
	GPUMemUsedMB  int64    `json:"gpu_mem_used_mb"`
	GPUMemTotalMB int64    `json:"gpu_mem_total_mb"`
	CPUUtil       float64  `json:"cpu_util"`
	RAMUsedMB     int64    `json:"ram_used_mb"`
	RAMTotalMB    int64    `json:"ram_total_mb"`
	DiskUsedGB    float64  `json:"disk_used_gb"`
	DiskFreeGB    float64  `json:"disk_free_gb"`
	ActiveRuns    int      `json:"active_runs"`
	Uptime        string   `json:"uptime"`
}

// RunInfo holds the state of a single simulation run.
type RunInfo struct {
	ID          string   `json:"id"`
	Status      RunState `json:"status"`
	SimTime     float64  `json:"sim_time"`
	TotalTime   float64  `json:"total_time"`
	WallSecs    float64  `json:"wall_seconds"`
	LastDiag    string   `json:"last_diag,omitempty"`
	Error       string   `json:"error,omitempty"`
	OutputFiles []string `json:"output_files,omitempty"`
}

// DownloadInfo holds the state of a file download.
type DownloadInfo struct {
	ID         string        `json:"id"`
	Status     DownloadState `json:"status"`
	RemotePath string        `json:"remote_path"`
	LocalPath  string        `json:"local_path"`
	BytesTotal int64         `json:"bytes_total"`
	BytesDone  int64         `json:"bytes_done"`
	Error      string        `json:"error,omitempty"`
}

// BuildResult holds the outcome of a compilation.
type BuildResult struct {
	Status string `json:"status"` // "built", "cached", "failed"
	Binary string `json:"binary"`
	Cached bool   `json:"cached"`
	Error  string `json:"error,omitempty"`
}

// FileInfo describes a file on disk.
type FileInfo struct {
	Path    string `json:"path"`
	Size    int64  `json:"size"`
	ModTime string `json:"mod_time"`
}

// Executor is the interface for both local and remote execution environments.
type Executor interface {
	Type() ExecType
	Setup(ctx context.Context) error
	Teardown(ctx context.Context) error
	Status(ctx context.Context) (*ExecutorStatus, error)

	Build(ctx context.Context, sources []string, cmd string) (*BuildResult, error)

	Run(ctx context.Context, config string, id string) error
	RunStatus(id string) *RunInfo
	RunCancel(id string) error

	Upload(ctx context.Context, localPath, remotePath string) error
	Download(ctx context.Context, remotePath, localPath string) (string, error)
	DownloadStatus(id string) *DownloadInfo
	ListFiles(ctx context.Context, pattern string) ([]FileInfo, error)

	Exec(ctx context.Context, cmd string, timeout time.Duration) (string, error)
}
