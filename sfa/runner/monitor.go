package main

import (
	"context"
	"encoding/json"
	"log"
	"os"
	"path/filepath"
	"sync"
	"time"
)

// StatusState holds the latest snapshot written to disk.
type StatusState struct {
	Executor  *ExecutorStatus        `json:"executor,omitempty"`
	Runs      map[string]*RunInfo    `json:"runs,omitempty"`
	Downloads map[string]*DownloadInfo `json:"downloads,omitempty"`
	UpdatedAt string                 `json:"updated_at"`
}

// Monitor writes periodic status snapshots and manages background tasks.
type Monitor struct {
	server  *MCPServer
	dir     string
	mu      sync.Mutex
	latest  StatusState
}

func NewMonitor(server *MCPServer, dir string) *Monitor {
	return &Monitor{
		server: server,
		dir:    dir,
	}
}

// Run starts the monitor loop. It blocks until ctx is cancelled.
func (m *Monitor) Run(ctx context.Context) {
	if err := os.MkdirAll(m.dir, 0755); err != nil {
		log.Printf("monitor: mkdir %s: %v", m.dir, err)
		return
	}

	ticker := time.NewTicker(5 * time.Second)
	defer ticker.Stop()

	for {
		select {
		case <-ctx.Done():
			return
		case <-ticker.C:
			m.snapshot(ctx)
		}
	}
}

func (m *Monitor) snapshot(ctx context.Context) {
	state := StatusState{
		Runs:      make(map[string]*RunInfo),
		Downloads: make(map[string]*DownloadInfo),
		UpdatedAt: time.Now().Format(time.RFC3339),
	}

	ex, err := m.server.getExecutor()
	if err == nil {
		status, err := ex.Status(ctx)
		if err == nil {
			state.Executor = status
		}

		// Collect run states.
		collectRuns(ex, state.Runs)

		// Collect download states.
		collectDownloads(ex, state.Downloads)
	}

	data, err := json.MarshalIndent(state, "", "  ")
	if err != nil {
		log.Printf("monitor: marshal: %v", err)
		return
	}

	path := filepath.Join(m.dir, "status.json")
	if err := os.WriteFile(path, data, 0644); err != nil {
		log.Printf("monitor: write %s: %v", path, err)
	}

	m.mu.Lock()
	m.latest = state
	m.mu.Unlock()
}

// collectRuns gathers run info from the executor's internal maps.
// Since the Executor interface doesn't expose iteration, we check known IDs.
// In practice, the MCP handlers track IDs and this is a best-effort snapshot.
func collectRuns(ex Executor, runs map[string]*RunInfo) {
	switch e := ex.(type) {
	case *LocalExecutor:
		e.runs.Range(func(key, value any) bool {
			id := key.(string)
			lr := value.(*localRun)
			lr.mu.Lock()
			info := lr.info
			lr.mu.Unlock()
			runs[id] = &info
			return true
		})
	case *RemoteExecutor:
		e.runs.Range(func(key, value any) bool {
			id := key.(string)
			rr := value.(*remoteRun)
			rr.mu.Lock()
			info := rr.info
			rr.mu.Unlock()
			runs[id] = &info
			return true
		})
	}
}

func collectDownloads(ex Executor, downloads map[string]*DownloadInfo) {
	switch e := ex.(type) {
	case *LocalExecutor:
		e.downloads.Range(func(key, value any) bool {
			id := key.(string)
			di := value.(*DownloadInfo)
			downloads[id] = di
			return true
		})
	case *RemoteExecutor:
		e.downloads.Range(func(key, value any) bool {
			id := key.(string)
			di := value.(*DownloadInfo)
			downloads[id] = di
			return true
		})
	}
}
