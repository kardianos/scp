package main

import (
	"context"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"sync"
	"time"
)

// StatusState holds the latest snapshot written to disk.
type StatusState struct {
	Instances map[string]*InstanceStatusState `json:"instances,omitempty"`
	UpdatedAt string                          `json:"updated_at"`
}

// InstanceStatusState holds status for a single named instance.
type InstanceStatusState struct {
	Executor  *ExecutorStatus          `json:"executor,omitempty"`
	Runs      map[string]*RunInfo      `json:"runs,omitempty"`
	Downloads map[string]*DownloadInfo `json:"downloads,omitempty"`
}

// ProgressNotification is the payload for notifications/sim_progress.
type ProgressNotification struct {
	Instance  string   `json:"instance"`
	ID        string   `json:"id"`
	Status    RunState `json:"status"`
	SimTime   float64  `json:"sim_time"`
	TotalTime float64  `json:"total_time"`
	WallSecs  float64  `json:"wall_seconds"`
	LastDiag  string   `json:"last_diag,omitempty"`
	Error     string   `json:"error,omitempty"`
	Done      bool     `json:"done"`
}

// notifyState tracks per-run notification timing.
type notifyState struct {
	lastNotify time.Time
	doneSent   bool
}

// notifyKey combines instance name and run ID for notification tracking.
type notifyKey struct {
	instance string
	runID    string
}

// Monitor writes periodic status snapshots and manages background tasks.
type Monitor struct {
	server  *MCPServer
	dir     string
	mu      sync.Mutex
	latest  StatusState
	notify  map[notifyKey]*notifyState // per-run notification tracking
	wakeup  chan struct{}              // signal for immediate notification check
	state   *RunnerState              // persistent state (shared with server)
}

func NewMonitor(server *MCPServer, dir string) *Monitor {
	return &Monitor{
		server: server,
		dir:    dir,
		notify: make(map[notifyKey]*notifyState),
		wakeup: make(chan struct{}, 4),
	}
}

// Wakeup triggers an immediate notification check (non-blocking).
func (m *Monitor) Wakeup() {
	select {
	case m.wakeup <- struct{}{}:
	default:
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
		case <-m.wakeup:
			m.snapshot(ctx)
		}
	}
}

func (m *Monitor) snapshot(ctx context.Context) {
	state := StatusState{
		Instances: make(map[string]*InstanceStatusState),
		UpdatedAt: time.Now().Format(time.RFC3339),
	}

	// Iterate over all executors.
	allExecs := m.server.getAllExecutors()
	for name, ex := range allExecs {
		instState := &InstanceStatusState{
			Runs:      make(map[string]*RunInfo),
			Downloads: make(map[string]*DownloadInfo),
		}

		status, err := ex.Status(ctx)
		if err == nil {
			instState.Executor = status
		}

		// Collect run states and send notifications.
		collectRuns(ex, instState.Runs)
		m.sendRunNotifications(name, ex, instState.Runs)

		// Collect download states.
		collectDownloads(ex, instState.Downloads)

		state.Instances[name] = instState
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

// sendRunNotifications checks each run and emits progress notifications as needed.
func (m *Monitor) sendRunNotifications(instanceName string, ex Executor, runs map[string]*RunInfo) {
	now := time.Now()

	for id, info := range runs {
		interval := getNotifyInterval(ex, id)
		if interval <= 0 {
			continue
		}

		nk := notifyKey{instance: instanceName, runID: id}
		ns, ok := m.notify[nk]
		if !ok {
			ns = &notifyState{}
			m.notify[nk] = ns
		}

		done := info.Status == RunComplete || info.Status == RunFailed || info.Status == RunCancelled

		if done {
			if !ns.doneSent {
				pn := ProgressNotification{
					Instance:  instanceName,
					ID:        id,
					Status:    info.Status,
					SimTime:   info.SimTime,
					TotalTime: info.TotalTime,
					WallSecs:  info.WallSecs,
					LastDiag:  info.LastDiag,
					Error:     info.Error,
					Done:      true,
				}
				m.server.sendNotification("notifications/sim_progress", pn)
				// Also send standard MCP log message so clients surface it.
				level := "info"
				msg := fmt.Sprintf("[%s/%s] %s — sim_time=%.1f/%.1f wall=%.1fs",
					instanceName, id, info.Status, info.SimTime, info.TotalTime, info.WallSecs)
				if info.Status == RunFailed {
					level = "error"
					msg = fmt.Sprintf("[%s/%s] FAILED (%.1fs): %s", instanceName, id, info.WallSecs, info.Error)
				} else if info.Status == RunCancelled {
					level = "warning"
					msg = fmt.Sprintf("[%s/%s] cancelled at sim_time=%.1f/%.1f (%.1fs)",
						instanceName, id, info.SimTime, info.TotalTime, info.WallSecs)
				}
				m.server.sendLogMessage(level, msg)

				// Send channel event so the agent can react immediately.
				switch info.Status {
				case RunComplete:
					m.server.sendChannelEvent("run_complete", map[string]string{
						"instance": instanceName,
						"run_id":   id,
						"sim_time":  fmt.Sprintf("%.1f", info.SimTime),
						"wall_secs": fmt.Sprintf("%.1f", info.WallSecs),
					}, msg)
				case RunFailed:
					m.server.sendChannelEvent("run_failed", map[string]string{
						"instance": instanceName,
						"run_id":   id,
						"error":    info.Error,
					}, msg)
				case RunCancelled:
					m.server.sendChannelEvent("run_cancelled", map[string]string{
						"instance": instanceName,
						"run_id":   id,
					}, msg)
				}

				ns.doneSent = true
				// Persist run completion to state file.
				if m.state != nil {
					m.state.UpdateRunStatus(instanceName, id, string(info.Status))
				}
			}
			continue
		}

		// Running — check if enough time has elapsed.
		if now.Sub(ns.lastNotify) >= interval {
			pn := ProgressNotification{
				Instance:  instanceName,
				ID:        id,
				Status:    info.Status,
				SimTime:   info.SimTime,
				TotalTime: info.TotalTime,
				WallSecs:  info.WallSecs,
				LastDiag:  info.LastDiag,
				Done:      false,
			}
			m.server.sendNotification("notifications/sim_progress", pn)
			// Also send standard MCP log message for progress.
			pct := 0.0
			if info.TotalTime > 0 {
				pct = 100 * info.SimTime / info.TotalTime
			}
			m.server.sendLogMessage("info", fmt.Sprintf("[%s/%s] running — sim_time=%.1f/%.1f (%.0f%%) wall=%.1fs",
				instanceName, id, info.SimTime, info.TotalTime, pct, info.WallSecs))
			ns.lastNotify = now
		}
	}
}

// getNotifyInterval extracts the notify interval from the run's internal struct.
func getNotifyInterval(ex Executor, id string) time.Duration {
	switch e := ex.(type) {
	case *LocalExecutor:
		v, ok := e.runs.Load(id)
		if !ok {
			return 0
		}
		return v.(*localRun).notifyInterval
	case *RemoteExecutor:
		v, ok := e.runs.Load(id)
		if !ok {
			return 0
		}
		return v.(*remoteRun).notifyInterval
	}
	return 0
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
