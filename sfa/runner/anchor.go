package main

// anchor.go implements the pull/file-based completion signaling layer.
// MCP push notifications are unreliable in some client harnesses, so the
// runner mirrors run/instance/liveness state into ~/.scp-runner/state/ where
// external file-watchers can observe it:
//
//   ~/.scp-runner/state/heartbeat.log            appended every 30s (runner aliveness)
//   ~/.scp-runner/state/heartbeat.json           atomic rewrite of the same info
//   ~/.scp-runner/state/<instance>/instance.json remote reachability {reachable, last_contact, ...}
//   ~/.scp-runner/state/<instance>/<run>.status.json  per-run status, updated on every poll
//   ~/.scp-runner/state/<instance>/<run>.done    empty marker, atomic create on completion
//   ~/.scp-runner/state/<instance>/<run>.failed  empty marker, atomic create on failure/cancel

import (
	"context"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"
)

// runnerStateRoot returns ~/.scp-runner/state/, creating it if necessary.
func runnerStateRoot() (string, error) {
	dir, err := stateDir()
	if err != nil {
		return "", err
	}
	root := filepath.Join(dir, "state")
	if err := os.MkdirAll(root, 0755); err != nil {
		return "", fmt.Errorf("mkdir %s: %w", root, err)
	}
	return root, nil
}

// instanceStateDir returns ~/.scp-runner/state/<instance>/, creating it.
func instanceStateDir(instance string) (string, error) {
	root, err := runnerStateRoot()
	if err != nil {
		return "", err
	}
	dir := filepath.Join(root, instance)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return "", fmt.Errorf("mkdir %s: %w", dir, err)
	}
	return dir, nil
}

// writeJSONAtomic writes v as JSON to path via a tmp file + rename.
func writeJSONAtomic(path string, v any) error {
	data, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return fmt.Errorf("marshal: %w", err)
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0644); err != nil {
		return fmt.Errorf("write %s: %w", tmp, err)
	}
	if err := os.Rename(tmp, path); err != nil {
		return fmt.Errorf("rename: %w", err)
	}
	return nil
}

// RunStatusFile is the per-run JSON status written for external watchers.
type RunStatusFile struct {
	State       string  `json:"state"`
	Phase       string  `json:"phase,omitempty"`
	SimTime     float64 `json:"sim_time"`
	TotalTime   float64 `json:"total_time"`
	LastUpdate  string  `json:"last_update"`
	WallSeconds float64 `json:"wall_seconds"`
}

// writeRunStatusFile writes <instance>/<runID>.status.json (atomic).
func writeRunStatusFile(instance, runID string, info *RunInfo) {
	dir, err := instanceStateDir(instance)
	if err != nil {
		log.Printf("anchor: %v", err)
		return
	}
	sf := RunStatusFile{
		State:       string(info.Status),
		Phase:       info.Phase,
		SimTime:     info.SimTime,
		TotalTime:   info.TotalTime,
		LastUpdate:  time.Now().Format(time.RFC3339),
		WallSeconds: info.WallSecs,
	}
	if err := writeJSONAtomic(filepath.Join(dir, runID+".status.json"), &sf); err != nil {
		log.Printf("anchor: run status file: %v", err)
	}
}

// writeRunMarker atomically creates the empty <runID>.done or <runID>.failed
// marker file. This is what external file-watchers watch for. Idempotent:
// an already-existing marker is left alone.
func writeRunMarker(instance, runID string, status RunState) {
	dir, err := instanceStateDir(instance)
	if err != nil {
		log.Printf("anchor: %v", err)
		return
	}
	name := runID + ".failed"
	if status == RunComplete {
		name = runID + ".done"
	}
	f, err := os.OpenFile(filepath.Join(dir, name), os.O_CREATE|os.O_EXCL|os.O_WRONLY, 0644)
	if err != nil {
		if !os.IsExist(err) {
			log.Printf("anchor: marker %s: %v", name, err)
		}
		return
	}
	f.Close()
}

// clearRunMarkers removes stale .done/.failed markers and status file for a
// run id, so re-running the same id does not trip watchers immediately.
func clearRunMarkers(instance, runID string) {
	dir, err := instanceStateDir(instance)
	if err != nil {
		return
	}
	os.Remove(filepath.Join(dir, runID+".done"))
	os.Remove(filepath.Join(dir, runID+".failed"))
	os.Remove(filepath.Join(dir, runID+".status.json"))
}

// InstanceLiveness is written to <instance>/instance.json for external watchers.
type InstanceLiveness struct {
	Reachable   bool   `json:"reachable"`
	Degraded    bool   `json:"degraded"`
	LastContact string `json:"last_contact,omitempty"`
	Host        string `json:"host,omitempty"`
	Port        int    `json:"port,omitempty"`
	MachineID   int    `json:"machine_id,omitempty"`
}

// writeInstanceFile writes <instance>/instance.json (atomic).
func writeInstanceFile(instance string, il *InstanceLiveness) {
	dir, err := instanceStateDir(instance)
	if err != nil {
		log.Printf("anchor: %v", err)
		return
	}
	if err := writeJSONAtomic(filepath.Join(dir, "instance.json"), il); err != nil {
		log.Printf("anchor: instance file: %v", err)
	}
}

// heartbeatJSON is the payload of ~/.scp-runner/state/heartbeat.json.
type heartbeatJSON struct {
	Time       string   `json:"time"`
	Instances  []string `json:"instances"`
	ActiveRuns int      `json:"active_runs"`
}

// writeHeartbeat appends one line to heartbeat.log and atomically rewrites
// heartbeat.json. External watchers can detect a dead runner by file age.
func writeHeartbeat(instances []string, activeRuns int) {
	root, err := runnerStateRoot()
	if err != nil {
		log.Printf("anchor: %v", err)
		return
	}
	now := time.Now().Format(time.RFC3339)

	line := fmt.Sprintf("%s instances=%s active_runs=%d\n", now, strings.Join(instances, ","), activeRuns)
	f, err := os.OpenFile(filepath.Join(root, "heartbeat.log"), os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0644)
	if err == nil {
		if _, werr := f.WriteString(line); werr != nil {
			log.Printf("anchor: heartbeat log: %v", werr)
		}
		f.Close()
	} else {
		log.Printf("anchor: heartbeat log: %v", err)
	}

	hb := heartbeatJSON{Time: now, Instances: instances, ActiveRuns: activeRuns}
	if hb.Instances == nil {
		hb.Instances = []string{}
	}
	if err := writeJSONAtomic(filepath.Join(root, "heartbeat.json"), &hb); err != nil {
		log.Printf("anchor: heartbeat json: %v", err)
	}
}

// runOnCompleteCommand executes a user-supplied on_complete hook via bash -c
// with SCP_RUN_ID / SCP_STATE / SCP_INSTANCE in the environment.
func runOnCompleteCommand(cmdStr, instance, runID string, state RunState) {
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Minute)
	defer cancel()
	cmd := exec.CommandContext(ctx, "bash", "-c", cmdStr)
	cmd.Env = append(os.Environ(),
		"SCP_RUN_ID="+runID,
		"SCP_STATE="+string(state),
		"SCP_INSTANCE="+instance,
	)
	out, err := cmd.CombinedOutput()
	if err != nil {
		log.Printf("anchor: on_complete for %s/%s failed: %v: %s", instance, runID, err, strings.TrimSpace(string(out)))
	} else {
		log.Printf("anchor: on_complete for %s/%s ok", instance, runID)
	}
}
