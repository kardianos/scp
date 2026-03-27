package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sync"
)

// RunnerState is the top-level persistent state written to disk.
type RunnerState struct {
	Instance *InstanceState          `json:"instance,omitempty"`
	Binary   string                  `json:"binary,omitempty"`
	Runs     map[string]*PersistRun  `json:"runs,omitempty"`
	mu       sync.Mutex
}

// InstanceState records the remote instance connection info.
type InstanceState struct {
	ID      int    `json:"id"`
	Host    string `json:"host"`
	Port    int    `json:"port"`
	GPUName string `json:"gpu_name"`
}

// PersistRun records the persistent state of a single simulation run.
type PersistRun struct {
	Status       string   `json:"status"`
	AutoDownload string   `json:"auto_download,omitempty"`
	OutputFiles  []string `json:"output_files,omitempty"`
	RemoteDir    string   `json:"remote_dir,omitempty"`
}

// stateDir returns ~/.scp-runner/, creating it if necessary.
func stateDir() (string, error) {
	home, err := os.UserHomeDir()
	if err != nil {
		return "", fmt.Errorf("home dir: %w", err)
	}
	dir := filepath.Join(home, ".scp-runner")
	if err := os.MkdirAll(dir, 0755); err != nil {
		return "", fmt.Errorf("mkdir %s: %w", dir, err)
	}
	return dir, nil
}

// statePath returns the full path to the state file.
func statePath() (string, error) {
	dir, err := stateDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, "state.json"), nil
}

// LoadState reads the state file from disk. Returns an empty state if the file
// does not exist. Returns an error only on actual I/O or parse failures.
func LoadState() (*RunnerState, error) {
	path, err := statePath()
	if err != nil {
		return &RunnerState{Runs: make(map[string]*PersistRun)}, nil
	}

	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return &RunnerState{Runs: make(map[string]*PersistRun)}, nil
		}
		return nil, fmt.Errorf("read state: %w", err)
	}

	var s RunnerState
	if err := json.Unmarshal(data, &s); err != nil {
		return nil, fmt.Errorf("parse state: %w", err)
	}
	if s.Runs == nil {
		s.Runs = make(map[string]*PersistRun)
	}
	return &s, nil
}

// SaveState writes the state to disk atomically (write .tmp then rename).
func SaveState(s *RunnerState) error {
	s.mu.Lock()
	defer s.mu.Unlock()
	return saveStateLocked(s)
}

// saveStateLocked writes state without acquiring the mutex (caller holds it).
func saveStateLocked(s *RunnerState) error {
	path, err := statePath()
	if err != nil {
		return err
	}

	data, err := json.MarshalIndent(s, "", "  ")
	if err != nil {
		return fmt.Errorf("marshal state: %w", err)
	}

	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0644); err != nil {
		return fmt.Errorf("write tmp: %w", err)
	}
	if err := os.Rename(tmp, path); err != nil {
		return fmt.Errorf("rename: %w", err)
	}
	return nil
}

// SetInstance updates the instance info and saves.
func (s *RunnerState) SetInstance(inst *InstanceState) {
	s.mu.Lock()
	s.Instance = inst
	saveStateLocked(s)
	s.mu.Unlock()
}

// SetBinary updates the binary path and saves.
func (s *RunnerState) SetBinary(binary string) {
	s.mu.Lock()
	s.Binary = binary
	saveStateLocked(s)
	s.mu.Unlock()
}

// SetRun updates (or creates) a run entry and saves.
func (s *RunnerState) SetRun(id string, run *PersistRun) {
	s.mu.Lock()
	s.Runs[id] = run
	saveStateLocked(s)
	s.mu.Unlock()
}

// UpdateRunStatus updates only the status field of a run and saves.
func (s *RunnerState) UpdateRunStatus(id string, status string) {
	s.mu.Lock()
	if r, ok := s.Runs[id]; ok {
		r.Status = status
	}
	saveStateLocked(s)
	s.mu.Unlock()
}

// ClearInstance removes instance info and saves.
func (s *RunnerState) ClearInstance() {
	s.mu.Lock()
	s.Instance = nil
	s.Binary = ""
	s.Runs = make(map[string]*PersistRun)
	saveStateLocked(s)
	s.mu.Unlock()
}
