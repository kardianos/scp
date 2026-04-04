package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sync"
)

// RunnerState is the top-level persistent state written to disk.
// It holds a map of named instances, each with their own binary, runs, etc.
type RunnerState struct {
	Instances map[string]*PersistedInstance `json:"instances,omitempty"`
	mu        sync.Mutex
}

// PersistedInstance records the persistent state of a single named instance.
type PersistedInstance struct {
	Connection *InstanceState         `json:"connection,omitempty"`
	Binary     string                 `json:"binary,omitempty"`
	Runs       map[string]*PersistRun `json:"runs,omitempty"`
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
		return newEmptyState(), nil
	}

	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return newEmptyState(), nil
		}
		return nil, fmt.Errorf("read state: %w", err)
	}

	var s RunnerState
	if err := json.Unmarshal(data, &s); err != nil {
		return nil, fmt.Errorf("parse state: %w", err)
	}
	if s.Instances == nil {
		s.Instances = make(map[string]*PersistedInstance)
	}
	// Ensure each instance has initialized maps.
	for _, inst := range s.Instances {
		if inst.Runs == nil {
			inst.Runs = make(map[string]*PersistRun)
		}
	}
	return &s, nil
}

func newEmptyState() *RunnerState {
	return &RunnerState{
		Instances: make(map[string]*PersistedInstance),
	}
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

// GetInstance returns the persisted instance state for a given name, or nil.
func (s *RunnerState) GetInstance(name string) *PersistedInstance {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.Instances[name]
}

// SetInstance updates (or creates) the instance connection info for a name and saves.
func (s *RunnerState) SetInstance(name string, inst *InstanceState) {
	s.mu.Lock()
	pi := s.getOrCreateInstance(name)
	pi.Connection = inst
	saveStateLocked(s)
	s.mu.Unlock()
}

// SetBinary updates the binary path for a named instance and saves.
func (s *RunnerState) SetBinary(name string, binary string) {
	s.mu.Lock()
	pi := s.getOrCreateInstance(name)
	pi.Binary = binary
	saveStateLocked(s)
	s.mu.Unlock()
}

// SetRun updates (or creates) a run entry for a named instance and saves.
func (s *RunnerState) SetRun(name string, id string, run *PersistRun) {
	s.mu.Lock()
	pi := s.getOrCreateInstance(name)
	pi.Runs[id] = run
	saveStateLocked(s)
	s.mu.Unlock()
}

// UpdateRunStatus updates only the status field of a run for a named instance and saves.
func (s *RunnerState) UpdateRunStatus(name string, id string, status string) {
	s.mu.Lock()
	if pi, ok := s.Instances[name]; ok {
		if r, ok := pi.Runs[id]; ok {
			r.Status = status
		}
	}
	saveStateLocked(s)
	s.mu.Unlock()
}

// ClearInstance removes instance info for a named instance and saves.
func (s *RunnerState) ClearInstance(name string) {
	s.mu.Lock()
	delete(s.Instances, name)
	saveStateLocked(s)
	s.mu.Unlock()
}

// InstanceNames returns the names of all persisted instances.
func (s *RunnerState) InstanceNames() []string {
	s.mu.Lock()
	defer s.mu.Unlock()
	names := make([]string, 0, len(s.Instances))
	for n := range s.Instances {
		names = append(names, n)
	}
	return names
}

// getOrCreateInstance returns the PersistedInstance for name, creating it if needed.
// Caller must hold s.mu.
func (s *RunnerState) getOrCreateInstance(name string) *PersistedInstance {
	pi, ok := s.Instances[name]
	if !ok {
		pi = &PersistedInstance{Runs: make(map[string]*PersistRun)}
		s.Instances[name] = pi
	}
	return pi
}
