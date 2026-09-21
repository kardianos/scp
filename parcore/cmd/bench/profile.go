package main

import (
	"os"
	"runtime/pprof"
)

// startProfile begins a CPU profile into path; the returned func stops it.
func startProfile(path string) func() {
	f, err := os.Create(path)
	if err != nil {
		return nil
	}
	_ = pprof.StartCPUProfile(f)
	return func() {
		pprof.StopCPUProfile()
		_ = f.Close()
	}
}
