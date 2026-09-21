//go:build ignore

// One-shot partial scorer for in-flight sustain FCS files.
//   go run theory/probe/partial_score.go theory/runs/sustain_T3000_...
package main

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("usage: go run partial_score.go <run_dir>")
		os.Exit(2)
	}
	dir := os.Args[1]
	v93 := filepath.Dir(filepath.Dir(dir)) // .../theory/runs/X -> .../theory; need v93
	// dir = v93/theory/runs/NAME
	v93 = filepath.Clean(filepath.Join(dir, "../../.."))
	if filepath.Base(filepath.Dir(dir)) == "runs" {
		v93 = filepath.Clean(filepath.Join(dir, "../.."))
	}
	fcsdump := filepath.Join(v93, "fcsdump")
	matches, _ := filepath.Glob(filepath.Join(dir, "*.fcs"))
	fmt.Printf("# partial score %s n=%d\n", dir, len(matches))
	fmt.Println("arm t_last Em EmMax nHot slope_proxy")
	for _, fcs := range matches {
		name := strings.TrimSuffix(filepath.Base(fcs), ".fcs")
		cmd := exec.Command(fcsdump, "-mode", "stats", fcs)
		out, err := cmd.Output()
		if err != nil {
			fmt.Printf("%s ERR %v\n", name, err)
			continue
		}
		var ts, ems []float64
		var lastEm, lastMax float64
		var lastHot int
		for _, line := range strings.Split(string(out), "\n") {
			if !strings.HasPrefix(line, "CELL") {
				continue
			}
			f := strings.Fields(line)
			if len(f) < 8 {
				continue
			}
			t, _ := strconv.ParseFloat(f[1], 64)
			em, _ := strconv.ParseFloat(f[3], 64)
			emx, _ := strconv.ParseFloat(f[6], 64)
			nh, _ := strconv.Atoi(f[7])
			ts = append(ts, t)
			ems = append(ems, em)
			lastEm, lastMax, lastHot = em, emx, nh
		}
		slope := 0.0
		if n := len(ts); n >= 2 {
			// last two frames
			dt := ts[n-1] - ts[n-2]
			if dt > 0 {
				slope = (ems[n-1] - ems[n-2]) / dt
			}
		}
		fmt.Printf("%s t=%.0f Em=%.1f max=%.2f hot=%d dEm/dt=%+.4f (nfr=%d)\n",
			name, ts[len(ts)-1], lastEm, lastMax, lastHot, slope, len(ts))
	}
}
