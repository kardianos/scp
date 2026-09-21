// Bench: parcore parallel core vs the serial baseline across sizes and worker
// counts. Reports wall time, steps/s, speedup vs 1-worker and vs serial, cell
// throughput (cells/sec/proc), and an optional phase-timer breakdown.
//
// Run:
//
//	GOWORK=off go run ./cmd/bench -n 64000 -steps 20 -workers 1,2,4,8,16
//	GOWORK=off go run ./cmd/bench -n 8000,64000,512000 -workers 1,2,4,8 -phases
package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
	"unsafe"

	"scp/parcore"
	"scp/parcore/example"
)

type row struct {
	N, NC, Slots, Workers int
	Ms, StepsPerSec       float64
	Speedup1W, SpeedupSer float64
	CellsPerSecPerProc    float64
	Boundary              int
}

func main() {
	nFlag := flag.String("n", "8000", "approx cell count(s), comma-separated (each rounded to n³)")
	steps := flag.Int("steps", 30, "timed steps per (size, workers) point")
	warmup := flag.Int("warmup", 3, "untimed warmup steps")
	workers := flag.String("workers", "1,2,4,8,16", "comma-separated worker counts")
	skipSerial := flag.Bool("skip-serial", false, "skip the serial baseline")
	phases := flag.Bool("phases", false, "print phase-timer breakdown per size")
	csvPath := flag.String("csv", "", "optional CSV output path")
	cpuprof := flag.String("cpuprofile", "", "write CPU profile to file (first size/worker only)")
	flag.Parse()

	sizes := parseInts(*nFlag)
	wlist := parseInts(*workers)
	if len(sizes) == 0 || len(wlist) == 0 {
		fmt.Fprintln(os.Stderr, "need at least one -n and one -workers value")
		os.Exit(2)
	}
	flush := func() { _ = os.Stdout.Sync() }

	var c0 example.Cell
	fmt.Printf("sizeof Cell=%d Slot=%d (single AoS slab; no pointers)\n",
		unsafe.Sizeof(c0), unsafe.Sizeof(example.Slot{}))
	fmt.Printf("GOMAXPROCS(default)=%d  NumCPU=%d  steps=%d warmup=%d\n",
		runtime.GOMAXPROCS(0), runtime.NumCPU(), *steps, *warmup)
	fmt.Printf("sizes=%v  workers=%v\n\n", sizes, wlist)
	flush()

	var rows []row

	for _, nApprox := range sizes {
		law := example.DefaultLaw()
		side := icbrt(nApprox)
		if side < 2 {
			side = 2
		}
		law.L = float64(side) * 1.0

		w := example.BuildWorld(law, nApprox)
		nc := len(w.Cells)
		nslots := len(w.SI)
		fmt.Printf("======== n≈%d  NC=%d  side=%d  slots=%d  L=%.1f ========\n",
			nApprox, nc, side, nslots, law.L)
		flush()

		var ser time.Duration
		if !*skipSerial {
			s := example.NewSerial(w, law)
			for i := 0; i < *warmup; i++ {
				s.Step()
			}
			t0 := time.Now()
			for i := 0; i < *steps; i++ {
				s.Step()
			}
			ser = time.Since(t0)
			fmt.Printf("serial          %10.1f ms  %8.1f steps/s\n",
				ser.Seconds()*1e3, float64(*steps)/ser.Seconds())
			flush()
		}

		var t1w time.Duration
		for _, nw := range wlist {
			if nw < 1 {
				continue
			}
			ncpu := runtime.NumCPU()
			gmp := nw
			if gmp > ncpu {
				gmp = ncpu
			}
			if gmp < 1 {
				gmp = 1
			}
			prev := runtime.GOMAXPROCS(gmp)

			var profStop func()
			if *cpuprof != "" && nApprox == sizes[0] && nw == wlist[0] {
				profStop = startProfile(*cpuprof)
			}
			d := &example.Driver{Law: law}
			cl := parcore.FromWorld[example.Cell, example.Slot](w, d, nw)
			bnd := cl.NBoundary()
			if nw <= 8 || nw == wlist[0] {
				fmt.Printf("  workers=%d GOMAXPROCS=%d boundary_edges=%d", nw, gmp, bnd)
				for i := 0; i < nw && i < 4; i++ {
					sh := cl.Shard(i)
					fmt.Printf("  sh%d:own=%d/gh=%d", i, sh.NOwn, sh.NGhost)
				}
				fmt.Println()
			}

			if *phases && nw == wlist[0] {
				parcore.EnablePhaseTimers(true)
			}
			for i := 0; i < *warmup; i++ {
				cl.Step()
			}
			t1 := time.Now()
			for i := 0; i < *steps; i++ {
				cl.Step()
			}
			el := time.Since(t1)
			if *phases && nw == wlist[0] {
				parcore.PhaseTimerReport(os.Stdout, *steps, d.PhaseNames())
				parcore.EnablePhaseTimers(false)
			}
			cl.Stop()
			if profStop != nil {
				profStop()
				profStop = nil
			}
			runtime.GOMAXPROCS(prev)

			if nw == 1 || t1w == 0 {
				t1w = el
			}
			sp1 := t1w.Seconds() / el.Seconds()
			spSer := 0.0
			if ser > 0 {
				spSer = ser.Seconds() / el.Seconds()
			}
			cellsPerSecPerProc := float64(nc) * float64(*steps) / el.Seconds() / float64(gmp)
			fmt.Printf("par workers=%-3d %10.1f ms  %8.1f steps/s  vs1w=%.2fx",
				nw, el.Seconds()*1e3, float64(*steps)/el.Seconds(), sp1)
			if ser > 0 {
				fmt.Printf("  vsSer=%.2fx", spSer)
			}
			fmt.Printf("  %.2f Mcell/s/proc\n", cellsPerSecPerProc/1e6)
			flush()

			rows = append(rows, row{
				N: nApprox, NC: nc, Slots: nslots, Workers: nw,
				Ms: el.Seconds() * 1e3, StepsPerSec: float64(*steps) / el.Seconds(),
				Speedup1W: sp1, SpeedupSer: spSer, CellsPerSecPerProc: cellsPerSecPerProc,
				Boundary: bnd,
			})
		}
		fmt.Println()
	}

	printMatrix("steps/s", rows, sizes, wlist, func(r row) float64 { return r.StepsPerSec })
	printMatrix("speedup vs 1w", rows, sizes, wlist, func(r row) float64 { return r.Speedup1W })
	printMatrix("Mcell/s/proc", rows, sizes, wlist, func(r row) float64 { return r.CellsPerSecPerProc / 1e6 })

	// efficiency = speedup/workers
	fmt.Println("======== efficiency (speedup/workers) ========")
	seen := map[int]bool{}
	for _, r := range rows {
		if seen[r.NC] {
			continue
		}
		seen[r.NC] = true
		fmt.Printf("NC=%d:", r.NC)
		for _, nw := range wlist {
			hit := findRow(rows, r.NC, nw)
			if hit == nil {
				continue
			}
			eff := hit.Speedup1W / float64(nw) * 100
			fmt.Printf("  %d:%.0f%%", nw, eff)
		}
		fmt.Println()
	}

	if *csvPath != "" {
		writeCSV(*csvPath, rows)
	}
}

func printMatrix(title string, rows []row, sizes, wlist []int, cell func(row) float64) {
	fmt.Printf("======== SUMMARY: %s ========\n", title)
	fmt.Printf("%10s", "NC\\W")
	for _, nw := range wlist {
		fmt.Printf(" %10d", nw)
	}
	fmt.Println()
	seen := map[int]bool{}
	for _, r := range rows {
		if seen[r.NC] {
			continue
		}
		seen[r.NC] = true
		fmt.Printf("%10d", r.NC)
		for _, nw := range wlist {
			hit := findRow(rows, r.NC, nw)
			if hit == nil {
				fmt.Printf(" %10s", "—")
			} else if title == "speedup vs 1w" {
				fmt.Printf(" %9.2fx", cell(*hit))
			} else {
				fmt.Printf(" %10.2f", cell(*hit))
			}
		}
		fmt.Println()
	}
}

func findRow(rows []row, nc, nw int) *row {
	for i := range rows {
		if rows[i].NC == nc && rows[i].Workers == nw {
			return &rows[i]
		}
	}
	return nil
}

func writeCSV(path string, rows []row) {
	f, err := os.Create(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "csv: %v\n", err)
		os.Exit(1)
	}
	fmt.Fprintln(f, "n_approx,nc,slots,workers,ms,steps_per_s,speedup_vs_1w,speedup_vs_serial,mcell_per_s_per_proc,boundary_edges")
	for _, r := range rows {
		fmt.Fprintf(f, "%d,%d,%d,%d,%.3f,%.3f,%.4f,%.4f,%.4f,%d\n",
			r.N, r.NC, r.Slots, r.Workers, r.Ms, r.StepsPerSec, r.Speedup1W, r.SpeedupSer,
			r.CellsPerSecPerProc/1e6, r.Boundary)
	}
	_ = f.Close()
	fmt.Printf("\ncsv → %s\n", path)
}

func parseInts(s string) []int {
	var out []int
	for _, p := range strings.Split(s, ",") {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil || n < 1 {
			fmt.Fprintf(os.Stderr, "skip bad int %q\n", p)
			continue
		}
		out = append(out, n)
	}
	return out
}

func icbrt(n int) int {
	x := float64(n)
	s := 1
	for s*s*s < int(x+0.5) {
		s++
	}
	return s
}
