// Bench: serial core vs owned-shard Cluster across sizes and worker counts.
//
// Default sweep: workers 1,2,4,8,16,32 at each -n size. Reports wall time,
// steps/s, speedup vs 1-worker (and vs serial core when workers includes 1).
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

	"scp/pfc/core"
	"scp/pfc/par"
)

func main() {
	nFlag := flag.String("n", "8000", "approx cell count(s), comma-separated (each rounded to n³ lattice)")
	steps := flag.Int("steps", 40, "timed steps per (size, workers) point")
	warmup := flag.Int("warmup", 2, "untimed warmup steps before each timed run")
	workers := flag.String("workers", "1,2,4,8,16,32", "comma-separated worker counts")
	L := flag.Float64("L", 0, "box size (0 = auto from n per size)")
	skipSerial := flag.Bool("skip-serial", false, "skip serial core.World baseline")
	csvPath := flag.String("csv", "", "optional CSV output path")
	flag.Parse()

	sizes := parseInts(*nFlag)
	wlist := parseInts(*workers)
	if len(sizes) == 0 || len(wlist) == 0 {
		fmt.Fprintln(os.Stderr, "need at least one -n and one -workers value")
		os.Exit(2)
	}

	// Line-buffer progress when stdout is a file/pipe (remote nohup runs).
	flush := func() {
		_ = os.Stdout.Sync()
	}

	var lo core.CellLo
	var hi core.CellHi
	fmt.Printf("sizeof CellLo=%d CellHi=%d (split AoS)\n", unsafe.Sizeof(lo), unsafe.Sizeof(hi))
	fmt.Printf("GOMAXPROCS(default)=%d  NumCPU=%d  steps=%d warmup=%d\n",
		runtime.GOMAXPROCS(0), runtime.NumCPU(), *steps, *warmup)
	fmt.Printf("sizes=%v  workers=%v\n\n", sizes, wlist)
	flush()

	type row struct {
		N, NC, Slots, Workers int
		Ms, StepsPerSec       float64
		Speedup1W, SpeedupSer float64
		Boundary              int
	}
	var rows []row

	for _, nApprox := range sizes {
		law := core.DefaultLaw()
		side := int(cbrt(float64(nApprox)) + 0.5)
		if side < 2 {
			side = 2
		}
		if *L > 0 {
			law.L = *L
		} else {
			law.L = float64(side) * 1.0
		}

		// Build once per size for topology stats; rebuild fresh per timed run.
		probe := core.BuildLattice(law, nApprox)
		nc := probe.NC()
		nslots := len(probe.Slots)
		fmt.Printf("======== n≈%d  NC=%d  side=%d  slots=%d  L=%.1f ========\n",
			nApprox, nc, probe.Side, nslots, law.L)
		flush()

		var ser time.Duration
		if !*skipSerial {
			w := core.BuildLattice(law, nApprox)
			for i := 0; i < *warmup; i++ {
				w.Step()
			}
			t0 := time.Now()
			for i := 0; i < *steps; i++ {
				w.Step()
			}
			ser = time.Since(t0)
			fmt.Printf("serial          %10.1f ms  %8.1f steps/s\n",
				ser.Seconds()*1e3, float64(*steps)/ser.Seconds())
			flush()
		}

		var t1w time.Duration // baseline for parallel speedup
		for _, nw := range wlist {
			if nw < 1 {
				continue
			}
			// Pin GOMAXPROCS to worker count so 1-thread is truly single-threaded.
			// Cap at NumCPU so we don't oversubscribe beyond hardware.
			ncpu := runtime.NumCPU()
			gmp := nw
			if gmp > ncpu {
				gmp = ncpu
			}
			if gmp < 1 {
				gmp = 1
			}
			prev := runtime.GOMAXPROCS(gmp)

			src := core.BuildLattice(law, nApprox)
			cl := par.FromWorld(src, nw)
			bnd := cl.NBoundary()
			if nw <= 8 || nw == wlist[0] {
				fmt.Printf("  workers=%d GOMAXPROCS=%d boundary_edges=%d", nw, gmp, bnd)
				for i := 0; i < nw && i < 4; i++ {
					sh := cl.Shard(i)
					fmt.Printf("  sh%d:own=%d/gh=%d", i, sh.NOwn, sh.NGhost)
				}
				fmt.Println()
			}

			for i := 0; i < *warmup; i++ {
				cl.Step()
			}
			t1 := time.Now()
			for i := 0; i < *steps; i++ {
				cl.Step()
			}
			el := time.Since(t1)
			cl.Stop()
			runtime.GOMAXPROCS(prev)

			if nw == 1 || t1w == 0 {
				t1w = el
			}
			sp1 := t1w.Seconds() / el.Seconds()
			spSer := 0.0
			if ser > 0 {
				spSer = ser.Seconds() / el.Seconds()
			}
			fmt.Printf("par workers=%-3d %10.1f ms  %8.1f steps/s  vs1w=%.2fx",
				nw, el.Seconds()*1e3, float64(*steps)/el.Seconds(), sp1)
			if ser > 0 {
				fmt.Printf("  vsSer=%.2fx", spSer)
			}
			fmt.Println()
			flush()

			rows = append(rows, row{
				N: nApprox, NC: nc, Slots: nslots, Workers: nw,
				Ms: el.Seconds() * 1e3, StepsPerSec: float64(*steps) / el.Seconds(),
				Speedup1W: sp1, SpeedupSer: spSer, Boundary: bnd,
			})
		}
		fmt.Println()
	}

	// Summary matrix: steps/s by (size × workers)
	fmt.Println("======== SUMMARY: steps/s ========")
	fmt.Printf("%10s", "NC\\W")
	for _, nw := range wlist {
		fmt.Printf(" %10d", nw)
	}
	fmt.Println()
	// group by NC
	seen := map[int]bool{}
	for _, r := range rows {
		if seen[r.NC] {
			continue
		}
		seen[r.NC] = true
		fmt.Printf("%10d", r.NC)
		for _, nw := range wlist {
			var hit *row
			for i := range rows {
				if rows[i].NC == r.NC && rows[i].Workers == nw {
					hit = &rows[i]
					break
				}
			}
			if hit == nil {
				fmt.Printf(" %10s", "—")
			} else {
				fmt.Printf(" %10.1f", hit.StepsPerSec)
			}
		}
		fmt.Println()
	}

	fmt.Println("======== SUMMARY: speedup vs 1 worker ========")
	fmt.Printf("%10s", "NC\\W")
	for _, nw := range wlist {
		fmt.Printf(" %10d", nw)
	}
	fmt.Println()
	seen = map[int]bool{}
	for _, r := range rows {
		if seen[r.NC] {
			continue
		}
		seen[r.NC] = true
		fmt.Printf("%10d", r.NC)
		for _, nw := range wlist {
			var hit *row
			for i := range rows {
				if rows[i].NC == r.NC && rows[i].Workers == nw {
					hit = &rows[i]
					break
				}
			}
			if hit == nil {
				fmt.Printf(" %10s", "—")
			} else {
				fmt.Printf(" %9.2fx", hit.Speedup1W)
			}
		}
		fmt.Println()
	}

	// Cap-out note: first worker count where efficiency vs ideal < 50%
	fmt.Println("======== scaling notes (efficiency = speedup/workers) ========")
	seen = map[int]bool{}
	for _, r := range rows {
		if seen[r.NC] {
			continue
		}
		seen[r.NC] = true
		fmt.Printf("NC=%d:", r.NC)
		for _, nw := range wlist {
			var hit *row
			for i := range rows {
				if rows[i].NC == r.NC && rows[i].Workers == nw {
					hit = &rows[i]
					break
				}
			}
			if hit == nil || nw < 1 {
				continue
			}
			eff := hit.Speedup1W / float64(nw) * 100
			fmt.Printf("  %d:%.0f%%", nw, eff)
		}
		fmt.Println()
	}

	if *csvPath != "" {
		f, err := os.Create(*csvPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "csv: %v\n", err)
			os.Exit(1)
		}
		fmt.Fprintln(f, "n_approx,nc,slots,workers,ms,steps_per_s,speedup_vs_1w,speedup_vs_serial,boundary_edges")
		for _, r := range rows {
			fmt.Fprintf(f, "%d,%d,%d,%d,%.3f,%.3f,%.4f,%.4f,%d\n",
				r.N, r.NC, r.Slots, r.Workers, r.Ms, r.StepsPerSec, r.Speedup1W, r.SpeedupSer, r.Boundary)
		}
		_ = f.Close()
		fmt.Printf("\ncsv → %s\n", *csvPath)
	}
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

func cbrt(x float64) float64 {
	return mathCbrt(x)
}
