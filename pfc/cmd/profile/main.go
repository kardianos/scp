// CPU + heap profile and per-phase timers for a large owned-shard run.
package main

import (
	"flag"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"runtime/pprof"
	"time"

	"scp/pfc/core"
	"scp/pfc/par"
)

func main() {
	n := flag.Int("n", 216000, "approx cells")
	steps := flag.Int("steps", 30, "timed steps")
	workers := flag.Int("workers", 8, "workers")
	cpuprof := flag.String("cpuprofile", "/tmp/pfc_cpu.pprof", "cpu profile path")
	memprof := flag.String("memprofile", "/tmp/pfc_mem.pprof", "heap profile path")
	flag.Parse()

	law := core.DefaultLaw()
	side := int(cbrt(float64(*n)) + 0.5)
	if side < 2 {
		side = 2
	}
	law.L = float64(side)

	src := core.BuildLattice(law, *n)
	fmt.Printf("build NC=%d slots=%d side=%d workers=%d\n", src.NC(), len(src.Slots), src.Side, *workers)
	cl := par.FromWorld(src, *workers)
	defer cl.Stop()
	fmt.Printf("boundary_edges=%d\n", cl.NBoundary())
	for i := 0; i < *workers && i < 4; i++ {
		sh := cl.Shard(i)
		fmt.Printf("  shard%d own=%d ghost=%d slots=%d interior=%d boundary=%d\n",
			i, sh.NOwn, sh.NGhost, len(sh.Slots), len(sh.Interior), len(sh.Boundary))
	}

	for i := 0; i < 3; i++ {
		cl.Step()
	}
	runtime.GC()

	par.EnablePhaseTimers(true)
	t0 := time.Now()
	for i := 0; i < *steps; i++ {
		cl.Step()
	}
	wall := time.Since(t0)
	par.EnablePhaseTimers(false)
	fmt.Printf("\nNC=%d workers=%d steps=%d wall=%.1fms (%.1f steps/s)\n",
		src.NC(), *workers, *steps, wall.Seconds()*1e3, float64(*steps)/wall.Seconds())
	par.PrintPhaseTimers(os.Stdout, *steps)

	// CPU profile
	cf, err := os.Create(*cpuprof)
	if err != nil {
		panic(err)
	}
	if err := pprof.StartCPUProfile(cf); err != nil {
		panic(err)
	}
	for i := 0; i < *steps; i++ {
		cl.Step()
	}
	pprof.StopCPUProfile()
	_ = cf.Close()
	fmt.Printf("\ncpu profile → %s\n", *cpuprof)

	for i := 0; i < 5; i++ {
		cl.Step()
	}
	runtime.GC()
	mf, err := os.Create(*memprof)
	if err != nil {
		panic(err)
	}
	if err := pprof.WriteHeapProfile(mf); err != nil {
		panic(err)
	}
	_ = mf.Close()
	fmt.Printf("heap profile → %s\n", *memprof)

	fmt.Println("\n=== go tool pprof -top (cpu) ===")
	runPprof("-top", "-nodecount=25", *cpuprof)
	fmt.Println("\n=== go tool pprof -top inuse_space (heap) ===")
	runPprof("-top", "-sample_index=inuse_space", "-nodecount=20", *memprof)
	fmt.Println("\n=== go tool pprof -top alloc_space (heap) ===")
	runPprof("-top", "-sample_index=alloc_space", "-nodecount=15", *memprof)
}

func runPprof(args ...string) {
	cmd := exec.Command("go", append([]string{"tool", "pprof"}, args...)...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	_ = cmd.Run()
}

func cbrt(x float64) float64 {
	z := x
	if z < 1 {
		z = 1
	}
	for i := 0; i < 16; i++ {
		z = (2*z + x/(z*z)) / 3
	}
	return z
}
