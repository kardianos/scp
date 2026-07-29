// PRESTRESS fleet orchestrator — bounded concurrency, no pgrep.
//
// Builds and runs: go run .  OR  go build -o fleet && ./fleet
//
//   fleet [-wave 1] [-max 2] [-threads 4] [-dry]
//
// Skips jobs whose run log already contains a terminal "done <name>" line
// from run_net.py (or # RESULT conservation). Never matches its own argv
// via process greps — only tracks child Cmds it started.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"time"
)

const root = "/home/d/code/scp"

type job struct {
	Name      string
	Net       string // absolute path
	T         float64
	Threads   int
	SnapEvery int // 0 => run_net default
	NoRender  bool
	Extra     []string
}

func wave1(threads int) []job {
	cand := filepath.Join(root, "v89/prestress/candidates")
	morph := filepath.Join(root, "v89/prestress/morpho/nets")
	return []job{
		// Phase 0 — spectra (short); may already be done
		{Name: "w1_c8_spectra", Net: filepath.Join(cand, "c8_ring12.net"), T: 200, Threads: threads, SnapEvery: 250},
		// Phase 1 — cube showdown
		{Name: "w1_c2_cube150", Net: filepath.Join(cand, "c2_cube150.net"), T: 3000, Threads: threads},
		{Name: "w1_c2_cube150_ctrl", Net: filepath.Join(cand, "c2_cube150_ctrl.net"), T: 3000, Threads: threads},
		// Phase 2 — validation + ladder proxy
		{Name: "w1_c8_ring12", Net: filepath.Join(cand, "c8_ring12.net"), T: 3000, Threads: threads},
		{Name: "w1_ring8_m3", Net: filepath.Join(morph, "ring8_m3.net"), T: 3000, Threads: threads},
		// Phase 3 — structural bets
		{Name: "w1_c5_tube6", Net: filepath.Join(cand, "c5_tube6.net"), T: 3000, Threads: threads},
		{Name: "w1_c4_ring12_chords", Net: filepath.Join(cand, "c4_ring12_chords.net"), T: 3000, Threads: threads},
	}
}

func runsDir() string {
	return filepath.Join(root, "v89/prestress/runs")
}

func jobDone(name string) bool {
	logPath := filepath.Join(runsDir(), name+".log")
	f, err := os.Open(logPath)
	if err != nil {
		return false
	}
	defer f.Close()
	// Prefer run_net's final "done <name>" line; also accept conservation RESULT.
	sc := bufio.NewScanner(f)
	// large log lines (diag rows) can be long
	buf := make([]byte, 0, 64*1024)
	sc.Buffer(buf, 4*1024*1024)
	doneLine := "done " + name
	for sc.Scan() {
		line := sc.Text()
		if line == doneLine || strings.HasPrefix(line, doneLine+"\t") {
			return true
		}
		if strings.HasPrefix(line, "# RESULT conservation ") {
			// full sim finished final_report
			// keep scanning for "done" but this is enough
			return true
		}
	}
	return false
}

func runJob(j job) error {
	if _, err := os.Stat(j.Net); err != nil {
		return fmt.Errorf("%s: net missing: %w", j.Name, err)
	}
	py := filepath.Join(root, "v89/prestress/run_net.py")
	args := []string{
		py,
		"--name", j.Name,
		"--net", j.Net,
		"--T", fmt.Sprintf("%g", j.T),
		"--threads", fmt.Sprintf("%d", j.Threads),
	}
	if j.SnapEvery > 0 {
		args = append(args, "--snap_every", fmt.Sprintf("%d", j.SnapEvery))
	}
	if j.NoRender {
		args = append(args, "--no-render")
	}
	for _, e := range j.Extra {
		args = append(args, "--extra", e)
	}

	launchLog := filepath.Join(runsDir(), j.Name+".launch.log")
	_ = os.MkdirAll(runsDir(), 0o755)
	lf, err := os.Create(launchLog)
	if err != nil {
		return err
	}
	defer lf.Close()

	cmd := exec.Command("python3", args...)
	cmd.Dir = root
	cmd.Stdout = lf
	cmd.Stderr = lf
	// Do not inherit a polluted process group hunt; plain Wait is enough.
	fmt.Fprintf(os.Stderr, "[fleet] START %s T=%g threads=%d net=%s\n", j.Name, j.T, j.Threads, j.Net)
	t0 := time.Now()
	if err := cmd.Run(); err != nil {
		fmt.Fprintf(os.Stderr, "[fleet] FAIL %s after %s: %v\n", j.Name, time.Since(t0).Round(time.Second), err)
		return fmt.Errorf("%s: %w", j.Name, err)
	}
	fmt.Fprintf(os.Stderr, "[fleet] DONE  %s in %s\n", j.Name, time.Since(t0).Round(time.Second))
	return nil
}

func main() {
	wave := flag.Int("wave", 1, "wave number (only 1 implemented)")
	maxN := flag.Int("max", 2, "max concurrent jobs")
	threads := flag.Int("threads", 4, "OMP threads per job")
	dry := flag.Bool("dry", false, "print plan only")
	force := flag.Bool("force", false, "re-run even if done")
	flag.Parse()

	if *wave != 1 {
		fmt.Fprintf(os.Stderr, "only -wave 1 is defined\n")
		os.Exit(2)
	}
	if *maxN < 1 {
		*maxN = 1
	}

	jobs := wave1(*threads)
	var pending []job
	for _, j := range jobs {
		if !*force && jobDone(j.Name) {
			fmt.Fprintf(os.Stderr, "[fleet] SKIP %s (already done)\n", j.Name)
			continue
		}
		pending = append(pending, j)
	}
	fmt.Fprintf(os.Stderr, "[fleet] wave=%d max=%d threads=%d pending=%d/%d\n",
		*wave, *maxN, *threads, len(pending), len(jobs))

	if *dry {
		for _, j := range pending {
			fmt.Printf("%s  T=%g  %s\n", j.Name, j.T, j.Net)
		}
		return
	}
	if len(pending) == 0 {
		fmt.Fprintf(os.Stderr, "[fleet] nothing to do\n")
		return
	}

	sem := make(chan struct{}, *maxN)
	var wg sync.WaitGroup
	errCh := make(chan error, len(pending))

	for _, j := range pending {
		j := j
		wg.Add(1)
		sem <- struct{}{}
		go func() {
			defer wg.Done()
			defer func() { <-sem }()
			if err := runJob(j); err != nil {
				errCh <- err
			}
		}()
	}
	wg.Wait()
	close(errCh)

	var nfail int
	for err := range errCh {
		fmt.Fprintf(os.Stderr, "[fleet] error: %v\n", err)
		nfail++
	}
	if nfail > 0 {
		fmt.Fprintf(os.Stderr, "[fleet] finished with %d failure(s)\n", nfail)
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "[fleet] WAVE1 complete\n")
}
