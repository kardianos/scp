// PRESTRESS fleet orchestrator — bounded concurrency, no pgrep.
//
//	cd v89/prestress && go build -o fleet . && ./fleet -wave all
//
//	fleet [-wave 1|2|3|4|all] [-max 2] [-threads 4] [-dry] [-force]
//
// Skips jobs whose run log already has "# RESULT conservation" or
// "done <name>". Tracks only child processes it starts.
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
	Net       string
	T         float64
	Threads   int
	SnapEvery int
	NoRender  bool
	Extra     []string
	Wave      int
}

func paths() (cand, morph string) {
	return filepath.Join(root, "v89/prestress/candidates"),
		filepath.Join(root, "v89/prestress/morpho/nets")
}

func wave1(threads int) []job {
	cand, morph := paths()
	return []job{
		{Wave: 1, Name: "w1_c8_spectra", Net: filepath.Join(cand, "c8_ring12.net"), T: 200, Threads: threads, SnapEvery: 250},
		{Wave: 1, Name: "w1_c2_cube150", Net: filepath.Join(cand, "c2_cube150.net"), T: 3000, Threads: threads},
		{Wave: 1, Name: "w1_c2_cube150_ctrl", Net: filepath.Join(cand, "c2_cube150_ctrl.net"), T: 3000, Threads: threads},
		{Wave: 1, Name: "w1_c8_ring12", Net: filepath.Join(cand, "c8_ring12.net"), T: 3000, Threads: threads},
		{Wave: 1, Name: "w1_ring8_m3", Net: filepath.Join(morph, "ring8_m3.net"), T: 3000, Threads: threads},
		{Wave: 1, Name: "w1_c5_tube6", Net: filepath.Join(cand, "c5_tube6.net"), T: 3000, Threads: threads},
		{Wave: 1, Name: "w1_c4_ring12_chords", Net: filepath.Join(cand, "c4_ring12_chords.net"), T: 3000, Threads: threads},
	}
}

// Wave 2 — plasticity + topology under κ (pin hunt).
// κ*=1.0 provisional (anneal probe); mid-window 0.5 also tested.
func wave2(threads int) []job {
	cand, morph := paths()
	k1 := []string{"kappa_plast=1.0"}
	k05 := []string{"kappa_plast=0.5"}
	k1h := []string{"kappa_plast=1.0", "tau_harden=50"}
	return []job{
		// PLAST-1 A/B: cube ±κ (κ=0 already w1; re-run named w2 for plast logs)
		{Wave: 2, Name: "w2_c2_plast1_k1", Net: filepath.Join(cand, "c2_cube150.net"), T: 3000, Threads: threads, Extra: k1},
		{Wave: 2, Name: "w2_c2_plast1_k05", Net: filepath.Join(cand, "c2_cube150.net"), T: 3000, Threads: threads, Extra: k05},
		{Wave: 2, Name: "w2_c2_plast1_ctrl_k1", Net: filepath.Join(cand, "c2_cube150_ctrl.net"), T: 3000, Threads: threads, Extra: k1},
		// PLAST-3: parasites + plasticity (self-seal)
		{Wave: 2, Name: "w2_c1_plast3_k1", Net: filepath.Join(cand, "c1_cube125.net"), T: 3000, Threads: threads, Extra: k1},
		{Wave: 2, Name: "w2_c1_frozen", Net: filepath.Join(cand, "c1_cube125.net"), T: 3000, Threads: threads},
		// Topo under plasticity (best W1 survivors + tube)
		{Wave: 2, Name: "w2_ring8_m3_k1", Net: filepath.Join(morph, "ring8_m3.net"), T: 3000, Threads: threads, Extra: k1},
		{Wave: 2, Name: "w2_c8_k1", Net: filepath.Join(cand, "c8_ring12.net"), T: 3000, Threads: threads, Extra: k1},
		{Wave: 2, Name: "w2_tube_k1", Net: filepath.Join(cand, "c5_tube6.net"), T: 3000, Threads: threads, Extra: k1},
		{Wave: 2, Name: "w2_chords_k1", Net: filepath.Join(cand, "c4_ring12_chords.net"), T: 3000, Threads: threads, Extra: k1},
		// hex prism topo (formfind strained) ± plast
		{Wave: 2, Name: "w2_hex_frozen", Net: filepath.Join(cand, "c3_hexprism.net"), T: 3000, Threads: threads},
		{Wave: 2, Name: "w2_hex_k1", Net: filepath.Join(cand, "c3_hexprism.net"), T: 3000, Threads: threads, Extra: k1},
		// PLAST-2: harden after anneal path (τ_harden)
		{Wave: 2, Name: "w2_c2_harden", Net: filepath.Join(cand, "c2_cube150.net"), T: 3000, Threads: threads, Extra: k1h},
		{Wave: 2, Name: "w2_ring8_harden", Net: filepath.Join(morph, "ring8_m3.net"), T: 3000, Threads: threads, Extra: k1h},
	}
}

// Wave 3 — discriminators / package / negatives / free-search
func wave3(threads int) []job {
	cand, morph := paths()
	return []job{
		{Wave: 3, Name: "w3_c8_flight233", Net: filepath.Join(cand, "c8_ring12_flight233.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_ring9_m3", Net: filepath.Join(morph, "ring9_m3.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_ring10_m4", Net: filepath.Join(morph, "ring10_m4.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_ring12_m5", Net: filepath.Join(morph, "ring12_m5.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_ring12_m6", Net: filepath.Join(morph, "ring12_m6.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_hopf", Net: filepath.Join(morph, "hopf12x12_m5.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_mobius", Net: filepath.Join(morph, "mobius6_d1.5.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_octahedron", Net: filepath.Join(morph, "octahedron_d1.5.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_free0", Net: filepath.Join(morph, "free_0.net"), T: 3000, Threads: threads},
		{Wave: 3, Name: "w3_free1", Net: filepath.Join(morph, "free_1.net"), T: 3000, Threads: threads},
		// negatives under plast too
		{Wave: 3, Name: "w3_mobius_k1", Net: filepath.Join(morph, "mobius6_d1.5.net"), T: 3000, Threads: threads, Extra: []string{"kappa_plast=1.0"}},
		{Wave: 3, Name: "w3_hopf_k1", Net: filepath.Join(morph, "hopf12x12_m5.net"), T: 3000, Threads: threads, Extra: []string{"kappa_plast=1.0"}},
		// torus / truncocta infeasible-class (still run as negative)
		{Wave: 3, Name: "w3_torus3x8", Net: filepath.Join(cand, "c6_torus3x8.net"), T: 2000, Threads: threads},
		{Wave: 3, Name: "w3_truncocta", Net: filepath.Join(cand, "c7_truncocta.net"), T: 2000, Threads: threads},
	}
}

// Wave 4 — mass-spectrum gates (single foam; multi-seed full P19 needs new foams)
func wave4(threads int) []job {
	cand, morph := paths()
	kh := []string{"kappa_plast=1.0", "tau_harden=50"}
	return []job{
		// P19-lite: hardened vs frozen on same foam (cluster prep; not full C3)
		{Wave: 4, Name: "w4_p19_cube_harden", Net: filepath.Join(cand, "c2_cube150.net"), T: 4000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p19_ring8_harden", Net: filepath.Join(morph, "ring8_m3.net"), T: 4000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p19_c8_harden", Net: filepath.Join(cand, "c8_ring12.net"), T: 4000, Threads: threads, Extra: kh},
		// P20 m-family ladder (coexistence of m-classes — separate runs, compare masses)
		{Wave: 4, Name: "w4_p20_ring12_m5", Net: filepath.Join(morph, "ring12_m5.net"), T: 3000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p20_ring12_m6", Net: filepath.Join(morph, "ring12_m6.net"), T: 3000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p20_ring8_m3", Net: filepath.Join(morph, "ring8_m3.net"), T: 3000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p20_ring9_m3", Net: filepath.Join(morph, "ring9_m3.net"), T: 3000, Threads: threads, Extra: kh},
		{Wave: 4, Name: "w4_p20_ring10_m4", Net: filepath.Join(morph, "ring10_m4.net"), T: 3000, Threads: threads, Extra: kh},
		// longer tube under harden (exception re-test)
		{Wave: 4, Name: "w4_tube_harden", Net: filepath.Join(cand, "c5_tube6.net"), T: 4000, Threads: threads, Extra: kh},
	}
}

func allJobs(threads int) []job {
	var j []job
	j = append(j, wave1(threads)...)
	j = append(j, wave2(threads)...)
	j = append(j, wave3(threads)...)
	j = append(j, wave4(threads)...)
	return j
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
	sc := bufio.NewScanner(f)
	buf := make([]byte, 0, 64*1024)
	sc.Buffer(buf, 4*1024*1024)
	doneLine := "done " + name
	for sc.Scan() {
		line := sc.Text()
		if line == doneLine || strings.HasPrefix(line, doneLine+"\t") {
			return true
		}
		if strings.HasPrefix(line, "# RESULT conservation ") {
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

	_ = os.MkdirAll(runsDir(), 0o755)
	launchLog := filepath.Join(runsDir(), j.Name+".launch.log")
	lf, err := os.Create(launchLog)
	if err != nil {
		return err
	}
	defer lf.Close()

	cmd := exec.Command("python3", args...)
	cmd.Dir = root
	cmd.Stdout = lf
	cmd.Stderr = lf
	fmt.Fprintf(os.Stderr, "[fleet] START w%d %s T=%g extra=%v\n", j.Wave, j.Name, j.T, j.Extra)
	t0 := time.Now()
	if err := cmd.Run(); err != nil {
		fmt.Fprintf(os.Stderr, "[fleet] FAIL %s after %s: %v\n", j.Name, time.Since(t0).Round(time.Second), err)
		return fmt.Errorf("%s: %w", j.Name, err)
	}
	fmt.Fprintf(os.Stderr, "[fleet] DONE  %s in %s\n", j.Name, time.Since(t0).Round(time.Second))
	return nil
}

func main() {
	wave := flag.String("wave", "all", "1|2|3|4|all|2-4")
	maxN := flag.Int("max", 2, "max concurrent jobs")
	threads := flag.Int("threads", 4, "OMP threads per job")
	dry := flag.Bool("dry", false, "print plan only")
	force := flag.Bool("force", false, "re-run even if done")
	flag.Parse()

	if *maxN < 1 {
		*maxN = 1
	}

	var jobs []job
	switch *wave {
	case "1":
		jobs = wave1(*threads)
	case "2":
		jobs = wave2(*threads)
	case "3":
		jobs = wave3(*threads)
	case "4":
		jobs = wave4(*threads)
	case "2-4":
		jobs = append(wave2(*threads), wave3(*threads)...)
		jobs = append(jobs, wave4(*threads)...)
	case "all":
		jobs = allJobs(*threads)
	default:
		fmt.Fprintf(os.Stderr, "unknown -wave %q (use 1|2|3|4|all|2-4)\n", *wave)
		os.Exit(2)
	}

	var pending []job
	for _, j := range jobs {
		if !*force && jobDone(j.Name) {
			fmt.Fprintf(os.Stderr, "[fleet] SKIP %s (already done)\n", j.Name)
			continue
		}
		pending = append(pending, j)
	}
	fmt.Fprintf(os.Stderr, "[fleet] wave=%s max=%d threads=%d pending=%d/%d\n",
		*wave, *maxN, *threads, len(pending), len(jobs))

	if *dry {
		for _, j := range pending {
			fmt.Printf("w%d  %s  T=%g  extra=%v\n  %s\n", j.Wave, j.Name, j.T, j.Extra, j.Net)
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

	nfail := 0
	for err := range errCh {
		fmt.Fprintf(os.Stderr, "[fleet] error: %v\n", err)
		nfail++
	}
	if nfail > 0 {
		fmt.Fprintf(os.Stderr, "[fleet] finished with %d failure(s)\n", nfail)
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "[fleet] complete wave=%s\n", *wave)
}
