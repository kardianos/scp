// battery — the v90 battery harness (Go replaces the python harness).
//
// Carries the measured v89 free-cell claims (FREECELL.md §9–§10) as gated
// bars against the C kernel of record (./freecell), and cross-checks the
// experimental Go kernel (./fabrun) against it. Kernel-agnostic: both
// binaries share the config surface and the RESULT-line protocol.
//
// Ratchet rule (v89/battery/README.md, carried): every kernel or law
// modification runs the FULL battery before commit; passing experiments
// join the suite and gate future modifications; bars encode claims —
// sharpen by measurement, never soften to pass.
//
// Usage: from v90/:  make all && ./battery [-only name] [-kernel c|go|both]
package main

import (
	"crypto/sha256"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
)

// ------------------------------------------------------------------
// law purity: the header the kernels must print (laws_V2g VERBATIM).
// Byte-drift in these lines = the table changed = automatic RED.
// ------------------------------------------------------------------

var lawHeader = []string{
	"# laws = laws_V2g VERBATIM (defaults): C=1 w1=1.65 w2=2.9 q_detune=1.2",
	"# gamma_res=0.25 gamma_res_m=0.1 p_gate=8 lock_floor=0.005 k_dep=1.2 k_dep_m=2 cap=2.5",
	"# e_s0=1 es_floor=0.05 e_cond=0 f_conv=0.25 f_evap=0.5 s_pull=0.15",
	"# kappa_lock=1 kappa_align=0.5 s_k=0.06 s_disp=0.3 sigma_tumble=0.01",
	"# comb_limit=6 rough_k=0.35 gamma_rough=0.5 mob_sym=1 mob_floor=0.004 field_J=1.8",
	"# quant_A0=1.15 quant_mode=2 (A0eff=1.15)",
}

// ------------------------------------------------------------------
// harness plumbing
// ------------------------------------------------------------------

type runSpec struct {
	kernel string // "c" or "go"
	label  string // log file stem
	args   []string
}

type barResult struct {
	exp, bar string
	value    string
	bound    string
	pass     bool
	gating   bool
}

var (
	binC    = flag.String("bin-c", "./freecell", "C kernel of record")
	binGo   = flag.String("bin-go", "./fabrun", "experimental Go kernel")
	only    = flag.String("only", "", "run only this experiment")
	kernels = flag.String("kernel", "both", "c|go|both")
	jobs    = flag.Int("j", 8, "parallel kernel runs")
	runsDir = flag.String("runs", "runs", "log output directory")
)

var (
	logsMu sync.Mutex
	logs   = map[string]string{} // label -> log content
)

func runOne(sp runSpec) error {
	bin := *binC
	if sp.kernel == "go" {
		bin = *binGo
	}
	cmd := exec.Command(bin, sp.args...)
	out, err := cmd.CombinedOutput()
	if err != nil {
		return fmt.Errorf("%s (%s): %v\n%s", sp.label, bin, err, string(out))
	}
	path := filepath.Join(*runsDir, sp.label+".log")
	if werr := os.WriteFile(path, out, 0o644); werr != nil {
		return werr
	}
	logsMu.Lock()
	logs[sp.label] = string(out)
	logsMu.Unlock()
	return nil
}

func get(label string) string { return logs[label] }

func exf(log, pattern string) (float64, bool) {
	re := regexp.MustCompile(pattern)
	m := re.FindStringSubmatch(log)
	if m == nil {
		return 0, false
	}
	v, err := strconv.ParseFloat(m[1], 64)
	return v, err == nil
}

func countStr(log, sub string) int { return strings.Count(log, sub) }

// first diag row: t=0 line, extract Em_tag (column after the second |...|)
var diagRe = regexp.MustCompile(`(?m)^t=\s*\S+ [^|]*\| [^|]*\| [^|]*\| [^|]*\|\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)`)

func firstEmTag(log string) (float64, bool) {
	m := diagRe.FindStringSubmatch(log)
	if m == nil {
		return 0, false
	}
	v, err := strconv.ParseFloat(m[1], 64)
	return v, err == nil
}

// per-EDGE gg values from the final report
var edgeRe = regexp.MustCompile(`(?m)^# EDGE \d+ \S+ d=\S+ dev=\S+ gg=(\S+) gpl=`)

func edgeGGs(log string) []float64 {
	var out []float64
	for _, m := range edgeRe.FindAllStringSubmatch(log, -1) {
		if v, err := strconv.ParseFloat(m[1], 64); err == nil {
			out = append(out, v)
		}
	}
	return out
}

type reporter struct {
	rows []barResult
}

func (r *reporter) add(exp, bar string, pass bool, value, bound string, gating bool) {
	r.rows = append(r.rows, barResult{exp, bar, value, bound, pass, gating})
}

func (r *reporter) rangeBar(exp, bar, label, pattern string, lo, hi float64, gating bool) {
	v, ok := exf(get(label), pattern)
	if !ok {
		r.add(exp, bar, false, "PARSE-FAIL", fmt.Sprintf("[%g,%g]", lo, hi), gating)
		return
	}
	r.add(exp, bar, v >= lo && v <= hi, fmt.Sprintf("%.6g", v), fmt.Sprintf("[%g,%g]", lo, hi), gating)
}

// ------------------------------------------------------------------
// the suite
// ------------------------------------------------------------------

type experiment struct {
	name  string
	specs []runSpec
	eval  func(r *reporter)
}

func suite() []experiment {
	cw := func(label string, args ...string) runSpec { return runSpec{"c", label, args} }
	gw := func(label string, args ...string) runSpec { return runSpec{"go", label, args} }
	drift := `# RESULT drift_rel (\S+)`
	wantC := *kernels == "c" || *kernels == "both"
	wantGo := *kernels == "go" || *kernels == "both"

	var exps []experiment

	// FB1 conserve+determinism (e1/F2-class): live bath, FP-floor ledger,
	// churn without densification. v89: drift 0.000e+00, z_live 16.4.
	sp := []runSpec{}
	if wantC {
		sp = append(sp, cw("conserve_c", "exp=bath", "T=40"), cw("conserve_c2", "exp=bath", "T=40"))
	}
	if wantGo {
		sp = append(sp, gw("conserve_go", "exp=bath", "T=40"), gw("conserve_go2", "exp=bath", "T=40"))
	}
	exps = append(exps, experiment{"conserve", sp, func(r *reporter) {
		for _, k := range []string{"c", "go"} {
			if (k == "c" && !wantC) || (k == "go" && !wantGo) {
				continue
			}
			l1, l2 := "conserve_"+k, "conserve_"+k+"2"
			v, ok := exf(get(l1), drift)
			r.add("conserve", k+" |drift| <= 1e-13", ok && v >= -1e-13 && v <= 1e-13,
				fmt.Sprintf("%.3g", v), "1e-13", true)
			r.rangeBar("conserve", k+" z_live final", l1, `# RESULT geometry phi_nom=\S+ z_live=(\S+)`, 8, 30, true)
			r.rangeBar("conserve", k+" churn births", l1, `# RESULT births (\S+)`, 100, 1e9, true)
			h1 := sha256.Sum256([]byte(get(l1)))
			h2 := sha256.Sum256([]byte(get(l2)))
			r.add("conserve", k+" byte-determinism", h1 == h2, "sha256", "equal", true)
			// law purity: header lines byte-exact
			pure := true
			for _, want := range lawHeader {
				if !strings.Contains(get(l1), want+"\n") {
					pure = false
				}
			}
			r.add("conserve", k+" law purity (laws_V2g header)", pure, "header", "byte-exact", true)
		}
		if wantC && wantGo {
			gc := strings.SplitN(get("conserve_c"), "\n", 2)
			gg := strings.SplitN(get("conserve_go"), "\n", 2)
			r.add("conserve", "cross C==Go (after line 1)", len(gc) == 2 && len(gg) == 2 && gc[1] == gg[1],
				"full log", "byte-equal", true)
		}
	}})

	// FB2 the bond, both sides + controls (§9.2). v89 measured:
	// up 1.538698 / dn 1.538697 (Δ=1e-6), off_live +0.0329 (=v_rep/K_b),
	// nobond control slides to contact: off +0.1319.
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("pair_up_c", "exp=pair", "bath=0", "T=60", "pair_doff=0.05"),
			cw("pair_dn_c", "exp=pair", "bath=0", "T=60", "pair_doff=-0.05"),
			cw("pair_nobond_c", "exp=pair", "bath=0", "T=60", "pair_doff=0.05", "kappa_bond=0"))
	}
	if wantGo {
		sp = append(sp, gw("pair_up_go", "exp=pair", "bath=0", "T=60", "pair_doff=0.05"))
	}
	exps = append(exps, experiment{"pair", sp, func(r *reporter) {
		if wantC {
			dU, ok1 := exf(get("pair_up_c"), `# RESULT pair d_final=(\S+)`)
			dD, ok2 := exf(get("pair_dn_c"), `# RESULT pair d_final=(\S+)`)
			r.add("pair", "both sides converge |d_up-d_dn| <= 5e-5", ok1 && ok2 && abs(dU-dD) <= 5e-5,
				fmt.Sprintf("%.2g", abs(dU-dD)), "5e-5", true)
			r.rangeBar("pair", "off_live = pressure/K_b (up)", "pair_up_c", `off_live=\+?(-?\S+?) Em_i`, 0.02, 0.05, true)
			r.rangeBar("pair", "off_live (dn)", "pair_dn_c", `off_live=\+?(-?\S+?) Em_i`, 0.02, 0.05, true)
			r.rangeBar("pair", "nobond control slides to contact", "pair_nobond_c", `off_live=\+?(-?\S+?) Em_i`, 0.10, 1.0, true)
			if wantGo {
				dG, ok3 := exf(get("pair_up_go"), `# RESULT pair d_final=(\S+)`)
				r.add("pair", "cross |d_go-d_c| <= 1e-4", ok3 && abs(dG-dU) <= 1e-4,
					fmt.Sprintf("%.2g", abs(dG-dU)), "1e-4", true)
			}
		}
	}})

	// FB3 species by parity (§9.3): ring6 rung-locked, ring5 π-frustrated.
	// v89: ring6 edge_dev_mean 0.0350, gg >= 0.97; ring5 0.1319, gg ~ 0.
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("ring6_c", "exp=ring", "ring_n=6", "bath=0", "T=300"),
			cw("ring5_c", "exp=ring", "ring_n=5", "bath=0", "T=300"))
	}
	exps = append(exps, experiment{"ring", sp, func(r *reporter) {
		if !wantC {
			return
		}
		r.rangeBar("ring", "ring6 edge_dev_mean <= 0.05", "ring6_c", `# RESULT truss edge_dev_mean=(\S+)`, 0, 0.05, true)
		ggs := edgeGGs(get("ring6_c"))
		mn := 1.0
		for _, g := range ggs {
			if g < mn {
				mn = g
			}
		}
		r.add("ring", "ring6 all 6 edges alive, min gg >= 0.9", len(ggs) == 6 && mn >= 0.9,
			fmt.Sprintf("n=%d min=%.3f", len(ggs), mn), ">=0.9", true)
		r.rangeBar("ring", "ring5 pi-frustrated: edge_dev_mean >= 0.10", "ring5_c", `# RESULT truss edge_dev_mean=(\S+)`, 0.10, 1.0, true)
		// v89-measured signature (freecell_ring5_cop.log): the odd cycle
		// distributes the pi defect — 3 of 5 gates dead (gg <= 0.01), the
		// other two partial (0.59/0.61); ALL edges slide to contact.
		ggs5 := edgeGGs(get("ring5_c"))
		dead := 0
		mn5 := 1.0
		for _, g := range ggs5 {
			if g <= 0.01 {
				dead++
			}
			if g < mn5 {
				mn5 = g
			}
		}
		r.add("ring", "ring5 pi defect: >=3 of 5 gates dead", len(ggs5) == 5 && dead >= 3 && mn5 <= 0.01,
			fmt.Sprintf("dead=%d min=%.4f", dead, mn5), ">=3 dead", true)
	}})

	// FB4 Tier-A pulse (§9.6): field packet on the live substrate.
	// v89: v/C = 0.5740 (frozen 0.605), conservation at floor.
	sp = []runSpec{}
	// diag_every=50 is part of the apparatus: the v89 measurement
	// (v/C=0.5740, n=6) fits t in [5,10] before the front saturates the
	// box; the coarser default cadence fits the saturating tail instead.
	if wantC {
		sp = append(sp, cw("pulse_c", "exp=pulse", "L=24", "T=10", "amp=0.25", "sigma=2.5", "kx=1.1", "diag_every=50"))
	}
	if wantGo {
		sp = append(sp, gw("pulse_go", "exp=pulse", "L=24", "T=10", "amp=0.25", "sigma=2.5", "kx=1.1", "diag_every=50"))
	}
	exps = append(exps, experiment{"pulse", sp, func(r *reporter) {
		if wantC {
			r.rangeBar("pulse", "v_over_C in the measured band", "pulse_c", `v_over_C=(\S+)`, 0.50, 0.65, true)
			v, ok := exf(get("pulse_c"), drift)
			r.add("pulse", "|drift| <= 1e-13", ok && abs(v) <= 1e-13, fmt.Sprintf("%.3g", v), "1e-13", true)
		}
		if wantC && wantGo {
			vc, _ := exf(get("pulse_c"), `v_over_C=(\S+)`)
			vg, ok := exf(get("pulse_go"), `v_over_C=(\S+)`)
			r.add("pulse", "cross |v_go-v_c| <= 5e-3", ok && abs(vg-vc) <= 5e-3,
				fmt.Sprintf("%.2g", abs(vg-vc)), "5e-3", true)
		}
	}})

	// FB5 blob at density (§9.5a): localisation holds on the live
	// substrate. v89 L16 T=160: ret 0.8052, conn 1.000, rms 3.11.
	sp = []runSpec{}
	if wantC {
		sp = append(sp, cw("blob_c", "exp=blob", "L=16", "T=160"))
	}
	exps = append(exps, experiment{"blob", sp, func(r *reporter) {
		if !wantC {
			return
		}
		em0, ok0 := firstEmTag(get("blob_c"))
		emf, okf := exf(get("blob_c"), `# RESULT blob Em_tag=(\S+)`)
		ret := 0.0
		if ok0 && okf && em0 > 0 {
			ret = emf / em0
		}
		r.add("blob", "retention >= 0.78 at T=160", ok0 && okf && ret >= 0.78,
			fmt.Sprintf("%.4f", ret), ">=0.78", true)
		r.rangeBar("blob", "connected: conn = 1", "blob_c", `conn=(\S+) z_tag`, 0.999, 1.001, true)
		r.rangeBar("blob", "no spreading: rms <= 3.6", "blob_c", `rms=(\S+) conn`, 0, 3.6, true)
		v, ok := exf(get("blob_c"), drift)
		r.add("blob", "|drift| <= 1e-13", ok && abs(v) <= 1e-13, fmt.Sprintf("%.3g", v), "1e-13", true)
		r.rangeBar("blob", "no densification: z_live <= 30", "blob_c", `# RESULT geometry phi_nom=\S+ z_live=(\S+)`, 0, 30, true)
	}})

	// FB6 FCQ (§10): UUD equilateral by the ladder identity, closed;
	// UDD open chain (D-D edge beyond contact, channel dead).
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("uud_c", "exp=tri", "tri_xU=0.28", "bath=0", "T=100"),
			cw("udd_c", "exp=tri", "tri_kind=1", "tri_xU=0.28", "bath=0", "T=300"))
	}
	exps = append(exps, experiment{"fcq", sp, func(r *reporter) {
		if !wantC {
			return
		}
		dUU, ok1 := exf(get("uud_c"), `d\*_UU=(\S+)`)
		dUD, ok2 := exf(get("uud_c"), `d\*_UD\(m=2\)=(\S+)`)
		r.add("fcq", "ladder identity d*_UU == d*_UD at x_D", ok1 && ok2 && abs(dUU-dUD) <= 1e-6,
			fmt.Sprintf("%.2g", abs(dUU-dUD)), "1e-6", true)
		r.add("fcq", "UUD closed: 0 dead channels", countStr(get("uud_c"), "CHANNEL DEAD") == 0,
			fmt.Sprintf("%d", countStr(get("uud_c"), "CHANNEL DEAD")), "0", true)
		r.add("fcq", "UDD D-D edge beyond contact (seed print)", strings.Contains(get("udd_c"), "OPEN (no contact)"),
			"seed", "OPEN", true)
		r.add("fcq", "UDD open chain: exactly 1 dead channel", countStr(get("udd_c"), "# EDGE 1 ") == 1 &&
			countStr(get("udd_c"), "CHANNEL DEAD") == 1,
			fmt.Sprintf("%d", countStr(get("udd_c"), "CHANNEL DEAD")), "1", true)
	}})

	return exps
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// ------------------------------------------------------------------

func main() {
	flag.Parse()
	if err := os.MkdirAll(*runsDir, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}
	exps := suite()

	// collect runs (filtered), execute with a worker pool
	var specs []runSpec
	for _, e := range exps {
		if *only != "" && e.name != *only {
			continue
		}
		specs = append(specs, e.specs...)
	}
	sem := make(chan struct{}, *jobs)
	var wg sync.WaitGroup
	errs := make(chan error, len(specs))
	for _, sp := range specs {
		wg.Add(1)
		go func(sp runSpec) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()
			fmt.Printf("# run %-14s %v\n", sp.label, sp.args)
			if err := runOne(sp); err != nil {
				errs <- err
			}
		}(sp)
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		fmt.Fprintf(os.Stderr, "RUN FAILED: %v\n", err)
		os.Exit(2)
	}

	// evaluate
	rep := &reporter{}
	for _, e := range exps {
		if *only != "" && e.name != *only {
			continue
		}
		e.eval(rep)
	}

	fmt.Printf("\n%-10s | %-45s | %-16s | %-10s | %s\n", "EXPERIMENT", "BAR", "MEASURED", "BOUND", "VERDICT")
	fmt.Println(strings.Repeat("-", 100))
	nfail := 0
	for _, b := range rep.rows {
		verdict := "PASS"
		if !b.pass {
			verdict = "FAIL"
			nfail++
		}
		fmt.Printf("%-10s | %-45s | %-16s | %-10s | %s\n", b.exp, b.bar, b.value, b.bound, verdict)
	}
	fmt.Println(strings.Repeat("-", 100))
	if nfail == 0 {
		fmt.Printf("ALL GREEN (%d bars)\n", len(rep.rows))
	} else {
		fmt.Printf("RED: %d/%d bars failed\n", nfail, len(rep.rows))
		os.Exit(1)
	}
}
