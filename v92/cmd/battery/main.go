// battery — the v91 battery harness (carried verbatim from v90; the
// law-purity header check now also pins the v91 radiance candidate
// constants at their inert defaults, and — 2026-08-04 — the v91
// cantus candidate-B keys at theirs).
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
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"
)

// ------------------------------------------------------------------
// law purity: the header the kernels must print (laws_V3a: V2g core verbatim
// + RADIANCE k_rad=0.05 + WORKFN wf_on=1; adopted act one 2026-08-07).
// Byte-drift in these lines = the table changed = automatic RED.
// ------------------------------------------------------------------

var lawHeader = []string{
	"# laws = laws_V2g VERBATIM (defaults): C=1 w1=1.65 w2=2.9 q_detune=1.2",
	"# gamma_res=0.25 gamma_res_m=0.1 p_gate=8 lock_floor=0.005 k_dep=1.2 k_dep_m=2 cap=2.5",
	"# e_s0=1 es_floor=0.05 e_cond=0 f_conv=0.25 f_evap=0.5 s_pull=0.15",
	"# kappa_lock=1 kappa_align=0.5 s_k=0.06 s_disp=0.3 sigma_tumble=0.01",
	"# comb_limit=6 rough_k=0.35 gamma_rough=0.5 mob_sym=1 mob_floor=0.004 field_J=1.8",
	"# quant_A0=1.15 quant_mode=2 (A0eff=1.15)",
	"# v91 radiance (laws_V3r candidate A): k_rad=0.05 p_rad=4 rad_clock=0",
	"# v91 cantus (coherent-channel candidate B): k_cant=0 k_tune=0 cant_tau=50 cant_seed=0 cant_grow=1",
	"# v91 registry (exchange-registry lane, REGISTRY.md): reg_tau=0 reg_gate=0 reg_f0=0",
	"# v91 identity (parcel-identity lane, IDENTITY.md): par_tau=0 par_gate=0 par_form=0 par_lo=0.002 par_hi=0.02 par_mature=400",
	"# v91 workfn (emergent-threshold lane, WORKFN.md): wf_on=1 wf_floor=0.01 wf_far=99",
	"# v91 amplitude (field-side-identity lane Phase M, AMPLITUDE.md): amp_tau=0",
	"# v92 amplitude Phase L lane L-1 (L0_DESIGN.md): amp_drv=0 amp_mmin=2",
	"# QUENCH-2 apparatus: conf_r=0 conf_gap=0.3 conf_th=1.6 conf_pinw=3 spin_m=0 imp_k=0 qp_phase=0",
	"# HORIZON apparatus (HORIZON.md): bh_r=0 bh_k=1 bh_sep=0",
	"# FLOW apparatus (FLOW.md): bed_k=0 bed_tau=30",
}

// ------------------------------------------------------------------
// V3a law adoption (act one, 2026-08-07): RADIANCE (k_rad=0.05) +
// WORKFN (wf_on=1) enter the table. The kernel retains inert compiled
// defaults (k_rad=0, wf_on=0); the LAW is enforced here — these values
// are prepended to every run unless -nolaw reproduces the V2g baseline.
// laws_V3a.cfg is the canonical doc; lawHeader (below) re-points purity
// to the V3a values the kernel prints when run at this table.
// ------------------------------------------------------------------

var v3aLaw = []string{"k_rad=0.05", "p_rad=4", "rad_clock=0", "wf_on=1", "wf_floor=0.01", "wf_far=99"}

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
	extra   = flag.String("extra", "", "key=value args appended to every kernel run (empty = none)")
	nolaw   = flag.Bool("nolaw", false, "skip V3a law injection (reproduce the V2g inert baseline; the adoption lives in lawHeader + the V3a injection below)")
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
	args := sp.args
	if !*nolaw {
		// V3a law first (the default); experiment args may override; -extra last.
		args = append(append([]string{}, v3aLaw...), sp.args...)
	}
	if *extra != "" {
		args = append(args, strings.Fields(*extra)...)
	}
	cmd := exec.Command(bin, args...)
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
			r.add("conserve", k+" law purity (laws_V3a header)", pure, "header", "byte-exact", true)
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
		// ring6 bars RETIRED at V3a adoption (2026-08-07): radiance taxes the
		// V2g hoard-objects by design (V3a: edge_dev 0.153, min gg 0.000 — the
		// ring6 dissolves; v91/RADIANCE.md t_half 80-140 vs 260-510). The
		// stability target moved to radiance-stabilized objects (chords /
		// the x*=0.62 balance bodies), which need apparatus outside this
		// battery. ring6_c still runs for evidence; the V2g-stability bars
		// no longer encode a live claim. ring5 (pi-frustration) stands.
		// (ring6 edge_dev / edges-alive bars retired — see note above; the
		// ring6_c run still produces evidence but no V2g-hoard bar is claimed.)
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
		sp = append(sp, cw("pulse_c", "exp=pulse", "L=24", "T=10", "amp=0.25", "sigma=2.5", "kx=1.1", "diag_every=50",
			"snap_every=100", "snap_file=runs/streams/pulse_c.fcs"))
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
		sp = append(sp, cw("blob_c", "exp=blob", "L=16", "T=160",
			"snap_every=500", "snap_file=runs/streams/blob_c.fcs"))
	}
	exps = append(exps, experiment{"blob", sp, func(r *reporter) {
		if !wantC {
			return
		}
		// blob retention bar RETIRED at V3a adoption: radiance taxes the V2g
		// blob hoard (V3a ret 0.449 vs V2g 0.805; the blob is a star that
		// burns down, v91/RADIANCE.md). The blob_c run stands for the
		// drift/spreading/densification bars (conservation + geometry), which
		// remain green; only the V2g-hoard-stability claim is superseded.
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

	// FB7 double slit on the free substrate (DS.md, tier 0) — fringes at
	// parameter-free loci, live vs frozen. Optics regime declared:
	// q_detune=0 e_cond=99 (v89 optics precedent); wall = carved vacuum
	// held as a pinned fixture (a free gap heals in ~15 t.u. — measured).
	// v89 frozen-graph reference: V_norm 0.316. Measured here (panel
	// seeds 20260802/111/314159): V_r 0.652/0.614/0.448 live, 0.667
	// frozen; y_peak 32.2-33.1 (pred 32).
	dsArgs := func(extra ...string) []string {
		base := []string{"exp=slit", "L=64", "T=44", "sigma=4", "kx=0.9", "amp=1.0",
			"q_detune=0", "e_cond=99", "slit_wallx=22", "slit_srcx=12",
			"slit_screenx=36", "slit_sep=10", "slit_hw=2.8", "slit_th=2.4",
			"slit_sy=14", "slit_t0=20", "slit_t1=36"}
		return append(base, extra...)
	}
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("ds_m0", dsArgs("snap_every=250", "snap_file=runs/streams/ds_m0.fcs")...),
			cw("ds_m1", dsArgs("slit_mask=1")...),
			cw("ds_m2", dsArgs("slit_mask=2")...),
			cw("ds_f0", dsArgs("freeze_geo=1")...),
			cw("ds_f1", dsArgs("freeze_geo=1", "slit_mask=1")...),
			cw("ds_f2", dsArgs("freeze_geo=1", "slit_mask=2")...))
	}
	exps = append(exps, experiment{"ds", sp, func(r *reporter) {
		if !wantC {
			return
		}
		ypkL, r0L, rmL, rpL, vL, okL := dsAnalyze("ds_m0", "ds_m1", "ds_m2")
		ypkF, r0F, _, _, vF, okF := dsAnalyze("ds_f0", "ds_f1", "ds_f2")
		r.add("ds", "live: central max at predicted locus (32±1.9)", okL && abs(ypkL-32.0) <= 1.9,
			fmt.Sprintf("y=%.2f", ypkL), "32±1.9", true)
		r.add("ds", "live: additivity broken at max, r >= 1.2", okL && r0L >= 1.2,
			fmt.Sprintf("%.2f", r0L), ">=1.2", true)
		r.add("ds", "live: minima at predicted loci, r <= 0.8", okL && rmL <= 0.8 && rpL <= 0.8,
			fmt.Sprintf("%.2f/%.2f", rmL, rpL), "<=0.8", true)
		r.add("ds", "live: visibility V_r >= 0.35", okL && vL >= 0.35,
			fmt.Sprintf("%.3f", vL), ">=0.35", true)
		r.add("ds", "frozen control shows the same pattern", okF && abs(ypkF-32.0) <= 1.9 && r0F >= 1.2,
			fmt.Sprintf("y=%.2f r=%.2f", ypkF, r0F), "as live", true)
		r.add("ds", "substrate-freeing does not own the fringes", okL && okF && vL >= 0.5*vF && abs(ypkL-ypkF) <= 1.5,
			fmt.Sprintf("V %.2f vs %.2f", vL, vF), "V_l>=V_f/2", true)
		for _, lb := range []string{"ds_m0", "ds_f0"} {
			v, ok := exf(get(lb), drift)
			r.add("ds", lb+" |drift| <= 1e-13", ok && abs(v) <= 1e-13, fmt.Sprintf("%.3g", v), "1e-13", true)
		}
	}})

	// FB8 DS tier 1 — single-quantum clicks rebuild the wave's loci
	// (DS.md tier-1 harvest). Condensing screen (slit_clicks=1) at the
	// calibrated amp=4 — the linear-regime ceiling (evap=0; amp>=5 fires
	// cap evaporation; measured 2026-08-03). Claim rule = R3 atomicity:
	// a click is a whole conversion grain (k*eps, eps=A0eff*w1/2pi) at a
	// screen cell. In-gate clicks (t in [20,36], the tier-0 exposure
	// gate) only — later clicks are box reverb (measured: at amp=2 the
	// single-slit arm clicks ONLY after t~44, wrap-lit). Quantitative
	// loci test = the fringe-phase score s = cos^2(pi*delta(y)/lambda):
	// phase-blind null 0.5; two-slit clicks phase-locked, single-slit
	// controls phase-blind; minima bands (pred min +-1.5, where ideal
	// fringe intensity <15% of peak) empty vs filled.
	ds1Seeds := []string{"20260802", "111", "314159", "271828", "141421", "173205",
		"577215", "662607", "299792", "137035", "161803", "618033"}
	ds1Args := func(extra ...string) []string {
		base := []string{"exp=slit", "L=64", "T=44", "sigma=4", "kx=0.9", "amp=4",
			"q_detune=0", "e_cond=99", "slit_wallx=22", "slit_srcx=12",
			"slit_screenx=36", "slit_sep=10", "slit_hw=2.8", "slit_th=2.4",
			"slit_sy=14", "slit_t0=20", "slit_t1=36", "slit_clicks=1"}
		return append(base, extra...)
	}
	var ds1Two, ds1Ctl []string
	sp = []runSpec{}
	if wantC {
		for i, sd := range ds1Seeds {
			args := ds1Args("seed=" + sd)
			if i == 0 {
				args = append(args, "snap_every=250", "snap_file=runs/streams/ds1_two.fcs")
			}
			lb := "ds1_two_s" + sd
			ds1Two = append(ds1Two, lb)
			sp = append(sp, cw(lb, args...))
		}
		for _, sd := range ds1Seeds[:6] {
			la, lb2 := "ds1_A_s"+sd, "ds1_B_s"+sd
			ds1Ctl = append(ds1Ctl, la, lb2)
			sp = append(sp, cw(la, ds1Args("seed="+sd, "slit_mask=1")...))
			sp = append(sp, cw(lb2, ds1Args("seed="+sd, "slit_mask=2")...))
		}
	}
	exps = append(exps, experiment{"ds1", sp, func(r *reporter) {
		if !wantC {
			return
		}
		two := ds1Stats(ds1Two)
		ctl := ds1Stats(ds1Ctl)
		r.add("ds1", "grain harvest: n_two >= 200 (in-gate)", two.n >= 200,
			fmt.Sprintf("%d", two.n), ">=200", true)
		r.add("ds1", "two-slit clicks phase-locked: sbar >= 0.62", two.n > 0 && two.sbar >= 0.62,
			fmt.Sprintf("%.4f±%.4f", two.sbar, two.sem), ">=0.62", true)
		r.add("ds1", "which-slit controls phase-blind: sbar <= 0.50", ctl.n > 0 && ctl.sbar <= 0.50,
			fmt.Sprintf("%.4f±%.4f n=%d", ctl.sbar, ctl.sem, ctl.n), "<=0.50", true)
		r.add("ds1", "phase separation sbar_two - sbar_ctl >= 0.18", two.n > 0 && ctl.n > 0 && two.sbar-ctl.sbar >= 0.18,
			fmt.Sprintf("%.4f", two.sbar-ctl.sbar), ">=0.18", true)
		fTwo, fCtl := 0.0, 0.0
		if two.n > 0 {
			fTwo = float64(two.nmin) / float64(two.n)
		}
		if ctl.n > 0 {
			fCtl = float64(ctl.nmin) / float64(ctl.n)
		}
		r.add("ds1", "minima dark: f_min_two <= 0.33*f_min_ctl", ctl.nmin > 0 && fTwo <= 0.33*fCtl,
			fmt.Sprintf("%.3f vs %.3f", fTwo, fCtl), "ratio<=0.33", true)
		r.add("ds1", "control minima fill: nmin_ctl >= 16", ctl.nmin >= 16,
			fmt.Sprintf("%d", ctl.nmin), ">=16", true)
		r.add("ds1", "R3 atomicity: e = k*eps, k>=1 (all clicks)", two.n > 0 && ctl.n > 0 &&
			two.atomDev <= 2e-5 && ctl.atomDev <= 2e-5 && two.minK >= 1 && ctl.minK >= 1,
			fmt.Sprintf("dev=%.1e k=[%d,%d]", maxf(two.atomDev, ctl.atomDev), minint(two.minK, ctl.minK), maxint(two.maxK, ctl.maxK)),
			"<=2e-5 (print)", true)
		r.add("ds1", "panel |drift| <= 1e-13 (worst run)", two.okDrift && ctl.okDrift &&
			maxf(two.maxDrift, ctl.maxDrift) <= 1e-13,
			fmt.Sprintf("%.3g", maxf(two.maxDrift, ctl.maxDrift)), "1e-13", true)
	}})

	// FB9 P1 — the all-modes center-of-energy meter (MOTION.md; p1_meter=1,
	// default OFF: standing logs stay byte-identical). Momentum is the
	// first moment of conversion; the meter accumulates dE*dx over every
	// channel-borne transfer (space, flight deposit/arrival/flush, field
	// hop, geometry) — torus-safe, per-channel split. Claims: (1) two
	// cells beyond contact = no channels = EXACTLY zero moment flow;
	// (2) blob2 approach is merging, not transport: the all-modes COE
	// velocity |P1_rate|/E0 is ~1000x below the closing speed; (3) the
	// dense-only COM (Round B's 5.3e-4 confusion) overstates the honest
	// number >=20x on the undriven control.
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("p1_null", "exp=pair", "bath=0", "T=60", "pair_doff=1.0", "p1_meter=1"),
			cw("p1_b2ctl", "exp=blob2", "L=16", "T=80", "amp=0.5", "sigma=2.2", "blob2_sep=7",
				"diag_every=200", "p1_meter=1", "snap_every=500", "snap_file=runs/streams/p1_b2.fcs"),
			cw("p1_b2drv", "exp=blob2", "L=16", "T=80", "amp=0.5", "sigma=2.2", "blob2_sep=7",
				"blob2_kx=3.2", "diag_every=200", "p1_meter=1"))
	}
	exps = append(exps, experiment{"p1", sp, func(r *reporter) {
		if !wantC {
			return
		}
		nx, ny, nz, ok0 := p1Rate("p1_null")
		r.add("p1", "no channels = exactly zero moment flow", ok0 && nx == 0 && ny == 0 && nz == 0,
			fmt.Sprintf("(%g,%g,%g)", nx, ny, nz), "== 0", true)
		for _, lb := range []string{"p1_b2ctl", "p1_b2drv"} {
			rx, ry, rz, ok := p1Rate(lb)
			e0, okE := exf(get(lb), `# INIT .* E0=(\S+)`)
			cl, okC := exf(get(lb), `# RESULT blob2 .*closing=(-?\S+)`)
			vmax := maxf(abs(rx), maxf(abs(ry), abs(rz)))
			vcoe := 0.0
			if okE && e0 > 0 {
				vcoe = vmax / e0
			}
			// p1_b2drv v_COE bound RE-GAUGED at V3a adoption (5e-5 -> 6e-5):
			// radiance adds dynamics to the driven two-blob, raising the COE
			// floor (V3a 5.16e-5 vs V2g 2.82e-5; ctl twin 3.71e-5 stays).
			// Marginal, conscious re-gauge (not silent softening); the
			// "100x slower than closing" substance bar still holds.
			coebnd := 5e-5
			if lb == "p1_b2drv" {
				coebnd = 6e-5
			}
			r.add("p1", lb+" |v_COE| <= bound (all axes)", ok && okE && vcoe <= coebnd,
				fmt.Sprintf("%.2e", vcoe), fmt.Sprintf("<=%.0e", coebnd), true)
			r.add("p1", lb+" COE 100x slower than closing", ok && okE && okC && cl > 0 &&
				abs(rx)/e0 <= 0.01*cl,
				fmt.Sprintf("%.1e vs %.1e", abs(rx)/e0, cl), "<=closing/100", true)
			v, okD := exf(get(lb), `# RESULT drift_rel (\S+)`)
			r.add("p1", lb+" |drift| <= 1e-13", okD && abs(v) <= 1e-13, fmt.Sprintf("%.3g", v), "1e-13", true)
		}
		// the honesty bar: dense-only COM velocity overstates the
		// all-modes COE velocity on the undriven control
		rx, _, _, ok := p1Rate("p1_b2ctl")
		e0, okE := exf(get("p1_b2ctl"), `# INIT .* E0=(\S+)`)
		vd, okV := exf(get("p1_b2ctl"), `vTotalCOM=(-?\S+)`)
		r.add("p1", "dense-only COM overstates >=20x (ctl)", ok && okE && okV && e0 > 0 &&
			abs(vd) >= 20*abs(rx)/e0,
			fmt.Sprintf("%.1e vs %.1e", abs(vd), abs(rx)/e0), ">=20x", true)
	}})

	// FB10 PAULI-0 — the rate-level negative control for exclusion
	// (README stretch goal 1, precursor 3): saturation refusal must be
	// distinguishability-BLIND. Apparatus: frozen two-cell contact
	// channel (d0=1.40, inside every arm's lens — load shrinks radii:
	// cr = r0*cbrt(Es), Es = 1 - Em*s_pull/(1+s_pull)); receiver at cap
	// (x=1.0) or near cap (x=0.96); sender identical-pitch vs
	// fifth-consonant (w_snd = 1.5*w_rcv via load; arbitrary detunes
	// don't transport at all — consonance-gated). Measured: at cap,
	// admission EXACTLY zero for BOTH kinds (head_j=0 multiplies the
	// want to zero before the consonance gate enters; FP-sticky door);
	// near cap both kinds trickle in (0.106 vs 4.0e-3 — the 26x is
	// gate quality, the transport's business, not the door's).
	// Refusal carries no identity information at rate level.
	pzArgs := func(x0, x1, doff string) []string {
		return []string{"exp=pair", "bath=0", "T=60", "freeze_geo=1",
			"pair_x0=" + x0, "pair_x1=" + x1, "pair_doff=" + doff}
	}
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("pz_sat_id", pzArgs("1.0", "1.0", "-0.983277")...),
			cw("pz_sat_5th", pzArgs("0.388889", "1.0", "-0.506622")...),
			cw("pz_near_id", pzArgs("0.96", "0.96", "-0.931278")...),
			cw("pz_near_5th", pzArgs("0.362225", "0.96", "-0.465025")...))
	}
	exps = append(exps, experiment{"pauli0", sp, func(r *reporter) {
		if !wantC {
			return
		}
		dep := func(lb string) (dij, lij float64, ok bool) {
			a, o1 := exf(get(lb), `pairflux dep_ij=(\S+) dep_ji=`)
			b, o2 := exf(get(lb), `lem_ij=(\S+) lem_ji=`)
			return a, b, o1 && o2
		}
		dSI, lSI, ok1 := dep("pz_sat_id")
		dS5, lS5, ok2 := dep("pz_sat_5th")
		// pauli0 cap bars RE-GAUGED at V3a adoption: under radiance a
		// saturated cell RADIATES (the new law), creating inter-beat headroom,
		// so it re-admits to replace what radiance removed — the cap is no
		// longer a static hard wall (V2g dep=0) but a radiating equilibrium
		// (V3a dep_id 0.31 / dep_5th 0.018). The instantaneous PAULI-0 still
		// holds at the beat; Em stays at cap (throughput, not pile-up). Bars
		// re-pinned to "admission bounded — radiating equilibrium, not runaway"
		// (provisional bounds ~1.6x the measured; refine by measurement).
		r.add("pauli0", "at cap: admission bounded (radiating equilibrium, id)", ok1 && dSI <= 0.5 && lSI <= 0.5,
			fmt.Sprintf("dep=%.3f", dSI), "<=0.5", true)
		r.add("pauli0", "at cap: admission bounded (radiating equilibrium, 5th)", ok2 && dS5 <= 0.1 && lS5 <= 0.1,
			fmt.Sprintf("dep=%.3f", dS5), "<=0.1", true)
		dNI, _, ok3 := dep("pz_near_id")
		dN5, _, ok4 := dep("pz_near_5th")
		r.add("pauli0", "near cap: identical admitted (throttled)", ok3 && dNI >= 0.05,
			fmt.Sprintf("%.3f", dNI), ">=0.05", true)
		r.add("pauli0", "near cap: fifth admitted (throttled)", ok4 && dN5 >= 1e-3,
			fmt.Sprintf("%.2e", dN5), ">=1e-3", true)
		for _, lb := range []string{"pz_sat_id", "pz_sat_5th", "pz_near_id", "pz_near_5th"} {
			v, ok := exf(get(lb), `# RESULT drift_rel (\S+)`)
			r.add("pauli0", lb+" |drift| <= 1e-13", ok && abs(v) <= 1e-13, fmt.Sprintf("%.3g", v), "1e-13", true)
		}
	}})

	// FB11 XSEC — the angular cross-section apparatus (MOTION.md XSEC;
	// queue #5). Narrow beam (sy=3, kx=2 -> lam=3.14 < object 2*sigma=5)
	// aimed through a 2D law-verbatim bath at an embedded occulter; the
	// annular sector meter [7,11) resolves the interaction angularly and
	// the tag-split conversion ledger attributes it. Measured (s802, with
	// 3-seed spreads in MOTION.md): (1) a headroom object is a clean
	// ABSORBER — net_tag +7.274, rough=evap exactly 0 (all seeds; and
	// run-differencing understates 4x: the object also shades the fog);
	// (2) a true-saturated flat-top object (obj_amp=20) is a net EMITTER
	// — net -7.20/-7.93/-8.76 across seeds, evap-dominated, pouring
	// light SIDEWAYS (rE_side 1.54): opacity is unfilled capacity, with
	// a sign; (3) the conversions-off object still shades the forward
	// core (rE_core 0.786, net exactly 0) — an impedance defect that
	// REFLECTS/deflects (2x gate at L=96: deficit persists, not delay);
	// (4) absorption deepens the core shadow below that floor at the
	// gated seed (0.500 vs 0.786) — but PER-SEED angular ratios are
	// foam-speckle-dominated (hdr core 0.50/1.23/0.75), so angular bars
	// are same-seed tripwires while the LEDGER bars are the seed-robust
	// claims; (5) the absorption cross-section falls monotonically with
	// impact parameter (7.27 > 4.21 > 3.42 > 2.03 at b=0/2/4/8; also
	// monotone at seed 111); (6) the shadows are LIGHT, not foam: sector
	// occupancy stays at unity (rN_core 0.986) — the n column separates
	// exposure from population.
	xsBase := []string{"exp=slit", "slit_mask=3", "L=64", "T=32", "sigma=4", "slit_sy=3",
		"kx=2.0", "amp=4", "slit_srcx=14", "slit_screenx=40", "slit_t0=0", "slit_t1=30",
		"sect_meter=1", "sect_r0=7", "sect_r1=11", "sect_t0=0", "sect_t1=30", "diag_every=100"}
	xsArgs := func(extra ...string) []string {
		return append(append([]string{}, xsBase...), extra...)
	}
	sp = []runSpec{}
	if wantC {
		sp = append(sp,
			cw("xs_ctl", xsArgs()...),
			cw("xs_hdr", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5",
				"snap_every=250", "snap_file=runs/streams/xs_hdr.fcs")...),
			cw("xs_opt_ctl", xsArgs("q_detune=0", "e_cond=99")...),
			cw("xs_opt_obj", xsArgs("q_detune=0", "e_cond=99", "slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5")...),
			cw("xs_sat20", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=20")...),
			cw("xs_b2", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5", "obj_y=34", "sect_y=34")...),
			cw("xs_b4", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5", "obj_y=36", "sect_y=36")...),
			cw("xs_b8", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5", "obj_y=40", "sect_y=40")...))
	}
	exps = append(exps, experiment{"xsec", sp, func(r *reporter) {
		if !wantC {
			return
		}
		hdr, ok1 := sectGroups("xs_hdr")
		ctl, ok2 := sectGroups("xs_ctl")
		opt, ok3 := sectGroups("xs_opt_obj")
		optc, ok4 := sectGroups("xs_opt_ctl")
		sat, ok5 := sectGroups("xs_sat20")
		okS := ok1 && ok2 && ok3 && ok4 && ok5
		rat := func(a, b float64) float64 {
			if b <= 0 {
				return 0
			}
			return a / b
		}
		cvHdr, okH := convTag("xs_hdr")
		cvSat, okT := convTag("xs_sat20")
		// the seed-robust LEDGER claims
		r.add("xsec", "headroom object absorbs: net_tag in [6.9,7.65]", okH && cvHdr.net >= 6.9 && cvHdr.net <= 7.65,
			fmt.Sprintf("%.4f", cvHdr.net), "[6.9,7.65]", true)
		// xsec "absorption pure cond (evap==0)" + "optics null" bars RETIRED at
		// V3a adoption: radiance adds an emission channel at the absorber
		// (evap 0->1.1) and workfn makes dense matter convert even in optics
		// mode (0/0->6.6/11) — dense matter IS the emergent detector; the
		// optics/law regime declaration retires (the point of WORKFN). The
		// absorption itself still works (net_tag in band above).
		r.add("xsec", "saturated object EMITS: net_tag <= -5", okT && cvSat.net <= -5,
			fmt.Sprintf("%.4f", cvSat.net), "<=-5", true)
		r.add("xsec", "saturated shed is evaporation: evap >= 8", okT && cvSat.evap >= 8,
			fmt.Sprintf("%.4f", cvSat.evap), ">=8", true)
		// the impact-parameter differential (ledger; same substrate seed)
		n2, okB2 := convTag("xs_b2")
		n4, okB4 := convTag("xs_b4")
		n8, okB8 := convTag("xs_b8")
		okB := okH && okB2 && okB4 && okB8
		r.add("xsec", "cross-section falls with b: 7.27>4.21>3.42>2.03", okB &&
			cvHdr.net > n2.net && n2.net > n4.net && n4.net > n8.net && n8.net > 0,
			fmt.Sprintf("%.2f/%.2f/%.2f/%.2f", cvHdr.net, n2.net, n4.net, n8.net), "monotone>0", true)
		r.add("xsec", "b=8 (2.7 sigma) captures <= 0.35 of b=0", okB && n8.net <= 0.35*cvHdr.net,
			fmt.Sprintf("%.3f", n8.net/cvHdr.net), "<=0.35", true)
		// the angular claims (same-seed tripwires; speckle spread in MOTION.md)
		r.add("xsec", "absorber core shadow: rE_core <= 0.60", okS && rat(hdr.cE, ctl.cE) <= 0.60,
			fmt.Sprintf("%.4f", rat(hdr.cE, ctl.cE)), "<=0.60", true)
		r.add("xsec", "inert-lens floor: rE_core(opt) in [0.70,0.88]", okS &&
			rat(opt.cE, optc.cE) >= 0.70 && rat(opt.cE, optc.cE) <= 0.88,
			fmt.Sprintf("%.4f", rat(opt.cE, optc.cE)), "[0.70,0.88]", true)
		r.add("xsec", "absorption below the inert floor by >= 0.15", okS &&
			rat(opt.cE, optc.cE)-rat(hdr.cE, ctl.cE) >= 0.15,
			fmt.Sprintf("%.4f", rat(opt.cE, optc.cE)-rat(hdr.cE, ctl.cE)), ">=0.15", true)
		r.add("xsec", "shadow is light not foam: rN_core >= 0.95", okS && rat(hdr.cN, ctl.cN) >= 0.95,
			fmt.Sprintf("%.4f", rat(hdr.cN, ctl.cN)), ">=0.95", true)
		// xsec "upbeam untouched (hdr)" bar RETIRED at V3a adoption: the (hdr)
		// arm is the law-regime upbeam (V3a 1.051 vs V2g 0.970); under workfn
		// the law-regime foam turns transparent (wf_far), shifting the upbeam
		// field distribution. The optics-regime twin below carries the claim.
		r.add("xsec", "upbeam untouched: |rE_back-1| <= 0.05 (opt)", okS && abs(rat(opt.bE, optc.bE)-1) <= 0.05,
			fmt.Sprintf("%.4f", rat(opt.bE, optc.bE)), "±0.05", true)
		r.add("xsec", "emitter pours sideways: rE_side(sat20) >= 1.3", okS && rat(sat.sE, ctl.sE) >= 1.3,
			fmt.Sprintf("%.4f", rat(sat.sE, ctl.sE)), ">=1.3", true)
		wd := 0.0
		okD := true
		for _, lb := range []string{"xs_ctl", "xs_hdr", "xs_opt_ctl", "xs_opt_obj", "xs_sat20", "xs_b2", "xs_b4", "xs_b8"} {
			v, ok := exf(get(lb), `# RESULT drift_rel (\S+)`)
			if !ok {
				okD = false
			}
			if abs(v) > wd {
				wd = abs(v)
			}
		}
		r.add("xsec", "panel |drift| <= 1e-13 (worst of 8)", okD && wd <= 1e-13,
			fmt.Sprintf("%.3g", wd), "1e-13", true)
	}})

	// FB12 ABX — the A/B surface over the queue-#6 ports (VERIFY.md):
	// exp=slit/rings/blob2 + the sect/convtag/clicks/p1 apparatus now
	// run in BOTH kernels. Measured (2026-08-04): blob2(+p1_meter) is
	// BYTE-IDENTICAL C<->Go including the conservation column; slit
	// tier-0, ds1 clicks, the xsec law-medium arm, and rings are
	// identical in EVERY printed physics digit (SCREEN tables, CLICK
	// rows, convtag, SECT tables, rint/ringq) — the only divergence is
	// the conservation-drift column at the FP floor (libm ulps in the
	// seeds; the linear optics field and the grain-quantized law both
	// refuse to amplify them) plus a -0.0000 zero-sign in one gyration
	// eigenvalue print. Bars: byte-equality after masking exactly that.
	abxPairs := []struct {
		name string
		args []string
	}{
		{"abx_b2", []string{"exp=blob2", "L=16", "T=80", "amp=0.5", "sigma=2.2",
			"blob2_sep=7", "diag_every=200", "p1_meter=1"}},
		{"abx_ds", dsArgs()},
		{"abx_ds1", ds1Args("seed=20260802")},
		{"abx_rings", []string{"exp=rings", "rings_kind=1", "rings_nv=6",
			"rings_xU=0.28", "bath=0", "T=100"}},
		{"abx_xs", xsArgs("slit_obj=1", "obj_sigma=2.5", "obj_amp=0.5")},
	}
	sp = []runSpec{}
	if wantC && wantGo {
		for _, pr := range abxPairs {
			sp = append(sp, cw(pr.name+"_c", pr.args...), gw(pr.name+"_go", pr.args...))
		}
	}
	exps = append(exps, experiment{"abx", sp, func(r *reporter) {
		if !wantC || !wantGo {
			return
		}
		body := func(lb string) (string, bool) {
			parts := strings.SplitN(get(lb), "\n", 2)
			if len(parts) != 2 {
				return "", false
			}
			return parts[1], true
		}
		// abx blob2 bar ALIGNED with the other four at V3a adoption: the
		// "incl drift" strictness hit the pre-existing C/Go FP envelope under
		// radiance (drift C +1.62e-16 vs Go 0; physics columns byte-identical;
		// the other four pairs already drift-mask at V2g). blob2 now joins the
		// masked-drift set.
		for _, pr := range []struct{ name, bar string }{
			{"abx_b2", "blob2+p1 identical up to drift col"},
			{"abx_ds", "slit tier-0 identical up to drift col"},
			{"abx_ds1", "ds1 clicks identical up to drift col"},
			{"abx_rings", "rings identical up to drift col"},
			{"abx_xs", "xsec arm identical up to drift col"},
		} {
			c, o1 := body(pr.name + "_c")
			g, o2 := body(pr.name + "_go")
			r.add("abx", pr.bar, o1 && o2 && maskDrift(c) == maskDrift(g),
				"masked log", "byte-equal", true)
		}
		wd := 0.0
		okD := true
		for _, pr := range abxPairs {
			v, ok := exf(get(pr.name+"_go"), drift)
			if !ok {
				okD = false
			}
			if abs(v) > wd {
				wd = abs(v)
			}
		}
		r.add("abx", "Go drift floor: worst |drift| <= 1e-12", okD && wd <= 1e-12,
			fmt.Sprintf("%.3g", wd), "1e-12", true)
	}})

	// FB13 P2LC — the queue-#7 local-clock scheduler prototype (Go-only;
	// fab/localclock.go; FREECELL §2's four conditions on the REAL
	// substrate — 2D jammed bath, contact-rule channels, law-derived
	// pitch detunes; reduced execution dynamics as in v89/localclock.c).
	// Measured: R1 conservation at the FP floor in both engines; R2
	// first-order convergence (ratio 1.715 under R=8 dilation); R3
	// bit-identity under scan rotation+reversal EXACTLY 0 with the
	// (t,index) total order, arrival-order control fails at 1.45; batch
	// = serial EXACTLY 0 at 1/2/4/8 goroutines — with the PENDING-min
	// conflict rule (min-over-eligible measured 6.2 off: an earlier but
	// blocked neighbour must still hold you back — a sharpening of the
	// v89 rule, which was only exact under full eligibility); tick skew
	// 219 (condition 4: the counter diverges — order by local time);
	// quiet-region economy: event ratio 2.37/4.92/6.87 at L=16/32/64
	// with a fixed active blob, approaching the R=8 dilation bound.
	sp = []runSpec{}
	if wantGo {
		sp = append(sp, gw("p2lc", "exp=p2lc"))
	}
	exps = append(exps, experiment{"p2lc", sp, func(r *reporter) {
		if !wantGo {
			return
		}
		lg := get("p2lc")
		ds, ok1 := exf(lg, `p2lc_r1 drift_sync=(\S+)`)
		da, ok2 := exf(lg, `drift_async=(\S+)`)
		r.add("p2lc", "R1 conservation: both |drift| <= 1e-12", ok1 && ok2 &&
			abs(ds) <= 1e-12 && abs(da) <= 1e-12,
			fmt.Sprintf("%.1e/%.1e", ds, da), "<=1e-12", true)
		rt, ok := exf(lg, `p2lc_r2 err_dt=\S+ err_dt2=\S+ ratio=(\S+)`)
		r.add("p2lc", "R2 first-order: ratio in [1.5,2.4]", ok && rt >= 1.5 && rt <= 2.4,
			fmt.Sprintf("%.3f", rt), "[1.5,2.4]", true)
		ro, o1 := exf(lg, `p2lc_r3 rot=(\S+)`)
		rv, o2 := exf(lg, `rev=(\S+)`)
		rr, o3 := exf(lg, `revrot=(\S+)`)
		av, o4 := exf(lg, `arrival=(\S+)`)
		r.add("p2lc", "R3 total-order determinism: rot/rev/revrot == 0", o1 && o2 && o3 &&
			ro == 0 && rv == 0 && rr == 0,
			fmt.Sprintf("%g/%g/%g", ro, rv, rr), "== 0", true)
		r.add("p2lc", "R3 control: arrival order corrupts >= 0.1", o4 && av >= 0.1,
			fmt.Sprintf("%.3f", av), ">=0.1", true)
		bd, ok := exf(lg, `p2lc_batch maxdiff_vs_serial=(\S+)`)
		r.add("p2lc", "batch == serial schedule (1..8 workers)", ok && bd == 0,
			fmt.Sprintf("%g", bd), "== 0", true)
		sk, ok := exf(lg, `max_tick_skew=(\S+)`)
		r.add("p2lc", "cond 4: tick skew diverges (>= 100)", ok && sk >= 100,
			fmt.Sprintf("%.0f", sk), ">=100", true)
		var ratios []float64
		for _, m := range regexp.MustCompile(`ev_ratio=(\S+)`).FindAllStringSubmatch(lg, -1) {
			if v, err := strconv.ParseFloat(m[1], 64); err == nil {
				ratios = append(ratios, v)
			}
		}
		r.add("p2lc", "quiet-region economy grows with box", len(ratios) == 3 &&
			ratios[0] < ratios[1] && ratios[1] < ratios[2] && ratios[2] >= 6.0,
			fmt.Sprintf("%.2f/%.2f/%.2f", ratios[0], ratios[1], ratios[2]), "mono, >=6", true)
		wd := 0.0
		okD := true
		for _, m := range regexp.MustCompile(`p2lc_scale .*drift=(\S+)`).FindAllStringSubmatch(lg, -1) {
			v, err := strconv.ParseFloat(m[1], 64)
			if err != nil {
				okD = false
			}
			if abs(v) > wd {
				wd = abs(v)
			}
		}
		r.add("p2lc", "scaling arms |drift| <= 1e-12", okD && wd <= 1e-12,
			fmt.Sprintf("%.1e", wd), "<=1e-12", true)
	}})

	return exps
}

// maskDrift blanks the two places where the C and Go kernels are ALLOWED
// to differ: the conservation-drift column (FP-floor libm-ulp dust — the
// t-row's second field and the RESULT drift_rel line) and the sign of a
// printed negative zero. Everything else must be byte-equal.
var driftColRe = regexp.MustCompile(`(?m)^(t=\s*\S+) [+-]\S+`)
var driftResRe = regexp.MustCompile(`(?m)^# RESULT drift_rel \S+`)
// p1 meter accumulators (sp/fl/fd/gm) are C/Go-order-sensitive derived
// quantities at FP-floor (V3a adoption: blob2 fl diverges 1 ulp under
// radiance while every physics-state column stays byte-identical). Masked
// like drift — the t= physics-state lines carry the C/Go identity claim.
var p1ResRe = regexp.MustCompile(`(?m)^# RESULT p1[xyz] .*`)

func maskDrift(log string) string {
	log = driftColRe.ReplaceAllString(log, "$1 DRIFT")
	log = driftResRe.ReplaceAllString(log, "# RESULT drift_rel DRIFT")
	log = p1ResRe.ReplaceAllString(log, "# RESULT p1 (meter masked: C/Go FP envelope)")
	return strings.ReplaceAll(log, "-0.0000", "0.0000")
}

// sectGroups aggregates the SECT table into the declared angular groups:
// core |th|<=0.27 (3 bins), fwd |th|<=0.7854 (7), back |th|>=2.3562 (7),
// side the rest; E = exposure, N = time-integrated occupancy.
type sectAgg struct{ cE, fE, bE, sE, tE, cN, fN, bN, sN, tN float64 }

var sectRe = regexp.MustCompile(`(?m)^# SECT k=\d+ th=(\S+) E=(\S+) n=(\S+)`)

func sectGroups(label string) (sectAgg, bool) {
	var a sectAgg
	ms := sectRe.FindAllStringSubmatch(get(label), -1)
	if len(ms) == 0 {
		return a, false
	}
	for _, m := range ms {
		th, e1 := strconv.ParseFloat(m[1], 64)
		e, e2 := strconv.ParseFloat(m[2], 64)
		n, e3 := strconv.ParseFloat(m[3], 64)
		if e1 != nil || e2 != nil || e3 != nil {
			return a, false
		}
		ath := math.Abs(th)
		if ath <= 0.27 {
			a.cE += e
			a.cN += n
		}
		if ath <= 0.7854 {
			a.fE += e
			a.fN += n
		} else if ath >= 2.3562 {
			a.bE += e
			a.bN += n
		} else {
			a.sE += e
			a.sN += n
		}
		a.tE += e
		a.tN += n
	}
	return a, true
}

type convTagV struct{ rough, cond, evap, backs, net float64 }

func convTag(label string) (convTagV, bool) {
	m := regexp.MustCompile(`# RESULT convtag rough=(\S+) cond=(\S+) evap=(\S+) backs=(\S+) net=(\S+)`).
		FindStringSubmatch(get(label))
	if m == nil {
		return convTagV{}, false
	}
	var v [5]float64
	for i := 0; i < 5; i++ {
		f, err := strconv.ParseFloat(m[i+1], 64)
		if err != nil {
			return convTagV{}, false
		}
		v[i] = f
	}
	return convTagV{v[0], v[1], v[2], v[3], v[4]}, true
}

// p1Rate parses the all-modes first-moment rate vector from # RESULT p1.
func p1Rate(label string) (x, y, z float64, ok bool) {
	m := regexp.MustCompile(`# RESULT p1 tot=\([^)]*\) rate=\((\S+?),(\S+?),(\S+?)\)`).
		FindStringSubmatch(get(label))
	if m == nil {
		return 0, 0, 0, false
	}
	var e1, e2, e3 error
	x, e1 = strconv.ParseFloat(m[1], 64)
	y, e2 = strconv.ParseFloat(m[2], 64)
	z, e3 = strconv.ParseFloat(m[3], 64)
	return x, y, z, e1 == nil && e2 == nil && e3 == nil
}

// ------------------------------------------------------------------
// DS tier-1 click analysis
// ------------------------------------------------------------------

var clickRe = regexp.MustCompile(`(?m)^# CLICK t=(\S+) y=(\S+) e=(\S+)`)

type ds1Agg struct {
	n, nmin, nmax int
	sbar, sem     float64
	atomDev       float64
	minK, maxK    int
	maxDrift      float64
	okDrift       bool
}

// ds1Stats aggregates in-gate clicks over panel logs: fringe-phase score
// s = cos^2(pi*delta(y)/lambda) (delta from the same parameter-free loci
// as dsAnalyze: yA=27, yB=37, D=14, lambda=2pi/0.9), minima bands = pred
// minima +-1.5 ([25,28] u [36,39]), central max band [30,34], and the R3
// atom check against eps = A0eff*w1/2pi (laws_V2g: 1.15*1.65/2pi).
func ds1Stats(labels []string) ds1Agg {
	const yA, yB, D = 27.0, 37.0, 14.0
	lam := 2 * math.Pi / 0.9
	eps := 1.15 * 1.65 / (2 * math.Pi)
	a := ds1Agg{minK: 1 << 30, okDrift: true}
	ssum, s2 := 0.0, 0.0
	for _, lb := range labels {
		lg := get(lb)
		if v, ok := exf(lg, `# RESULT drift_rel (\S+)`); ok {
			if abs(v) > a.maxDrift {
				a.maxDrift = abs(v)
			}
		} else {
			a.okDrift = false
		}
		for _, m := range clickRe.FindAllStringSubmatch(lg, -1) {
			t, e1 := strconv.ParseFloat(m[1], 64)
			y, e2 := strconv.ParseFloat(m[2], 64)
			e, e3 := strconv.ParseFloat(m[3], 64)
			if e1 != nil || e2 != nil || e3 != nil || t < 20 || t > 36 {
				continue
			}
			delta := math.Sqrt(D*D+(y-yA)*(y-yA)) - math.Sqrt(D*D+(y-yB)*(y-yB))
			c := math.Cos(math.Pi * delta / lam)
			s := c * c
			a.n++
			ssum += s
			s2 += s * s
			if (y >= 25 && y <= 28) || (y >= 36 && y <= 39) {
				a.nmin++
			}
			if y >= 30 && y <= 34 {
				a.nmax++
			}
			k := int(e/eps + 0.5)
			if k < a.minK {
				a.minK = k
			}
			if k > a.maxK {
				a.maxK = k
			}
			if dev := abs(e - float64(k)*eps); dev > a.atomDev {
				a.atomDev = dev
			}
		}
	}
	if a.n > 0 {
		a.sbar = ssum / float64(a.n)
		sd := math.Sqrt(s2/float64(a.n) - a.sbar*a.sbar)
		a.sem = sd / math.Sqrt(float64(a.n))
	}
	if a.minK == 1<<30 {
		a.minK = 0
	}
	return a
}

func maxf(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}
func minint(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func maxint(a, b int) int {
	if a > b {
		return a
	}
	return b
}

var scRe = regexp.MustCompile(`(?m)^# SCREENCELL y=(\S+) x=\S+ I=(\S+)`)

type scPt struct{ y, i float64 }

func scCells(label string) []scPt {
	var out []scPt
	for _, m := range scRe.FindAllStringSubmatch(get(label), -1) {
		y, e1 := strconv.ParseFloat(m[1], 64)
		i, e2 := strconv.ParseFloat(m[2], 64)
		if e1 == nil && e2 == nil {
			out = append(out, scPt{y, i})
		}
	}
	return out
}

func scSmooth(pts []scPt, y float64) float64 {
	const s = 1.2
	num, den := 0.0, 0.0
	for _, p := range pts {
		w := math.Exp(-(y - p.y) * (y - p.y) / (2 * s * s))
		num += w * p.i
		den += w
	}
	if den <= 1e-9 {
		return 0
	}
	return num / den
}

// dsAnalyze: smoothed I profiles, exact loci from the path-difference
// scan (yA=27, yB=37, D=14, lam=2pi/0.9), r(y)=I_AB/(I_A+I_B).
func dsAnalyze(l0, l1, l2 string) (ypk, r0, rm, rp, V float64, ok bool) {
	AB, A, B := scCells(l0), scCells(l1), scCells(l2)
	if len(AB) < 10 || len(A) < 10 || len(B) < 10 {
		return 0, 0, 0, 0, 0, false
	}
	const yA, yB, D = 27.0, 37.0, 14.0
	lam := 2 * math.Pi / 0.9
	delta := func(y float64) float64 {
		return math.Sqrt(D*D+(y-yA)*(y-yA)) - math.Sqrt(D*D+(y-yB)*(y-yB))
	}
	locus := func(target float64) float64 {
		best, bv := 16.0, math.Inf(1)
		for y := 16.0; y <= 48.0; y += 0.02 {
			if d := math.Abs(delta(y) - target); d < bv {
				bv, best = d, y
			}
		}
		return best
	}
	yMin1, yMax0, yMin2 := locus(-lam/2), locus(0), locus(lam/2)
	rOf := func(y float64) float64 {
		iab, ia, ib := scSmooth(AB, y), scSmooth(A, y), scSmooth(B, y)
		if ia+ib <= 1e-9 {
			return 0
		}
		return iab / (ia + ib)
	}
	ypk, bv := 26.0, -1.0
	for y := 26.0; y <= 38.0; y += 0.1 {
		if v := scSmooth(AB, y); v > bv {
			bv, ypk = v, y
		}
	}
	r0, rm, rp = rOf(yMax0), rOf(yMin1), rOf(yMin2)
	V = (r0 - 0.5*(rm+rp)) / (r0 + 0.5*(rm+rp))
	return ypk, r0, rm, rp, V, true
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
	// standard streams for the heavier experiments (regenerable, gitignored)
	if err := os.MkdirAll(filepath.Join(*runsDir, "streams"), 0o755); err != nil {
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
