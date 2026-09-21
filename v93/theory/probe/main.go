// Command probe: parallel fabrun/freecell campaigns for shear genesis,
// progressive (Bohmian-method) meters, and long-horizon SUSTAIN isolation.
// Pure Go — no Python. Lives under v93/theory.
//
//	go build -o theory/bin/probe ./theory/probe
//	./theory/bin/probe -suite sustain -kernel fabrun -T 3000 -j 4
package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

type ArmResult struct {
	Arm     Arm         `json:"arm"`
	Log     LogResult   `json:"log"`
	Prog    ProgMeters  `json:"prog"`
	Series  []EmSample  `json:"series,omitempty"` // Em(t) from fcs stats
	Seconds float64     `json:"seconds"`
	Err     string      `json:"err,omitempty"`
}

type EmSample struct {
	T, EmTot, EmMax, SumEe float64
	NHot                   int
}

func main() {
	j := flag.Int("j", max(2, runtime.NumCPU()/4), "parallel sim workers")
	seed := flag.Int("seed", 20260802, "RNG seed")
	root := flag.String("root", "", "v93 root")
	out := flag.String("out", "", "output dir")
	only := flag.String("only", "", "family filter: shear,prog,ctrl,sust,genesis")
	nameFilter := flag.String("name", "", "substring filter on arm name")
	suite := flag.String("suite", "default", "default | sustain | r1 | mono | p2")
	kernel := flag.String("kernel", "fabrun", "fabrun (Go) or freecell (C)")
	T := flag.Float64("T", 3000, "horizon for sustain/r1/mono/p2 suite")
	skipFCS := flag.Bool("no-fcs-meters", false, "skip final progressive meters")
	snapEvery := flag.Int("snap-every", 5000, "fcs snapshot period in steps (dt=0.02 → 5000=100 t.u.)")
	flag.Parse()

	v93 := *root
	if v93 == "" {
		v93 = findV93()
	}
	bin := filepath.Join(v93, *kernel)
	if _, err := os.Stat(bin); err != nil {
		// fallback
		alt := filepath.Join(v93, "freecell")
		if *kernel == "fabrun" {
			if _, e2 := os.Stat(alt); e2 == nil {
				fmt.Fprintf(os.Stderr, "warn: fabrun missing, using freecell\n")
				bin = alt
			}
		}
		if _, err := os.Stat(bin); err != nil {
			fmt.Fprintf(os.Stderr, "kernel not found: %s\n", bin)
			os.Exit(1)
		}
	}
	fcsdump := filepath.Join(v93, "fcsdump")

	runDir := *out
	if runDir == "" {
		tag := *suite
		runDir = filepath.Join(v93, "theory", "runs",
			fmt.Sprintf("%s_%s", tag, time.Now().Format("20060102_150405")))
	}
	if err := os.MkdirAll(runDir, 0o755); err != nil {
		fatal(err)
	}

	var arms []Arm
	switch *suite {
	case "sustain":
		arms = buildSustainArms(*seed, *T)
	case "r1":
		if *T == 3000 {
			*T = 10000 // r1 default endurance horizon
		}
		arms = buildR1Arms(*seed, *T)
	case "mono":
		if *T == 3000 {
			*T = 2000 // mono default shorter scan
		}
		arms = buildMonoArms(*seed, *T)
	case "p2":
		if *T == 3000 {
			*T = 600 // P2 thermometer: ≥500 t.u. per bar
		}
		arms = buildP2Arms(*seed, *T)
	default:
		arms = buildArms(*seed)
	}
	if *only != "" {
		want := map[string]bool{}
		for _, f := range strings.Split(*only, ",") {
			want[strings.TrimSpace(f)] = true
		}
		var f []Arm
		for _, a := range arms {
			if want[a.Family] {
				f = append(f, a)
			}
		}
		arms = f
	}
	if *nameFilter != "" {
		parts := strings.Split(*nameFilter, "|")
		var f []Arm
		for _, a := range arms {
			for _, p := range parts {
				p = strings.TrimSpace(p)
				if p != "" && strings.Contains(a.Name, p) {
					f = append(f, a)
					break
				}
			}
		}
		arms = f
	}

	// meta
	meta := fmt.Sprintf("suite=%s kernel=%s T=%g seed=%d j=%d snap_every=%d n=%d\n",
		*suite, filepath.Base(bin), *T, *seed, *j, *snapEvery, len(arms))
	_ = os.WriteFile(filepath.Join(runDir, "META.txt"), []byte(meta), 0o644)

	fmt.Printf("# theory/probe\n# %s", meta)
	fmt.Printf("# out=%s\n", runDir)

	jobs := make(chan Arm, len(arms))
	var wg sync.WaitGroup
	var mu sync.Mutex
	var results []ArmResult

	worker := func() {
		defer wg.Done()
		for arm := range jobs {
			res := runArm(bin, fcsdump, runDir, arm, !*skipFCS, *snapEvery)
			mu.Lock()
			results = append(results, res)
			status := "ok"
			if res.Err != "" {
				status = "ERR"
			}
			emEnd := res.Prog.EmTot
			if len(res.Series) > 0 {
				emEnd = res.Series[len(res.Series)-1].EmTot
			}
			slope := sustainSlope(res.Series)
			fmt.Printf("[%s] %-22s cond=%.4g Em=%.4g slope=%+.3g/1k  fd=%+.3g Jx=%+.3g  %s (%.0fs)\n",
				status, arm.Name, res.Log.Cond, emEnd, slope*1000, res.Log.Fd, res.Prog.Jx,
				arm.Family, res.Seconds)
			// incremental save
			writeAll(runDir, results, *seed)
			mu.Unlock()
		}
	}

	nW := *j
	if nW < 1 {
		nW = 1
	}
	if nW > len(arms) && len(arms) > 0 {
		nW = len(arms)
	}
	wg.Add(nW)
	for i := 0; i < nW; i++ {
		go worker()
	}
	for _, a := range arms {
		jobs <- a
	}
	close(jobs)
	wg.Wait()

	sort.Slice(results, func(i, j int) bool {
		return results[i].Arm.Name < results[j].Arm.Name
	})
	writeAll(runDir, results, *seed)
	writeSustainReport(filepath.Join(runDir, "SUSTAIN_REPORT.md"), results, *T, *seed, bin)
	fmt.Printf("\n# done → %s\n", runDir)
	printSustainVerdict(results)
}

func runArm(bin, fcsdump, runDir string, arm Arm, doMeters bool, snapEvery int) ArmResult {
	res := ArmResult{Arm: arm}
	logPath := filepath.Join(runDir, arm.Name+".log")
	args := append([]string{}, arm.Args...)
	var fcs string
	if arm.Snap {
		fcs = fcsPath(runDir, arm.Name)
		args = append(args,
			fmt.Sprintf("snap_every=%d", snapEvery),
			"snap_file="+fcs,
		)
	}
	t0 := time.Now()
	cmd := exec.Command(bin, args...)
	cmd.Dir = filepath.Dir(bin)
	// stream log live (long T=3000 arms)
	logf, err := os.Create(logPath)
	if err != nil {
		res.Err = err.Error()
		return res
	}
	cmd.Stdout = logf
	cmd.Stderr = logf
	err = cmd.Run()
	_ = logf.Close()
	res.Seconds = time.Since(t0).Seconds()
	if err != nil {
		res.Err = err.Error()
	}
	if lr, e := parseLog(logPath); e == nil {
		res.Log = lr
	} else if res.Err == "" {
		res.Err = e.Error()
	}
	if arm.Snap && fcs != "" {
		res.Series = fcsEmSeries(fcsdump, fcs)
		if doMeters {
			res.Prog = measureFCS(fcsdump, fcs, 16)
		} else if len(res.Series) > 0 {
			last := res.Series[len(res.Series)-1]
			res.Prog.EmTot = last.EmTot
			res.Prog.EmMax = last.EmMax
			res.Prog.EeTot = last.SumEe
			res.Prog.OK = true
		}
	}
	// write series csv
	if len(res.Series) > 0 {
		var b strings.Builder
		b.WriteString("t,Em_tot,Em_max,Ee,nHot\n")
		for _, s := range res.Series {
			b.WriteString(fmt.Sprintf("%.2f,%.6g,%.6g,%.6g,%d\n",
				s.T, s.EmTot, s.EmMax, s.SumEe, s.NHot))
		}
		_ = os.WriteFile(filepath.Join(runDir, arm.Name+".em.csv"), []byte(b.String()), 0o644)
	}
	return res
}

// fcsEmSeries: parse fcsdump -mode stats for Em(t).
func fcsEmSeries(fcsdump, fcs string) []EmSample {
	if _, err := os.Stat(fcs); err != nil {
		return nil
	}
	cmd := exec.Command(fcsdump, "-mode", "stats", fcs)
	out, err := cmd.Output()
	if err != nil {
		return nil
	}
	var ser []EmSample
	sc := bufio.NewScanner(strings.NewReader(string(out)))
	for sc.Scan() {
		line := sc.Text()
		if !strings.HasPrefix(line, "CELL") {
			continue
		}
		// CELL t N sumEm meanEm varEm maxEm nHot sumEe ...
		f := strings.Fields(line)
		if len(f) < 9 {
			continue
		}
		t, _ := strconv.ParseFloat(f[1], 64)
		em, _ := strconv.ParseFloat(f[3], 64)
		emx, _ := strconv.ParseFloat(f[6], 64)
		nh, _ := strconv.Atoi(f[7])
		ee, _ := strconv.ParseFloat(f[8], 64)
		ser = append(ser, EmSample{T: t, EmTot: em, EmMax: emx, SumEe: ee, NHot: nh})
	}
	return ser
}

// sustainSlope: linear fit dEm/dt on last half of series (per unit t).
func sustainSlope(ser []EmSample) float64 {
	if len(ser) < 4 {
		return 0
	}
	half := len(ser) / 2
	s := ser[half:]
	var sx, sy, sxx, sxy float64
	n := float64(len(s))
	for _, p := range s {
		sx += p.T
		sy += p.EmTot
		sxx += p.T * p.T
		sxy += p.T * p.EmTot
	}
	den := n*sxx - sx*sx
	if den == 0 {
		return 0
	}
	return (n*sxy - sx*sy) / den
}

func writeAll(runDir string, results []ArmResult, seed int) {
	sort.Slice(results, func(i, j int) bool {
		return results[i].Arm.Name < results[j].Arm.Name
	})
	js, _ := json.MarshalIndent(results, "", "  ")
	_ = os.WriteFile(filepath.Join(runDir, "results.json"), js, 0o644)
	var b strings.Builder
	b.WriteString("name,family,cond,evap,rad,fd,fl,Em_end,Em_max,Em_mid,slope,Jx,R2f,Ee,seconds,err\n")
	for _, r := range results {
		emEnd, emMax, emMid := r.Prog.EmTot, r.Prog.EmMax, 0.0
		if n := len(r.Series); n > 0 {
			emEnd = r.Series[n-1].EmTot
			emMax = r.Series[n-1].EmMax
			emMid = r.Series[n/2].EmTot
		}
		b.WriteString(fmt.Sprintf("%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.1f,%s\n",
			r.Arm.Name, r.Arm.Family, r.Log.Cond, r.Log.Evap, r.Log.Rad, r.Log.Fd, r.Log.Fl,
			emEnd, emMax, emMid, sustainSlope(r.Series), r.Prog.Jx, r.Prog.R2f, r.Prog.EeTot,
			r.Seconds, strings.ReplaceAll(r.Err, ",", ";")))
	}
	_ = os.WriteFile(filepath.Join(runDir, "results.csv"), []byte(b.String()), 0o644)
}

func writeSustainReport(path string, results []ArmResult, T float64, seed int, bin string) {
	var b strings.Builder
	b.WriteString("# SUSTAIN REPORT\n\n")
	b.WriteString(fmt.Sprintf("T=%g seed=%d kernel=`%s` n=%d\n\n", T, seed, bin, len(results)))
	b.WriteString("## Criterion\n\n")
	b.WriteString("A condensate **sustains** if over the last half of [0,T]:\n")
	b.WriteString("- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;\n")
	b.WriteString("- not pure melt (Em→0 / PR→N);\n")
	b.WriteString("- preferably net book `cond - evap - rad` near zero if doors are open.\n\n")
	b.WriteString("| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |\n")
	b.WriteString("|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
	for _, r := range results {
		emEnd, emMax := r.Prog.EmTot, r.Prog.EmMax
		if n := len(r.Series); n > 0 {
			emEnd = r.Series[n-1].EmTot
			emMax = r.Series[n-1].EmMax
		}
		sl := sustainSlope(r.Series)
		ver := classifySustain(r)
		b.WriteString(fmt.Sprintf("| %s | %.1f | %.2f | %+.2f | %.1f | %.2f | %.2f | %+.1f | **%s** |\n",
			r.Arm.Name, emEnd, emMax, sl*1000, r.Log.Cond, r.Log.Evap, r.Log.Rad, r.Log.Fd, ver))
	}
	b.WriteString("\n## Isolation map\n\n")
	b.WriteString("| switch | arms |\n|---|---|\n")
	b.WriteString("| transport uni vs add | S00 vs S01 |\n")
	b.WriteString("| gate | S00 vs S02 |\n")
	b.WriteString("| q_detune trap | S03 / S00 / S04 |\n")
	b.WriteString("| radiance outtake | S00 vs S05; S04 vs S06 |\n")
	b.WriteString("| field fuel (noise) | S00 vs S07; S05 vs S08 |\n")
	b.WriteString("| workfn door | S00 vs S10; S07 vs S11 |\n")
	b.WriteString("| metabolic (open door+fuel+tax) | S12 |\n")
	b.WriteString("| progressive dense phase | S14, S15, S16 |\n")
	b.WriteString("| genesis long shear | G01 vs G02 |\n")
	b.WriteString("\n## Bohmian-method aspects\n\n")
	b.WriteString("See SUSTAIN_BOHMIAN.md / report body: progressive meters (J, fd, support), ")
	b.WriteString("not pilot particles. S14–S16 vary phase progression on dense amplitude.\n")
	_ = os.WriteFile(path, []byte(b.String()), 0o644)
}

func classifySustain(r ArmResult) string {
	if len(r.Series) < 4 {
		if r.Prog.EmTot > 10 {
			return "incomplete-series"
		}
		return "no-data"
	}
	n := len(r.Series)
	em0 := r.Series[0].EmTot
	emH := r.Series[n/2].EmTot
	em1 := r.Series[n-1].EmTot
	sl := sustainSlope(r.Series)
	// melt: ends low
	if em1 < 5 && em0 > 20 {
		return "MELT"
	}
	if em1 < 5 && em0 < 5 && r.Log.Cond < 1 {
		return "empty"
	}
	// runaway growth last half
	if sl > 0.05 { // Em units per t.u.
		return "GROW"
	}
	if sl < -0.05 {
		return "DRAIN"
	}
	// hold: last half slope small and Em substantial
	if em1 > 20 && abs(sl) <= 0.05 {
		// distinguish plateau vs slow
		if em1 > 0.5*emH && em1 < 2*emH {
			return "HOLD"
		}
		return "quasi-hold"
	}
	if em1 > 20 {
		return "alive"
	}
	return "weak"
}

func printSustainVerdict(results []ArmResult) {
	fmt.Println("\n# SUSTAIN VERDICT")
	counts := map[string]int{}
	for _, r := range results {
		counts[classifySustain(r)]++
	}
	for k, v := range counts {
		fmt.Printf("  %s: %d\n", k, v)
	}
	// list HOLD
	for _, r := range results {
		if classifySustain(r) == "HOLD" {
			em := r.Prog.EmTot
			if n := len(r.Series); n > 0 {
				em = r.Series[n-1].EmTot
			}
			fmt.Printf("  HOLD %s Em=%.1f slope=%+.4g\n", r.Arm.Name, em, sustainSlope(r.Series))
		}
	}
}

func findV93() string {
	wd, _ := os.Getwd()
	for d := wd; d != "/" && d != "."; d = filepath.Dir(d) {
		if _, err := os.Stat(filepath.Join(d, "fabrun")); err == nil {
			return d
		}
		if _, err := os.Stat(filepath.Join(d, "freecell")); err == nil {
			return d
		}
	}
	return wd
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, err)
	os.Exit(1)
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
