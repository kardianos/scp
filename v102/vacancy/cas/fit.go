// fit.go — v102 vacancy step (../PREREG.md): the gate, the tables and the
// fits, read from the outputs in out/.  Usage: go run fit.go t0|q3|q2|q4
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

const J0 = 0.0095861869

// ---------- parsing ----------

type Run struct {
	E      []float64 // rate units
	Res    []float64
	S2, I2 map[int]float64
	Lines  []string
	D      int
}

func readRun(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{S2: map[int]float64{}, I2: map[int]float64{}}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	for sc.Scan() {
		ln := sc.Text()
		r.Lines = append(r.Lines, ln)
		fs := strings.Fields(ln)
		if len(fs) >= 3 && fs[0] == "E" {
			v, _ := strconv.ParseFloat(fs[2], 64)
			r.E = append(r.E, v)
			if k := idx(fs, "resid"); k >= 0 {
				x, _ := strconv.ParseFloat(fs[k+1], 64)
				r.Res = append(r.Res, x)
			}
		}
		if len(fs) >= 8 && fs[0] == "state" && fs[2] == "E" {
			v, _ := strconv.Atoi(fs[1])
			s2, _ := strconv.ParseFloat(fs[6], 64)
			i2, _ := strconv.ParseFloat(fs[8], 64)
			r.S2[v], r.I2[v] = s2, i2
		}
		if strings.HasPrefix(ln, "# vx") {
			for _, w := range fs {
				if strings.HasPrefix(w, "D=") {
					r.D, _ = strconv.Atoi(w[2:])
				}
			}
		}
	}
	if len(r.E) == 0 {
		return nil
	}
	return r
}

func idx(fs []string, k string) int {
	for i, f := range fs {
		if f == k {
			return i
		}
	}
	return -1
}

// the numbers after a line prefix such as "state 0 corr SzSz"
func series(r *Run, prefix string) []float64 {
	for _, ln := range r.Lines {
		if strings.HasPrefix(ln, prefix+" ") {
			var out []float64
			for _, w := range strings.Fields(ln[len(prefix):]) {
				v, err := strconv.ParseFloat(w, 64)
				if err == nil {
					out = append(out, v)
				}
			}
			return out
		}
	}
	return nil
}

func lineSeries(path, prefix string) []float64 {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	for sc.Scan() {
		ln := sc.Text()
		if strings.HasPrefix(ln, prefix+" ") {
			var out []float64
			for _, w := range strings.Fields(ln[len(prefix):]) {
				if v, err := strconv.ParseFloat(w, 64); err == nil {
					out = append(out, v)
				}
			}
			return out
		}
	}
	return nil
}

func firstLine(path string) string {
	f, err := os.Open(path)
	if err != nil {
		return ""
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	if sc.Scan() {
		return sc.Text()
	}
	return ""
}

// ---------- fits ----------

// least-squares slope of y on x
func slope(x, y []float64) float64 {
	n := float64(len(x))
	var sx, sy, sxx, sxy float64
	for i := range x {
		sx += x[i]
		sy += y[i]
		sxx += x[i] * x[i]
		sxy += x[i] * y[i]
	}
	return (n*sxy - sx*sy) / (n*sxx - sx*sx)
}

// pure exponential: -1/slope of log|a_d|, d in [lo, hi] (skipping |a| < 1e-12)
func decayExp(a []float64, lo, hi int) float64 {
	var x, y []float64
	for d := lo; d <= hi && d < len(a); d++ {
		if math.Abs(a[d]) > 1e-12 {
			x = append(x, float64(d))
			y = append(y, math.Log(math.Abs(a[d])))
		}
	}
	if len(x) < 2 {
		return math.NaN()
	}
	return -1 / slope(x, y)
}

// Ornstein-Zernike: -1/slope of log(|C(r)| sqrt r), r in [lo, hi]
func decayOZ(c []float64, lo, hi int) float64 {
	var x, y []float64
	for r := lo; r <= hi && r < len(c); r++ {
		if math.Abs(c[r]) > 1e-14 {
			x = append(x, float64(r))
			y = append(y, math.Log(math.Abs(c[r])*math.Sqrt(float64(r))))
		}
	}
	if len(x) < 2 {
		return math.NaN()
	}
	return -1 / slope(x, y)
}

// ---------- T0 ----------

func t0() {
	type chk struct {
		file string
		what string
		got  func(r *Run) float64
		want float64
		dig  float64 // half a unit in the last printed digit
	}
	J := func(i int) func(r *Run) float64 { return func(r *Run) float64 { return r.E[i] / J0 } }
	rel := func(i int) func(r *Run) float64 { return func(r *Run) float64 { return (r.E[i] - r.E[0]) / J0 } }
	raw := func(i int) func(r *Run) float64 { return func(r *Run) float64 { return r.E[i] } }
	var e0n8 float64
	if r := readRun("out/t0_ao_32n8_e0.out"); r != nil {
		e0n8 = r.E[0]
	}
	maxBand := func(lim float64) func(r *Run) float64 {
		return func(r *Run) float64 {
			m := 0.0
			for _, e := range r.E {
				x := (e - r.E[0]) / J0
				if x < lim && x > m {
					m = x
				}
			}
			return m
		}
	}
	cs := []chk{
		{"out/t0_lo_W6_t0.out", "LO W n=6 mobile t=0: E0 (J0)", J(0), -20.51293, 5e-6},
		{"out/t0_lo_W6_t0_31.out", "LO W n=6 t=0: lowest (3/2,1/2) - E0, via sector (3,1) (J0)", func(r *Run) float64 { return r.E[0]/J0 - (-20.51292595) }, 1.9917, 5e-5},
		{"out/t0_lo_W6_t0.1.out", "LO W n=6 mobile t=0.1: E0 (J0)", J(0), -20.73466, 5e-6},
		{"out/t0_lo_W6_t0.1.out", "LO W n=6 t=0.1: band width (J0)", maxBand(1.0), 0.38794, 5e-6},
		{"out/t0_lo_32n8_t0.out", "LO 3+2 n=8 mobile t=0: E0 (J0)", J(0), -19.66294, 5e-6},
		{"out/t0_lo_32n8_t0.1.out", "LO 3+2 n=8 mobile t=0.1: E0 (J0)", J(0), -19.72715, 5e-6},
		{"out/t0_lo_32n8_t0.1.out", "LO 3+2 n=8 t=0.1: band width (J0)", maxBand(0.5), 0.12258, 5e-6},
		{"out/t0_ao_32n4.out", "AO 3+2 n=4 sector 0: level 0 (rate)", raw(0), -0.091490, 5e-7},
		{"out/t0_ao_32n4.out", "AO 3+2 n=4 sector 0: level 1", raw(1), -0.081985, 5e-7},
		{"out/t0_ao_32n4.out", "AO 3+2 n=4 sector 0: level 2", raw(2), -0.063837, 5e-7},
		{"out/t0_ao_W4.out", "AO W n=4 sector (0,0): level 0", raw(0), -0.127331, 5e-7},
		{"out/t0_ao_W4.out", "AO W n=4 sector (0,0): level 1", raw(1), -0.108256, 5e-7},
		{"out/t0_ao_32n8_rip.out", "AO 3+2 n=8 sum c<=3: Delta_8 (rate)", func(r *Run) float64 { return r.E[0] - e0n8 }, 0.005549, 5e-7},
		{"out/t0_ao_32n8_sing.out", "AO 3+2 n=8: lowest S=0 excitation (J0)", rel(1), 0.024014 / J0, 5e-7 / J0},
	}
	_ = J
	all := true
	for _, c := range cs {
		r := readRun(c.file)
		if r == nil {
			fmt.Printf("  MISSING %s\n", c.file)
			all = false
			continue
		}
		g := c.got(r)
		ok := math.Abs(g-c.want) <= c.dig
		if !ok {
			all = false
		}
		fmt.Printf("  %-62s got %+.6f  recorded %+.6f  %s\n", c.what, g, c.want, map[bool]string{true: "ok", false: "MISS"}[ok])
	}
	if r := readRun("out/t0_ao_32n8_sing.out"); r != nil {
		rip := readRun("out/t0_ao_32n8_rip.out")
		R := (r.E[1] - r.E[0]) / (2 * (rip.E[0] - e0n8))
		ok := math.Abs(R-2.1639) <= 5e-5
		if !ok {
			all = false
		}
		fmt.Printf("  %-62s got %.5f  recorded 2.1639  %s\n", "AO 3+2 n=8: ratio to 2 Delta_8", R, map[bool]string{true: "ok", false: "MISS"}[ok])
	}
	fmt.Printf("T0: %s\n", map[bool]string{true: "PASS (every recorded value reproduced to its printed digits)", false: "FAIL"}[all])
}

// ---------- Q3 ----------

type q3pt struct {
	n                    int
	E0, Del, Es, R, Eany float64
	s2, i2               float64
	ripS2, ripI2         float64
	ok                   bool
}

func q3cell(model string, M int, geom string) []q3pt {
	var pts []q3pt
	files, _ := filepath.Glob(fmt.Sprintf("out/q3_%s_M%d_%s_n*_e0.out", model, M, geom))
	for _, f := range files {
		base := strings.TrimSuffix(f, "_e0.out")
		var n int
		fmt.Sscanf(base[strings.LastIndex(base, "_n")+2:], "%d", &n)
		e0, rp, sg := readRun(f), readRun(base+"_rip.out"), readRun(base+"_sing.out")
		if e0 == nil || rp == nil || sg == nil || len(sg.E) < 2 {
			continue
		}
		p := q3pt{n: n, E0: e0.E[0]}
		p.Del = rp.E[0] - p.E0
		p.Es = sg.E[1] - p.E0
		p.R = p.Es / (2 * p.Del)
		if len(e0.E) > 1 {
			p.Eany = e0.E[1] - p.E0
		}
		p.s2, p.i2 = sg.S2[1], sg.I2[1]
		p.ripS2, p.ripI2 = rp.S2[0], rp.I2[0]
		p.ok = math.Abs(sg.E[0]-p.E0) < 1e-9 && math.Abs(p.s2) < 1e-6 && math.Abs(p.i2) < 1e-6
		pts = append(pts, p)
	}
	sort.Slice(pts, func(a, b int) bool { return pts[a].n < pts[b].n })
	return pts
}

func q3() {
	fmt.Println("# Q3: R = eps_s / (2 Delta_n); eps_s = lowest S=0 (W: (S,I)=(0,0)) excitation, Delta_n = lowest ripple; energies in J0")
	verdicts := map[string]string{}
	for _, M := range []int{2, 4} {
		for _, model := range []string{"lo", "ao"} {
			for _, geom := range []string{"ring", "open"} {
				pts := q3cell(model, M, geom)
				fab := map[int]string{2: "3+2", 4: "W"}[M]
				key := fmt.Sprintf("%s %s %s", fab, strings.ToUpper(model), geom)
				fmt.Printf("  %s:\n", key)
				for _, p := range pts {
					sz := fmt.Sprintf("n=%d", p.n)
					if geom == "open" {
						sz = fmt.Sprintf("L_bulk=%d", p.n-2)
					}
					fmt.Printf("    %-10s Delta %.5f  eps_s %.5f  R %.4f  (lowest excitation of the zero sector, any (S,I): %.5f; ripple <S2> %.3f <I2> %.3f; eps_s state <S2> %.1e <I2> %.1e %s)\n",
						sz, p.Del/J0, p.Es/J0, p.R, p.Eany/J0, p.ripS2, p.ripI2, p.s2, p.i2, map[bool]string{true: "", false: "CHECK"}[p.ok])
				}
				v := "unclear"
				if len(pts) >= 2 {
					a, b := pts[len(pts)-2], pts[len(pts)-1]
					if b.R >= 1 {
						v = "no line"
					} else if a.R < 1 && (1-b.R) >= 0.5*(1-a.R) {
						v = "line"
					}
				}
				verdicts[key] = v
				fmt.Printf("    verdict: %s\n", v)
			}
		}
	}
	nl, nw := 0, 0
	for k, v := range verdicts {
		if strings.HasPrefix(k, "W ") {
			nw++
			if v == "line" {
				nl++
			}
		}
	}
	// the ripple beside a particle: W LO ring n=8 with one static vacancy, labelled levels of sector (1,1)
	if v, b0, br := readRun("out/q3v_lo_M4_n8_static.out"), readRun("out/q3_lo_M4_ring_n8_e0.out"), readRun("out/q3_lo_M4_ring_n8_rip.out"); v != nil && b0 != nil && br != nil {
		del := br.E[0] - b0.E[0]
		fmt.Printf("  ripple beside a vacancy (LO W ring n=8, static vacancy, sector (1,1)); threshold = vacancy ground + Delta_8 = +%.5f J0\n", del/J0)
		for i := range v.E {
			s2, i2 := v.S2[i], v.I2[i]
			fmt.Printf("    level %2d  %+.5f J0  <S2> %.4f  <I2> %.4f  %s\n", i, (v.E[i]-v.E[0])/J0, s2, i2,
				map[bool]string{true: "below threshold", false: ""}[i > 0 && v.E[i]-v.E[0] < del-1e-9])
		}
		if vr := readRun("out/q3v_lo_M4_n8_static_rip.out"); vr != nil {
			fmt.Printf("    lowest (3/2,3/2) (sector (3,3)): %+.5f J0\n", (vr.E[0]-v.E[0])/J0)
		}
	}
	end := "0 (no W cell is a line: the spin-0, isospin-0 ripple falls apart into two ripples)"
	if nl == nw && nw > 0 {
		end = "H (every W cell is a line)"
	} else if nl > 0 {
		end = "h (some W cell is a line)"
	}
	fmt.Printf("Q3 ending: %s\n", end)
}

// ---------- Q2 ----------

func q2() {
	fmt.Println("# Q2: the vacancy's size against the background's correlation length")
	// W background, long chains
	for _, L := range []int{24, 32} {
		f := fmt.Sprintf("out/q2_W_nohole_%d.out", L)
		c := lineSeries(f, "SzSz")
		ci := lineSeries(f, "IzIz")
		if c == nil {
			fmt.Printf("  W background L_bulk=%d: missing\n", L)
			continue
		}
		hi := 12
		if L == 24 {
			hi = 8
		}
		fmt.Printf("  W background L_bulk=%d (%s)\n", L, strings.Fields(firstLine(f))[4])
		fmt.Printf("    <Sz Sz>(r) r=0..14:")
		for r := 0; r < len(c) && r <= 14; r++ {
			fmt.Printf(" %+.3e", c[r])
		}
		fmt.Println()
		fmt.Printf("    xi_bg (OZ, r in [3,%d]) Sz %.4f  Iz %.4f;  pure exponential, same window: Sz %.4f\n", hi, decayOZ(c, 3, hi), decayOZ(ci, 3, hi), decayExp(c, 3, hi))
		fmt.Printf("    OZ window scan (report): ")
		for _, w := range [][2]int{{2, 8}, {3, 10}, {4, 12}, {5, 14}} {
			fmt.Printf(" [%d,%d] %.4f", w[0], w[1], decayOZ(c, w[0], w[1]))
		}
		fmt.Println()
	}
	// W mobile vacancy: the registered cell runs R12 ('score'); the first run (p11 = zero, pre-R12) is a report
	for _, mf := range []string{"out/q2_W_mob_24_t0.1_R12.out", "out/q2_W_mob_24_t0.1.out"} {
		wmob(mf)
	}
	// 3+2 all-order static vacancy at g/Delta 0.05, 0.2
	q2rest()
}

func wmob(mf string) {
	if nh, C := lineSeries(mf, "nh"), lineSeries(mf, "C"); nh != nil && C != nil {
		n := len(nh)
		cc := make([]float64, n)
		for d := 0; d < n; d++ {
			w := 0.0
			for i := 0; i < n; i++ {
				k := 0.0
				if i+d < n {
					k++
				}
				if i-d >= 0 && d > 0 {
					k++
				}
				w += nh[i] * k
			}
			cc[d] = C[d] / math.Max(w, 1e-12)
		}
		L := 24
		fmt.Printf("  W mobile vacancy t=0.1 L_bulk=24 (%s): conditional profile", strings.Fields(firstLine(mf))[7])
		for d := 0; d < 8; d++ {
			fmt.Printf(" %+.4f", cc[d])
		}
		fmt.Printf("; decay length (d=2..%d) %.4f; raw C(d) %.4f\n", (L+1)/2-3, decayExp(cc, 2, (L+1)/2-3), decayExp(C, 2, (L+1)/2-3))
	} else {
		fmt.Printf("  W mobile vacancy (%s): missing\n", mf)
	}
}

func q2rest() {
	// 3+2 all-order static vacancy at g/Delta 0.05, 0.2
	for _, g := range []string{"0.05", "0.2"} {
		f := "out/q2_32_ao_hole_40_g" + g + ".out"
		p := lineSeries(f, "profile")
		if p == nil {
			fmt.Printf("  3+2 AO vacancy g/D=%s: missing\n", g)
			continue
		}
		Lh := 20
		h := 1 + Lh
		a := make([]float64, h+1)
		a[0] = p[h]
		for d := 1; d <= h; d++ {
			a[d] = p[h-d] + p[h+d]
		}
		fmt.Printf("  3+2 AO static vacancy g/D=%s L_bulk=40: decay length (d=2..%d) %.4f\n", g, Lh-8, decayExp(a, 2, Lh-8))
	}
	// ring scale: the vacancy profile's decay against the background correlation's decay, same window [1, n/2]
	fmt.Println("  ring scale (window d, r in [1, n/2]):")
	for _, c := range []struct {
		M, n int
	}{{2, 8}, {4, 6}} {
		for _, kind := range []string{"static", "mobile"} {
			var rat [2]float64
			for k, model := range []string{"lo", "ao"} {
				bg := readRun(fmt.Sprintf("out/q3_%s_M%d_ring_n%d_e0.out", model, c.M, c.n))
				vr := readQ2r(model, c.M, c.n, kind)
				if bg == nil || vr == nil {
					fmt.Printf("    %s M=%d n=%d %s: missing\n", model, c.M, c.n, kind)
					continue
				}
				cz := series(bg, "state 0 corr SzSz")
				var a []float64
				if kind == "static" {
					p := series(vr, "state 0 profile Sz")
					a = make([]float64, c.n/2+1)
					for d := 0; d <= c.n/2; d++ {
						a[d] = (p[(0+d)%c.n] + p[(c.n-d)%c.n]) / 2
					}
				} else {
					a = series(vr, "state 0 cond C(d)")
				}
				lv, xb := decayExp(a, 1, c.n/2), decayExp(cz, 1, c.n/2)
				rat[k] = lv / xb
				fmt.Printf("    %s %s n=%d %-6s vacancy decay %.4f  background decay %.4f  ratio %.4f\n", strings.ToUpper(model), map[int]string{2: "3+2", 4: "W"}[c.M], c.n, kind, lv, xb, rat[k])
			}
			if rat[0] != 0 && rat[1] != 0 {
				ok := math.Abs(rat[1]/rat[0]-1) <= 0.10
				fmt.Printf("      AO/LO ratio %.4f -> %s\n", rat[1]/rat[0], map[bool]string{true: "keeps", false: "does not keep"}[ok])
			}
		}
	}
	// supplement (not registered): the profiles themselves, AO against LO, each divided by its d = 1 value
	fmt.Println("  supplement, ring scale: vacancy profile a(d)/a(1) and background C(r)/C(1), LO | AO")
	for _, c := range []struct{ M, n int }{{2, 8}, {4, 6}} {
		for _, kind := range []string{"static", "mobile"} {
			var prof [2][]float64
			for k, model := range []string{"lo", "ao"} {
				vr := readQ2r(model, c.M, c.n, kind)
				if vr == nil {
					continue
				}
				if kind == "static" {
					p := series(vr, "state 0 profile Sz")
					for d := 0; d <= c.n/2; d++ {
						prof[k] = append(prof[k], (p[d%c.n]+p[(c.n-d)%c.n])/2)
					}
				} else {
					prof[k] = series(vr, "state 0 cond C(d)")
				}
			}
			if prof[0] == nil || prof[1] == nil {
				continue
			}
			fmt.Printf("    %s n=%d %-6s vacancy:", map[int]string{2: "3+2", 4: "W"}[c.M], c.n, kind)
			mx := 0.0
			for d := 1; d <= c.n/2; d++ {
				a, b := prof[0][d]/prof[0][1], prof[1][d]/prof[1][1]
				fmt.Printf("  d=%d %+.4f|%+.4f", d, a, b)
				mx = math.Max(mx, math.Abs(b-a))
			}
			fmt.Printf("   max |AO-LO| %.4f\n", mx)
		}
		var bg [2][]float64
		for k, model := range []string{"lo", "ao"} {
			if r := readRun(fmt.Sprintf("out/q3_%s_M%d_ring_n%d_e0.out", model, c.M, c.n)); r != nil {
				bg[k] = series(r, "state 0 corr SzSz")
			}
		}
		if bg[0] != nil && bg[1] != nil {
			fmt.Printf("    %s n=%d background:     ", map[int]string{2: "3+2", 4: "W"}[c.M], c.n)
			for r := 1; r <= c.n/2; r++ {
				fmt.Printf("  r=%d %+.4f|%+.4f", r, bg[0][r]/bg[0][1], bg[1][r]/bg[1][1])
			}
			fmt.Println()
		}
	}
}

// ---------- Q4 ----------

func q4() {
	fmt.Println("# Q4: |CS> = prod_l A^+_l |0> (one converted space quantum per link)")
	files, _ := filepath.Glob("out/q3_*_e0.out")
	more, _ := filepath.Glob("out/q2r_*_static.out")
	files = append(files, more...)
	sort.Strings(files)
	for _, f := range files {
		r := readRun(f)
		if r == nil {
			continue
		}
		for _, ln := range r.Lines {
			if strings.Contains(ln, " cs mode") {
				fmt.Printf("  %-40s %s\n", strings.TrimSuffix(filepath.Base(f), ".out"), strings.TrimPrefix(ln, "state 0 "))
			}
		}
		if c := series(r, "state 0 link c"); c != nil {
			s := 0.0
			for _, x := range c {
				s += x
			}
			N := series(r, "state 0 profile N")
			fmt.Printf("  %-40s space per link: mean %.5f  max %.5f;  <N> %v\n", "", s/float64(len(c)), maxf(c), fmtv(N))
		}
	}
}

func maxf(a []float64) float64 {
	m := math.Inf(-1)
	for _, x := range a {
		m = math.Max(m, x)
	}
	return m
}
func fmtv(a []float64) string {
	var s []string
	for _, x := range a {
		s = append(s, fmt.Sprintf("%.4f", x))
	}
	return strings.Join(s, " ")
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("usage: fit t0|q3|q2|q4")
		return
	}
	switch os.Args[1] {
	case "t0":
		t0()
	case "q3":
		q3()
	case "q2":
		q2()
	case "q4":
		q4()
	case "q3s":
		q3s()
	}
}

// the ring-scale vacancy run; the AO 3+2 mobile cell runs at sum c <= 3 (the full closure did not fit)
func readQ2r(model string, M, n int, kind string) *Run {
	if r := readRun(fmt.Sprintf("out/q2r_%s_M%d_n%d_%s.out", model, M, n, kind)); r != nil {
		return r
	}
	return readRun(fmt.Sprintf("out/q2r_%s_M%d_n%d_%s_c3.out", model, M, n, kind))
}

// the labelled supplement: the next size in two W cells
func q3s() {
	for _, b := range []string{"q3s_ao_M4_ring_n8", "q3s_lo_M4_open_n9"} {
		e0, rp, sg := readRun("out/"+b+"_e0.out"), readRun("out/"+b+"_rip.out"), readRun("out/"+b+"_sing.out")
		if e0 == nil || rp == nil || sg == nil || len(sg.E) < 2 {
			fmt.Printf("  %s: missing\n", b)
			continue
		}
		del, es := rp.E[0]-e0.E[0], sg.E[1]-e0.E[0]
		fmt.Printf("  %s (D=%d): Delta %.5f  eps_s %.5f  R %.4f  (eps_s state <S2> %.1e <I2> %.1e; E0 match %.1e)\n",
			b, e0.D, del/J0, es/J0, es/(2*del), sg.S2[1], sg.I2[1], sg.E[0]-e0.E[0])
	}
}
