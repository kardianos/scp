// fit.go — v103 lane c1: pair potentials, fits, signs, anemia and the
// registered ending (../PREREG.md §5) from the chain outputs in out/.
//
//	go run fit.go
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

var rules = []string{"none", "a", "b", "c", "ab", "ac", "bc", "abc"}

type Bin struct {
	Shell int
	R     float64
	Opp   bool
	Ideal float64
	Cnt   []float64
	G     []float64
}

type Run struct {
	Name, Lat      string
	N, Dim, Bip, D int
	RB             bool
	Bins           []*Bin
	MGrid          []float64
	Anemia         map[float64][]float64
	Ledger         string
	Acc            string
}

func atof(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }

var hdr = regexp.MustCompile(`lat=(\S+) .* N=(\d+) .* dim=(\d+) bipartite=(\d+) rules=\S+ \(a=\d b=(\d) c=\d\) CAP=\d q=\d D=(\d+)`)

func parse(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{Name: path, Anemia: map[float64][]float64{}}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	done := false
	for sc.Scan() {
		l := sc.Text()
		switch {
		case strings.HasPrefix(l, "# chain"):
			m := hdr.FindStringSubmatch(l)
			if m == nil {
				return nil
			}
			r.Lat = m[1]
			r.N, _ = strconv.Atoi(m[2])
			r.Dim, _ = strconv.Atoi(m[3])
			r.Bip, _ = strconv.Atoi(m[4])
			r.RB = m[5] == "1"
			r.D, _ = strconv.Atoi(m[6])
		case strings.HasPrefix(l, "MGRID"):
			for _, x := range strings.Fields(l)[1:] {
				r.MGrid = append(r.MGrid, atof(x))
			}
		case strings.HasPrefix(l, "SHELL"):
			fs := strings.Fields(l)
			b := &Bin{}
			b.Shell, _ = strconv.Atoi(fs[1])
			b.R = atof(strings.TrimPrefix(fs[2], "r="))
			b.Opp = fs[3] == "type=opp"
			b.Ideal = atof(strings.TrimPrefix(fs[4], "ideal="))
			i := 6
			for ; fs[i] != "G"; i++ {
				b.Cnt = append(b.Cnt, atof(fs[i]))
			}
			for _, x := range fs[i+1:] {
				b.G = append(b.G, atof(x))
			}
			r.Bins = append(r.Bins, b)
		case strings.HasPrefix(l, "ANEMIA"):
			fs := strings.Fields(l)
			rr := atof(strings.TrimPrefix(fs[2], "r="))
			for _, x := range fs[3:] {
				r.Anemia[rr] = append(r.Anemia[rr], atof(x))
			}
		case strings.HasPrefix(l, "LEDGER"):
			r.Ledger = l
		case strings.HasPrefix(l, "ACC"):
			r.Acc = l
		case strings.HasPrefix(l, "CHECK"):
			done = true
		}
	}
	if !done {
		return nil
	}
	return r
}

// allowed: the bins the rule set can visit (b on and bipartite: opposite only for D=2)
func allowed(r *Run, b *Bin, typ int) bool {
	if typ >= 0 && (typ == 1) != b.Opp {
		return false
	}
	if r.D == 2 && r.RB && r.Bip == 1 && !b.Opp {
		return false
	}
	return b.Ideal > 0
}

type VB struct {
	R, V, SE float64
	G        []float64
	Opp      bool
}

// potential per bin: V = -ln(P/Pideal), per seed, mean and se
func potential(r *Run, typ int) []VB {
	ns := len(r.Bins[0].Cnt)
	var bins []*Bin
	for _, b := range r.Bins {
		if allowed(r, b, typ) {
			bins = append(bins, b)
		}
	}
	idealTot := 0.0
	for _, b := range bins {
		idealTot += b.Ideal
	}
	tot := make([]float64, ns)
	for _, b := range bins {
		for s := 0; s < ns; s++ {
			tot[s] += b.Cnt[s]
		}
	}
	var out []VB
	for _, b := range bins {
		vs := make([]float64, 0, ns)
		ok := true
		for s := 0; s < ns; s++ {
			if b.Cnt[s] <= 0 || tot[s] <= 0 {
				ok = false
				break
			}
			vs = append(vs, -math.Log((b.Cnt[s]/tot[s])/(b.Ideal/idealTot)))
		}
		if !ok {
			continue
		}
		m, se := meanSE(vs)
		out = append(out, VB{R: b.R, V: m, SE: se, G: b.G, Opp: b.Opp})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].R < out[j].R })
	return out
}

func meanSE(v []float64) (float64, float64) {
	n := float64(len(v))
	m := 0.0
	for _, x := range v {
		m += x
	}
	m /= n
	s := 0.0
	for _, x := range v {
		s += (x - m) * (x - m)
	}
	return m, math.Sqrt(s / (n - 1) / n)
}

// weighted fit V = alpha + beta*x; returns chi2, alpha, beta, beta_se
func wfit(vb []VB, x func(VB) float64, withBeta bool) (float64, float64, float64, float64) {
	var S, Sx, Sy, Sxx, Sxy float64
	for _, p := range vb {
		w := 1 / (p.SE * p.SE)
		xv := x(p)
		S += w
		Sx += w * xv
		Sy += w * p.V
		Sxx += w * xv * xv
		Sxy += w * xv * p.V
	}
	var a, b, bse float64
	if withBeta {
		det := S*Sxx - Sx*Sx
		a = (Sxx*Sy - Sx*Sxy) / det
		b = (S*Sxy - Sx*Sy) / det
		bse = math.Sqrt(S / det)
	} else {
		a = Sy / S
	}
	chi := 0.0
	for _, p := range vb {
		f := a + b*x(p)
		chi += (p.V - f) * (p.V - f) / (p.SE * p.SE)
	}
	return chi, a, b, bse
}

type FitRes struct {
	NBins                   int
	Chi0, Chi1, Chi2        float64
	Beta, BetaSE, MBest, Xi float64
	Long                    bool
	Verdict                 string
}

func fitE2(r *Run) FitRes {
	vb := potential(r, -1)
	var use []VB
	for _, p := range vb {
		if p.R >= 2-1e-9 && p.SE > 0 {
			use = append(use, p)
		}
	}
	var fr FitRes
	fr.NBins = len(use)
	if len(use) < 4 {
		fr.Verdict = "too few bins"
		return fr
	}
	fr.Chi0, _, _, _ = wfit(use, func(p VB) float64 { return 0 }, false)
	fr.Chi1, _, fr.Beta, fr.BetaSE = wfit(use, func(p VB) float64 { return p.G[0] }, true)
	fr.Chi2 = fr.Chi1
	fr.MBest = 0
	for m := 1; m < len(r.MGrid); m++ {
		mm := m
		c, _, _, _ := wfit(use, func(p VB) float64 { return p.G[mm] }, true)
		if c < fr.Chi2 {
			fr.Chi2 = c
			fr.MBest = r.MGrid[m]
		}
	}
	fr.Xi = math.Inf(1)
	if fr.MBest > 0 {
		fr.Xi = 1 / fr.MBest
	}
	L := math.Pow(float64(r.N), 1/float64(r.Dim)) // cells per side (approximate for basis lattices)
	dof := float64(len(use) - 2)
	screenedWins := fr.Chi1-fr.Chi2 >= 9 && fr.Xi < L/4
	fr.Long = fr.Chi0-fr.Chi1 >= 25 && fr.Chi1/dof <= 3 && !screenedWins
	switch {
	case fr.Long:
		fr.Verdict = "LONG (Coulomb on this lattice)"
	case screenedWins:
		fr.Verdict = fmt.Sprintf("short (screened, xi=%.2f)", fr.Xi)
	case fr.Chi0-fr.Chi1 < 25:
		fr.Verdict = "short (no significant G term)"
	default:
		fr.Verdict = fmt.Sprintf("short (F1 chi2/dof=%.1f > 3)", fr.Chi1/dof)
	}
	return fr
}

// E4: g at the two nearest allowed shells, per type
type Sign struct {
	Same, Opp [][3]float64 // r, g, se
	Verdict   string
}

func gtype(r *Run, typ int) [][3]float64 {
	ns := len(r.Bins[0].Cnt)
	var bins []*Bin
	for _, b := range r.Bins {
		if b.Ideal > 0 && (b.Opp == (typ == 1)) {
			bins = append(bins, b)
		}
	}
	idealTot := 0.0
	tot := make([]float64, ns)
	for _, b := range bins {
		idealTot += b.Ideal
		for s := 0; s < ns; s++ {
			tot[s] += b.Cnt[s]
		}
	}
	if tot[0] == 0 {
		return nil
	}
	sort.Slice(bins, func(i, j int) bool { return bins[i].R < bins[j].R })
	var out [][3]float64
	for _, b := range bins {
		gs := make([]float64, ns)
		for s := 0; s < ns; s++ {
			gs[s] = (b.Cnt[s] / tot[s]) / (b.Ideal / idealTot)
		}
		m, se := meanSE(gs)
		out = append(out, [3]float64{b.R, m, se})
		if len(out) == 2 {
			break
		}
	}
	return out
}

func signE4(r *Run) Sign {
	var s Sign
	s.Same = gtype(r, 0)
	if r.Bip == 1 {
		s.Opp = gtype(r, 1)
	}
	sig := func(g [3]float64, above bool) bool {
		if g[2] == 0 {
			return false
		}
		z := (g[1] - 1) / g[2]
		if above {
			return z >= 5
		}
		return z <= -5
	}
	all := func(gs [][3]float64, above bool) bool {
		if len(gs) == 0 {
			return false
		}
		for _, g := range gs {
			if !sig(g, above) {
				return false
			}
		}
		return true
	}
	switch {
	case r.Bip == 1 && all(s.Opp, true) && all(s.Same, false):
		s.Verdict = "signed"
	case all(s.Same, true) && (r.Bip == 0 || all(s.Opp, true)):
		s.Verdict = "universal attraction"
	case all(s.Same, false) && (r.Bip == 0 || all(s.Opp, false)):
		s.Verdict = "universal repulsion"
	default:
		s.Verdict = "none"
	}
	return s
}

func get(dim, r string, D int) *Run {
	return parse(fmt.Sprintf("out/c1_%s_%s_E%d.out", dim, r, D))
}

func main() {
	fmt.Printf("## c1 E2: pair potential fits (r >= 2)\n\n| rules | d | lattice | bins | chi2 F0 / F1 / F2 | beta (se) | best m | verdict |\n|---|---|---|---|---|---|---|---|\n")
	long := map[string]map[string]bool{}
	for _, r := range rules {
		long[r] = map[string]bool{}
		for _, d := range []string{"2d", "3d"} {
			run := get(d, r, 2)
			if run == nil {
				fmt.Printf("| %s | %s | (missing) | | | | | |\n", r, d)
				continue
			}
			fr := fitE2(run)
			long[r][d] = fr.Long
			fmt.Printf("| %s | %s | %s | %d | %.1f / %.1f / %.1f | %.4g (%.2g) | %g | %s |\n", r, d, run.Lat, fr.NBins, fr.Chi0, fr.Chi1, fr.Chi2, fr.Beta, fr.BetaSE, fr.MBest, fr.Verdict)
		}
	}
	fmt.Printf("\n## c1 E2: V(r) at selected shells (3-D)\n\n")
	for _, r := range rules {
		run := get("3d", r, 2)
		if run == nil {
			continue
		}
		vb := potential(run, -1)
		fmt.Printf("- %s (%s):", r, run.Lat)
		for i, p := range vb {
			if i < 6 || i%8 == 0 {
				fmt.Printf(" r=%.2f V=%.4f±%.4f;", p.R, p.V, p.SE)
			}
		}
		fmt.Println()
	}
	fmt.Printf("\n## c1 E4: sign at the two nearest allowed shells\n\n| rules | d | g_same (r: g ± se) | g_opp (r: g ± se) | verdict |\n|---|---|---|---|---|\n")
	sign := map[string]map[string]string{}
	for _, r := range rules {
		sign[r] = map[string]string{}
		for _, d := range []string{"2d", "3d"} {
			run := get(d, r, 4)
			if run == nil {
				fmt.Printf("| %s | %s | (missing) | | |\n", r, d)
				continue
			}
			s := signE4(run)
			sign[r][d] = s.Verdict
			f := func(gs [][3]float64) string {
				var b strings.Builder
				for _, g := range gs {
					fmt.Fprintf(&b, "%.2f: %.4f ± %.4f; ", g[0], g[1], g[2])
				}
				return b.String()
			}
			fmt.Printf("| %s | %s | %s | %s | %s |\n", r, d, f(s.Same), f(s.Opp), s.Verdict)
		}
	}
	fmt.Printf("\n## c1 1-D control (ring, E2)\n\n")
	for _, r := range []string{"c", "ac", "bc", "abc"} {
		run := get("1d", r, 2)
		if run == nil {
			continue
		}
		vb := potential(run, -1)
		ns := len(run.Bins[0].Cnt)
		tot, near, mean := 0.0, 0.0, 0.0
		for _, b := range run.Bins {
			for s := 0; s < ns; s++ {
				tot += b.Cnt[s]
				mean += b.R * b.Cnt[s]
				if b.R <= 1+1e-9 {
					near += b.Cnt[s]
				}
			}
		}
		fmt.Printf("- %s: mean separation %.2f (ring N=128; independent ~32), P(r<=1)=%.4f, bins visited %d\n", r, mean/tot, near/tot, len(vb))
	}
	fmt.Printf("\n## c1 anemia (b off): unbonded quanta per cell vs distance to nearest void (3-D, E2)\n\n")
	for _, r := range []string{"none", "a", "c", "ac"} {
		run := get("3d", r, 2)
		if run == nil {
			continue
		}
		var rs []float64
		for k := range run.Anemia {
			rs = append(rs, k)
		}
		sort.Float64s(rs)
		fmt.Printf("- %s (%s):", r, run.Lat)
		for i, k := range rs {
			if i < 8 {
				m, se := meanSE(run.Anemia[k])
				fmt.Printf(" r=%.2f %.4f±%.4f;", k, m, se)
			}
		}
		fmt.Println()
	}
	fmt.Printf("\n## checks\n\n")
	for _, r := range rules {
		for _, d := range []string{"1d", "2d", "3d"} {
			for _, D := range []int{2, 4} {
				if run := get(d, r, D); run != nil {
					fmt.Printf("- %s %s E%d: %s; %s\n", r, d, D, run.Ledger, run.Acc)
				}
			}
		}
	}
	fmt.Printf("\n## Ending (PREREG §5)\n\n")
	ending := "S"
	anyLong := false
	for _, r := range rules {
		if long[r]["2d"] && long[r]["3d"] {
			anyLong = true
			fmt.Printf("- %s: long range in 2-D and 3-D; E4 sign 3-D: %s, 2-D: %s\n", r, sign[r]["3d"], sign[r]["2d"])
			if sign[r]["3d"] == "universal attraction" {
				ending = "G"
			}
		}
	}
	if ending != "G" && anyLong {
		ending = "E"
		for _, r := range rules {
			if long[r]["2d"] && long[r]["3d"] && sign[r]["3d"] != "signed" {
				ending = "E? (a long-range set is not signed: read the table)"
			}
		}
	}
	fmt.Printf("\n**Ending: %s**\n", ending)
}
