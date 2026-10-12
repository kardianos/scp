// tab.go — v103 gate g1: ledger, bars and the registered ending
// (../PREREG.md §7) from the geo3d outputs in out/.
//
//	go run tab.go > out/tab.out
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Sample map[string]float64

type ProfAcc struct {
	k                   int // samples accumulated
	n, sqd, mis, dm, gf map[int]float64
}

func newProf() *ProfAcc {
	return &ProfAcc{n: map[int]float64{}, sqd: map[int]float64{}, mis: map[int]float64{}, dm: map[int]float64{}, gf: map[int]float64{}}
}

type Sector struct {
	Samples []Sample
	Prof    *ProfAcc
	Check   int // -1 absent, 0 fail, 1 ok
	Made    string
	GreenOK bool
}

type Run struct {
	Sec   map[string]*Sector
	Rigid Sample
	Flat  float64
	Done  bool
}

func kv(fields []string) Sample {
	s := Sample{}
	for _, f := range fields {
		p := strings.SplitN(f, "=", 2)
		if len(p) != 2 {
			continue
		}
		v, err := strconv.ParseFloat(strings.TrimPrefix(p[1], "inf:"), 64)
		if err != nil {
			continue
		}
		s[p[0]] = v
	}
	return s
}

func secName(fields []string) string {
	for _, f := range fields {
		if strings.HasPrefix(f, "sec=") {
			return strings.TrimPrefix(f, "sec=")
		}
	}
	return ""
}

func parse(path string) *Run {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	r := &Run{Sec: map[string]*Sector{}}
	get := func(s string) *Sector {
		if r.Sec[s] == nil {
			r.Sec[s] = &Sector{Prof: newProf(), Check: -1, GreenOK: true}
		}
		return r.Sec[s]
	}
	cur := "" // sector being made (GAP/EXCESS lines precede their sector)
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	for sc.Scan() {
		l := sc.Text()
		fs := strings.Fields(l)
		if len(fs) == 0 {
			continue
		}
		switch fs[0] {
		case "FLAT":
			r.Flat = kv(fs[1:])["E"]
		case "RIGID":
			r.Rigid = kv(fs[1:])
		case "GAP":
			cur = "G"
			get(cur).Made = l
		case "EXCESS":
			cur = "X"
			get(cur).Made = l
		case "SAMPLE":
			s := get(secName(fs))
			s.Samples = append(s.Samples, kv(fs[1:]))
		case "CENTRE":
			if kv(fs[1:])["greenit"] < 0 {
				get(secName(fs)).GreenOK = false
			}
		case "PROF":
			s := get(secName(fs))
			v := kv(fs[1:])
			rr := int(v["r"])
			s.Prof.n[rr] += v["n"]
			s.Prof.sqd[rr] += v["sqd"]
			s.Prof.mis[rr] += v["mis"]
			s.Prof.dm[rr] += v["dm"] * v["n"]
			s.Prof.gf[rr] += v["gf"] * v["n"]
		case "CHECK":
			s := get(secName(fs))
			s.Check = int(kv(fs[1:])["ok"])
		case "DONE":
			r.Done = true
		}
	}
	_ = cur
	for _, s := range r.Sec {
		s.Prof.k = len(s.Samples)
	}
	return r
}

func mean(ss []Sample, key string, from, to int) float64 {
	t, c := 0.0, 0
	for i, s := range ss {
		if i < from || i >= to {
			continue
		}
		if v, ok := s[key]; ok && !math.IsNaN(v) {
			t += v
			c++
		}
	}
	if c == 0 {
		return math.NaN()
	}
	return t / float64(c)
}

func all(ss []Sample, f func(Sample) bool) bool {
	if len(ss) == 0 {
		return false
	}
	for _, s := range ss {
		if !f(s) {
			return false
		}
	}
	return true
}

var Ls = []int{6, 8, 10, 12, 16}

func name(a, b, c int) string {
	s := ""
	if a == 1 {
		s += "a"
	}
	if b == 1 {
		s += "b"
	}
	if c == 1 {
		s += "c"
	}
	if s == "" {
		s = "none"
	}
	return s
}

type key struct{ a, b, c int }

func f4(v float64) string {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return "—"
	}
	return fmt.Sprintf("%.4f", v)
}

func yn(v bool) string {
	if v {
		return "✔"
	}
	return "✘"
}

// two-size readings, 12 -> 16 (PREREG §7)
func dims(r12, r16 *Run, sec string) (dH, ds float64) {
	if r12 == nil || r16 == nil || r12.Sec[sec] == nil || r16.Sec[sec] == nil {
		return math.NaN(), math.NaN()
	}
	s12, s16 := r12.Sec[sec].Samples, r16.Sec[sec].Samples
	lr := math.Log(4096.0 / 1728.0)
	dH = lr / math.Log(mean(s16, "meand", 0, 99)/mean(s12, "meand", 0, 99))
	ds = -2 * lr / math.Log(mean(s16, "lam2", 0, 99)/mean(s12, "lam2", 0, 99))
	return
}

// the medium reading at one size: one piece, stationary, Lanczos converged
func steady(r *Run, sec string) (ok bool, dT, dO float64) {
	if r == nil || r.Sec[sec] == nil || len(r.Sec[sec].Samples) != 20 {
		return false, math.NaN(), math.NaN()
	}
	ss := r.Sec[sec].Samples
	dT = mean(ss, "psiT", 10, 20) - mean(ss, "psiT", 0, 10)
	dO = mean(ss, "psiO", 10, 20) - mean(ss, "psiO", 0, 10)
	one := all(ss, func(s Sample) bool { return s["pieces"] == 1 && s["lamok"] == 1 })
	return one && math.Abs(dT) <= 0.01 && math.Abs(dO) <= 0.01, dT, dO
}

func in(v, lo, hi float64) bool { return !math.IsNaN(v) && v >= lo && v <= hi }

// a readable, persistent defect at one size (PREREG §7, P1 and P2)
func readable(r *Run, sec string) bool {
	if r == nil || r.Sec["M"] == nil || r.Sec[sec] == nil || len(r.Sec[sec].Samples) != 20 {
		return false
	}
	if !(mean(r.Sec["M"].Samples, "nmis", 0, 99) < 1) {
		return false
	}
	return all(r.Sec[sec].Samples, func(s Sample) bool { return s["ncl"] == 1 && s["dspread"] <= 6 })
}

type Foot struct {
	rdef        int
	dmeand, dsq float64
	dm, dgf     map[int]float64
	slope       float64
	sign        string
	rmax        int
}

func footprint(r *Run, sec string, L int) *Foot {
	if r == nil || r.Sec[sec] == nil || r.Sec["M"] == nil || r.Sec[sec].Prof.k == 0 {
		return nil
	}
	s, m := r.Sec[sec], r.Sec["M"]
	k := float64(s.Prof.k)
	ft := &Foot{rdef: -1, dm: map[int]float64{}, dgf: map[int]float64{}}
	mref := mean(m.Samples, "meand", 0, 99)
	ft.dmeand = mean(s.Samples, "meand", 0, 99) - mref
	ft.dsq = mean(s.Samples, "sqdef", 0, 99) - mean(m.Samples, "sqdef", 0, 99)
	rs := []int{}
	for rr := range s.Prof.n {
		rs = append(rs, rr)
	}
	sort.Ints(rs)
	for _, rr := range rs {
		if s.Prof.mis[rr]/k >= 0.5 && rr > ft.rdef {
			ft.rdef = rr
		}
		if n := s.Prof.n[rr]; n > 0 {
			ft.dm[rr] = s.Prof.dm[rr]/n - mref
			if mn := m.Prof.n[rr]; mn > 0 {
				ft.dgf[rr] = s.Prof.gf[rr]/n - m.Prof.gf[rr]/mn
			} else {
				ft.dgf[rr] = math.NaN()
			}
		}
		if rr > ft.rmax {
			ft.rmax = rr
		}
	}
	// distance footprint: slope of ln|dm| against ln r over 2 <= r <= L/2, if one sign
	xs, ys := []float64{}, []float64{}
	pos, neg := 0, 0
	for rr := 2; rr <= L/2; rr++ {
		v, ok := ft.dm[rr]
		if !ok || v == 0 || math.IsNaN(v) {
			continue
		}
		if v > 0 {
			pos++
		} else {
			neg++
		}
		xs = append(xs, math.Log(float64(rr)))
		ys = append(ys, math.Log(math.Abs(v)))
	}
	ft.slope = math.NaN()
	switch {
	case pos > 0 && neg > 0:
		ft.sign = "mixed"
	case pos > 0:
		ft.sign = "+ (more where)"
	case neg > 0:
		ft.sign = "− (less where)"
	default:
		ft.sign = "0"
	}
	if len(xs) >= 2 && !(pos > 0 && neg > 0) {
		mx, my := 0.0, 0.0
		for i := range xs {
			mx += xs[i]
			my += ys[i]
		}
		mx /= float64(len(xs))
		my /= float64(len(xs))
		sxy, sxx := 0.0, 0.0
		for i := range xs {
			sxy += (xs[i] - mx) * (ys[i] - my)
			sxx += (xs[i] - mx) * (xs[i] - mx)
		}
		if sxx > 0 {
			ft.slope = -sxy / sxx // dm ~ r^(-p)
		}
	}
	return ft
}

func main() {
	runs := map[key]map[int]*Run{}
	order := []key{}
	for m := 0; m < 8; m++ {
		k := key{m >> 2 & 1, m >> 1 & 1, m & 1}
		order = append(order, k)
		runs[k] = map[int]*Run{}
		for _, L := range Ls {
			runs[k][L] = parse(fmt.Sprintf("out/g1_a%db%dc%d_L%d.out", k.a, k.b, k.c, L))
		}
	}
	missing := 0
	for _, k := range order {
		for _, L := range Ls {
			if r := runs[k][L]; r == nil || !r.Done {
				missing++
			}
		}
	}

	fmt.Printf("## g1 the medium from the flat start (M), L = 12 | L = 16\n\n")
	fmt.Printf("| set | cheapest flip from flat (t) | ψ_T | ψ_O | ⟨d⟩ | λ2 | d_H | d_s | steady (ΔψT, ΔψO) | holds |\n|---|---|---|---|---|---|---|---|---|---|\n")
	holds := map[key]bool{}
	asm := map[key]bool{}
	for _, k := range order {
		r12, r16 := runs[k][12], runs[k][16]
		cell := func(f func(*Run) string) string {
			out := []string{}
			for _, r := range []*Run{r12, r16} {
				if r == nil {
					out = append(out, "(missing)")
					continue
				}
				out = append(out, f(r))
			}
			return strings.Join(out, " \\| ")
		}
		dH, ds := dims(r12, r16, "M")
		ok12, dT12, dO12 := steady(r12, "M")
		ok16, dT16, dO16 := steady(r16, "M")
		h := in(dH, 2.7, 3.3) && in(ds, 2.7, 3.3) && ok12 && ok16
		holds[k] = h
		fmt.Printf("| %s | %s | %s | %s | %s | %s | %s | %s | %+.3g, %+.3g \\| %+.3g, %+.3g | %s |\n", name(k.a, k.b, k.c),
			cell(func(r *Run) string { return fmt.Sprintf("%.4g", r.Rigid["min_dE"]) }),
			cell(func(r *Run) string { return f4(mean(r.Sec["M"].Samples, "psiT", 0, 99)) }),
			cell(func(r *Run) string { return f4(mean(r.Sec["M"].Samples, "psiO", 0, 99)) }),
			cell(func(r *Run) string { return f4(mean(r.Sec["M"].Samples, "meand", 0, 99)) }),
			cell(func(r *Run) string { return f4(mean(r.Sec["M"].Samples, "lam2", 0, 99)) }),
			f4(dH), f4(ds), dT12, dO12, dT16, dO16, yn(h))
	}

	fmt.Printf("\n## assembly (A: melted, then cooled to t), L = 12 | L = 16 — reported\n\n")
	fmt.Printf("| set | the melt: d_H, d_s | ψ_T | ψ_O | triangles per quantum | d_H | d_s | assembles |\n|---|---|---|---|---|---|---|---|\n")
	for _, k := range order {
		r12, r16 := runs[k][12], runs[k][16]
		cell := func(key string) string {
			out := []string{}
			for _, r := range []*Run{r12, r16} {
				if r == nil || r.Sec["A"] == nil {
					out = append(out, "(missing)")
					continue
				}
				out = append(out, f4(mean(r.Sec["A"].Samples, key, 0, 99)))
			}
			return strings.Join(out, " \\| ")
		}
		mH, mS := dims(r12, r16, "A0")
		dH, ds := dims(r12, r16, "A")
		ok12, _, _ := steady(r12, "A")
		ok16, _, _ := steady(r16, "A")
		a := in(dH, 2.7, 3.3) && in(ds, 2.7, 3.3) && ok12 && ok16
		asm[k] = a
		fmt.Printf("| %s | %s, %s | %s | %s | %s | %s | %s | %s |\n", name(k.a, k.b, k.c), f4(mH), f4(mS), cell("psiT"), cell("psiO"), cell("tri"), f4(dH), f4(ds), yn(a))
	}

	fmt.Printf("\n## size-ladder fits, M (all five sizes; reported)\n\n| set | d_H (ln N vs ln ⟨d⟩) | d_s (ln N vs ln λ2) |\n|---|---|---|\n")
	for _, k := range order {
		xs, yd, yl := []float64{}, []float64{}, []float64{}
		for _, L := range Ls {
			r := runs[k][L]
			if r == nil || r.Sec["M"] == nil {
				continue
			}
			xs = append(xs, math.Log(float64(L*L*L)))
			yd = append(yd, math.Log(mean(r.Sec["M"].Samples, "meand", 0, 99)))
			yl = append(yl, math.Log(mean(r.Sec["M"].Samples, "lam2", 0, 99)))
		}
		slope := func(ys []float64) float64 {
			if len(xs) < 2 {
				return math.NaN()
			}
			mx, my := 0.0, 0.0
			for i := range xs {
				mx += xs[i]
				my += ys[i]
			}
			mx /= float64(len(xs))
			my /= float64(len(xs))
			sxy, sxx := 0.0, 0.0
			for i := range xs {
				sxy += (xs[i] - mx) * (ys[i] - my)
				sxx += (xs[i] - mx) * (xs[i] - mx)
			}
			return sxy / sxx
		}
		fmt.Printf("| %s | %s | %s |\n", name(k.a, k.b, k.c), f4(1/slope(yd)), f4(-2/slope(yl)))
	}

	fmt.Printf("\n## one gap (G) and one excess (X), L = 12 | L = 16\n\n")
	fmt.Printf("| set | sector | made | clusters (min–max) | spread (max) | readable and persistent | δ⟨d⟩ | δ(square deficit) | defect radius r_def | distance footprint: sign, p in Δm ∝ r^−p |\n|---|---|---|---|---|---|---|---|---|---|\n")
	P := map[key]bool{}
	PX := map[key]bool{}
	cls := map[key]string{}
	for _, k := range order {
		for _, sec := range []string{"G", "X"} {
			cells := make([][]string, 8)
			for _, L := range []int{12, 16} {
				r := runs[k][L]
				if r == nil || r.Sec[sec] == nil {
					for i := range cells {
						cells[i] = append(cells[i], "(missing)")
					}
					continue
				}
				s := r.Sec[sec]
				mn, mxv, sp := 1e9, -1e9, -1e9
				for _, x := range s.Samples {
					mn = math.Min(mn, x["ncl"])
					mxv = math.Max(mxv, x["ncl"])
					sp = math.Max(sp, x["dspread"])
				}
				made := "—"
				if f := strings.Fields(s.Made); len(f) >= 3 {
					made = strings.Join(f[1:3], " ")
				}
				cells[0] = append(cells[0], made)
				cells[1] = append(cells[1], fmt.Sprintf("%.0f–%.0f", mn, mxv))
				cells[2] = append(cells[2], fmt.Sprintf("%.0f", sp))
				cells[3] = append(cells[3], yn(readable(r, sec)))
				ft := footprint(r, sec, L)
				if ft == nil {
					for i := 4; i < 8; i++ {
						cells[i] = append(cells[i], "—")
					}
					continue
				}
				cells[4] = append(cells[4], fmt.Sprintf("%+.5f", ft.dmeand))
				cells[5] = append(cells[5], fmt.Sprintf("%+.2f", ft.dsq))
				cells[6] = append(cells[6], fmt.Sprintf("%d", ft.rdef))
				cells[7] = append(cells[7], fmt.Sprintf("%s, %s", ft.sign, f4(ft.slope)))
			}
			pr := readable(runs[k][12], sec) && readable(runs[k][16], sec)
			if sec == "G" {
				P[k] = pr
				f12, f16 := footprint(runs[k][12], "G", 12), footprint(runs[k][16], "G", 16)
				switch {
				case f12 == nil || f16 == nil:
					cls[k] = "—"
				case f12.rdef <= 3 && f16.rdef <= 3 && f16.rdef <= f12.rdef+1:
					cls[k] = "contact"
				default:
					cls[k] = "extended"
				}
			} else {
				PX[k] = pr
			}
			row := []string{}
			for _, c := range cells {
				row = append(row, strings.Join(c, " \\| "))
			}
			fmt.Printf("| %s | %s | %s |\n", name(k.a, k.b, k.c), sec, strings.Join(row, " | "))
		}
	}

	fmt.Printf("\n## distance and field footprints of the gap, L = 16 (reported): Δm(r) = m̄(r) − m_ref, Δgf(r) = G_gap(r) − G_ref(r)\n\n")
	for _, k := range order {
		ft := footprint(runs[k][16], "G", 16)
		if ft == nil {
			continue
		}
		rs := []int{}
		for rr := range ft.dm {
			rs = append(rs, rr)
		}
		sort.Ints(rs)
		parts := []string{}
		for _, rr := range rs {
			parts = append(parts, fmt.Sprintf("%d:%+.4f/%+.2e", rr, ft.dm[rr], ft.dgf[rr]))
		}
		fmt.Printf("- %s: %s\n", name(k.a, k.b, k.c), strings.Join(parts, " "))
	}

	// the ending (PREREG §7)
	fmt.Printf("\n## bars and the ending\n\n| set | holds | gap readable | phase (M, L = 16) | footprint | excess readable | assembles |\n|---|---|---|---|---|---|---|\n")
	gs, hs := []string{}, []string{}
	for _, k := range order {
		phase := "—"
		if r := runs[k][16]; r != nil && r.Sec["M"] != nil {
			pt, po := mean(r.Sec["M"].Samples, "psiT", 0, 99), mean(r.Sec["M"].Samples, "psiO", 0, 99)
			switch {
			case pt >= 0.99 && po >= 0.99:
				phase = "crystal"
			case po >= 0.99:
				phase = "nematic"
			default:
				phase = "disordered"
			}
		}
		n := name(k.a, k.b, k.c)
		fmt.Printf("| %s | %s | %s | %s | %s | %s | %s |\n", n, yn(holds[k]), yn(P[k]), phase, cls[k], yn(PX[k]), yn(asm[k]))
		if holds[k] && P[k] {
			gs = append(gs, n+" ("+phase+", "+cls[k]+")")
		} else if holds[k] {
			hs = append(hs, n)
		}
	}
	ending := "F (no set holds three dimensions)"
	switch {
	case len(gs) > 0:
		ending = "G — the gate passes: " + strings.Join(gs, "; ")
	case len(hs) > 0:
		ending = "H — holds without a readable gap: " + strings.Join(hs, ", ")
	}
	fmt.Printf("\n**Ending: %s**\n", ending)

	fmt.Printf("\n## checks\n\n")
	bad, flat, green := 0, 0, 0
	for _, k := range order {
		for _, L := range Ls {
			r := runs[k][L]
			if r == nil {
				continue
			}
			if r.Flat != 0 {
				flat++
			}
			for sn, s := range r.Sec {
				if s.Check != 1 && sn != "A0" {
					bad++
					fmt.Printf("- CHECK not ok: %s L=%d sector %s (%d)\n", name(k.a, k.b, k.c), L, sn, s.Check)
				}
				if !s.GreenOK {
					green++
				}
				for _, x := range s.Samples {
					if x["lamok"] != 1 && x["pieces"] == 1 {
						fmt.Printf("- Lanczos not converged: %s L=%d sector %s k=%.0f\n", name(k.a, k.b, k.c), L, sn, x["k"])
					}
				}
			}
		}
	}
	fmt.Printf("- missing or unfinished runs: %d\n- sectors failing CHECK: %d\n- flat start with nonzero energy: %d\n- Green's-function solves not converged: %d\n", missing, bad, flat, green)
}
