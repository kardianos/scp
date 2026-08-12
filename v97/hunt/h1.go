package hunt

import (
	"fmt"
	"math"
	"strings"
)

// H1 window and clause floors. Pre-registered against the H0 streams
// so Face C / HOLD / rings / leaks cannot sneak through as dipoles.
const (
	h1MinSpan   = 100.0 // t.u.
	h1Stretch   = 1000.0
	h1MinFrames = 4

	// Φ live vs inventory. Face C was ⟨Φ⟩/Em ~ 0.0025; HOLD ~ 7e-5.
	h1PhiOverEm = 0.015
	h1PhiAbs    = 0.05
	// Spatial split. Rings were 0; Face C ~1.1 with tiny Φ.
	h1SepAbs = 0.8
	// Polarization: intake and outtake on opposite halves.
	h1Polar = 0.25
	// Mass parked: |Em_end/Em_start − 1| and |Σ dEm|/Em.
	h1MassDrift = 0.20
	// J dipole is a dipole, not a loop, and does not slosh.
	h1JdipOverCirc = 0.25
	h1DirStable    = 0.80
	// Es cycle: throughput above the leak floor, and not tracking dM.
	h1EsThruOverDM = 2.0
	h1EsCorrMax    = 0.40
)

type H1Window struct {
	I0, I1  int
	T0, T1  float64
	Span    float64
	NPass   int
	Clauses []Check
	OK      bool
}

type H1Report struct {
	Path      string
	N         int
	T0, T1    float64
	Span      float64
	Best      H1Window
	Inhabited bool
	Stretch   bool
	Notes     []string
}

func (r H1Report) String() string {
	var b strings.Builder
	ver := "EMPTY"
	if r.Inhabited {
		ver = "INHABITED"
	}
	fmt.Fprintf(&b, "H1 %s  frames=%d  span=%.1f  %s\n", ver, r.N, r.Span, r.Path)
	if r.N == 0 {
		fmt.Fprintf(&b, "  no frames\n")
		return b.String()
	}
	fmt.Fprintf(&b, "  stream t=[%.1f, %.1f]  best window t=[%.1f, %.1f] Δt=%.1f  clauses %d/5\n",
		r.T0, r.T1, r.Best.T0, r.Best.T1, r.Best.Span, r.Best.NPass)
	for _, c := range r.Best.Clauses {
		mark := "FAIL"
		if c.Skip {
			mark = "SKIP"
		} else if c.OK {
			mark = "PASS"
		}
		fmt.Fprintf(&b, "  [%s] %s — %s\n", mark, c.Name, c.Detail)
	}
	if r.Stretch {
		fmt.Fprintf(&b, "  stretch: a passing window reached ≥ %.0f t.u.\n", h1Stretch)
	}
	for _, n := range r.Notes {
		fmt.Fprintf(&b, "  note: %s\n", n)
	}
	return b.String()
}

func GradeH1(path string, rows []Row) H1Report {
	rep := H1Report{Path: path, N: len(rows)}
	if len(rows) == 0 {
		rep.Notes = []string{"no frames"}
		return rep
	}
	rep.T0, rep.T1 = rows[0].T, rows[len(rows)-1].T
	rep.Span = rep.T1 - rep.T0
	if rep.Span < h1MinSpan {
		// still score the whole stream so the clause table is visible
		w := scoreWindow(rows, 0, len(rows)-1)
		rep.Best = w
		rep.Notes = append(rep.Notes, fmt.Sprintf("stream shorter than %.0f t.u. — cannot inhabit", h1MinSpan))
		return rep
	}
	best := H1Window{NPass: -1}
	var bestPass H1Window
	havePass := false
	for i := 0; i < len(rows); i++ {
		for j := len(rows) - 1; j >= i+h1MinFrames-1; j-- {
			span := rows[j].T - rows[i].T
			if span < h1MinSpan {
				break
			}
			w := scoreWindow(rows, i, j)
			if !havePass || betterWindow(w, best) {
				best = w
			}
			if w.OK && (!havePass || w.Span > bestPass.Span) {
				bestPass = w
				havePass = true
			}
		}
	}
	if havePass {
		rep.Best = bestPass
		rep.Inhabited = true
		rep.Stretch = bestPass.Span >= h1Stretch
	} else {
		rep.Best = best
	}
	return rep
}

func betterWindow(a, b H1Window) bool {
	if a.NPass != b.NPass {
		return a.NPass > b.NPass
	}
	if a.OK != b.OK {
		return a.OK
	}
	return a.Span > b.Span
}

func scoreWindow(rows []Row, i0, i1 int) H1Window {
	w := H1Window{I0: i0, I1: i1, T0: rows[i0].T, T1: rows[i1].T}
	w.Span = w.T1 - w.T0
	sl := rows[i0 : i1+1]
	w.Clauses = []Check{
		clausePhi(sl),
		clauseSep(sl),
		clauseMass(sl),
		clauseJDip(sl),
		clauseEsCycle(sl),
	}
	w.NPass = 0
	w.OK = true
	for _, c := range w.Clauses {
		if c.OK {
			w.NPass++
		} else if !c.Skip {
			w.OK = false
		}
	}
	if w.NPass < 5 {
		w.OK = false
	}
	return w
}

func liveRows(rows []Row) []Row {
	out := make([]Row, 0, len(rows))
	for _, r := range rows {
		if r.N > 0 {
			out = append(out, r)
		}
	}
	return out
}

func clausePhi(rows []Row) Check {
	rs := liveRows(rows)
	if len(rs) == 0 {
		return Check{Name: "Φ > 0", OK: false, Detail: "no live motif"}
	}
	var phi, em, ii, oo []float64
	for _, r := range rs {
		phi = append(phi, r.Phi)
		em = append(em, r.Em)
		ii = append(ii, r.I)
		oo = append(oo, r.O)
	}
	emBar := mean(em)
	if emBar <= 0 {
		emBar = 1
	}
	p, iBar, oBar := mean(phi), mean(ii), mean(oo)
	ok := p >= h1PhiAbs && p/emBar >= h1PhiOverEm && iBar >= h1PhiAbs/2 && oBar >= h1PhiAbs/2
	return Check{Name: "Φ > 0 (live throughput)", OK: ok,
		Detail: fmt.Sprintf("⟨Φ⟩=%.4g  ⟨I⟩=%.4g  ⟨O⟩=%.4g  Φ/Em=%.4g", p, iBar, oBar, p/emBar)}
}

func clauseSep(rows []Row) Check {
	rs := liveRows(rows)
	if len(rs) == 0 {
		return Check{Name: "I/O spatially separated", OK: false, Detail: "no live motif"}
	}
	var sep, pol []float64
	for _, r := range rs {
		sep = append(sep, r.Sep)
		// opposite-half traffic: intake fore & outtake aft, or the reverse
		a := r.Ifore + r.Oaft
		b := r.Iaft + r.Ofore
		den := a + b
		if den > 1e-12 {
			pol = append(pol, math.Abs(a-b)/den)
		}
	}
	sepBar := meanAbs(sep)
	polBar := mean(pol)
	ok := sepBar >= h1SepAbs && polBar >= h1Polar
	return Check{Name: "I/O spatially separated", OK: ok,
		Detail: fmt.Sprintf("⟨|sep|⟩=%.3f  ⟨polar⟩=%.3f", sepBar, polBar)}
}

func clauseMass(rows []Row) Check {
	rs := liveRows(rows)
	if len(rs) < 2 {
		return Check{Name: "dM/dt ≈ 0", OK: false, Detail: "too few live frames"}
	}
	e0, e1 := rs[0].Em, rs[len(rs)-1].Em
	emBar := mean(func() []float64 {
		xs := make([]float64, len(rs))
		for i, r := range rs {
			xs[i] = r.Em
		}
		return xs
	}())
	if emBar <= 0 {
		emBar = 1
	}
	// prefer box inventory when the motif is a giant foam (n huge)
	dBox := math.Abs(rs[len(rs)-1].EmTot - rs[0].EmTot)
	box0 := rs[0].EmTot
	if box0 <= 0 {
		box0 = 1
	}
	relMotif := math.Abs(e1-e0) / emBar
	relBox := dBox / box0
	// parked if either the motif or the box holds (giant motifs ≈ box)
	ok := relMotif <= h1MassDrift
	return Check{Name: "dM/dt ≈ 0 (parked)", OK: ok,
		Detail: fmt.Sprintf("Em %.4g→%.4g  |Δ|/⟨Em⟩=%.3f  |ΔEmTot|/EmTot0=%.3f", e0, e1, relMotif, relBox)}
}

func clauseJDip(rows []Row) Check {
	rs := liveRows(rows)
	if len(rs) < 3 {
		return Check{Name: "dipole of J steady", OK: false, Detail: "too few live frames"}
	}
	var jdip, circ, jdx, jdy []float64
	for _, r := range rs {
		jdip = append(jdip, r.Jdip)
		circ = append(circ, r.Circ)
		jdx = append(jdx, r.JdipX)
		jdy = append(jdy, r.JdipY)
	}
	jd, ci := meanAbs(jdip), meanAbs(circ)
	dir := math.Max(fracSameSign(jdx), fracSameSign(jdy))
	// not a loop: dipole moment competes with circulation
	notLoop := jd >= h1JdipOverCirc*math.Max(ci, 1e-9)
	ok := jd > 1e-3 && notLoop && dir >= h1DirStable
	return Check{Name: "dipole of J steady", OK: ok,
		Detail: fmt.Sprintf("⟨|Jdip|⟩=%.4g  ⟨|Circ|⟩=%.4g  dir-stable=%.2f", jd, ci, dir)}
}

func clauseEsCycle(rows []Row) Check {
	rs := liveRows(rows)
	if len(rs) < 3 {
		return Check{Name: "Es throughput at constant M", OK: false, Detail: "too few live frames"}
	}
	var thru, net, dM []float64
	for _, r := range rs {
		thru = append(thru, r.CoreEsThru)
		net = append(net, r.CoreEsNet)
		dM = append(dM, r.DEmTot)
	}
	th, dm := meanAbs(thru), meanAbs(dM)
	c := pearson(net, dM)
	cAbs := pearson(absSlice(thru), absSlice(dM))
	tracks := (!math.IsNaN(c) && math.Abs(c) >= h1EsCorrMax) ||
		(!math.IsNaN(cAbs) && cAbs >= h1EsCorrMax)
	aboveLeak := th >= h1EsThruOverDM*math.Max(dm, 1e-6)
	ok := th > 0.05 && aboveLeak && !tracks
	return Check{Name: "Es throughput at constant M", OK: ok,
		Detail: fmt.Sprintf("⟨coreThru⟩=%.4g  ⟨|dEmTot|⟩=%.4g  corr(net,dM)=%.3f  |Thru|~|dM|=%.3f  tracks=%v",
			th, dm, c, cAbs, tracks)}
}
