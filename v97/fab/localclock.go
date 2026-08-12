package fab

// localclock.go — P2 PROTOTYPE (queue #7): the FREECELL §2 local-clock
// scheduler carried to the v90 FREE-CELL SUBSTRATE, in Go first (README
// P2/P4: "the Go kernel prototypes the localclock batch schedule").
//
// WHAT THIS IS. v89/localclock.c proved the execution half on a
// synthetic degree-7 graph. This prototype re-proves the four measured
// conditions ON THE REAL SUBSTRATE — the actual 2D dart-throw bath, jam
// settle, contact-rule channels, and law-derived pitches (w2e detuned
// by a seeded blob load) — and adds the two things v90 needs from P2:
//   * goroutine BATCH execution, bit-identical to the serial async
//     schedule at any worker count (the batch IS the schedule, wide);
//   * the QUIET-REGION ECONOMY, measured: activity-derived dilation
//     (active cells step at dt, quiet bath cells at R*dt) makes the
//     event cost scale with the ACTIVE region, not the box — the
//     scaling table is the acceptance artifact (README P2).
//
// WHAT THIS IS NOT. The dynamics here is the REDUCED execution skeleton
// (phase advance + paired antisymmetric channel transfers), exactly as
// in v89/localclock.c — NOT the 12-pass law. The full-law async engine
// is the production P2 build (C, after this prototype); its open items
// are recorded in NEXT.md. No law constant is touched; exp=p2lc is a
// Go-only apparatus experiment.
//
// THE FOUR CONDITIONS (FREECELL §2, each a measured requirement):
//   1. order events by a total function of state (t, index) — never
//      arrival (the arrival-order arm is kept as the failing control);
//   2. bound skew in LOCAL TIME (the lookahead), never tick count;
//   3. channels own their transfers and their own clock;
//   4. never order by the tick counter — track the skew it would need
//      to hold and assert K < M/2 (M=256 here).

import (
	"math"
	"sync"
	"time"
)

const lcMaxDeg = 24

type lcFab struct {
	N      int
	nb     [][]int32 // neighbours (cells)
	ce     [][]int32 // incident edge ids
	th     []float64
	pub    []float64
	w      []float64 // pitch DETUNE from the bath (0 in vacuum — quiet)
	E      []float64
	t, h   []float64
	tau    []int64
	ne     int
	ei, ej []int32
	et, eh []float64
	kc, ge float64
}

func (f *lcFab) clone() *lcFab {
	g := &lcFab{N: f.N, nb: f.nb, ce: f.ce, ne: f.ne, ei: f.ei, ej: f.ej,
		kc: f.kc, ge: f.ge}
	g.th = append([]float64(nil), f.th...)
	g.pub = append([]float64(nil), f.pub...)
	g.w = f.w
	g.E = append([]float64(nil), f.E...)
	g.h = f.h
	g.eh = f.eh
	g.t = make([]float64, f.N)
	g.tau = make([]int64, f.N)
	g.et = make([]float64, f.ne)
	return g
}

// lcBuild constructs the reduced-model state ON the real substrate:
// 2D bath at box L, jam-settled, channels from the live contact rule,
// blob load at centre -> per-cell law pitch w2e -> detune w. Active =
// tagged blob cells (h=dt); quiet bath cells take R*dt steps.
func (s *Sim) lcBuild(L float64, dil float64, kc, ge, dt float64) (*lcFab, int) {
	P := &s.P
	savedL := P.L
	P.L = L
	defer func() { P.L = savedL }()

	V := L * L * L
	nmax := int(0.8*V/(P.Dmin*P.Dmin*P.Dmin)) + 4096
	bx := make([]float64, nmax)
	by := make([]float64, nmax)
	bz := make([]float64, nmax)
	nb := s.buildBath2d(bx, by, bz, nmax)
	s.allocAll(nb)
	for i := 0; i < nb; i++ {
		s.px[i], s.py[i], s.pz[i] = bx[i], by[i], bz[i]
		s.cr0[i] = P.R0 * (1.0 + P.Rjit*(2.0*s.frand()-1.0))
		s.cr[i] = s.cr0[i]
		s.n1x[i], s.n1y[i], s.n1z[i] = s.randUnit()
		s.n2x[i], s.n2y[i], s.n2z[i] = s.randUnit()
		_ = s.frand()
		s.th2[i] = s.frand() * TwoPi
		s.cbeta[i] = s.frand() * TwoPi
		s.Es[i] = P.Es0
	}
	s.jamSettle()
	for i := 0; i < s.NC; i++ {
		s.cr[i] = s.cr0[i] // Es = e_s0 everywhere at this stage
	}
	s.topoRefresh()

	// blob load at centre (the ACTIVE region): Em -> x -> w2e detune
	cx0, cy0, cz0 := 0.5*L, 0.5*L, 0.5*L
	nact := 0
	f := &lcFab{N: s.NC, kc: kc, ge: ge}
	f.th = make([]float64, f.N)
	f.pub = make([]float64, f.N)
	f.w = make([]float64, f.N)
	f.E = make([]float64, f.N)
	f.t = make([]float64, f.N)
	f.h = make([]float64, f.N)
	f.tau = make([]int64, f.N)
	f.nb = make([][]int32, f.N)
	f.ce = make([][]int32, f.N)
	for i := 0; i < s.NC; i++ {
		dx := s.wr(s.px[i] - cx0)
		dy := s.wr(s.py[i] - cy0)
		dz := s.wr(s.pz[i] - cz0)
		add := P.Amp * math.Exp(-(dx*dx+dy*dy+dz*dz)/(2.0*P.Sigma*P.Sigma))
		em := 0.0
		active := false
		if add >= 1e-4 {
			em = add
			if add > 0.05*P.Amp {
				active = true
				nact++
			}
		}
		x := em / P.Cap
		w2e := P.W2 / (1.0 + P.QDetune*x)
		f.th[i] = s.th2[i]
		f.pub[i] = f.th[i]
		f.w[i] = w2e - P.W2 // detune from the bath: EXACTLY 0 in vacuum
		f.E[i] = 1.0 + em
		if active {
			f.h[i] = dt
		} else {
			f.h[i] = dil * dt
		}
	}
	// channels from the live slot table
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree {
			continue
		}
		i, j := s.sli[sl], s.slj[sl]
		if len(f.nb[i]) >= lcMaxDeg || len(f.nb[j]) >= lcMaxDeg {
			continue
		}
		e := int32(f.ne)
		f.ei = append(f.ei, int32(i))
		f.ej = append(f.ej, int32(j))
		f.eh = append(f.eh, 0.5*(f.h[i]+f.h[j]))
		f.et = append(f.et, 0)
		f.nb[i] = append(f.nb[i], int32(j))
		f.nb[j] = append(f.nb[j], int32(i))
		f.ce[i] = append(f.ce[i], e)
		f.ce[j] = append(f.ce[j], e)
		f.ne++
	}
	return f, nact
}

func (f *lcFab) dtheta(i int) float64 {
	sm := 0.0
	for _, j := range f.nb[i] {
		sm += math.Sin(f.pub[j] - f.th[i])
	}
	return f.w[i] + f.kc*sm
}

func (f *lcFab) advance(i int) {
	hi := f.h[i]
	f.th[i] += hi * f.dtheta(i)
	f.pub[i] = f.th[i]
	f.t[i] += hi
	f.tau[i]++
}

func (f *lcFab) advanceEdge(e int) {
	i, j := f.ei[e], f.ej[e]
	q := f.eh[e] * f.ge * (f.E[i] - f.E[j])
	f.E[i] -= q
	f.E[j] += q
	f.et[e] += f.eh[e]
}

func (f *lcFab) totalE() float64 {
	sm := 0.0
	for i := 0; i < f.N; i++ {
		sm += f.E[i]
	}
	return sm
}

// runSync: the global-tick reference integrator at step dtv.
func (f *lcFab) runSync(T, dtv float64) int64 {
	N := f.N
	nth := make([]float64, N)
	dE := make([]float64, N)
	steps := int(T/dtv + 0.5)
	for st := 0; st < steps; st++ {
		for i := 0; i < N; i++ {
			nth[i] = f.th[i] + dtv*f.dtheta(i)
		}
		for i := range dE {
			dE[i] = 0
		}
		for e := 0; e < f.ne; e++ {
			i, j := f.ei[e], f.ej[e]
			q := dtv * f.ge * (f.E[i] - f.E[j])
			dE[i] -= q
			dE[j] += q
		}
		for i := 0; i < N; i++ {
			f.th[i] = nth[i]
			f.E[i] += dE[i]
			f.pub[i] = f.th[i]
		}
	}
	return int64(steps) * int64(N+f.ne)
}

// runLocal: the serial async schedule. sel 0 = canonical (t, index)
// order; 1 = arrival order (the measured FAILURE mode, kept as control).
// rot/dir vary the scan origin/direction: under sel=0 they must not
// change one bit. Returns events executed; maxSkew = largest neighbour
// tick difference seen (condition 4's quantity).
func (f *lcFab) runLocal(T, look float64, sel, rot, dir int) (ev int64, maxSkew int64) {
	N, NE := f.N, f.ne
	NT := N + NE
	for {
		pick := -1
		bt := 0.0
		for q := 0; q < NT; q++ {
			u := (rot + q) % NT
			if dir != 0 {
				u = (rot + NT - 1 - q) % NT
			}
			var tu, lim float64
			lim = 1e300
			if u < N {
				if f.t[u] >= T {
					continue
				}
				tu = f.t[u]
				for _, j := range f.nb[u] {
					if f.t[j] < lim {
						lim = f.t[j]
					}
				}
				if tu+f.h[u] > lim+look {
					continue
				}
			} else {
				e := u - N
				if f.et[e] >= T {
					continue
				}
				tu = f.et[e]
				a, b := f.t[f.ei[e]], f.t[f.ej[e]]
				lim = a
				if b < a {
					lim = b
				}
				if tu+f.eh[e] > lim+look {
					continue
				}
			}
			if sel == 1 {
				pick = u
				break
			}
			if pick < 0 || tu < bt || (tu == bt && u < pick) {
				pick = u
				bt = tu
			}
		}
		if pick < 0 {
			break
		}
		if pick < N {
			f.advance(pick)
			for _, j := range f.nb[pick] {
				dd := f.tau[pick] - f.tau[j]
				if dd < 0 {
					dd = -dd
				}
				if dd > maxSkew {
					maxSkew = dd
				}
			}
		} else {
			f.advanceEdge(pick - N)
		}
		ev++
	}
	return ev, maxSkew
}

func pfor(n, workers int, fn func(lo, hi int)) {
	if workers <= 1 {
		fn(0, n)
		return
	}
	var wg sync.WaitGroup
	chunk := (n + workers - 1) / workers
	for w := 0; w < workers; w++ {
		lo := w * chunk
		hi := lo + chunk
		if hi > n {
			hi = n
		}
		if lo >= hi {
			continue
		}
		wg.Add(1)
		go func(lo, hi int) { defer wg.Done(); fn(lo, hi) }(lo, hi)
	}
	wg.Wait()
}

// runBatch: the batch-parallel execution of the SAME schedule — an
// eligible event joins the batch iff it is the minimum (t, index) over
// its conflict neighbourhood (incident channels PLUS adjacent cells:
// a cell event reads its neighbours' published phase, so adjacent cell
// events race even though they write disjoint memory — the v89 bug).
// Bit-identical to runLocal(sel=0) at any worker count, by measurement.
func (f *lcFab) runBatch(T, look float64, workers int) (ev int64, rounds int64, width float64) {
	N, NE := f.N, f.ne
	NT := N + NE
	elig := make([]bool, NT)
	sel := make([]bool, NT)
	evt := make([]float64, NT)
	acc := 0.0
	for {
		pfor(NT, workers, func(lo, hi int) {
			for u := lo; u < hi; u++ {
				elig[u] = false
				var tu, lim float64
				lim = 1e300
				if u < N {
					if f.t[u] >= T {
						continue
					}
					tu = f.t[u]
					for _, j := range f.nb[u] {
						if f.t[j] < lim {
							lim = f.t[j]
						}
					}
					if tu+f.h[u] > lim+look {
						continue
					}
				} else {
					e := u - N
					if f.et[e] >= T {
						continue
					}
					tu = f.et[e]
					a, b := f.t[f.ei[e]], f.t[f.ej[e]]
					lim = a
					if b < a {
						lim = b
					}
					if tu+f.eh[e] > lim+look {
						continue
					}
				}
				elig[u] = true
				evt[u] = tu
			}
		})
		pfor(NT, workers, func(lo, hi int) {
			for u := lo; u < hi; u++ {
				sel[u] = false
				if !elig[u] {
					continue
				}
				win := true
				var cells [2]int32
				nc := 0
				if u < N {
					cells[0] = int32(u)
					nc = 1
				} else {
					cells[0] = f.ei[u-N]
					cells[1] = f.ej[u-N]
					nc = 2
				}
				// schedule equivalence needs min over PENDING events
				// (eligible or not): an earlier-but-blocked neighbour
				// must still hold this event back, exactly as the
				// serial (t, index) order would. Min-over-eligible
				// measured 6.2 vs serial; min-over-pending is exact.
				pendT := func(v int) (float64, bool) {
					if v < N {
						if f.t[v] >= T {
							return 0, false
						}
						return f.t[v], true
					}
					e := v - N
					if f.et[e] >= T {
						return 0, false
					}
					return f.et[e], true
				}
				beats := func(v int) bool {
					tv, pend := pendT(v)
					return pend && (tv < evt[u] || (tv == evt[u] && v < u))
				}
				for c := 0; c < nc && win; c++ {
					i := int(cells[c])
					if i != u && beats(i) {
						win = false
					}
					for _, jv := range f.nb[i] {
						v := int(jv)
						if v != u && beats(v) {
							win = false
							break
						}
					}
					if !win {
						break
					}
					for _, ev2 := range f.ce[i] {
						v := N + int(ev2)
						if v != u && beats(v) {
							win = false
							break
						}
					}
				}
				sel[u] = win
			}
		})
		cnt := 0
		for u := 0; u < NT; u++ {
			if sel[u] {
				cnt++
			}
		}
		if cnt == 0 {
			break
		}
		pfor(NT, workers, func(lo, hi int) {
			for u := lo; u < hi; u++ {
				if !sel[u] {
					continue
				}
				if u < N {
					f.advance(u)
				} else {
					f.advanceEdge(u - N)
				}
			}
		})
		acc += float64(cnt)
		ev += int64(cnt)
		rounds++
	}
	if rounds > 0 {
		width = acc / float64(rounds)
	}
	return ev, rounds, width
}

func lcMaxdiffTh(a, b *lcFab) float64 {
	m := 0.0
	for i := 0; i < a.N; i++ {
		d := math.Abs(a.th[i] - b.th[i])
		if d > m {
			m = d
		}
	}
	return m
}

// RunP2LC — the P2 readiness experiment (battery `p2lc`).
func (s *Sim) RunP2LC() {
	P := &s.P
	const (
		// Kc=0.2: strong enough that arrival-order scheduling visibly
		// corrupts the solution (the R3 control), weak enough that the
		// Kuramoto sector stays single-basin so R2 measures integrator
		// order, not attractor selection (Kc=0.8 measured ratio 1.02 —
		// multistable phase-locking, not a scheduler property).
		kc   = 0.2
		ge   = 0.05
		dil  = 8.0 // quiet-cell dilation R
		T    = 5.0
		refD = 64.0
	)
	dt := P.Dt
	look := 1.25 * dil * dt
	fprintf(s.Out, "# freecell-go — v90 P2 local-clock prototype (exp=p2lc)\n")
	fprintf(s.Out, "# P2LC reduced dynamics: dth=w_detune + Kc*sum sin(pub_j-th); paired antisymmetric E transfers\n")
	fprintf(s.Out, "# P2LC params: Kc=%g gE=%g dt=%g T=%g dilation R=%g look=%g blob amp=%g sigma=%g seed=%d\n",
		kc, ge, dt, T, dil, look, P.Amp, P.Sigma, P.Seed)

	// ---- criteria substrate (L=16) ----
	s.rng = P.Seed
	if s.rng == 0 {
		s.rng = 1
	}
	for w := 0; w < 8; w++ {
		s.xrand()
	}
	f0, nact := s.lcBuild(16, dil, kc, ge, dt)
	degMax := 0
	for i := 0; i < f0.N; i++ {
		if len(f0.nb[i]) > degMax {
			degMax = len(f0.nb[i])
		}
	}
	fprintf(s.Out, "# P2LC substrate L=16: N=%d ne=%d deg_max=%d active=%d quiet=%d\n",
		f0.N, f0.ne, degMax, nact, f0.N-nact)

	// R1 conservation
	fs := f0.clone()
	fs.runSync(T, dt)
	fa := f0.clone()
	evA, skew := fa.runLocal(T, look, 0, 0, 0)
	e0 := f0.totalE()
	fprintf(s.Out, "# RESULT p2lc_r1 drift_sync=%.3e drift_async=%.3e (E0=%.6f, events=%d)\n",
		(fs.totalE()-e0)/e0, (fa.totalE()-e0)/e0, e0, evA)

	// R2 convergence: async at dt and dt/2 vs a dt/refD sync reference
	fr := f0.clone()
	fr.runSync(T, dt/refD)
	e1 := lcMaxdiffTh(fa, fr)
	fh := f0.clone()
	// halve every clock (the dilation pattern is preserved)
	fh.h = append([]float64(nil), f0.h...)
	fh.eh = append([]float64(nil), f0.eh...)
	for i := range fh.h {
		fh.h[i] *= 0.5
	}
	for e := range fh.eh {
		fh.eh[e] *= 0.5
	}
	fh.runLocal(T, 0.5*look, 0, 0, 0)
	e2 := lcMaxdiffTh(fh, fr)
	fprintf(s.Out, "# RESULT p2lc_r2 err_dt=%.6f err_dt2=%.6f ratio=%.3f (first order: ~2)\n",
		e1, e2, e1/e2)

	// R3 determinism: rotated/reversed scans must not change a bit;
	// arrival order (sel=1) is the failing control
	f1 := f0.clone()
	f1.runLocal(T, look, 0, 17, 0)
	f2 := f0.clone()
	f2.runLocal(T, look, 0, 0, 1)
	f3 := f0.clone()
	f3.runLocal(T, look, 0, 17, 1)
	f4 := f0.clone()
	f4.runLocal(T, look, 1, 0, 0)
	fprintf(s.Out, "# RESULT p2lc_r3 rot=%.3e rev=%.3e revrot=%.3e arrival=%.3e\n",
		lcMaxdiffTh(f1, fa), lcMaxdiffTh(f2, fa), lcMaxdiffTh(f3, fa), lcMaxdiffTh(f4, fa))

	// condition 4: the tick counter diverges BY CONSTRUCTION (that
	// divergence IS the dilation) — record the skew a wrapped counter
	// would have to hold; any cyclic use needs M > 2*skew or it fails
	// silently. This prototype orders by local TIME, never by tau.
	fprintf(s.Out, "# RESULT p2lc_skew max_tick_skew=%d (a mod-M byte needs M > %d; ordering uses local time, never tau)\n",
		skew, 2*skew)

	// batch: bit-identical to the serial schedule at every worker count
	wmax := 0.0
	bmax := 0.0
	var rounds int64
	for _, wk := range []int{1, 2, 4, 8} {
		fb := f0.clone()
		_, r, wd := fb.runBatch(T, look, wk)
		d := lcMaxdiffTh(fb, fa)
		if d > bmax {
			bmax = d
		}
		if wd > wmax {
			wmax = wd
		}
		rounds = r
	}
	fprintf(s.Out, "# RESULT p2lc_batch maxdiff_vs_serial=%.3e workers=1,2,4,8 width=%.1f rounds=%d\n",
		bmax, wmax, rounds)

	// ---- the quiet-region economy: the scaling table ----
	// Fixed active blob, growing box: sync pays per cell; the local
	// clocks pay ~1/R on the quiet bath. Events are the deterministic
	// cost metric; ms are recorded for the wall-clock record.
	fprintf(s.Out, "# P2LC scaling: fixed blob (amp=%g sigma=%g), growing quiet bath\n", P.Amp, P.Sigma)
	for _, Lv := range []float64{16, 32, 64} {
		s.rng = P.Seed + uint64(Lv)
		for w := 0; w < 8; w++ {
			s.xrand()
		}
		fL, na := s.lcBuild(Lv, dil, kc, ge, dt)
		fsy := fL.clone()
		t0 := time.Now()
		evS := fsy.runSync(T, dt)
		msS := float64(time.Since(t0).Microseconds()) / 1000.0
		fba := fL.clone()
		t0 = time.Now()
		evB, _, _ := fba.runBatch(T, look, 8)
		msB := float64(time.Since(t0).Microseconds()) / 1000.0
		drift := (fba.totalE() - fL.totalE()) / fL.totalE()
		fprintf(s.Out, "# RESULT p2lc_scale L=%g N=%d ne=%d active=%d ev_sync=%d ev_async=%d ev_ratio=%.2f ms_sync=%.1f ms_async=%.1f drift=%.1e\n",
			Lv, fL.N, fL.ne, na, evS, evB, float64(evS)/float64(evB), msS, msB, drift)
	}
	fprintf(s.Out, "# RESULT p2lc done\n")
}
