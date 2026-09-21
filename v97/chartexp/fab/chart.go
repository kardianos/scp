package fab

// v97/chartexp — forked inner-chart probe. Byte-inert at Chart==0.
//
// A chart is a small contact Laplacian (a cell's star, or a small
// live-Em motif). We keep 3 shells + a residual bucket (3+1) and
// optionally feed the occupancy-weighted pitch back into the door
// (grain) or into w2e (mix). This is an experiment fork, not a face.

import "math"

const chartNMax = 24
const chartKeep = 3

func (s *Sim) chartAlloc(nc int) {
	s.chW = make([]float64, nc)
	s.chN0 = make([]float64, nc)
	s.chN1 = make([]float64, nc)
	s.chN2 = make([]float64, nc)
	s.chRes = make([]float64, nc)
	s.chO0 = make([]float64, nc)
	s.chO1 = make([]float64, nc)
	s.chO2 = make([]float64, nc)
}

// chartStep rebuilds 3+1 occupations.
//   chart=1  diagnostic only
//   chart=2  clean grain: rebuild here; door ε uses chW; w2e UNTOUCHED
//   chart=3  mix: overwrite w2e after pass-1 (hops+beat+door)
//   chart=4  mode door: rebuild here; ε = A0|ω0−w1e|/2π; price residual
func (s *Sim) chartStep(when string) {
	if s.P.Chart <= 0 {
		return
	}
	if when == "door" && s.P.Chart != 2 && s.P.Chart != 4 {
		return
	}
	if when == "freq" && s.P.Chart != 3 {
		return
	}
	s.chartRebuild()
	if when == "freq" && s.P.Chart == 3 {
		for i := 0; i < s.NC; i++ {
			if s.chW[i] > 1e-12 {
				s.w2e[i] = s.chW[i]
			}
		}
	}
}

// chartWDense is the dense-side pitch for door atoms. Table w2e unless
// clean grain (weff) or mode door (dominant-shell ω).
func (s *Sim) chartWDense(i int) float64 {
	if i < 0 || i >= s.NC {
		return s.w2e[i]
	}
	switch s.P.Chart {
	case 2:
		if s.chW[i] > 1e-12 {
			return s.chW[i]
		}
	case 4:
		if s.chO0[i] > 1e-12 {
			return s.chO0[i]
		}
	}
	return s.w2e[i]
}

// chartScale prices residual occupation (chart=4). 1 if the book is
// clean; 0 if there is nothing occupied (refuse).
func (s *Sim) chartScale(i int) float64 {
	if s.P.Chart != 4 || i < 0 || i >= s.NC {
		return 1
	}
	kept := s.chN0[i] + s.chN1[i]
	tot := kept + s.chN2[i] + s.chRes[i]
	if tot < 1e-12 || kept < 1e-9 {
		return 0
	}
	return kept / tot
}

// chartGapEps is the L3 promotion atom A0|ω_a−ω_b|/2π between the
// occupied shell and the field pitch. ok=false → caller uses table ε.
func (s *Sim) chartGapEps(i int, wField float64) (eps float64, ok bool) {
	if s.P.Chart != 4 || i < 0 || i >= s.NC || s.chO0[i] <= 1e-12 {
		return 0, false
	}
	gap := math.Abs(s.chO0[i] - wField)
	if gap < 1e-12 {
		return 0, false
	}
	return s.A0eff * gap / TwoPi, true
}

func (s *Sim) chartDoorPair(i int, wField, wDense float64) (eSrc, eDst, scale float64) {
	scale = s.chartScale(i)
	wd := wDense
	if s.P.Chart == 2 || s.P.Chart == 4 {
		wd = s.chartWDense(i)
	}
	eSrc = s.A0eff * wField / TwoPi
	eDst = s.A0eff * wd / TwoPi
	if eps, ok := s.chartGapEps(i, wField); ok {
		eSrc, eDst = eps, eps
	}
	return eSrc, eDst, scale
}

func (s *Sim) chartRebuild() {
	for i := 0; i < s.NC; i++ {
		s.chW[i] = 0
		s.chN0[i], s.chN1[i], s.chN2[i], s.chRes[i] = 0, 0, 0, 0
		s.chO0[i], s.chO1[i], s.chO2[i] = 0, 0, 0
	}
	scope := s.P.ChartScope
	if scope == "" {
		scope = "star"
	}
	if scope == "motif" {
		s.chartMotifs()
		return
	}
	for i := 0; i < s.NC; i++ {
		if s.Em[i] <= 0 {
			continue
		}
		s.chartStar(i)
	}
}

func (s *Sim) chartStar(i int) {
	idx := make([]int, 0, 8)
	idx = append(idx, i)
	seen := map[int]bool{i: true}
	for q := s.cls[i]; q < s.cls[i+1]; q++ {
		sl := s.clidx[q]
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		j, _, _, _ := s.slotFromI(sl, i)
		if j < 0 || j >= s.NC || seen[j] {
			continue
		}
		seen[j] = true
		idx = append(idx, j)
		if len(idx) >= chartNMax {
			break
		}
	}
	s.chartFill(idx)
}

func (s *Sim) chartMotifs() {
	parent := make([]int, s.NC)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(a int) int {
		for parent[a] != a {
			parent[a] = parent[parent[a]]
			a = parent[a]
		}
		return a
	}
	union := func(a, b int) {
		ra, rb := find(a), find(b)
		if ra != rb {
			parent[rb] = ra
		}
	}
	live := func(i int) bool { return i >= 0 && i < s.NC && s.Em[i] > 0.05 }
	for sl := 0; sl < s.NSLOT; sl++ {
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		a, b := s.sli[sl], s.slj[sl]
		if live(a) && live(b) {
			union(a, b)
		}
	}
	groups := map[int][]int{}
	for i := 0; i < s.NC; i++ {
		if !live(i) {
			continue
		}
		r := find(i)
		groups[r] = append(groups[r], i)
	}
	for _, g := range groups {
		if len(g) < 2 {
			continue
		}
		if len(g) > chartNMax {
			for _, i := range g {
				s.chartStar(i)
			}
			continue
		}
		s.chartFill(g)
	}
}

func (s *Sim) chartFill(idx []int) {
	n := len(idx)
	if n < 2 {
		return
	}
	loc := make(map[int]int, n)
	for a, id := range idx {
		loc[id] = a
	}
	L := make([]float64, n*n)
	acct := s.P.ChartAcct
	if acct == "" {
		acct = "energy"
	}
	for a, i := range idx {
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sst[sl] == sFree || s.sA[sl] <= 0 {
				continue
			}
			j, _, _, _ := s.slotFromI(sl, i)
			b, ok := loc[j]
			if !ok || b <= a {
				continue
			}
			w := 1.0
			switch acct {
			case "mass":
				w = math.Sqrt(math.Max(s.Em[i], 0) * math.Max(s.Em[j], 0))
				if w < 1e-12 {
					w = 1e-12
				}
			case "area":
				w = s.sA[sl]
				if w < 1e-12 {
					w = 1e-12
				}
			}
			L[a*n+a] += w
			L[b*n+b] += w
			L[a*n+b] -= w
			L[b*n+a] -= w
		}
	}
	eval, evec, ok := jacobiSym(L, n)
	if !ok {
		return
	}
	// project ψ_m = √Em (cos θ, sin θ) — occupation is |⟨φ|ψ⟩|²
	// using the real/imag pair as a 2-vector dotted per mode.
	cr := make([]float64, n)
	ci := make([]float64, n)
	var etot float64
	for a, i := range idx {
		amp := math.Sqrt(math.Max(s.Em[i], 0))
		cr[a] = amp * lutCos(s.th2[i])
		ci[a] = amp * lutSin(s.th2[i])
		etot += s.Em[i]
	}
	type shell struct{ ω, occ float64 }
	sh := make([]shell, n)
	for k := 0; k < n; k++ {
		var pr, pi float64
		for a := 0; a < n; a++ {
			v := evec[a*n+k]
			pr += v * cr[a]
			pi += v * ci[a]
		}
		ω := 0.0
		if eval[k] > 0 {
			ω = math.Sqrt(eval[k])
		}
		sh[k] = shell{ω: ω, occ: pr*pr + pi*pi}
	}
	// keep the 3 most occupied shells (not the 3 lowest ω — composition,
	// not a ground-state truncation).
	for a := 0; a < n; a++ {
		for b := a + 1; b < n; b++ {
			if sh[b].occ > sh[a].occ {
				sh[a], sh[b] = sh[b], sh[a]
			}
		}
	}
	var kept, wsum float64
	o0, o1, o2 := 0.0, 0.0, 0.0
	n0, n1, n2 := 0.0, 0.0, 0.0
	if n > 0 {
		n0, o0 = sh[0].occ, sh[0].ω
		kept += n0
		wsum += n0 * o0
	}
	if n > 1 {
		n1, o1 = sh[1].occ, sh[1].ω
		kept += n1
		wsum += n1 * o1
	}
	if n > 2 {
		n2, o2 = sh[2].occ, sh[2].ω
		kept += n2
		wsum += n2 * o2
	}
	res := etot - kept
	if res < 0 {
		res = 0
	}
	weff := 0.0
	if kept > 1e-12 {
		weff = wsum / kept
	}
	for _, i := range idx {
		s.chW[i] = weff
		s.chN0[i], s.chN1[i], s.chN2[i] = n0, n1, n2
		s.chRes[i] = res
		s.chO0[i], s.chO1[i], s.chO2[i] = o0, o1, o2
	}
}

func (s *Sim) chartDiag() {
	if s.P.Chart <= 0 {
		return
	}
	var nLive, nOn int
	var wbar, n0, n1, n2, res, dw float64
	for i := 0; i < s.NC; i++ {
		if s.Em[i] <= 0 {
			continue
		}
		nLive++
		if s.chW[i] > 0 {
			nOn++
			wbar += s.chW[i]
			n0 += s.chN0[i]
			n1 += s.chN1[i]
			n2 += s.chN2[i]
			res += s.chRes[i]
			d := s.chW[i] - s.w2e[i]
			if d < 0 {
				d = -d
			}
			dw += d
		}
	}
	if nOn == 0 {
		fprintf(s.Out, "# CHART t=%.2f live=%d on=0 W=%.3f\n", s.simT, nLive, s.chartWind())
		return
	}
	inv := 1.0 / float64(nOn)
	o0, o1, o2 := 0.0, 0.0, 0.0
	for i := 0; i < s.NC; i++ {
		if s.chW[i] > 0 {
			o0, o1, o2 = s.chO0[i], s.chO1[i], s.chO2[i]
			break
		}
	}
	fprintf(s.Out, "# CHART t=%.2f live=%d on=%d weff=%.5f n0=%.4f n1=%.4f n2=%.4f res=%.4f |weff-w2e|=%.4f o=(%.3f,%.3f,%.3f) W=%.4f\n",
		s.simT, nLive, nOn, wbar*inv, n0*inv, n1*inv, n2*inv, res*inv, dw*inv,
		o0, o1, o2, s.chartWind())
}

func (s *Sim) chartWind() float64 {
	if s.P.Exp != "ring" || s.NC < 3 || s.NC > 16 {
		return 0
	}
	var acc float64
	for i := 0; i < s.NC; i++ {
		j := (i + 1) % s.NC
		acc += wrapPi(s.th2[j] - s.th2[i])
	}
	return acc / TwoPi
}

func (s *Sim) chartResult() {
	if s.P.Chart <= 0 {
		return
	}
	var nOn int
	var wbar, n0, n1, n2, res float64
	for i := 0; i < s.NC; i++ {
		if s.chW[i] > 0 {
			nOn++
			wbar += s.chW[i]
			n0 += s.chN0[i]
			n1 += s.chN1[i]
			n2 += s.chN2[i]
			res += s.chRes[i]
		}
	}
	if nOn == 0 {
		fprintf(s.Out, "# RESULT chart on=0\n")
		return
	}
	inv := 1.0 / float64(nOn)
	fprintf(s.Out, "# RESULT chart on=%d weff=%.6g n0=%.6g n1=%.6g n2=%.6g res=%.6g\n",
		nOn, wbar*inv, n0*inv, n1*inv, n2*inv, res*inv)
	if s.P.Exp == "ring" && s.NC >= 3 && s.NC <= 16 {
		// circulation of th2 around roster order (ring seed is cyclic)
		var acc float64
		for i := 0; i < s.NC; i++ {
			j := (i + 1) % s.NC
			acc += wrapPi(s.th2[j] - s.th2[i])
		}
		fprintf(s.Out, "# RESULT chart_wind %.6g\n", acc/TwoPi)
	}
}

// jacobiSym: symmetric eigen-decomp. eval[k], evec column k.
func jacobiSym(A []float64, n int) (eval []float64, evec []float64, ok bool) {
	if n < 1 || len(A) < n*n {
		return nil, nil, false
	}
	a := append([]float64(nil), A...)
	v := make([]float64, n*n)
	for i := 0; i < n; i++ {
		v[i*n+i] = 1
	}
	const maxIt = 32 * 24
	for it := 0; it < maxIt; it++ {
		p, q := 0, 1
		best := 0.0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				x := math.Abs(a[i*n+j])
				if x > best {
					best, p, q = x, i, j
				}
			}
		}
		if best < 1e-14 {
			break
		}
		app, aqq, apq := a[p*n+p], a[q*n+q], a[p*n+q]
		τ := (aqq - app) / (2 * apq)
		t := 1.0
		if τ != 0 {
			sgn := 1.0
			if τ < 0 {
				sgn = -1
			}
			t = sgn / (math.Abs(τ) + math.Sqrt(1+τ*τ))
		}
		c := 1 / math.Sqrt(1+t*t)
		s := t * c
		for k := 0; k < n; k++ {
			if k == p || k == q {
				continue
			}
			aik, aiq := a[k*n+p], a[k*n+q]
			a[k*n+p] = c*aik - s*aiq
			a[p*n+k] = a[k*n+p]
			a[k*n+q] = s*aik + c*aiq
			a[q*n+k] = a[k*n+q]
		}
		a[p*n+p] = c*c*app - 2*c*s*apq + s*s*aqq
		a[q*n+q] = s*s*app + 2*c*s*apq + c*c*aqq
		a[p*n+q], a[q*n+p] = 0, 0
		for k := 0; k < n; k++ {
			vip, viq := v[k*n+p], v[k*n+q]
			v[k*n+p] = c*vip - s*viq
			v[k*n+q] = s*vip + c*viq
		}
	}
	eval = make([]float64, n)
	for i := 0; i < n; i++ {
		eval[i] = a[i*n+i]
	}
	return eval, v, true
}
