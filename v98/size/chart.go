package size

import (
	"math"

	"scp/v97/hunt"
)

const chartNMax = 24

// LapSpec is the unweighted combinatorial Laplacian of a motif
// plus a 3+1 occupation of ψ = √Em exp(i θ) on that basis.
type LapSpec struct {
	N      int
	Evals  []float64 // sorted
	Omega  [3]float64
	Occ    [3]float64
	Res    float64
	OK     bool
	Match6 bool // |evals − {0,1,1,3,3,4}|_∞ ≤ 1e-6 (C6 book)
}

func motifLaplacian(fr hunt.Frame, cells []int) ([]float64, bool) {
	n := len(cells)
	if n < 1 || n > chartNMax {
		return nil, false
	}
	idx := make(map[int]int, n)
	for k, i := range cells {
		idx[i] = k
	}
	L := make([]float64, n*n)
	for _, ln := range fr.Links {
		a, okA := idx[ln.I]
		b, okB := idx[ln.J]
		if !okA || !okB || a == b {
			continue
		}
		L[a*n+a] += 1
		L[b*n+b] += 1
		L[a*n+b] -= 1
		L[b*n+a] -= 1
	}
	return L, true
}

func project31(fr hunt.Frame, cells []int, evals, evec []float64) (om, occ [3]float64, res float64) {
	n := len(cells)
	// group near-degenerate eigenvalues
	type shell struct {
		lam float64
		cs  []int
	}
	var sh []shell
	used := make([]bool, n)
	for i := 0; i < n; i++ {
		if used[i] {
			continue
		}
		s := shell{lam: evals[i], cs: []int{i}}
		used[i] = true
		for j := i + 1; j < n; j++ {
			if !used[j] && math.Abs(evals[j]-evals[i]) < 1e-6 {
				s.cs = append(s.cs, j)
				used[j] = true
			}
		}
		sh = append(sh, s)
	}
	// occupation of each eigenspace: sum_k |⟨v_k, ψ⟩|²
	psiR := make([]float64, n)
	psiI := make([]float64, n)
	var nrm float64
	for a, ci := range cells {
		c := fr.Cells[ci]
		amp := math.Sqrt(math.Max(c.Em, 0))
		psiR[a] = amp * math.Cos(c.Th2)
		psiI[a] = amp * math.Sin(c.Th2)
		nrm += c.Em
	}
	if nrm < 1e-18 {
		return
	}
	type scored struct {
		lam, n float64
	}
	sc := make([]scored, len(sh))
	for si, s := range sh {
		var o float64
		for _, k := range s.cs {
			var re, im float64
			for a := 0; a < n; a++ {
				v := evec[a*n+k]
				re += v * psiR[a]
				im += v * psiI[a]
			}
			o += re*re + im*im
		}
		sc[si] = scored{lam: s.lam, n: o / nrm}
	}
	for i := 0; i < len(sc); i++ {
		for j := i + 1; j < len(sc); j++ {
			if sc[j].n > sc[i].n {
				sc[i], sc[j] = sc[j], sc[i]
			}
		}
	}
	keep := 3
	if keep > len(sc) {
		keep = len(sc)
	}
	var kept float64
	for i := 0; i < keep; i++ {
		om[i] = math.Sqrt(math.Max(sc[i].lam, 0))
		occ[i] = sc[i].n
		kept += sc[i].n
	}
	res = 1 - kept
	if res < 0 {
		res = 0
	}
	return
}

func analyzeChart(fr hunt.Frame, cells []int) LapSpec {
	var out LapSpec
	L, ok := motifLaplacian(fr, cells)
	if !ok {
		return out
	}
	n := len(cells)
	eval, evec, jok := jacobiSym(L, n)
	if !jok {
		return out
	}
	sortFloats(eval)
	// evec columns are in original (unsorted) order; re-sort eval
	// already done. For projection we need evals aligned with
	// columns — recompute from the unsorted copy.
	evalRaw, evec, jok := jacobiSym(L, n)
	if !jok {
		return out
	}
	out.N = n
	out.Evals = eval
	out.Omega, out.Occ, out.Res = project31(fr, cells, evalRaw, evec)
	out.OK = true
	if n == 6 && len(eval) == 6 {
		want := []float64{0, 1, 1, 3, 3, 4}
		out.Match6 = true
		for i := range want {
			if math.Abs(eval[i]-want[i]) > 1e-6 {
				out.Match6 = false
				break
			}
		}
	}
	return out
}
