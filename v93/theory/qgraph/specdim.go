package qgraph

import (
	"math"
	"sort"
)

// ---------------------------------------------------------------------------
// Spectral dimension — the measurable that fixes the Stefan exponent.
//
// UNIFIED_LAW §4: p_rad = d_s + 1. The adopted p_rad = 4 is therefore the
// claim that the field sector has d_s = 3, and NOTHING ELSE. On a graph d_s is
// not a modelling choice, it is a measurement, and these are the two ways to
// take it.
// ---------------------------------------------------------------------------

// TorusSpectrum returns the exact Laplacian eigenvalues of the discrete
// d-torus of side L, without diagonalising anything:
//
//	lambda(k) = sum_i 2 (1 - cos(2 pi k_i / L)),  k_i = 0..L-1
//
// It exists so d_s can be checked against a graph whose answer is known.
func TorusSpectrum(d, L int) []float64 {
	n := 1
	for i := 0; i < d; i++ {
		n *= L
	}
	out := make([]float64, 0, n)
	k := make([]int, d)
	for idx := 0; idx < n; idx++ {
		r := idx
		lam := 0.0
		for i := 0; i < d; i++ {
			k[i] = r % L
			r /= L
			lam += 2 * (1 - math.Cos(2*math.Pi*float64(k[i])/float64(L)))
		}
		out = append(out, lam)
	}
	sort.Float64s(out)
	return out
}

// SpectralDimensionFromEigenvalues fits the Weyl staircase
// N(lambda) ~ lambda^(d_s/2) over the lowest frac of the non-zero spectrum by
// least squares in log-log, and returns d_s = 2 * slope.
//
// ACCURACY (measured against exact tori, rates_test.go): at a few thousand
// vertices this estimator carries a stable finite-size bias of about +0.25
// (d=1 -> 1.24, d=2 -> 2.20, d=3 -> 3.31). It separates d_s = 1, 2, 3 cleanly
// but must NOT be quoted to better than ~0.3.
//
// Use ReturnDimension for anything load-bearing: it needs no diagonalisation
// and reproduces the exact tori to ~0.02 (d=1 -> 0.99, d=2 -> 2.00,
// d=3 -> 3.02). That is the estimator to run on a free-cell contact graph.
func SpectralDimensionFromEigenvalues(lam []float64, frac float64) float64 {
	var nz []float64
	for _, l := range lam {
		if l > 1e-9 {
			nz = append(nz, l)
		}
	}
	if len(nz) < 8 {
		return math.NaN()
	}
	sort.Float64s(nz)
	k := max(int(frac*float64(len(nz))), 8)
	if k > len(nz) {
		k = len(nz)
	}
	var sx, sy, sxx, sxy, nn float64
	for i := 0; i < k; i++ {
		x := math.Log(nz[i])
		y := math.Log(float64(i + 1))
		sx, sy, sxx, sxy, nn = sx+x, sy+y, sxx+x*x, sxy+x*y, nn+1
	}
	den := nn*sxx - sx*sx
	if nn < 4 || math.Abs(den) < 1e-30 {
		return math.NaN()
	}
	return 2 * (nn*sxy - sx*sy) / den
}

// ReturnDimension estimates d_s from the lazy random-walk return probability
//
//	p(t) ~ t^(-d_s/2),
//
// measured exactly by evolving one delta per source vertex. This is the
// estimator to use on a real free-cell contact graph, where the graph is far
// too large to diagonalise: it costs O(|E| * tmax * sources) and needs no
// eigenvectors at all.
//
// The walk is made lazy (P -> (I+P)/2) so a bipartite contact graph does not
// alternate, which would otherwise fake a dimension.
func ReturnDimension(g *Graph, tmax int, sources []int) float64 {
	if tmax < 8 || len(sources) == 0 {
		return math.NaN()
	}
	deg := make([]float64, g.N)
	for i := 0; i < g.N; i++ {
		for j := range g.Adj[i] {
			deg[i] += g.W[i][g.Adj[i][j]]
		}
		if deg[i] == 0 {
			deg[i] = 1
		}
	}
	pt := make([]float64, tmax+1)
	for _, s := range sources {
		if s < 0 || s >= g.N {
			continue
		}
		v := make([]float64, g.N)
		nv := make([]float64, g.N)
		v[s] = 1
		pt[0] += v[s]
		for t := 1; t <= tmax; t++ {
			for i := range nv {
				nv[i] = 0.5 * v[i] // lazy half
			}
			for i := 0; i < g.N; i++ {
				if v[i] == 0 {
					continue
				}
				share := 0.5 * v[i] / deg[i]
				for _, j := range g.Adj[i] {
					nv[j] += share * g.W[i][j]
				}
			}
			v, nv = nv, v
			pt[t] += v[s]
		}
	}
	inv := 1 / float64(len(sources))
	for t := range pt {
		pt[t] *= inv
	}
	// fit over the second half of the window: early times see the lattice,
	// late times see the finite volume.
	lo, hi := tmax/4, tmax/2
	if hi <= lo+3 {
		lo, hi = 2, tmax
	}
	var sx, sy, sxx, sxy, nn float64
	for t := lo; t <= hi; t++ {
		if pt[t] <= 0 {
			continue
		}
		x := math.Log(float64(t))
		y := math.Log(pt[t])
		sx, sy, sxx, sxy, nn = sx+x, sy+y, sxx+x*x, sxy+x*y, nn+1
	}
	if nn < 4 {
		return math.NaN()
	}
	den := nn*sxx - sx*sx
	if math.Abs(den) < 1e-30 {
		return math.NaN()
	}
	slope := (nn*sxy - sx*sy) / den
	return -2 * slope
}
