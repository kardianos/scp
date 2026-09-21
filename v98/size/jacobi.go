package size

import "math"

// jacobiSym: symmetric eigen-decomp of n×n row-major A.
// eval[k], evec column k. Copied from v97/chartexp (meter only).
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

func sortFloats(a []float64) {
	for i := 0; i < len(a); i++ {
		for j := i + 1; j < len(a); j++ {
			if a[j] < a[i] {
				a[i], a[j] = a[j], a[i]
			}
		}
	}
}
