package fab

import (
	"math"
	"testing"
)

func TestJacobiC6Spectrum(t *testing.T) {
	n := 6
	L := make([]float64, n*n)
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		L[i*n+i] += 1
		L[j*n+j] += 1
		L[i*n+j] -= 1
		L[j*n+i] -= 1
	}
	eval, evec, ok := jacobiSym(L, n)
	if !ok {
		t.Fatal("jacobi failed")
	}
	// C6: 0, 1, 1, 3, 3, 4
	want := []float64{0, 1, 1, 3, 3, 4}
	got := append([]float64(nil), eval...)
	for a := 0; a < n; a++ {
		for b := a + 1; b < n; b++ {
			if got[b] < got[a] {
				got[a], got[b] = got[b], got[a]
			}
		}
	}
	for k := range want {
		if math.Abs(got[k]-want[k]) > 1e-8 {
			t.Fatalf("eval[%d]=%g want %g (raw %v)", k, got[k], want[k], eval)
		}
	}
	// orthonormal columns
	for k := 0; k < n; k++ {
		var nrm float64
		for i := 0; i < n; i++ {
			nrm += evec[i*n+k] * evec[i*n+k]
		}
		if math.Abs(nrm-1) > 1e-10 {
			t.Fatalf("col %d norm %g", k, nrm)
		}
	}
}

func TestJacobiPath2(t *testing.T) {
	// two vertices, one edge: evals 0, 2
	L := []float64{1, -1, -1, 1}
	eval, _, ok := jacobiSym(L, 2)
	if !ok {
		t.Fatal("jacobi failed")
	}
	lo, hi := eval[0], eval[1]
	if lo > hi {
		lo, hi = hi, lo
	}
	if math.Abs(lo) > 1e-10 || math.Abs(hi-2) > 1e-10 {
		t.Fatalf("path2 evals %v want 0,2", eval)
	}
}
