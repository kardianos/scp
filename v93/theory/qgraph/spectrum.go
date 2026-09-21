package qgraph

import (
	"math"
	"sort"
)

// Mode is one graph harmonic: an eigenpair of the Laplacian.
//
// Omega = sqrt(Lambda) is the pitch (wave-like convention: L is the discrete
// -nabla^2, so omega = sqrt(lambda) is the frequency of the mode). The action
// atom of the mode is Eps = A0*Omega/(2*pi), matching the programme's door.
type Mode struct {
	Index  int
	Lambda float64   // Laplacian eigenvalue
	Omega  float64   // sqrt(Lambda) — the mode's pitch
	Phi    []float64 // L2-normalised real eigenvector
}

// Eps returns the action atom of the mode, eps = A0*omega/(2*pi).
func (m Mode) Eps(a0 float64) float64 { return a0 * m.Omega / (2 * math.Pi) }

// Spectrum returns the Laplacian modes of g, sorted by eigenvalue.
func Spectrum(g *Graph) []Mode {
	vals, vecs := jacobiEigen(g.Laplacian())
	idx := make([]int, len(vals))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool { return vals[idx[a]] < vals[idx[b]] })

	modes := make([]Mode, len(vals))
	for r, i := range idx {
		lam := vals[i]
		if lam < 0 && lam > -1e-9 {
			lam = 0 // the Laplacian is PSD; clean the numerical zero
		}
		phi := make([]float64, len(vals))
		for k := range phi {
			phi[k] = vecs[k][i]
		}
		modes[r] = Mode{Index: r, Lambda: lam, Omega: math.Sqrt(math.Max(lam, 0)), Phi: phi}
	}
	return modes
}

func sq(x float64) float64 { return x * x }

// jacobiEigen diagonalises a real symmetric matrix by cyclic Jacobi rotations.
// It returns eigenvalues and a matrix whose COLUMNS are eigenvectors.
//
// Jacobi is chosen deliberately: it is exact-in-the-limit, dependency-free and
// orthogonality-preserving by construction, which is the same reason pass U
// uses Givens rotations for dense transport.
func jacobiEigen(a [][]float64) ([]float64, [][]float64) {
	n := len(a)
	A := make([][]float64, n)
	V := make([][]float64, n)
	for i := 0; i < n; i++ {
		A[i] = append([]float64(nil), a[i]...)
		V[i] = make([]float64, n)
		V[i][i] = 1
	}
	for sweep := 0; sweep < 200; sweep++ {
		off := 0.0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				off += sq(A[i][j])
			}
		}
		if off < 1e-26 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(A[p][q]) < 1e-18 {
					continue
				}
				theta := (A[q][q] - A[p][p]) / (2 * A[p][q])
				t := 1.0
				if theta >= 0 {
					t = 1 / (theta + math.Sqrt(theta*theta+1))
				} else {
					t = -1 / (-theta + math.Sqrt(theta*theta+1))
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := A[k][p], A[k][q]
					A[k][p] = c*akp - s*akq
					A[k][q] = s*akp + c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := A[p][k], A[q][k]
					A[p][k] = c*apk - s*aqk
					A[q][k] = s*apk + c*aqk
				}
				for k := 0; k < n; k++ {
					vkp, vkq := V[k][p], V[k][q]
					V[k][p] = c*vkp - s*vkq
					V[k][q] = s*vkp + c*vkq
				}
			}
		}
	}
	vals := make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = A[i][i]
	}
	return vals, V
}
