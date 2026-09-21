package qgraph

import (
	"math"
	"math/cmplx"
)

// ---------------------------------------------------------------------------
// Doors and selection rules.
//
// A conversion door is a LOCAL operator: it multiplies the state by a profile
// d_v on the vertices. Its matrix element between harmonics is
//
//	M[a][b] = sum_v conj(phi_a(v)) * d_v * phi_b(v)
//
// which on a ring is exactly the discrete Fourier transform of d at the
// winding difference. That single fact sorts every measured door result in the
// programme:
//
//	uniform door   d_v = 1              -> only Delta m = 0     (neutral light)
//	local door     d_v = delta(v,v0)    -> ALL Delta m allowed  (scrambles)
//	wound door     d_v = exp(2 pi i k v/N) -> exactly Delta m = k (pumps)
//
// So Face B's "door traffic destroys winding" is not a bug: a point-contact
// door is, by construction, the operator with flat support over every winding
// transfer. And "shear genesis" is not a metaphor: an organized field with
// phase wind k is the operator that moves circulation in quantized steps.
// ---------------------------------------------------------------------------

// DoorProfile is the per-vertex complex weight of a conversion door.
type DoorProfile []complex128

// UniformDoor couples every vertex equally: the neutral-radiation door.
func UniformDoor(n int) DoorProfile {
	d := make(DoorProfile, n)
	for i := range d {
		d[i] = 1
	}
	return d
}

// LocalDoor couples a single vertex: the point-contact door.
func LocalDoor(n, v int) DoorProfile {
	d := make(DoorProfile, n)
	if v >= 0 && v < n {
		d[v] = 1
	}
	return d
}

// WoundDoor couples with a phase wind of k turns around a ring of n vertices:
// the organized / sheared field.
func WoundDoor(n, k int) DoorProfile {
	d := make(DoorProfile, n)
	for j := range d {
		th := 2 * math.Pi * float64(k) * float64(j) / float64(n)
		d[j] = cmplx.Exp(complex(0, th))
	}
	return d
}

// ringHarmonic returns chi_m(j) = exp(2 pi i m j/N)/sqrt(N).
func ringHarmonic(n, m, j int) complex128 {
	th := 2 * math.Pi * float64(m) * float64(j) / float64(n)
	return cmplx.Exp(complex(0, th)) / complex(math.Sqrt(float64(n)), 0)
}

// RingMatrixElement returns <chi_a | D | chi_b> on C_n for a door profile.
func RingMatrixElement(n, a, b int, d DoorProfile) complex128 {
	var s complex128
	for j := 0; j < n && j < len(d); j++ {
		s += cmplx.Conj(ringHarmonic(n, a, j)) * d[j] * ringHarmonic(n, b, j)
	}
	return s
}

// AllowedTransfers returns the set of winding transfers Delta m = a - b for
// which the door has a non-negligible matrix element on C_n.
//
// This is the selection rule, derived rather than imposed.
func AllowedTransfers(n int, d DoorProfile, tol float64) []int {
	seen := map[int]bool{}
	var out []int
	for a := 0; a < n; a++ {
		for b := 0; b < n; b++ {
			if cmplx.Abs(RingMatrixElement(n, a, b, d)) > tol {
				dm := ((a-b)%n + n) % n
				if !seen[dm] {
					seen[dm] = true
					out = append(out, dm)
				}
			}
		}
	}
	// deterministic order
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j] < out[j-1]; j-- {
			out[j], out[j-1] = out[j-1], out[j]
		}
	}
	return out
}

// CirculationTransfer is the mean winding a door imparts per event, weighted
// by |matrix element|^2 over all transitions.
//
// A wound door returns exactly k. A uniform door returns 0 with zero variance
// (it cannot touch circulation). A local door returns 0 with LARGE variance —
// it moves winding every event, just with no preferred sign. Mean alone is
// therefore not the meter; the pair (mean, variance) is.
func CirculationTransfer(n int, d DoorProfile) (mean, variance float64) {
	var w, sum, sumsq float64
	for a := 0; a < n; a++ {
		for b := 0; b < n; b++ {
			p := cmplx.Abs(RingMatrixElement(n, a, b, d))
			p *= p
			if p < 1e-18 {
				continue
			}
			dm := ((a-b)%n + n) % n
			if 2*dm > n {
				dm -= n
			}
			w += p
			if 2*dm == n {
				// the zone-boundary transfer is +-n/2 with no defined sign
				// (the mode is self-conjugate): it contributes to the spread
				// but not to the mean.
				sumsq += p * float64(dm) * float64(dm)
				continue
			}
			sum += p * float64(dm)
			sumsq += p * float64(dm) * float64(dm)
		}
	}
	if w == 0 {
		return 0, 0
	}
	mean = sum / w
	variance = sumsq/w - mean*mean
	return mean, variance
}

// GraphMatrixElement is the general (non-ring) door element between two real
// Laplacian harmonics: M = sum_v phi_a(v) * Re(d_v) * phi_b(v).
func GraphMatrixElement(modes []Mode, a, b int, d DoorProfile) float64 {
	s := 0.0
	for v := 0; v < len(modes[a].Phi) && v < len(d); v++ {
		s += modes[a].Phi[v] * real(d[v]) * modes[b].Phi[v]
	}
	return s
}
