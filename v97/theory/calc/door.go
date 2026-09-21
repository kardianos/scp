// Package calc — v97 theory checks (PLAN step 3; symbolic-adjacent
// numerics, no kernel). Two results live here:
//
//  1. door.go — the door selection rule on C_N is the Fourier
//     transform of the door profile (qgraph L5 at closed form).
//  2. phibar.go — the corrected Φ-bar theorem (REVIEW_grok_next Q1).
package calc

import "math"

// Harmonics of the C_N cycle: chi_m(v) = exp(i 2π m v/N)/√N.
// For a DIAGONAL door profile d_v the matrix element between
// harmonics m (source) and m' (destination),
//
//	M(m'←m) = Σ_v conj(chi_m'(v)) d_v chi_m(v)
//	        = (1/N) Σ_v d_v exp(i 2π (m−m') v / N),
//
// depends only on Δm = m−m' and IS the DFT of d at Δm. Selectivity
// in Δm is therefore bounded by the uncertainty relation: a door
// supported on w adjacent vertices has |M| support of width ≳ N/w
// in Δm. A one-vertex or one-edge door cannot be Δm-selective.

// DoorAmp returns |M(m'←m)| for a diagonal (complex) door profile d.
func DoorAmp(d []complex128, m, mp int) float64 {
	n := len(d)
	var re, im float64
	for v := 0; v < n; v++ {
		ph := 2 * math.Pi * float64((m-mp)*v) / float64(n)
		c, s := math.Cos(ph), math.Sin(ph)
		dr, di := real(d[v]), imag(d[v])
		re += dr*c - di*s
		im += dr*s + di*c
	}
	return math.Hypot(re, im) / float64(n)
}

// Named door profiles on C_N.

// VertexDoor: d = δ_0 (the standing cell-local door).
func VertexDoor(n int) []complex128 {
	d := make([]complex128, n)
	d[0] = 1
	return d
}

// EdgeSymDoor: d = δ_0 + δ_1 (slot door, unsigned).
func EdgeSymDoor(n int) []complex128 {
	d := make([]complex128, n)
	d[0], d[1] = 1, 1
	return d
}

// EdgeDipoleDoor: d = δ_0 − δ_1 (oriented slot door, signed —
// the "promotion with a lever arm" read as a linear profile).
func EdgeDipoleDoor(n int) []complex128 {
	d := make([]complex128, n)
	d[0], d[1] = 1, -1
	return d
}

// WoundDoor: d_v = exp(i k 2π v / N) (the SHEAR-G organized-field
// door; delocalized over the whole cycle).
func WoundDoor(n, k int) []complex128 {
	d := make([]complex128, n)
	for v := 0; v < n; v++ {
		ph := 2 * math.Pi * float64(k*v) / float64(n)
		d[v] = complex(math.Cos(ph), math.Sin(ph))
	}
	return d
}

// Closed forms (proved by the DFT line above; asserted in tests):
//
//	vertex:       |M| = 1/N                    for every Δm  (flat)
//	edge sym:     |M| = (2/N)·|cos(π Δm / N)|  (mild low-pass)
//	edge dipole:  |M| = (2/N)·|sin(π Δm / N)|  (zero at Δm=0, RISING
//	              with |Δm| — anti-selective; elastic channel closed)
//	wound k:      |M| = δ_{Δm,k}               (exact pump)

// CurrentDoorAmp: the non-diagonal case — door coupling through the
// oriented link current operator J = i(|1⟩⟨0| − |0⟩⟨1|) on edge (0,1).
// ⟨chi_m'|J|chi_m⟩ = (i/N)(e^{i2πm/N} − e^{−i2πm'/N}), so
// |M| = (2/N)·|sin(π (m+m') / N)| — selective on the SUM m+m',
// not on Δm.
func CurrentDoorAmp(n, m, mp int) float64 {
	a1 := 2 * math.Pi * float64(m) / float64(n)
	a2 := -2 * math.Pi * float64(mp) / float64(n)
	re := math.Cos(a1) - math.Cos(a2)
	im := math.Sin(a1) - math.Sin(a2)
	return math.Hypot(re, im) / float64(n)
}
