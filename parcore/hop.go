package parcore

import "math"

// Givens2 is the blessed pairwise unitary primitive for dual-endpoint updates.
// It rotates the two-component amplitudes (a1i,a2i) and (a1j,a2j) by angle tau
// in a way that preserves the SUM OF SQUARES over the pair EXACTLY:
//
//	|a'i|² + |a'j|² == |ai|² + |aj|²   (a theorem of the construction)
//
// This is the v93/freecell pass-F (field) and pass-U (unitary dense channel)
// rotation (see v93/fab/step.go applyUhop and pass F). Drivers SHOULD compose
// dual-endpoint slot updates from it so the conservation contract holds by
// construction rather than by luck.
//
// The mapping is:
//
//	a1'i = cc*a1i + ss*a2j
//	a2'i = cc*a2i - ss*a1j
//	a1'j = cc*a1j + ss*a2i
//	a2'j = cc*a2j - ss*a1i
//
// (cc, ss) = (cos(tau), sin(tau)). At tau=0 this is the identity exactly.
func Givens2(a1i, a2i, a1j, a2j, tau float64) (n1i, n2i, n1j, n2j float64) {
	ss, cc := math.Sincos(tau)
	n1i = cc*a1i + ss*a2j
	n2i = cc*a2i - ss*a1j
	n1j = cc*a1j + ss*a2i
	n2j = cc*a2j - ss*a1i
	return
}

// Givens2Buf applies the same rotation in place on two cells of a slice, used
// by the boundary hop over the global publish buffer. get/set are the driver's
// accessors for the two amplitude components being rotated.
func Givens2Buf[C any](
	buf []C, i, j int, tau float64,
	get func(c *C) (a1, a2 float64),
	set func(c *C, a1, a2 float64),
) {
	if tau <= 0 {
		return
	}
	if tau > 0.5 {
		tau = 0.5
	}
	a1i, a2i := get(&buf[i])
	a1j, a2j := get(&buf[j])
	n1i, n2i, n1j, n2j := Givens2(a1i, a2i, a1j, a2j, tau)
	set(&buf[i], n1i, n2i)
	set(&buf[j], n1j, n2j)
}

// LensArea is the channel overlap area for two radii at separation d (the
// freecell lens). Drivers call it from their geom pass; it lives in the core
// because the geometry is shared across any free-cell-shaped driver.
func LensArea(d, ri, rj float64) float64 {
	if d >= ri+rj {
		return 0
	}
	t := d*d - rj*rj + ri*ri
	a2 := (4*d*d*ri*ri - t*t) / (4 * d * d)
	if a2 > 0 {
		return math.Pi * a2
	}
	if d < math.Abs(ri-rj) {
		rm := ri
		if rj < ri {
			rm = rj
		}
		return math.Pi * rm * rm
	}
	return 0
}
