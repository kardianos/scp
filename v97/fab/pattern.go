package fab

import "math"

// pattern_only (v99 Tact 1): mass is not a second hoppable fluid.
//   0 = legacy Em + dense hop (byte-inert default)
//   1 = Es + fa only; condensation is space↔pattern; no dense hop
//   2 = 1 plus a derived wrap-gate on field hops (triangles and 4-cycles)
//
// Ledger when on: Σ(Es + |fa|² + |χ|²). Em is absorbed into fa after
// seed and at the start of each step so leftover seed paths cannot
// mint a dense species.

func (s *Sim) patternOnly() bool { return s.P.PatternOnly != 0 }

func (s *Sim) spaceDef(i int) float64 {
	d := s.P.Es0 - s.Es[i]
	if d < 0 {
		return 0
	}
	return d
}

func (s *Sim) occupancyX(i int) float64 {
	P := &s.P
	var x float64
	if s.patternOnly() {
		x = (s.Ee[i] + s.flload[i]) / P.Cap
	} else {
		x = (s.Em[i] + s.flload[i]) / P.Cap
	}
	if P.DoorG > 0 && P.ChiLoad > 0 {
		x += (s.dc1[i]*s.dc1[i] + s.dc2[i]*s.dc2[i]) / P.Cap
	}
	return x
}

func (s *Sim) presentMatter(i int) bool {
	if s.patternOnly() {
		return s.spaceDef(i) >= s.P.WfFloor
	}
	return s.Em[i] >= s.P.WfFloor
}

func (s *Sim) holdings(i int) float64 {
	if s.patternOnly() {
		return s.Ee[i]
	}
	return s.Em[i]
}

// absorbEm folds any dense tank into fa and zeros Em.
// Called after seed and at the start of step when pattern_only>0.
func (s *Sim) absorbEm() {
	if !s.patternOnly() {
		return
	}
	for i := 0; i < s.NC; i++ {
		if s.Em[i] > 0 {
			s.fieldInject(i, s.Em[i])
		}
		s.Em[i] = 0
	}
}

// applyPatternCond: beat-gated conversion with no Em tank.
// Screen (scond): field → space (the click absorbs the beam).
// Bulk: space → pattern (Es down, energy into fa). d1 is the atom
// that sets the space demand; traveling fa is not a source tank.
func (s *Sim) applyPatternCond(i int, d1, dsp float64) {
	if d1 <= 0 {
		return
	}
	if s.scond[i] != 0 {
		if d1 > s.Ee[i] {
			d1 = s.Ee[i]
		}
		if d1 <= 0 {
			return
		}
		if s.Ee[i] > 0 {
			fac := math.Sqrt((s.Ee[i] - d1) / s.Ee[i])
			s.fa1[i] *= fac
			s.fa2[i] *= fac
		}
		s.Ee[i] -= d1
		s.Es[i] += d1
		return
	}
	avail := s.Es[i] - s.P.EsFloor
	if avail < 0 {
		avail = 0
	}
	if dsp > avail {
		dsp = avail
	}
	if dsp > 0 {
		s.Es[i] -= dsp
		s.fieldInject(i, dsp)
	}
}

// applyPatternEvap: over-cap pattern → space.
func (s *Sim) applyPatternEvap(i int, d2 float64) {
	if d2 > s.Ee[i] {
		d2 = s.Ee[i]
	}
	if d2 <= 0 {
		return
	}
	if s.Ee[i] > 0 {
		fac := math.Sqrt((s.Ee[i] - d2) / s.Ee[i])
		s.fa1[i] *= fac
		s.fa2[i] *= fac
	}
	s.Ee[i] -= d2
	s.Es[i] += d2
}

func (s *Sim) faPhase(i int) float64 {
	if s.fa1[i] == 0 && s.fa2[i] == 0 {
		return 0
	}
	return math.Atan2(s.fa2[i], s.fa1[i])
}

func (s *Sim) linked(i, j int) bool {
	if i < 0 || j < 0 || i == j {
		return false
	}
	for q := s.cls[i]; q < s.cls[i+1]; q++ {
		sl := s.clidx[q]
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		o := s.sli[sl] + s.slj[sl] - i
		if o == j {
			return true
		}
	}
	return false
}

func (s *Sim) faceWinding(ids ...int) float64 {
	acc := 0.0
	n := len(ids)
	if n < 3 {
		return 0
	}
	for a := 0; a < n; a++ {
		acc += wrapPi(s.faPhase(ids[(a+1)%n]) - s.faPhase(ids[a]))
	}
	return acc
}

func (s *Sim) facePresent(ids ...int) bool {
	for _, i := range ids {
		if !s.presentMatter(i) {
			return false
		}
	}
	return true
}

// wrapHopBlocked is the derived door: a 3-space field hop does not
// fire if the link sits on a space-deficit face whose fa-phase winds.
// Triangles and 4-cycles only; C_6 has girth 6 and is NOT gated here.
func (s *Sim) wrapHopBlocked(sl int) bool {
	if s.P.PatternOnly < 2 {
		return false
	}
	if s.sst[sl] == sFree || s.sA[sl] <= 0 {
		return false
	}
	i, j := s.sli[sl], s.slj[sl]
	if !s.presentMatter(i) || !s.presentMatter(j) {
		return false
	}
	const windMin = 3.0 // |Σ Δθ| ≳ π
	for q := s.cls[i]; q < s.cls[i+1]; q++ {
		s2 := s.clidx[q]
		if s.sst[s2] == sFree || s.sA[s2] <= 0 {
			continue
		}
		k := s.slj[s2]
		if s.sli[s2] != i {
			k = s.sli[s2]
		}
		if k == j || k == i {
			continue
		}
		if s.linked(j, k) && s.facePresent(i, j, k) {
			if math.Abs(s.faceWinding(i, j, k)) > windMin {
				return true
			}
		}
		// 4-cycles i-j-b-a
		for qj := s.cls[j]; qj < s.cls[j+1]; qj++ {
			s3 := s.clidx[qj]
			if s.sst[s3] == sFree || s.sA[s3] <= 0 {
				continue
			}
			b := s.slj[s3]
			if s.sli[s3] != j {
				b = s.sli[s3]
			}
			if b == i || b == j || b == k {
				continue
			}
			if s.linked(k, b) && s.facePresent(i, j, b, k) {
				if math.Abs(s.faceWinding(i, j, b, k)) > windMin {
					return true
				}
			}
		}
	}
	return false
}
