package fab

import "math"

// v97 oriented-door leads. Byte-inert at Promo==0.

func (s *Sim) linkJ(i, j int) float64 {
	if i < 0 || j < 0 || i >= s.NC || j >= s.NC {
		return 0
	}
	if s.Em[i] <= 0 || s.Em[j] <= 0 {
		return 0
	}
	ai, aj := math.Sqrt(s.Em[i]), math.Sqrt(s.Em[j])
	c1, s1 := lutCos(s.th2[i]), lutSin(s.th2[i])
	c2, s2 := lutCos(s.th2[j]), lutSin(s.th2[j])
	return 2 * (ai*c1*aj*s2 - ai*s1*aj*c2)
}

func (s *Sim) linkJfield(i, j int) float64 {
	return 2 * (s.fa1[i]*s.fa2[j] - s.fa2[i]*s.fa1[j])
}

func (s *Sim) slotFromI(sl, i int) (j int, ux, uy, uz float64) {
	if s.sli[sl] == i {
		return s.slj[sl], s.sux[sl], s.suy[sl], s.suz[sl]
	}
	return s.sli[sl], -s.sux[sl], -s.suy[sl], -s.suz[sl]
}

func (s *Sim) cellDir(i int) (x, y, z float64) {
	for q := s.cls[i]; q < s.cls[i+1]; q++ {
		sl := s.clidx[q]
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		j, ux, uy, uz := s.slotFromI(sl, i)
		J := s.linkJ(i, j)
		x += J * ux
		y += J * uy
		z += J * uz
	}
	if x*x+y*y+z*z < 1e-24 {
		for q := s.cls[i]; q < s.cls[i+1]; q++ {
			sl := s.clidx[q]
			if s.sst[sl] == sFree || s.sA[sl] <= 0 {
				continue
			}
			j, ux, uy, uz := s.slotFromI(sl, i)
			J := s.linkJfield(i, j)
			x += J * ux
			y += J * uy
			z += J * uz
		}
	}
	return
}

// promoPick returns the oriented neighbour, or i if none qualifies.
// fore=true → max û·J (ahead); fore=false → min û·J (behind).
func (s *Sim) promoPick(i int, fore bool) int {
	if s.P.Promo <= 0 {
		return i
	}
	jx, jy, jz := s.cellDir(i)
	if jx*jx+jy*jy+jz*jz < 1e-24 {
		return i
	}
	best := i
	bestDot := -1e300
	if !fore {
		bestDot = 1e300
	}
	p := s.P.PGate
	if p <= 0 {
		p = 8
	}
	for q := s.cls[i]; q < s.cls[i+1]; q++ {
		sl := s.clidx[q]
		if s.sst[sl] == sFree || s.sA[sl] <= 0 {
			continue
		}
		j, ux, uy, uz := s.slotFromI(sl, i)
		if j < 0 || j >= s.NC || j == i {
			continue
		}
		if s.P.Promo >= 3 {
			g := math.Abs(math.Cos(s.th2[i] - s.th2[j]))
			if math.Pow(g, p) < 0.05 {
				continue
			}
		}
		dot := ux*jx + uy*jy + uz*jz
		if fore && dot > bestDot && dot > 0 {
			bestDot, best = dot, j
		}
		if !fore && dot < bestDot && dot < 0 {
			bestDot, best = dot, j
		}
	}
	return best
}

func (s *Sim) applyPromoCond(i int, d1, dsp float64) {
	dest, esSrc := i, i
	if s.P.Promo > 0 {
		if j := s.promoPick(i, true); j != i {
			dest = j
			if s.P.Promo >= 2 {
				esSrc = j
			}
		}
	}
	avail := s.Es[esSrc] - s.P.EsFloor
	if avail < 0 {
		avail = 0
	}
	if dsp > avail {
		dsp = avail
	}
	fac := 0.0
	if s.Ee[i] > 0 {
		fac = math.Sqrt((s.Ee[i] - d1) / s.Ee[i])
	}
	s.fa1[i] *= fac
	s.fa2[i] *= fac
	s.Ee[i] -= d1
	s.Es[esSrc] -= dsp
	emOld := s.Em[dest]
	s.Em[dest] += d1 + dsp
	s.pI[dest] += d1
	if s.P.AmpDoor > 0 {
		aph := math.Atan2(s.fa2[i], s.fa1[i])
		if emOld < 0 {
			emOld = 0
		}
		ro, rc := math.Sqrt(emOld), math.Sqrt(s.P.AmpDoor*d1)
		s1 := ro*lutCos(s.th2[dest]) + rc*lutCos(aph)
		s2 := ro*lutSin(s.th2[dest]) + rc*lutSin(aph)
		s.th2[dest] = math.Mod(math.Atan2(s2, s1)+8.0*TwoPi, TwoPi)
	} else if s.P.QpPhase > 0 && s.Em[dest] > 0 {
		aph := math.Atan2(s.fa2[i], s.fa1[i])
		mix := s.P.QpPhase * d1 / s.Em[dest]
		s.th2[dest] = math.Mod(s.th2[dest]+mix*wrapPi(aph-s.th2[dest])+8.0*TwoPi, TwoPi)
	}
}

func (s *Sim) applyPromoOut(i int, d2, bs float64) {
	dest := i
	if s.P.Promo > 0 {
		if k := s.promoPick(i, false); k != i {
			dest = k
		}
	}
	s.Em[i] -= d2
	s.pO[dest] += d2
	s.fieldInject(dest, d2-bs)
	s.Es[dest] += bs
}
