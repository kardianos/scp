package core

import "math"

// Passes take Lo and/or Hi slices. Prefer Lo-only signatures when a pass
// does not touch scratch/planes so the Hi slab stays cold in cache.

// PassPitch updates radii and load-detuned pitches (Lo only).
func PassPitch(lo []CellLo, law Law, a, b int) {
	for i := a; i < b; i++ {
		c := &lo[i]
		ratio := c.Es / law.Es0
		if ratio < 0 {
			ratio = 0
		}
		c.Cr = c.Cr0 * math.Cbrt(ratio)
		x := (c.Em + c.Ee) / law.Cap
		det := 1.0 + law.QDetune*x
		c.W1e = law.W1 / det
		c.W2e = law.W2 / det
	}
}

// PassSpaceWants writes slot Swl from Lo energies.
func PassSpaceWants(lo []CellLo, slots []Slot, law Law, a, b int) {
	dt := law.Dt
	for s := a; s < b; s++ {
		sl := &slots[s]
		if sl.Alive == 0 || sl.A <= 0 {
			sl.Swl = 0
			continue
		}
		i, j := sl.I, sl.J
		pi := lo[i].Es + law.SDisp*(lo[i].Em+lo[i].Ee)
		pj := lo[j].Es + law.SDisp*(lo[j].Em+lo[j].Ee)
		dp := pi - pj
		if dp == 0 {
			sl.Swl = 0
			continue
		}
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sl.Swl = law.SK * dt * w * dp
	}
}

// PassSpaceScale fills sscl for owned indices using Lo.Es.
func PassSpaceScale(lo []CellLo, slots []Slot, cls, inc []int32, law Law, sscl []float64, a, b int) {
	for i := a; i < b; i++ {
		rq := 0.0
		for q := cls[i]; q < cls[i+1]; q++ {
			sl := &slots[inc[q]]
			f := sl.Swl
			if f == 0 {
				continue
			}
			if (f > 0 && sl.I == int32(i)) || (f < 0 && sl.J == int32(i)) {
				if f < 0 {
					rq += -f
				} else {
					rq += f
				}
			}
		}
		avail := lo[i].Es - law.EsFloor
		if avail <= 0 {
			sscl[i] = 0
		} else if rq > avail {
			sscl[i] = avail / rq
		} else {
			sscl[i] = 1
		}
	}
}

// PassSpaceApply writes EsNext on Hi from published sscl.
func PassSpaceApply(lo []CellLo, hi []CellHi, slots []Slot, cls, inc []int32, sscl []float64, a, b int) {
	for i := a; i < b; i++ {
		de := 0.0
		for q := cls[i]; q < cls[i+1]; q++ {
			sl := &slots[inc[q]]
			f := sl.Swl
			if f == 0 {
				continue
			}
			var src int32
			if f > 0 {
				src = sl.I
			} else {
				src = sl.J
			}
			mag := math.Abs(f) * sscl[src]
			if src == int32(i) {
				de -= mag
			} else {
				de += mag
			}
		}
		hi[i].EsNext = lo[i].Es + de
	}
}

// PassCommitSpace copies EsNext (Hi) → Es (Lo).
func PassCommitSpace(lo []CellLo, hi []CellHi, a, b int) {
	for i := a; i < b; i++ {
		lo[i].Es = hi[i].EsNext
	}
}

// PassFieldPrecess writes Fa*Next on Hi from Lo Fa + W1e.
func PassFieldPrecess(lo []CellLo, hi []CellHi, law Law, a, b int) {
	for i := a; i < b; i++ {
		c := &lo[i]
		h := &hi[i]
		ang := c.W1e * law.Dt
		cc, ss := math.Cos(ang), math.Sin(ang)
		a1, a2 := c.Fa1, c.Fa2
		h.Fa1Next = cc*a1 + ss*a2
		h.Fa2Next = -ss*a1 + cc*a2
	}
}

// PassFieldSum fills fsum from slot geometry (Lo radii already in slots).
func PassFieldSum(lo []CellLo, slots []Slot, cls, inc []int32, law Law, fsum []float64, a, b int) {
	_ = lo
	for i := a; i < b; i++ {
		sacc := 0.0
		for q := cls[i]; q < cls[i+1]; q++ {
			sl := &slots[inc[q]]
			if sl.A <= 0 {
				continue
			}
			sacc += (sl.A / law.Aref) * (law.Dref / sl.D)
		}
		fsum[i] = sacc
	}
}

// PassFieldHopWeight caches per-slot hop angles.
func PassFieldHopWeight(slots []Slot, fsum []float64, law Law, a, b int) {
	for s := a; s < b; s++ {
		sl := &slots[s]
		if sl.Alive == 0 || sl.A <= 0 {
			sl.WField = 0
			continue
		}
		i, j := int(sl.I), int(sl.J)
		si, sj := fsum[i], fsum[j]
		if si <= 1e-12 || sj <= 1e-12 {
			sl.WField = 0
			continue
		}
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sl.WField = law.FieldJ * w / math.Sqrt(si*sj) * law.Dt
	}
}

// ApplyFieldHopColor applies Givens on Hi.Fa*Next for one color class.
func ApplyFieldHopColor(hi []CellHi, slots []Slot, color []int, want, a, b int) {
	for s := a; s < b; s++ {
		if color[s] != want {
			continue
		}
		applyHopOne(hi, slots, s)
	}
}

// ApplyHopOne is the single-slot hop (used by color buckets).
func ApplyHopOne(hi []CellHi, slots []Slot, s int) {
	applyHopOne(hi, slots, s)
}

func applyHopOne(hi []CellHi, slots []Slot, s int) {
	sl := &slots[s]
	if sl.Alive == 0 || sl.WField == 0 {
		return
	}
	tau := sl.WField
	if tau > 0.5 {
		tau = 0.5
	}
	cc, ss := math.Cos(tau), math.Sin(tau)
	i, j := sl.I, sl.J
	a1i, a2i := hi[i].Fa1Next, hi[i].Fa2Next
	a1j, a2j := hi[j].Fa1Next, hi[j].Fa2Next
	hi[i].Fa1Next = cc*a1i + ss*a2j
	hi[i].Fa2Next = cc*a2i - ss*a1j
	hi[j].Fa1Next = cc*a1j + ss*a2i
	hi[j].Fa2Next = cc*a2j - ss*a1i
}

// PassCommitField: Fa <- FaNext, sync Ee (Lo+Hi).
func PassCommitField(lo []CellLo, hi []CellHi, a, b int) {
	for i := a; i < b; i++ {
		lo[i].Fa1 = hi[i].Fa1Next
		lo[i].Fa2 = hi[i].Fa2Next
		lo[i].SyncEe()
	}
}

// PassClock advances dense clocks (Lo only).
func PassClock(lo []CellLo, law Law, a, b int) {
	for i := a; i < b; i++ {
		c := &lo[i]
		c.Th2 = math.Mod(c.Th2+c.W2e*law.Dt, 2*math.Pi)
		c.Cbeta += (c.W1e - c.W2e) * law.Dt
		if c.Cbeta >= 2*math.Pi {
			c.Cbeta -= 2 * math.Pi
		} else if c.Cbeta <= -2*math.Pi {
			c.Cbeta += 2 * math.Pi
		}
	}
}

// PassForceGatherDisjoint accumulates into Hi.F* for one color.
func PassForceGatherDisjoint(lo []CellLo, hi []CellHi, slots []Slot, color []int, want int, law Law, a, b int) {
	for s := a; s < b; s++ {
		if color[s] != want {
			continue
		}
		ForceOne(lo, hi, slots, s, law)
	}
}

// ForceOne applies soft repulsion for one slot into Hi forces.
func ForceOne(lo []CellLo, hi []CellHi, slots []Slot, s int, law Law) {
	sl := &slots[s]
	if sl.Alive == 0 {
		return
	}
	i, j := sl.I, sl.J
	d0 := lo[i].Cr + lo[j].Cr
	if sl.D >= d0 || sl.D <= 1e-12 {
		return
	}
	mag := law.KRep * (d0 - sl.D)
	fx, fy, fz := mag*sl.Ux, mag*sl.Uy, mag*sl.Uz
	hi[i].Fx -= fx
	hi[i].Fy -= fy
	hi[i].Fz -= fz
	hi[j].Fx += fx
	hi[j].Fy += fy
	hi[j].Fz += fz
}

// PassMotion integrates positions from Hi forces into Lo, clears forces.
func PassMotion(lo []CellLo, hi []CellHi, law Law, a, b int) {
	dt := law.Dt * law.MobGeo
	L := law.L
	for i := a; i < b; i++ {
		c := &lo[i]
		h := &hi[i]
		c.X = fold(c.X+dt*h.Fx, L)
		c.Y = fold(c.Y+dt*h.Fy, L)
		c.Z = fold(c.Z+dt*h.Fz, L)
		h.Fx, h.Fy, h.Fz = 0, 0, 0
	}
}

func fold(x, L float64) float64 {
	x = math.Mod(x, L)
	if x < 0 {
		x += L
	}
	return x
}

// Step is the serial two-phase freecell-shaped step.
func (w *World) Step() {
	law := w.Law
	nc := len(w.Lo)
	ns := len(w.Slots)
	if nc == 0 {
		return
	}

	if len(w.sscl) < nc {
		w.sscl = make([]float64, nc)
		w.fsum = make([]float64, nc)
	}

	PassPitch(w.Lo, law, 0, nc)
	w.refreshGeom()

	if law.SK > 0 {
		PassSpaceWants(w.Lo, w.Slots, law, 0, ns)
		PassSpaceScale(w.Lo, w.Slots, w.Cls, w.Inc, law, w.sscl, 0, nc)
		PassSpaceApply(w.Lo, w.Hi, w.Slots, w.Cls, w.Inc, w.sscl, 0, nc)
		PassCommitSpace(w.Lo, w.Hi, 0, nc)
		PassPitch(w.Lo, law, 0, nc)
		w.refreshGeom()
	}

	PassFieldPrecess(w.Lo, w.Hi, law, 0, nc)
	PassFieldSum(w.Lo, w.Slots, w.Cls, w.Inc, law, w.fsum, 0, nc)
	PassFieldHopWeight(w.Slots, w.fsum, law, 0, ns)
	if w.Color == nil {
		w.Color, w.NColor = w.EdgeColor()
	}
	for c := 0; c < w.NColor; c++ {
		ApplyFieldHopColor(w.Hi, w.Slots, w.Color, c, 0, ns)
	}
	PassCommitField(w.Lo, w.Hi, 0, nc)

	for i := 0; i < nc; i++ {
		w.Hi[i].Fx, w.Hi[i].Fy, w.Hi[i].Fz = 0, 0, 0
	}
	for c := 0; c < w.NColor; c++ {
		PassForceGatherDisjoint(w.Lo, w.Hi, w.Slots, w.Color, c, law, 0, ns)
	}
	PassMotion(w.Lo, w.Hi, law, 0, nc)
	PassClock(w.Lo, law, 0, nc)
}
