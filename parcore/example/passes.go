package example

import (
	"math"

	"scp/parcore"
)

// PassPitch updates the live radius and load-detuned field frequency (owned
// cells only). Lo-only in pfc terms.
func PassPitch(law Law, c []Cell, a, b int) {
	for i := a; i < b; i++ {
		ci := &c[i]
		ratio := ci.Es / law.Es0
		if ratio < 0 {
			ratio = 0
		}
		ci.Cr = law.R0 * math.Cbrt(ratio)
		x := (ci.Em + ci.Ee) / law.Cap
		det := 1.0 + law.QDetune*x
		ci.W1e = law.W1 / det
	}
}

// PassGeom recomputes slot geometry (D, U, A) from current positions/radii and
// rebuilds the live-slot list. live is reset ([:0]) and repopulated. Mark a
// slot live when its lens area A > 0.
func PassGeom(law Law, c []Cell, slots []Slot, si, sj []int32, live *[]int32) {
	*live = (*live)[:0]
	for s := range slots {
		sl := &slots[s]
		ci := &c[si[s]]
		cj := &c[sj[s]]
		dx := parcore.Wrap(cj.X-ci.X, law.L)
		dy := parcore.Wrap(cj.Y-ci.Y, law.L)
		dz := parcore.Wrap(cj.Z-ci.Z, law.L)
		d2 := dx*dx + dy*dy + dz*dz
		if d2 < 1e-18 {
			d2 = 1e-18
		}
		rsum := ci.Cr + cj.Cr
		if d2 >= rsum*rsum {
			// no overlap: record distance/unit but A=0, not live
			d := math.Sqrt(d2)
			sl.D = d
			inv := 1.0 / d
			sl.Ux, sl.Uy, sl.Uz = dx*inv, dy*inv, dz*inv
			sl.A = 0
			continue
		}
		d := math.Sqrt(d2)
		sl.D = d
		inv := 1.0 / d
		sl.Ux, sl.Uy, sl.Uz = dx*inv, dy*inv, dz*inv
		sl.A = parcore.LensArea(d, ci.Cr, cj.Cr)
		if sl.A > 0 {
			*live = append(*live, int32(s))
		}
	}
}

// PassSpaceWants writes Swl from the Es pressure difference over live slots.
func PassSpaceWants(law Law, c []Cell, slots []Slot, si, sj []int32, live []int32) {
	dt := law.Dt
	for _, si32 := range live {
		s := int(si32)
		sl := &slots[s]
		if sl.A <= 0 {
			sl.Swl = 0
			continue
		}
		i, j := si[s], sj[s]
		pi := c[i].Es + law.SDisp*(c[i].Em+c[i].Ee)
		pj := c[j].Es + law.SDisp*(c[j].Em+c[j].Ee)
		dp := pi - pj
		if dp == 0 {
			sl.Swl = 0
			continue
		}
		w := (sl.A / law.Aref) * (law.Dref / sl.D)
		sl.Swl = law.SK * dt * w * dp
	}
}

// PassSpaceScale computes the per-cell outflow scale Sscl (clamped so Es stays
// above EsFloor). Operates over owned cells [0:nOwn). A slot's flow has its
// SOURCE at the higher-pressure end (si when f>0, sj when f<0); that source's
// outflow request accumulates |f|.
func PassSpaceScale(law Law, c []Cell, slots []Slot, si, sj []int32, cls, inc []int32, nOwn int) {
	for i := 0; i < nOwn; i++ {
		rq := 0.0
		for q := cls[i]; q < cls[i+1]; q++ {
			s := int(inc[q])
			sl := &slots[s]
			f := sl.Swl
			if f == 0 {
				continue
			}
			src := si[s]
			if f < 0 {
				src = sj[s]
			}
			if int(src) == i {
				rq += math.Abs(f)
			}
		}
		avail := c[i].Es - law.EsFloor
		if avail <= 0 {
			c[i].Sscl = 0
		} else if rq > avail {
			c[i].Sscl = avail / rq
		} else {
			c[i].Sscl = 1
		}
	}
}

// PassSpaceApply moves the scaled wants: source debits mag, dest credits mag.
// Reads ghost Sscl (refreshed by the mid-step halo exchange). Conserves Es
// exactly (each slot moves |f|*Sscl[src] from src to dst).
func PassSpaceApply(law Law, c []Cell, slots []Slot, si, sj []int32, cls, inc []int32, nOwn int) {
	_ = law
	for i := 0; i < nOwn; i++ {
		de := 0.0
		for q := cls[i]; q < cls[i+1]; q++ {
			s := int(inc[q])
			sl := &slots[s]
			f := sl.Swl
			if f == 0 {
				continue
			}
			src := si[s]
			if f < 0 {
				src = sj[s]
			}
			mag := math.Abs(f) * c[src].Sscl
			if int(src) == i {
				de -= mag
			} else {
				de += mag
			}
		}
		c[i].Es += de
	}
}

// PassFieldPrecess rotates the field amplitude by W1e*dt (owned cells, in
// place). Per-cell local; no cross-cell dependency.
func PassFieldPrecess(law Law, c []Cell, a, b int) {
	for i := a; i < b; i++ {
		ci := &c[i]
		ang := ci.W1e * law.Dt
		ss, cc := math.Sincos(ang)
		a1, a2 := ci.Fa1, ci.Fa2
		ci.Fa1 = cc*a1 + ss*a2
		ci.Fa2 = -ss*a1 + cc*a2
	}
}

// PassFieldSum accumulates the geometric weight sum per cell (fsum) over live
// slots. fsum must be sized to len(c) (owned+ghost) and is zeroed here.
func PassFieldSum(slots []Slot, si, sj []int32, fsum []float64, live []int32) {
	for i := range fsum {
		fsum[i] = 0
	}
	for _, si32 := range live {
		s := int(si32)
		sl := &slots[s]
		if sl.A <= 0 {
			continue
		}
		w := sl.A / sl.D
		fsum[si[s]] += w
		fsum[sj[s]] += w
	}
}

// PassHopWeight writes the per-slot hop angle into param (the core's
// SlotParam). tau = FieldJ * w / sqrt(fi*fj) * dt, clamped.
func PassHopWeight(law Law, slots []Slot, si, sj []int32, fsum []float64, param []float64, live []int32) {
	for _, si32 := range live {
		s := int(si32)
		sl := &slots[s]
		if sl.A <= 0 {
			param[s] = 0
			continue
		}
		fi, fj := fsum[si[s]], fsum[sj[s]]
		if fi <= 1e-12 || fj <= 1e-12 {
			param[s] = 0
			continue
		}
		w := sl.A / sl.D
		tau := law.FieldJ * w / math.Sqrt(fi*fj) * law.Dt
		if tau > 0.5 {
			tau = 0.5
		}
		param[s] = tau
	}
}

// ApplyHop applies one unitary Givens rotation on the field amplitude pair of
// slot s's endpoints. Uses parcore.Givens2 so Fa1²+Fa2² is conserved over the
// pair exactly. For interior colored rounds (disjoint endpoints within a
// color) and for the boundary round over the publish buffer.
func ApplyHop(c []Cell, si, sj []int32, param []float64, s int) {
	tau := param[s]
	if tau <= 0 {
		return
	}
	i, j := int(si[s]), int(sj[s])
	n1i, n2i, n1j, n2j := parcore.Givens2(c[i].Fa1, c[i].Fa2, c[j].Fa1, c[j].Fa2, tau)
	c[i].Fa1, c[i].Fa2 = n1i, n2i
	c[j].Fa1, c[j].Fa2 = n1j, n2j
}

// ApplyHopPub is the publish-buffer variant for cross-shard boundary hops.
func ApplyHopPub(pub []Cell, gi, gj int, tau float64) {
	if tau <= 0 {
		return
	}
	if tau > 0.5 {
		tau = 0.5
	}
	n1i, n2i, n1j, n2j := parcore.Givens2(pub[gi].Fa1, pub[gi].Fa2, pub[gj].Fa1, pub[gj].Fa2, tau)
	pub[gi].Fa1, pub[gi].Fa2 = n1i, n2i
	pub[gj].Fa1, pub[gj].Fa2 = n1j, n2j
}

// PassCommit syncs Ee = Fa1²+Fa2² (owned cells).
func PassCommit(c []Cell, a, b int) {
	for i := a; i < b; i++ {
		c[i].Ee = c[i].Fa1*c[i].Fa1 + c[i].Fa2*c[i].Fa2
	}
}
