package par

import (
	"math"

	"scp/pfc/core"
)

// Shard is one worker's permanent memory. Lo and Hi are separately
// allocated contiguous slabs that stay on this shard for the cluster
// lifetime. Layout:
//
//	Lo/Hi [0 : NOwn)           — cells this worker owns (global G0..G0+NOwn-1)
//	Lo/Hi [NOwn : NOwn+NGhost) — ghost copies of remote endpoints
//
// Slot I,J are indices into these local slabs (not global ids).
// Worker id == shard index; that mapping is fixed for the Cluster life.
type Shard struct {
	ID       int
	G0, NOwn int
	NGhost   int

	Lo []core.CellLo // exclusive make() for this worker
	Hi []core.CellHi // exclusive make() for this worker

	Slots []core.Slot
	Cls   []int32
	Inc   []int32

	// ghosts
	GhostGID []int32 // global id of local index NOwn+k
	g2l      map[int32]int32

	// interior vs boundary slot lists (local slot indices)
	Interior []int32
	Boundary []int32

	// owned local indices that touch at least one boundary edge
	// (only these need FaNext published/merged for boundary hops)
	BndTouch []int32

	// Publish: owned local indices that some other shard ghosts.
	// Face plane for z-slabs — much smaller than NOwn.
	Publish []int32

	// live slots: bitset + dense list (A>0 after last geom)
	LiveBits []uint64
	Live     []int32

	// active owned cells (endpoints of live slots) — dense list
	Active []int32
	actBits []uint64 // NOwn bits

	// per-shard scratch (sized to len(Lo))
	sscl []float64
	fsum []float64

	// local edge colors for interior slots only
	color   []int
	nColor  int
	buckets [][]int32 // color -> interior slot indices
}

func (sh *Shard) rebuildCSR() {
	n := sh.NOwn + sh.NGhost
	sh.Cls = make([]int32, n+1)
	for i := range sh.Slots {
		s := &sh.Slots[i]
		if s.Alive == 0 {
			continue
		}
		sh.Cls[s.I+1]++
		sh.Cls[s.J+1]++
	}
	for i := 0; i < n; i++ {
		sh.Cls[i+1] += sh.Cls[i]
	}
	sh.Inc = make([]int32, sh.Cls[n])
	fill := make([]int32, n)
	copy(fill, sh.Cls[:n])
	for si := range sh.Slots {
		s := &sh.Slots[si]
		if s.Alive == 0 {
			continue
		}
		sh.Inc[fill[s.I]] = int32(si)
		fill[s.I]++
		sh.Inc[fill[s.J]] = int32(si)
		fill[s.J]++
	}
}

func (sh *Shard) localToGlobal(li int) int {
	if li < sh.NOwn {
		return sh.G0 + li
	}
	return int(sh.GhostGID[li-sh.NOwn])
}

// ensureLiveCap resets Live/Active without allocating when pre-sized in FromWorld.
func (sh *Shard) ensureLiveCap(nSlot int) {
	nb := (nSlot + 63) / 64
	if len(sh.LiveBits) < nb {
		// first use only if FromWorld skipped (should not happen)
		sh.LiveBits = make([]uint64, nb)
		sh.Live = make([]int32, 0, nSlot)
	} else {
		sh.LiveBits = sh.LiveBits[:nb]
		clear(sh.LiveBits)
	}
	sh.Live = sh.Live[:0]
	na := (sh.NOwn + 63) / 64
	if len(sh.actBits) < na {
		sh.actBits = make([]uint64, na)
		sh.Active = make([]int32, 0, sh.NOwn)
	} else {
		sh.actBits = sh.actBits[:na]
		clear(sh.actBits)
	}
	sh.Active = sh.Active[:0]
}

func (sh *Shard) markLive(s int) {
	sh.LiveBits[s>>6] |= 1 << (uint(s) & 63)
	sh.Live = append(sh.Live, int32(s))
	// active owned endpoints
	sl := &sh.Slots[s]
	if int(sl.I) < sh.NOwn {
		i := int(sl.I)
		if sh.actBits[i>>6]&(1<<(uint(i)&63)) == 0 {
			sh.actBits[i>>6] |= 1 << (uint(i) & 63)
			sh.Active = append(sh.Active, sl.I)
		}
	}
	if int(sl.J) < sh.NOwn {
		j := int(sl.J)
		if sh.actBits[j>>6]&(1<<(uint(j)&63)) == 0 {
			sh.actBits[j>>6] |= 1 << (uint(j) & 63)
			sh.Active = append(sh.Active, sl.J)
		}
	}
}

func (sh *Shard) refreshGeom(L float64) {
	sh.ensureLiveCap(len(sh.Slots))
	for s := range sh.Slots {
		sl := &sh.Slots[s]
		if sl.Alive == 0 {
			sl.A = 0
			continue
		}
		ci := &sh.Lo[sl.I]
		cj := &sh.Lo[sl.J]
		dx := wrap(cj.X-ci.X, L)
		dy := wrap(cj.Y-ci.Y, L)
		dz := wrap(cj.Z-ci.Z, L)
		d2 := dx*dx + dy*dy + dz*dz
		if d2 < 1e-18 {
			d2 = 1e-18
		}
		// cheap reject before sqrt: no overlap possible
		rsum := ci.Cr + cj.Cr
		if d2 >= rsum*rsum {
			sl.D = math.Sqrt(d2)
			inv := 1.0 / sl.D
			sl.Ux, sl.Uy, sl.Uz = dx*inv, dy*inv, dz*inv
			sl.A = 0
			continue
		}
		d := math.Sqrt(d2)
		sl.D = d
		inv := 1.0 / d
		sl.Ux, sl.Uy, sl.Uz = dx*inv, dy*inv, dz*inv
		sl.A = lensArea(d, ci.Cr, cj.Cr)
		if sl.A > 0 {
			sh.markLive(s)
		}
	}
}

// refreshAreaOnly recomputes A from existing d and current Cr (no √).
func (sh *Shard) refreshAreaOnly() {
	sh.ensureLiveCap(len(sh.Slots))
	for s := range sh.Slots {
		sl := &sh.Slots[s]
		if sl.Alive == 0 {
			sl.A = 0
			continue
		}
		// reject without full lens formula
		rsum := sh.Lo[sl.I].Cr + sh.Lo[sl.J].Cr
		if sl.D >= rsum {
			sl.A = 0
			continue
		}
		sl.A = lensArea(sl.D, sh.Lo[sl.I].Cr, sh.Lo[sl.J].Cr)
		if sl.A > 0 {
			sh.markLive(s)
		}
	}
}

func (sh *Shard) slotLive(s int) bool {
	return sh.LiveBits[s>>6]&(1<<(uint(s)&63)) != 0
}

// pullGhostsFromPub copies published Lo for each ghost gid.
func (sh *Shard) pullGhostsFromPub(pub []core.CellLo) {
	for k, gid := range sh.GhostGID {
		sh.Lo[sh.NOwn+k] = pub[gid]
	}
}

// publishFace writes Publish-set cells into pub (face plane, not all owned).
func (sh *Shard) publishFace(pub []core.CellLo) {
	for _, li := range sh.Publish {
		i := int(li)
		pub[sh.G0+i] = sh.Lo[i]
	}
}

// publishFaceCr writes only Cr for the publish set into pubS (global id index).
func (sh *Shard) publishFaceCr(pubS []float64) {
	for _, li := range sh.Publish {
		i := int(li)
		pubS[sh.G0+i] = sh.Lo[i].Cr
	}
}

// pullGhostCr loads Cr for ghosts from pubS.
func (sh *Shard) pullGhostCr(pubS []float64) {
	for k, gid := range sh.GhostGID {
		sh.Lo[sh.NOwn+k].Cr = pubS[gid]
	}
}

func wrap(d, L float64) float64 {
	if d > 0.5*L {
		d -= L
	}
	if d < -0.5*L {
		d += L
	}
	return d
}

func lensArea(d, ri, rj float64) float64 {
	if d >= ri+rj {
		return 0
	}
	t := d*d - rj*rj + ri*ri
	a2 := (4*d*d*ri*ri - t*t) / (4 * d * d)
	if a2 > 0 {
		return math.Pi * a2
	}
	rm := ri
	if rj < ri {
		rm = rj
	}
	if d < math.Abs(ri-rj) {
		return math.Pi * rm * rm
	}
	return 0
}

// buildInteriorColors edge-colors only interior slots (both ends owned).
func (sh *Shard) buildInteriorColors() {
	// temporary World-like color on a filtered slot list is awkward;
	// greedy color on interior only using local cell indices 0..NOwn.
	nV := sh.NOwn
	nColor := 1
	sh.color = make([]int, len(sh.Slots))
	for i := range sh.color {
		sh.color[i] = -1
	}
	used := make([][]bool, nV)
	for i := range used {
		used[i] = make([]bool, 8)
	}
	for _, si := range sh.Interior {
		s := &sh.Slots[si]
		i, j := int(s.I), int(s.J)
		c := 0
		for {
			if c >= len(used[0]) {
				for k := range used {
					used[k] = append(used[k], false)
				}
			}
			if c >= nColor {
				nColor = c + 1
			}
			if !used[i][c] && !used[j][c] {
				break
			}
			c++
		}
		sh.color[si] = c
		used[i][c] = true
		used[j][c] = true
		if c+1 > nColor {
			nColor = c + 1
		}
	}
	sh.nColor = nColor
	sh.buckets = make([][]int32, nColor)
	for _, si := range sh.Interior {
		c := sh.color[si]
		if c >= 0 {
			sh.buckets[c] = append(sh.buckets[c], si)
		}
	}
}
