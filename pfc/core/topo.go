package core

import "math"

// World is a serial free-cell roster: two contiguous slabs Lo/Hi + slots.
// Topology is frozen after Build* (experiment v1).
type World struct {
	Law  Law
	Lo   []CellLo
	Hi   []CellHi
	Slots []Slot
	Cls  []int32
	Inc  []int32
	Color  []int
	NColor int
	// Side is the cubic lattice edge length when built by BuildLattice
	// (0 if unknown). Used by par for spatial z-slab ownership.
	Side int

	// Reused scratch for Step (zero alloc after first Step or Build).
	sscl []float64
	fsum []float64
}

// NC returns cell count (len(Lo) == len(Hi)).
func (w *World) NC() int { return len(w.Lo) }

// BuildLattice places cells on a cubic lattice (6-connected).
func BuildLattice(law Law, nApprox int) *World {
	if nApprox < 8 {
		nApprox = 8
	}
	side := int(math.Cbrt(float64(nApprox)) + 0.5)
	if side < 2 {
		side = 2
	}
	nc := side * side * side
	dx := law.L / float64(side)
	w := &World{
		Law:  law,
		Lo:   make([]CellLo, nc),
		Hi:   make([]CellHi, nc),
		Side: side,
	}
	idx := func(ix, iy, iz int) int32 {
		return int32((iz*side+iy)*side + ix)
	}
	for iz := 0; iz < side; iz++ {
		for iy := 0; iy < side; iy++ {
			for ix := 0; ix < side; ix++ {
				i := int(idx(ix, iy, iz))
				lo := &w.Lo[i]
				hi := &w.Hi[i]
				lo.X = (float64(ix) + 0.5) * dx
				lo.Y = (float64(iy) + 0.5) * dx
				lo.Z = (float64(iz) + 0.5) * dx
				lo.Cr0 = law.R0
				lo.Cr = law.R0
				lo.Es = law.Es0
				lo.Em = 0
				lo.W1e = law.W1
				lo.W2e = law.W2
				hi.N1x, hi.N1y, hi.N1z = 1, 0, 0
				hi.N2x, hi.N2y, hi.N2z = 0, 1, 0
				lo.Fa1 = 0.05 * math.Sin(0.3*lo.X+0.1*lo.Y)
				lo.Fa2 = 0.05 * math.Cos(0.2*lo.Z+0.1*lo.X)
				lo.SyncEe()
				lo.Th2 = 0.1 * lo.X
			}
		}
	}
	slots := make([]Slot, 0, 3*nc)
	add := func(a, b int32) {
		if a > b {
			a, b = b, a
		}
		slots = append(slots, Slot{I: a, J: b, Alive: 1})
	}
	for iz := 0; iz < side; iz++ {
		for iy := 0; iy < side; iy++ {
			for ix := 0; ix < side; ix++ {
				a := idx(ix, iy, iz)
				add(a, idx((ix+1)%side, iy, iz))
				add(a, idx(ix, (iy+1)%side, iz))
				add(a, idx(ix, iy, (iz+1)%side))
			}
		}
	}
	type pair struct{ i, j int32 }
	seen := make(map[pair]struct{}, len(slots))
	out := slots[:0]
	for _, s := range slots {
		p := pair{s.I, s.J}
		if _, ok := seen[p]; ok {
			continue
		}
		seen[p] = struct{}{}
		out = append(out, s)
	}
	w.Slots = out
	w.rebuildCSR()
	w.Color, w.NColor = w.EdgeColor()
	w.sscl = make([]float64, nc)
	w.fsum = make([]float64, nc)
	w.refreshGeom()
	return w
}

func (w *World) rebuildCSR() {
	nc := len(w.Lo)
	w.Cls = make([]int32, nc+1)
	for i := range w.Slots {
		s := &w.Slots[i]
		if s.Alive == 0 {
			continue
		}
		w.Cls[s.I+1]++
		w.Cls[s.J+1]++
	}
	for i := 0; i < nc; i++ {
		w.Cls[i+1] += w.Cls[i]
	}
	w.Inc = make([]int32, w.Cls[nc])
	fill := make([]int32, nc)
	copy(fill, w.Cls[:nc])
	for si := range w.Slots {
		s := &w.Slots[si]
		if s.Alive == 0 {
			continue
		}
		w.Inc[fill[s.I]] = int32(si)
		fill[s.I]++
		w.Inc[fill[s.J]] = int32(si)
		fill[s.J]++
	}
}

func (w *World) wrap(d float64) float64 {
	L := w.Law.L
	if d > 0.5*L {
		d -= L
	}
	if d < -0.5*L {
		d += L
	}
	return d
}

func (w *World) refreshGeom() {
	for i := range w.Slots {
		s := &w.Slots[i]
		if s.Alive == 0 {
			s.A = 0
			continue
		}
		ci := &w.Lo[s.I]
		cj := &w.Lo[s.J]
		dx := w.wrap(cj.X - ci.X)
		dy := w.wrap(cj.Y - ci.Y)
		dz := w.wrap(cj.Z - ci.Z)
		d := math.Sqrt(dx*dx + dy*dy + dz*dz)
		if d < 1e-9 {
			d = 1e-9
		}
		s.D = d
		s.Ux, s.Uy, s.Uz = dx/d, dy/d, dz/d
		s.A = lensArea(d, ci.Cr, cj.Cr)
	}
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
	if d < math.Abs(ri-rj) {
		rm := ri
		if rj < ri {
			rm = rj
		}
		return math.Pi * rm * rm
	}
	return 0
}

// TotalE returns Em+Ee+Es.
func (w *World) TotalE() float64 {
	var e float64
	for i := range w.Lo {
		c := &w.Lo[i]
		e += c.Em + c.Ee + c.Es
	}
	return e
}

// EdgeColor assigns colors so no two slots sharing a cell share a color.
// Uses a flat used[cell*stride+color] table — not nc separate []bool slices
// (those 1.7M+ tiny allocs dominate BuildLattice at large NC).
func (w *World) EdgeColor() (color []int, nColor int) {
	nc := len(w.Lo)
	deg := make([]int, nc)
	for i := range w.Slots {
		s := &w.Slots[i]
		if s.Alive == 0 {
			continue
		}
		deg[s.I]++
		deg[s.J]++
	}
	maxDeg := 0
	for _, d := range deg {
		if d > maxDeg {
			maxDeg = d
		}
	}
	nColor = maxDeg + 1
	if nColor < 1 {
		nColor = 1
	}
	// Lattice degree ≤ 6; keep a little headroom without re-alloc thrash.
	if nColor < 8 {
		nColor = 8
	}
	color = make([]int, len(w.Slots))
	for i := range color {
		color[i] = -1
	}
	stride := nColor
	used := make([]bool, nc*stride)
	for si := range w.Slots {
		s := &w.Slots[si]
		if s.Alive == 0 {
			continue
		}
		i, j := int(s.I), int(s.J)
		c := 0
		for c < nColor && (used[i*stride+c] || used[j*stride+c]) {
			c++
		}
		if c >= nColor {
			// Grow flat table (rare on lattices).
			newN := nColor + 4
			nu := make([]bool, nc*newN)
			for cell := 0; cell < nc; cell++ {
				copy(nu[cell*newN:cell*newN+nColor], used[cell*stride:cell*stride+nColor])
			}
			used = nu
			nColor = newN
			stride = nColor
		}
		color[si] = c
		used[i*stride+c] = true
		used[j*stride+c] = true
	}
	// Trim reported nColor to max color id + 1.
	maxC := -1
	for _, c := range color {
		if c > maxC {
			maxC = c
		}
	}
	if maxC+1 > 0 {
		nColor = maxC + 1
	}
	return color, nColor
}
