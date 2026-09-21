package parcore

import "math"

// Topology is the seam between frozen (v1) and dynamic roster graphs.
//
// v1 ships only FrozenTopology: FromWorld builds the CSR + edge coloring once
// from the World's endpoints, and Cluster.Step never mutates the graph. The
// free-cell programme requires that cells be born and die (no background), so
// the seam is named here explicitly: a future DynamicTopology would implement
// Refresh to do (a) spatial re-neighboring from current positions, (b)
// birth/death of cells into a freelist, (c) open-addressing pair→slot hash
// rebuild, (d) canonical CSR rebuild, (e) re-edge-coloring — all per shard,
// then a cross-shard ghost-set reconciliation at a barrier.
//
// The reason v1 freezes the graph is exactly pfc/README's pushback item 1:
// dynamic topo_refresh across shards (freelist + pair hash + cross-proc
// birth/death + canonical CSR rebuild) is the hard parallel-systems problem
// and is deferred on purpose. The FrozenTopology contract is that the roster
// is constant; geometry (positions, lens areas, hop weights) may still vary
// every step, which is what the driver's geom pass recomputes.
type Topology interface {
	// Dynamic reports whether Refresh can change the roster. FrozenTopology
	// always returns false; the core uses this to skip refresh work.
	Dynamic() bool
}

// frozen is the v1 topology policy: built once, never refreshed.
type frozen struct{}

func (frozen) Dynamic() bool { return false }

// BuildLattice places cells on a cubic lattice (6-connected, periodic) and
// returns a World ready for FromWorld. The cellSeed callback fills the driver-
// specific cell fields from the (ix,iy,iz) index and the lattice spacing dx;
// the slot values are zero-valued (the driver's geom pass fills A/D/U later).
//
// R0 is the rest radius (lens area is positive on a unit lattice when
// R0 > 0.5). Box is the periodic box length L.
func BuildLattice[C any, S any](
	ncApprox int, box, r0 float64,
	cellSeed func(ix, iy, iz int, dx float64) C,
) *World[C, S] {
	if ncApprox < 8 {
		ncApprox = 8
	}
	side := int(math.Cbrt(float64(ncApprox)) + 0.5)
	if side < 2 {
		side = 2
	}
	nc := side * side * side
	dx := box / float64(side)
	w := &World[C, S]{
		Cells: make([]C, nc),
		Slots: make([]S, 0, 3*nc),
		SI:    make([]int32, 0, 3*nc),
		SJ:    make([]int32, 0, 3*nc),
		Alive: make([]bool, 0, 3*nc),
		Side:  side,
		Box:   box,
	}
	idx := func(ix, iy, iz int) int32 {
		return int32((iz*side+iy)*side + ix)
	}
	for iz := 0; iz < side; iz++ {
		for iy := 0; iy < side; iy++ {
			for ix := 0; ix < side; ix++ {
				w.Cells[int(idx(ix, iy, iz))] = cellSeed(ix, iy, iz, dx)
			}
		}
	}
	// 6-connected periodic lattice, undirected edges (canonical i<j).
	add := func(a, b int32) {
		if a > b {
			a, b = b, a
		}
		w.SI = append(w.SI, a)
		w.SJ = append(w.SJ, b)
		w.Alive = append(w.Alive, true)
		var zero S
		w.Slots = append(w.Slots, zero)
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
	// dedup edges (a<<32|b keys) — lattice edges are already unique, but be safe.
	dedupEdges(w)
	return w
}

func dedupEdges[C any, S any](w *World[C, S]) {
	type pair struct{ a, b int32 }
	seen := make(map[pair]struct{}, len(w.SI))
	si, sj, al, sl := w.SI[:0], w.SJ[:0], w.Alive[:0], w.Slots[:0]
	for k := range w.SI {
		p := pair{w.SI[k], w.SJ[k]}
		if _, ok := seen[p]; ok {
			continue
		}
		seen[p] = struct{}{}
		si = append(si, w.SI[k])
		sj = append(sj, w.SJ[k])
		al = append(al, w.Alive[k])
		sl = append(sl, w.Slots[k])
	}
	w.SI, w.SJ, w.Alive, w.Slots = si, sj, al, sl
}

// buildCSR constructs the cell→incident-slot CSR (Cls/Inc) over local slots.
// Endpoints in si/sj are LOCAL cell indices in [0,n). Returns (cls, inc).
func buildCSR(si, sj []int32, alive []bool, n int) (cls, inc []int32) {
	return BuildCSRPublic(si, sj, alive, n)
}

// BuildCSRPublic is the exported CSR builder (drivers use it for serial
// baselines and custom topologies).
func BuildCSRPublic(si, sj []int32, alive []bool, n int) (cls, inc []int32) {
	cls = make([]int32, n+1)
	for s := range si {
		if !alive[s] {
			continue
		}
		cls[si[s]+1]++
		cls[sj[s]+1]++
	}
	for i := 0; i < n; i++ {
		cls[i+1] += cls[i]
	}
	inc = make([]int32, cls[n])
	fill := make([]int32, n)
	copy(fill, cls[:n])
	for s := range si {
		if !alive[s] {
			continue
		}
		inc[fill[si[s]]] = int32(s)
		fill[si[s]]++
		inc[fill[sj[s]]] = int32(s)
		fill[sj[s]]++
	}
	return cls, inc
}

// edgeColor colors slots so no two slots sharing a cell share a color. Used
// ONLY over interior slots (both endpoints owned) per shard, so each color's
// slot set is endpoint-disjoint and a dual-endpoint write is race-free.
//
// Flat used-table (nc*stride) — NOT one []bool per cell (the per-cell make was
// pfc's catastrophic-allocation bug at large NC; see pfc/README "Allocations").
func edgeColor(si, sj []int32, alive []bool, interior []int32, nCell int) (buckets [][]int32) {
	return EdgeColorPublic(si, sj, alive, interior, nCell)
}

// EdgeColorPublic is the exported edge-colorer (drivers use it for serial
// baselines). Returns color buckets; each bucket's slots are endpoint-disjoint.
func EdgeColorPublic(si, sj []int32, alive []bool, interior []int32, nCell int) (buckets [][]int32) {
	deg := make([]int, nCell)
	for _, s := range interior {
		deg[si[s]]++
		deg[sj[s]]++
	}
	maxDeg := 0
	for _, d := range deg {
		if d > maxDeg {
			maxDeg = d
		}
	}
	nColor := maxDeg + 1
	if nColor < 8 {
		nColor = 8 // lattice headroom without re-alloc thrash
	}
	color := make([]int, len(si))
	for i := range color {
		color[i] = -1
	}
	stride := nColor
	used := make([]bool, nCell*stride)
	for _, s := range interior {
		i, j := int(si[s]), int(sj[s])
		c := 0
		for c < nColor && (used[i*stride+c] || used[j*stride+c]) {
			c++
		}
		if c >= nColor {
			// grow (rare on lattices)
			newN := nColor + 4
			nu := make([]bool, nCell*newN)
			for cell := 0; cell < nCell; cell++ {
				copy(nu[cell*newN:cell*newN+nColor], used[cell*stride:cell*stride+nColor])
			}
			used = nu
			nColor = newN
			stride = nColor
		}
		color[s] = c
		used[i*stride+c] = true
		used[j*stride+c] = true
	}
	maxC := -1
	for _, s := range interior {
		if color[s] > maxC {
			maxC = color[s]
		}
	}
	if maxC < 0 {
		return nil
	}
	buckets = make([][]int32, maxC+1)
	for _, s := range interior {
		buckets[color[s]] = append(buckets[color[s]], s)
	}
	return buckets
}

// partitionCells returns [start,count) per worker. With side>0, assigns whole
// xy-planes (z-slabs) so the boundary is a single cut-plane (minimal ghost
// halo). Falls back to contiguous id ranges otherwise.
func partitionCells(nc, nWorkers, side int) (starts, counts []int) {
	starts = make([]int, nWorkers)
	counts = make([]int, nWorkers)
	if side > 0 && side*side*side == nc && nWorkers <= side {
		plane := side * side
		baseP := 0
		for p := 0; p < nWorkers; p++ {
			np := side / nWorkers
			if p < side%nWorkers {
				np++
			}
			starts[p] = baseP * plane
			counts[p] = np * plane
			baseP += np
		}
		return starts, counts
	}
	base := 0
	for p := 0; p < nWorkers; p++ {
		n := nc / nWorkers
		if p < nc%nWorkers {
			n++
		}
		starts[p] = base
		counts[p] = n
		base += n
	}
	return starts, counts
}

// Wrap folds a separation into the minimum image for a periodic box of side L.
func Wrap(d, L float64) float64 {
	if d > 0.5*L {
		d -= L
	}
	if d < -0.5*L {
		d += L
	}
	return d
}
