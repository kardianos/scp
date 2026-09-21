package parcore

// shard is one worker's permanent memory. Cells and Slots are separately
// allocated contiguous slices that stay on this shard for the cluster's life
// (distinct backing arrays per shard — see TestDistinctSlabs).
//
// Layout:
//
//	Cells [0:NOwn)           — cells this worker owns (global G0..G0+NOwn-1)
//	Cells [NOwn:NOwn+NGhost) — ghost copies of remote endpoints
//
// Slot SI/SJ are LOCAL cell indices. Worker id == shard index, fixed for life.
type shard[C any, S any] struct {
	id   int
	g0   int
	nOwn int

	cells []C
	slots []S

	// topology (local indices)
	si, sj   []int32
	cls, inc []int32
	colors   [][]int32 // interior color buckets
	live     []int32
	liveCap  int      // capacity Live was sized to (for zero-alloc reset)

	slotParam []float64

	// ghosts
	nGhost   int
	ghostGID []int32

	// interior vs boundary slot lists (local slot indices)
	interior  []int32
	boundary  []int32
	bndTouch  []int32 // owned local cells touching a boundary edge (publish for hops)
	publish   []int32 // owned local cells some other shard ghosts (the face)

	scratch any

	// the shard's public view (reused every Phase.Run to avoid alloc)
	view View[C, S]
}

// public returns a stable ShardView snapshot for diagnostics/tests.
func (s *shard[C, S]) public() ShardView[C, S] {
	return ShardView[C, S]{
		ID:     s.id,
		G0:     s.g0,
		NOwn:   s.nOwn,
		NGhost: s.nGhost,
		Cells:  s.cells,
		Slots:  s.slots,
	}
}

// resetLive clears the Live list in place (zero heap) for the next geom pass.
func (s *shard[C, S]) resetLive() {
	s.live = s.live[:0]
}

// FromWorld partitions w into nWorkers shards with exclusive Cells/Slots and
// builds the per-shard CSR, interior edge coloring, ghost halo, and the
// cross-shard boundary edge table. It then calls driver.Step() once to cache
// the phase graph and starts the procedural worker goroutines (if Workers>1).
func FromWorld[C any, S any](w *World[C, S], driver Driver[C, S], nWorkers int) *Cluster[C, S] {
	if nWorkers < 1 {
		nWorkers = 1
	}
	nc := w.NC()
	if nWorkers > nc {
		nWorkers = nc
	}
	c := &Cluster[C, S]{
		driver:   driver,
		Workers:  nWorkers,
		shards:   make([]shard[C, S], nWorkers),
		pub:      make([]C, nc),
		globalNC: nc,
		side:     w.Side,
		box:      w.Box,
	}

	// seed the global publish buffer so the first Pull sees real state.
	copy(c.pub, w.Cells)

	// ownership (z-slabs when Side known, else id ranges)
	owner := make([]int32, nc)
	starts, counts := partitionCells(nc, nWorkers, w.Side)
	for p := 0; p < nWorkers; p++ {
		base, n := starts[p], counts[p]
		sh := &c.shards[p]
		sh.id = p
		sh.g0 = base
		sh.nOwn = n
		sh.cells = make([]C, n)
		sh.slots = make([]S, 0, 3*n)
		copy(sh.cells, w.Cells[base:base+n])
		for i := 0; i < n; i++ {
			owner[base+i] = int32(p)
		}
	}

	// distribute slots: an edge lives on both endpoints' shards (ghosting the
	// remote end). Boundary edges (endpoints on different shards) feed the
	// global boundary table once (owned by the min-gid endpoint's shard).
	type gslot struct{ gi, gj int32 }
	per := make([][]gslot, nWorkers)
	ghostNeed := make([]map[int32]struct{}, nWorkers)
	for p := 0; p < nWorkers; p++ {
		ghostNeed[p] = make(map[int32]struct{})
	}
	for s := range w.SI {
		if !w.Alive[s] {
			continue
		}
		gi, gj := w.SI[s], w.SJ[s]
		pi, pj := owner[gi], owner[gj]
		add := func(p int32) {
			per[p] = append(per[p], gslot{gi, gj})
			if owner[gi] != p {
				ghostNeed[p][gi] = struct{}{}
			}
			if owner[gj] != p {
				ghostNeed[p][gj] = struct{}{}
			}
		}
		add(pi)
		if pj != pi {
			add(pj)
		}
	}

	for p := 0; p < nWorkers; p++ {
		sh := &c.shards[p]
		gs := sortedKeys(ghostNeed[p])
		sh.nGhost = len(gs)
		sh.ghostGID = gs
		// grow owned slice to hold ghosts
		sh.cells = append(sh.cells, make([]C, sh.nGhost)...)
		for k, gid := range gs {
			sh.cells[sh.nOwn+k] = w.Cells[gid]
		}

		// local slots with dedup
		g2l := make(map[int32]int32, sh.nOwn+sh.nGhost)
		for i := 0; i < sh.nOwn; i++ {
			g2l[int32(sh.g0+i)] = int32(i)
		}
		for k, gid := range gs {
			g2l[gid] = int32(sh.nOwn + k)
		}
		seen := make(map[[2]int32]struct{}, len(per[p]))
		for _, g := range per[p] {
			a, b := g.gi, g.gj
			if a > b {
				a, b = b, a
			}
			key := [2]int32{a, b}
			if _, ok := seen[key]; ok {
				continue
			}
			seen[key] = struct{}{}
			li, okI := g2l[g.gi]
			lj, okJ := g2l[g.gj]
			if !okI || !okJ {
				continue
			}
			var zero S
			sh.slots = append(sh.slots, zero)
			sh.si = append(sh.si, li)
			sh.sj = append(sh.sj, lj)
			idx := int32(len(sh.si)) - 1
			if int(li) < sh.nOwn && int(lj) < sh.nOwn {
				sh.interior = append(sh.interior, idx)
			} else {
				sh.boundary = append(sh.boundary, idx)
				if int(li) < sh.nOwn {
					sh.bndTouch = append(sh.bndTouch, li)
				}
				if int(lj) < sh.nOwn {
					sh.bndTouch = append(sh.bndTouch, lj)
				}
				// boundary table: one entry per undirected edge, owned by
				// the min-gid endpoint's shard.
				ga, gb := g.gi, g.gj
				if ga > gb {
					ga, gb = gb, ga
				}
				if owner[ga] == int32(p) {
					c.bndGi = append(c.bndGi, g.gi)
					c.bndGj = append(c.bndGj, g.gj)
					c.bndShard = append(c.bndShard, int32(p))
					c.bndLocal = append(c.bndLocal, idx)
				}
			}
		}

		nLoc := sh.nOwn + sh.nGhost
		sh.cls, sh.inc = buildCSR(sh.si, sh.sj, allAlive(len(sh.si)), nLoc)
		sh.colors = edgeColor(sh.si, sh.sj, allAlive(len(sh.si)), sh.interior, nLoc)
		sh.bndTouch = uniqInt32(sh.bndTouch)
		sh.slotParam = make([]float64, len(sh.si))
		sh.live = make([]int32, 0, len(sh.si))
		sh.liveCap = cap(sh.live)

		// per-shard driver scratch
		if sd, ok := driver.(ScratchDriver[C, S]); ok {
			sh.scratch = sd.InitShard(sh.nOwn, sh.nGhost, len(sh.si))
		}

		sh.view = View[C, S]{
			G0:        sh.g0,
			NOwn:      sh.nOwn,
			Cells:     sh.cells,
			Slots:     sh.slots,
			SI:        sh.si,
			SJ:        sh.sj,
			Cls:       sh.cls,
			Inc:       sh.inc,
			Colors:    sh.colors,
			Live:      sh.live,
			SlotParam: sh.slotParam,
			BndTouch:  sh.bndTouch,
			Scratch:   sh.scratch,
			NGhost:    sh.nGhost,
			box:       c.box,
		}
	}

	// publish set: owned cells that appear as ghosts on any shard (the face).
	needPub := make([][]int32, nWorkers)
	for p := range c.shards {
		for _, gid := range c.shards[p].ghostGID {
			op := int(owner[gid])
			local := int32(int(gid) - c.shards[op].g0)
			needPub[op] = append(needPub[op], local)
		}
	}
	for p := range c.shards {
		c.shards[p].publish = uniqInt32(needPub[p])
	}

	// wire each shard's exchange internals into its view.
	for p := range c.shards {
		sh := &c.shards[p]
		sh.view.pub = c.pub
		sh.view.publishSet = sh.publish
		sh.view.ghostGID = sh.ghostGID
		// publish the face once so the first Pull works.
		sh.view.Publish()
	}

	c.phases = driver.Step()
	c.bndTau = make([]float64, len(c.bndGi))

	if nWorkers > 1 {
		c.startWorkers()
	}
	return c
}

// allAlive is a shortcut for the frozen-topology case where every slot is
// alive (the lattice builder marks them so; a dynamic topology would supply a
// real alive mask).
func allAlive(n int) []bool {
	a := make([]bool, n)
	for i := range a {
		a[i] = true
	}
	return a
}

func uniqInt32(a []int32) []int32 {
	if len(a) == 0 {
		return a
	}
	// insertion sort (small face sets)
	for i := 1; i < len(a); i++ {
		v := a[i]
		j := i - 1
		for j >= 0 && a[j] > v {
			a[j+1] = a[j]
			j--
		}
		a[j+1] = v
	}
	w := 1
	for i := 1; i < len(a); i++ {
		if a[i] != a[w-1] {
			a[w] = a[i]
			w++
		}
	}
	return a[:w]
}

func sortedKeys(m map[int32]struct{}) []int32 {
	gs := make([]int32, 0, len(m))
	for gid := range m {
		gs = append(gs, gid)
	}
	for i := 1; i < len(gs); i++ {
		v := gs[i]
		j := i - 1
		for j >= 0 && gs[j] > v {
			gs[j+1] = gs[j]
			j--
		}
		gs[j+1] = v
	}
	return gs
}

// TotalEnergy reduces a per-cell energy across owned cells using the driver-
// supplied measure. Used by tests and the bench's conservation check.
func TotalEnergy[C any, S any](c *Cluster[C, S], measure func(C) float64) float64 {
	var e float64
	for p := range c.shards {
		sh := &c.shards[p]
		for i := 0; i < sh.nOwn; i++ {
			e += measure(sh.cells[i])
		}
	}
	return e
}
