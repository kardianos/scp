// Package par is the parallel harness.
//
// Each worker owns a permanent Shard with exclusive make()'d Lo/Hi slices.
// Ownership is a spatial z-slab when the world is a cubic lattice (fewer
// boundary faces than arbitrary id strips). Worker id never migrates.
//
// Ghost exchange is selective: Lo-only when geometry/space need remote
// positions/energies; FaNext-only (two float arrays) for boundary field hops.
package par

import (
	"sync"

	"scp/pfc/core"
)

// Cluster is a domain-decomposed simulation: one Shard per worker.
type Cluster struct {
	Law     core.Law
	Workers int
	shards  []Shard

	// publish buffers (global id indexing)
	pubLo  []core.CellLo
	pubFa1 []float64 // Fa1Next only — hot path for boundary hops
	pubFa2 []float64
	pubS   []float64 // space scale / fsum exchange
	owner  []int32

	// boundary edges owned by min-global-id endpoint, as global pairs +
	// a pointer back to a representative local slot for WField/geom.
	// Applied single-threaded after interior hops (no region mutex storm).
	bndGi []int32
	bndGj []int32
	// which shard holds the slot for WField (owner of min gid)
	bndShard []int32
	bndLocal []int32 // local slot index on that shard

	globalNC int
	side     int // lattice side, 0 if unknown

	// procedural workers (one goroutine per shard when Workers>1)
	wmu    sync.Mutex
	wcond  *sync.Cond
	wgen   uint64
	wdone  int
	wquit  bool
	wphase int
}

// FromWorld partitions w into nWorkers shards with exclusive Lo/Hi.
func FromWorld(w *core.World, nWorkers int) *Cluster {
	if nWorkers < 1 {
		nWorkers = 1
	}
	nc := w.NC()
	if nWorkers > nc {
		nWorkers = nc
	}
	cl := &Cluster{
		Law:      w.Law,
		Workers:  nWorkers,
		shards:   make([]Shard, nWorkers),
		pubLo:    make([]core.CellLo, nc),
		pubFa1:   make([]float64, nc),
		pubFa2:   make([]float64, nc),
		pubS:     make([]float64, nc),
		owner:    make([]int32, nc),
		globalNC: nc,
		side:     w.Side,
	}
	copy(cl.pubLo, w.Lo)

	// --- spatial ownership: z-slabs when Side known, else id ranges ---
	starts, counts := partitionCells(nc, nWorkers, w.Side)
	for p := 0; p < nWorkers; p++ {
		base, n := starts[p], counts[p]
		sh := &cl.shards[p]
		sh.ID = p
		sh.G0 = base
		sh.NOwn = n
		for i := 0; i < n; i++ {
			cl.owner[base+i] = int32(p)
		}
		sh.Lo = make([]core.CellLo, n)
		sh.Hi = make([]core.CellHi, n)
		copy(sh.Lo, w.Lo[base:base+n])
		copy(sh.Hi, w.Hi[base:base+n])
		sh.g2l = make(map[int32]int32, n*2)
		for i := 0; i < n; i++ {
			sh.g2l[int32(base+i)] = int32(i)
		}
	}

	// slots + ghosts
	ghostNeed := make([]map[int32]struct{}, nWorkers)
	type gslot struct{ gi, gj int32 }
	per := make([][]gslot, nWorkers)
	for p := 0; p < nWorkers; p++ {
		ghostNeed[p] = make(map[int32]struct{})
	}
	for si := range w.Slots {
		s := &w.Slots[si]
		if s.Alive == 0 {
			continue
		}
		gi, gj := s.I, s.J
		pi, pj := cl.owner[gi], cl.owner[gj]
		add := func(p int32) {
			per[p] = append(per[p], gslot{gi, gj})
			if cl.owner[gi] != p {
				ghostNeed[p][gi] = struct{}{}
			}
			if cl.owner[gj] != p {
				ghostNeed[p][gj] = struct{}{}
			}
		}
		add(pi)
		if pj != pi {
			add(pj)
		}
	}

	for p := 0; p < nWorkers; p++ {
		sh := &cl.shards[p]
		gs := sortedKeys(ghostNeed[p])
		sh.NGhost = len(gs)
		sh.GhostGID = gs
		sh.Lo = append(sh.Lo, make([]core.CellLo, sh.NGhost)...)
		sh.Hi = append(sh.Hi, make([]core.CellHi, sh.NGhost)...)
		for k, gid := range gs {
			li := int32(sh.NOwn + k)
			sh.g2l[gid] = li
			sh.Lo[li] = w.Lo[gid]
			sh.Hi[li] = w.Hi[gid]
		}

		sh.Slots = make([]core.Slot, 0, len(per[p]))
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
			li, okI := sh.g2l[g.gi]
			lj, okJ := sh.g2l[g.gj]
			if !okI || !okJ {
				continue
			}
			idx := int32(len(sh.Slots))
			sh.Slots = append(sh.Slots, core.Slot{I: li, J: lj, Alive: 1})
			if int(li) < sh.NOwn && int(lj) < sh.NOwn {
				sh.Interior = append(sh.Interior, idx)
			} else {
				sh.Boundary = append(sh.Boundary, idx)
				// owned endpoints of boundary edges need FaNext in pub
				if int(li) < sh.NOwn {
					sh.BndTouch = append(sh.BndTouch, li)
				}
				if int(lj) < sh.NOwn {
					sh.BndTouch = append(sh.BndTouch, lj)
				}
				// global boundary table: one entry per undirected edge
				ga, gb := g.gi, g.gj
				if ga > gb {
					ga, gb = gb, ga
				}
				if cl.owner[ga] == int32(p) {
					cl.bndGi = append(cl.bndGi, g.gi)
					cl.bndGj = append(cl.bndGj, g.gj)
					cl.bndShard = append(cl.bndShard, int32(p))
					cl.bndLocal = append(cl.bndLocal, idx)
				}
			}
		}
		sh.rebuildCSR()
		sh.buildInteriorColors()
		sh.BndTouch = uniqInt32(sh.BndTouch)
		nloc := sh.NOwn + sh.NGhost
		sh.sscl = make([]float64, nloc)
		sh.fsum = make([]float64, nloc)
		sh.refreshGeom(cl.Law.L)
	}

	// Publish set: owned cells that appear as ghosts on any shard (faces).
	needPub := make([][]int32, nWorkers)
	for p := range cl.shards {
		for _, gid := range cl.shards[p].GhostGID {
			op := int(cl.owner[gid])
			local := int32(int(gid) - cl.shards[op].G0)
			needPub[op] = append(needPub[op], local)
		}
	}
	for p := range cl.shards {
		cl.shards[p].Publish = uniqInt32(needPub[p])
	}

	// seed pubLo so first step can pull ghosts
	for p := range cl.shards {
		cl.shards[p].publishFace(cl.pubLo)
		// also publish all owned once at init (ghosts may need non-face if topology odd)
		for i := 0; i < cl.shards[p].NOwn; i++ {
			cl.pubLo[cl.shards[p].G0+i] = cl.shards[p].Lo[i]
		}
	}

	// Pre-size Live/Active so Step appends never allocate.
	for p := range cl.shards {
		sh := &cl.shards[p]
		nsl := len(sh.Slots)
		sh.Live = make([]int32, 0, nsl)
		sh.Active = make([]int32, 0, sh.NOwn)
		sh.LiveBits = make([]uint64, (nsl+63)/64)
		sh.actBits = make([]uint64, (sh.NOwn+63)/64)
	}

	if nWorkers > 1 {
		cl.startWorkers() // one go routine per shard; see workers.go
	}
	return cl
}

// partitionCells returns [start,count) per worker.
// With side>0, assigns whole xy-planes (z-slabs) for face-minimal cuts.
func partitionCells(nc, nWorkers, side int) (starts, counts []int) {
	starts = make([]int, nWorkers)
	counts = make([]int, nWorkers)
	if side > 0 && side*side*side == nc && nWorkers <= side {
		plane := side * side
		// distribute `side` planes across workers
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
	// fallback: contiguous id ranges
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

func uniqInt32(a []int32) []int32 {
	if len(a) == 0 {
		return a
	}
	// sort + unique
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

// Stop parks the per-shard goroutines.
func (c *Cluster) Stop() {
	c.stopWorkers()
}

// Shard returns worker id's permanent shard.
func (c *Cluster) Shard(id int) *Shard { return &c.shards[id] }

// NBoundary returns count of cross-shard edges (diagnostic).
func (c *Cluster) NBoundary() int { return len(c.bndGi) }

// TotalE sums owned Lo only.
func (c *Cluster) TotalE() float64 {
	var e float64
	for p := range c.shards {
		sh := &c.shards[p]
		for i := 0; i < sh.NOwn; i++ {
			lo := &sh.Lo[i]
			e += lo.Em + lo.Ee + lo.Es
		}
	}
	return e
}


