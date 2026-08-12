package hunt

import "math"

// Motif is a connected support of live Em (links, both ends live).
type Motif struct {
	ID         int
	Cells      []int
	Em         float64
	Ee         float64
	Es         float64
	Cx, Cy, Cz float64
	Rms        float64
	EmMin      float64
	EmMax      float64
	N          int
	NTag       int
}

type CensusCfg struct {
	EmFloor float64 // cell is live if Em > EmFloor (default 0.05)
	MinN    int     // drop components smaller than this (default 3)
	TagOnly bool    // restrict to tag>0.5 (painted ring / blob)
}

func (c CensusCfg) norm() CensusCfg {
	if c.EmFloor <= 0 {
		c.EmFloor = 0.05
	}
	if c.MinN <= 0 {
		c.MinN = 3
	}
	return c
}

func liveMask(fr Frame, floor float64) []bool {
	live := make([]bool, len(fr.Cells))
	for i, c := range fr.Cells {
		live[i] = c.Em > floor
	}
	return live
}

func unionFind(n int) (parent []int, find func(int) int, union func(int, int)) {
	parent = make([]int, n)
	rank := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	find = func(i int) int {
		for parent[i] != i {
			parent[i] = parent[parent[i]]
			i = parent[i]
		}
		return i
	}
	union = func(a, b int) {
		ra, rb := find(a), find(b)
		if ra == rb {
			return
		}
		if rank[ra] < rank[rb] {
			parent[ra] = rb
		} else if rank[ra] > rank[rb] {
			parent[rb] = ra
		} else {
			parent[rb] = ra
			rank[ra]++
		}
	}
	return parent, find, union
}

// Census connected live-Em components on the frame's link graph.
func Census(fr Frame, cfg CensusCfg) []Motif {
	cfg = cfg.norm()
	n := len(fr.Cells)
	live := liveMask(fr, cfg.EmFloor)
	if cfg.TagOnly {
		for i, c := range fr.Cells {
			if c.Tag <= 0.5 {
				live[i] = false
			}
		}
	}
	_, find, union := unionFind(n)
	for _, ln := range fr.Links {
		if ln.I < 0 || ln.J < 0 || ln.I >= n || ln.J >= n {
			continue
		}
		if live[ln.I] && live[ln.J] {
			union(ln.I, ln.J)
		}
	}
	groups := map[int][]int{}
	for i := 0; i < n; i++ {
		if !live[i] {
			continue
		}
		r := find(i)
		groups[r] = append(groups[r], i)
	}
	out := make([]Motif, 0, len(groups))
	for _, idxs := range groups {
		if len(idxs) < cfg.MinN {
			continue
		}
		out = append(out, motifOf(fr, idxs))
	}
	// stable-ish order: heaviest first
	for i := 0; i < len(out); i++ {
		for j := i + 1; j < len(out); j++ {
			if out[j].Em > out[i].Em {
				out[i], out[j] = out[j], out[i]
			}
		}
	}
	return out
}

func motifOf(fr Frame, idxs []int) Motif {
	m := Motif{Cells: append([]int(nil), idxs...), N: len(idxs), EmMin: math.Inf(1)}
	xs := make([]float64, len(idxs))
	ys := make([]float64, len(idxs))
	zs := make([]float64, len(idxs))
	ws := make([]float64, len(idxs))
	for k, i := range idxs {
		c := fr.Cells[i]
		m.Em += c.Em
		m.Ee += c.Ee
		m.Es += c.Es
		if c.Em < m.EmMin {
			m.EmMin = c.Em
		}
		if c.Em > m.EmMax {
			m.EmMax = c.Em
		}
		if c.Tag > 0.5 {
			m.NTag++
		}
		xs[k], ys[k], zs[k], ws[k] = c.X, c.Y, c.Z, c.Em
	}
	m.Cx, m.Cy, m.Cz, _ = centroid(xs, ys, zs, ws, fr.L)
	if m.Em > 0 {
		var r2 float64
		for k := range idxs {
			dx := wrap(xs[k]-m.Cx, fr.L)
			dy := wrap(ys[k]-m.Cy, fr.L)
			dz := wrap(zs[k]-m.Cz, fr.L)
			r2 += ws[k] * (dx*dx + dy*dy + dz*dz)
		}
		m.Rms = math.Sqrt(r2 / m.Em)
	}
	return m
}

// Track assigns stable IDs by wrapped-centroid proximity.
type Track struct {
	nextID int
	prev   []Motif
	MaxR   float64 // max match distance (default 4)
}

func (tr *Track) Assign(ms []Motif, L float64) []Motif {
	if tr.MaxR <= 0 {
		tr.MaxR = 4
	}
	used := make([]bool, len(tr.prev))
	for i := range ms {
		best := -1
		bestD := tr.MaxR
		for j := range tr.prev {
			if used[j] {
				continue
			}
			d := dist3(ms[i].Cx, ms[i].Cy, ms[i].Cz, tr.prev[j].Cx, tr.prev[j].Cy, tr.prev[j].Cz, L)
			if d < bestD {
				bestD = d
				best = j
			}
		}
		if best >= 0 {
			ms[i].ID = tr.prev[best].ID
			used[best] = true
		} else {
			tr.nextID++
			ms[i].ID = tr.nextID
		}
	}
	tr.prev = append([]Motif(nil), ms...)
	return ms
}

// Dominant is the heaviest motif, or a zero motif if none.
func Dominant(ms []Motif) Motif {
	if len(ms) == 0 {
		return Motif{}
	}
	return ms[0]
}
