package hunt

import "math"

// Dest-only book scope (REVIEW_grok_state confound A): the promo
// door books each event at the oriented link's far cell. `-dilate 1`
// recovered those books but imported every LINK neighbour — in a
// z≈16 foam that is tens of bath cells whose own box-fill books then
// pollute the motif's Φ (78% of dilated I on rg2c sat on the halo).
//
// The honest scope is the set of cells a motif door can actually
// have booked to: the fore/aft picks of promoPick, replayed here on
// the frame state, mirroring fab/promo.go (linkJ with matter-J and
// field-J fallback, best û·Ĵ with dot>0 fore / dot<0 aft).
//
// Approximations, named: picks are replayed from the `cur` snapshot
// while real events fired throughout the interval on evolving state;
// and a dest cell can also carry books from OTHER doors (bath cells
// firing into the same cell). The set is minimal (≤2 cells per motif
// cell), so the contamination is bounded by construction, unlike the
// dilate halo.

type adjEntry struct {
	j          int
	ux, uy, uz float64
}

func buildAdj(fr Frame) [][]adjEntry {
	adj := make([][]adjEntry, len(fr.Cells))
	L := fr.L
	for _, lk := range fr.Links {
		if lk.I < 0 || lk.J < 0 || lk.I >= len(fr.Cells) || lk.J >= len(fr.Cells) || lk.A <= 0 {
			continue
		}
		ci, cj := fr.Cells[lk.I], fr.Cells[lk.J]
		dx := wrap(cj.X-ci.X, L)
		dy := wrap(cj.Y-ci.Y, L)
		dz := wrap(cj.Z-ci.Z, L)
		d := math.Sqrt(dx*dx + dy*dy + dz*dz)
		if d <= 0 {
			continue
		}
		ux, uy, uz := dx/d, dy/d, dz/d
		adj[lk.I] = append(adj[lk.I], adjEntry{lk.J, ux, uy, uz})
		adj[lk.J] = append(adj[lk.J], adjEntry{lk.I, -ux, -uy, -uz})
	}
	return adj
}

// linkJmatter mirrors fab/promo.go linkJ: 2·√(Em_i·Em_j)·sin(θj−θi).
func linkJmatter(a, b Cell) float64 {
	if a.Em <= 0 || b.Em <= 0 {
		return 0
	}
	return 2 * math.Sqrt(a.Em) * math.Sqrt(b.Em) * math.Sin(b.Th2-a.Th2)
}

func linkJfld(a, b Cell) float64 {
	return 2 * (a.Fa1*b.Fa2 - a.Fa2*b.Fa1)
}

// destPicks replays promoPick(i, fore) and promoPick(i, aft) for one
// cell. Returns (-1,-1) components when no neighbour qualifies (the
// kernel then books at i itself, which the motif set already covers).
func destPicks(fr Frame, adj [][]adjEntry, i int) (fore, aft int) {
	fore, aft = -1, -1
	var jx, jy, jz float64
	for _, e := range adj[i] {
		J := linkJmatter(fr.Cells[i], fr.Cells[e.j])
		jx += J * e.ux
		jy += J * e.uy
		jz += J * e.uz
	}
	if jx*jx+jy*jy+jz*jz < 1e-24 {
		jx, jy, jz = 0, 0, 0
		for _, e := range adj[i] {
			J := linkJfld(fr.Cells[i], fr.Cells[e.j])
			jx += J * e.ux
			jy += J * e.uy
			jz += J * e.uz
		}
	}
	if jx*jx+jy*jy+jz*jz < 1e-24 {
		return
	}
	bestF, bestA := -1e300, 1e300
	for _, e := range adj[i] {
		dot := e.ux*jx + e.uy*jy + e.uz*jz
		if dot > 0 && dot > bestF {
			bestF, fore = dot, e.j
		}
		if dot < 0 && dot < bestA {
			bestA, aft = dot, e.j
		}
	}
	return
}

// DestSets returns the I-collection and O-collection cell sets for a
// motif under dest-only booking: the motif's own cells plus, per
// motif cell, its replayed fore pick (intake dest) and aft pick
// (outtake dest).
func DestSets(fr Frame, m Motif) (iSet, oSet []bool) {
	iSet = make([]bool, len(fr.Cells))
	oSet = make([]bool, len(fr.Cells))
	adj := buildAdj(fr)
	for _, i := range m.Cells {
		if i < 0 || i >= len(fr.Cells) {
			continue
		}
		iSet[i], oSet[i] = true, true
		fore, aft := destPicks(fr, adj, i)
		if fore >= 0 {
			iSet[fore] = true
		}
		if aft >= 0 {
			oSet[aft] = true
		}
	}
	return
}
