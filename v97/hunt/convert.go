package hunt

import "math"

// Convert is the frame-to-frame mode-crossing estimate.
// Intake at a cell: Em rose while Ee fell (condensation-shaped).
// Outtake: Em fell while Ee rose (radiance/evap-shaped).
// Dense hops and field hops contaminate this; H0 uses the *pattern*
// (near-zero on sealed/vacuum, not the absolute grain).

type Convert struct {
	I, O      float64
	Ifore     float64
	Ofore     float64
	Iaft      float64
	Oaft      float64
	Sep       float64 // (c_I − c_O) · ĵ   ; 0 = point book
	Phi       float64 // min(I, O)
	Matched   int
	Unmatched int
	// Exact: I/O came from the promo_i/promo_o ledger columns
	// (gross door traffic), not the mode-crossing heuristic.
	Exact bool
	// Ih/Oh: the heuristic estimate, always computed, so the
	// inflation of the mode-crossing meter stays measurable
	// (H0 `exact` class; REGRADE).
	Ih, Oh float64
}

// Match cells of b to a by wrapped position (roster indices are not identity).
func matchCells(a, b Frame, maxR float64) []int {
	if maxR <= 0 {
		maxR = 0.6
	}
	na, nb := len(a.Cells), len(b.Cells)
	out := make([]int, nb)
	for i := range out {
		out[i] = -1
	}
	if na == 0 || nb == 0 {
		return out
	}
	// spatial hash
	h := maxR
	type key struct{ x, y, z int }
	bins := map[key][]int{}
	bin := func(x, y, z float64) key {
		return key{int(math.Floor(x / h)), int(math.Floor(y / h)), int(math.Floor(z / h))}
	}
	L := b.L
	if L <= 0 {
		L = a.L
	}
	for i, c := range a.Cells {
		bins[bin(c.X, c.Y, c.Z)] = append(bins[bin(c.X, c.Y, c.Z)], i)
	}
	for j, c := range b.Cells {
		best := -1
		bestD := maxR
		b0 := bin(c.X, c.Y, c.Z)
		for dx := -1; dx <= 1; dx++ {
			for dy := -1; dy <= 1; dy++ {
				for dz := -1; dz <= 1; dz++ {
					for _, i := range bins[key{b0.x + dx, b0.y + dy, b0.z + dz}] {
						d := dist3(c.X, c.Y, c.Z, a.Cells[i].X, a.Cells[i].Y, a.Cells[i].Z, L)
						if d < bestD {
							bestD = d
							best = i
						}
					}
				}
			}
		}
		// also search wrapped image bins when near a face
		if best < 0 && L > 0 {
			for _, img := range [][3]float64{
				{c.X + L, c.Y, c.Z}, {c.X - L, c.Y, c.Z},
				{c.X, c.Y + L, c.Z}, {c.X, c.Y - L, c.Z},
				{c.X, c.Y, c.Z + L}, {c.X, c.Y, c.Z - L},
			} {
				k := bin(img[0], img[1], img[2])
				for _, i := range bins[k] {
					d := dist3(c.X, c.Y, c.Z, a.Cells[i].X, a.Cells[i].Y, a.Cells[i].Z, L)
					if d < bestD {
						bestD = d
						best = i
					}
				}
			}
		}
		out[j] = best
	}
	return out
}

func cellIO(prev, cur Cell) (I, O float64) {
	dEm := cur.Em - prev.Em
	dEe := cur.Ee - prev.Ee
	if dEm > 0 && dEe < 0 {
		I = math.Min(dEm, -dEe)
	}
	if dEm < 0 && dEe > 0 {
		O = math.Min(-dEm, dEe)
	}
	return I, O
}

// MeasureConvert estimates I/O on motif m (cells of `cur`) vs `prev`.
// Axis (jx,jy,jz) splits fore/aft; if near-zero, the motif's longest
// principal direction is used.
//
// dilate: expand the I/O collection set by that many rings of LINK
// neighbours (Em accounting is untouched). Rationale, measured on
// rg2a: promo>0 books each door event at the oriented link's FAR
// cell (fab/promo.go dest), so a machine's own deposits land one
// cell outside its census and ⟨O⟩ reads 0 while the box ledger
// carries the full rad total. Direction acknowledged: dilation can
// only raise motif I/O; it is scope correction (the link's far end
// is part of the machine's book by the promo semantics), not a
// softer bar — default remains 0 and both gradings are reported.
func MeasureConvert(prev, cur Frame, m Motif, jx, jy, jz float64, dilate int, dest bool) Convert {
	var out Convert
	if m.N == 0 || len(prev.Cells) == 0 {
		return out
	}
	inside := make([]bool, len(cur.Cells))
	for _, i := range m.Cells {
		if i >= 0 && i < len(inside) {
			inside[i] = true
		}
	}
	// Dest-only scope (dest.go): overrides dilate. I counted only on
	// the motif ∪ replayed fore picks; O only on motif ∪ aft picks.
	var iSet, oSet []bool
	if dest && cur.HasPromo {
		iSet, oSet = DestSets(cur, m)
		for k := range inside {
			inside[k] = iSet[k] || oSet[k]
		}
		dilate = 0
	}
	for d := 0; d < dilate; d++ {
		grow := make([]bool, len(inside))
		copy(grow, inside)
		for _, lk := range cur.Links {
			if lk.I >= 0 && lk.I < len(inside) && lk.J >= 0 && lk.J < len(inside) {
				if inside[lk.I] {
					grow[lk.J] = true
				}
				if inside[lk.J] {
					grow[lk.I] = true
				}
			}
		}
		inside = grow
	}
	jn := hypot3(jx, jy, jz)
	if jn < 1e-12 {
		jx, jy, jz = motifAxis(cur, m)
		jn = hypot3(jx, jy, jz)
	}
	if jn > 0 {
		jx, jy, jz = jx/jn, jy/jn, jz/jn
	}
	// Exact ledger vs heuristic is a SCHEMA property, not a per-frame
	// one: when the stream carries promo_i/promo_o, use them always,
	// zeros included. The old per-frame opportunistic switch fell back
	// to the hop-shaped heuristic on every quiet snap (atom grain),
	// inflating I/O and faking simultaneity (REVIEW_grok_next §F1).
	usePromo := cur.HasPromo
	out.Exact = usePromo
	mt := matchCells(prev, cur, 0.6)
	var ixs, iys, izs, iws []float64
	var oxs, oys, ozs, ows []float64
	L := cur.L
	for j, c := range cur.Cells {
		if !inside[j] {
			continue
		}
		var I, O float64
		pi := mt[j]
		var hI, hO float64
		if pi >= 0 {
			hI, hO = cellIO(prev.Cells[pi], c)
		}
		out.Ih += hI
		out.Oh += hO
		if usePromo {
			I, O = c.PromoI, c.PromoO
			if iSet != nil {
				if !iSet[j] {
					I = 0
				}
				if !oSet[j] {
					O = 0
				}
			}
			out.Matched++
		} else {
			if pi < 0 {
				out.Unmatched++
				continue
			}
			out.Matched++
			I, O = hI, hO
		}
		proj := wrap(c.X-m.Cx, L)*jx + wrap(c.Y-m.Cy, L)*jy + wrap(c.Z-m.Cz, L)*jz
		out.I += I
		out.O += O
		if proj >= 0 {
			out.Ifore += I
			out.Ofore += O
		} else {
			out.Iaft += I
			out.Oaft += O
		}
		if I > 0 {
			ixs = append(ixs, c.X)
			iys = append(iys, c.Y)
			izs = append(izs, c.Z)
			iws = append(iws, I)
		}
		if O > 0 {
			oxs = append(oxs, c.X)
			oys = append(oys, c.Y)
			ozs = append(ozs, c.Z)
			ows = append(ows, O)
		}
	}
	out.Phi = math.Min(out.I, out.O)
	if len(iws) > 0 && len(ows) > 0 && jn > 0 {
		cix, ciy, ciz, _ := centroid(ixs, iys, izs, iws, L)
		cox, coy, coz, _ := centroid(oxs, oys, ozs, ows, L)
		out.Sep = wrap(cix-cox, L)*jx + wrap(ciy-coy, L)*jy + wrap(ciz-coz, L)*jz
	}
	return out
}

func motifAxis(fr Frame, m Motif) (x, y, z float64) {
	// largest RMS axis via simple covariance (xy-preferred for rings).
	var xx, yy, zz float64
	L := fr.L
	for _, i := range m.Cells {
		c := fr.Cells[i]
		if c.Em <= 0 {
			continue
		}
		dx := wrap(c.X-m.Cx, L)
		dy := wrap(c.Y-m.Cy, L)
		dz := wrap(c.Z-m.Cz, L)
		xx += c.Em * dx * dx
		yy += c.Em * dy * dy
		zz += c.Em * dz * dz
	}
	if yy >= xx && yy >= zz {
		return 0, 1, 0
	}
	if zz >= xx {
		return 0, 0, 1
	}
	return 1, 0, 0
}
