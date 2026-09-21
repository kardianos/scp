package hunt

import (
	"math"
	"testing"
)

func ringFrame(n int, mWind int, em float64, L float64) Frame {
	fr := Frame{T: 0, L: L, Cells: make([]Cell, n), Links: make([]Link, n)}
	R := 2.0
	cx, cy := L/2, L/2
	for i := 0; i < n; i++ {
		phi := 2 * math.Pi * float64(i) / float64(n)
		fr.Cells[i] = Cell{
			X: cx + R*math.Cos(phi), Y: cy + R*math.Sin(phi), Z: L / 2,
			R: 0.8, Es: 1, Em: em, Th2: float64(mWind) * phi,
		}
		j := (i + 1) % n
		dx := fr.Cells[j].X - fr.Cells[i].X
		dy := fr.Cells[j].Y - fr.Cells[i].Y
		fr.Links[i] = Link{I: i, J: j, D: math.Hypot(dx, dy), A: 1}
	}
	return fr
}

func TestCensusRing(t *testing.T) {
	fr := ringFrame(6, 2, 0.8, 16)
	ms := Census(fr, CensusCfg{})
	if len(ms) != 1 {
		t.Fatalf("motifs=%d want 1", len(ms))
	}
	if ms[0].N != 6 {
		t.Fatalf("n=%d want 6", ms[0].N)
	}
	if math.Abs(ms[0].Em-4.8) > 1e-9 {
		t.Fatalf("Em=%g", ms[0].Em)
	}
}

func TestCensusTwoBlobs(t *testing.T) {
	fr := Frame{L: 16, Cells: []Cell{
		{X: 2, Y: 2, Z: 8, Em: 1}, {X: 2.4, Y: 2, Z: 8, Em: 1}, {X: 2, Y: 2.4, Z: 8, Em: 1},
		{X: 10, Y: 10, Z: 8, Em: 1}, {X: 10.4, Y: 10, Z: 8, Em: 1}, {X: 10, Y: 10.4, Z: 8, Em: 1},
	}, Links: []Link{
		{I: 0, J: 1, D: 0.4, A: 1}, {I: 0, J: 2, D: 0.4, A: 1}, {I: 1, J: 2, D: 0.5, A: 1},
		{I: 3, J: 4, D: 0.4, A: 1}, {I: 3, J: 5, D: 0.4, A: 1}, {I: 4, J: 5, D: 0.5, A: 1},
	}}
	ms := Census(fr, CensusCfg{})
	if len(ms) != 2 {
		t.Fatalf("motifs=%d want 2", len(ms))
	}
}

func TestWindingM2(t *testing.T) {
	fr := ringFrame(12, 2, 1, 16)
	ms := Census(fr, CensusCfg{})
	W, Rm, n := Winding(fr, ms[0])
	if n != 12 {
		t.Fatalf("n=%d", n)
	}
	if math.Abs(W-2) > 0.05 {
		t.Fatalf("W=%g want 2", W)
	}
	if Rm < 0.99 {
		t.Fatalf("Rm=%g", Rm)
	}
}

func TestLoopNeutral(t *testing.T) {
	fr := ringFrame(8, 2, 1, 16)
	ms := Census(fr, CensusCfg{})
	cur := MeasureCurrents(fr, ms[0], nil)
	if math.Abs(cur.Circ) < 1e-6 {
		t.Fatalf("expected circulation, got %g", cur.Circ)
	}
	if cur.JnetAbs() > 0.25*math.Abs(cur.Circ) {
		t.Fatalf("loop should be Jnet~0: Jnet=%g Circ=%g", cur.JnetAbs(), cur.Circ)
	}
	if cur.JdipAbs() > 0.4*math.Abs(cur.Circ) {
		t.Fatalf("loop should be Jdip~0: Jdip=%g Circ=%g", cur.JdipAbs(), cur.Circ)
	}
}

func TestTrackPersists(t *testing.T) {
	fr := ringFrame(6, 2, 0.8, 16)
	tr := &Track{}
	a := tr.Assign(Census(fr, CensusCfg{}), 16)
	// shift the ring
	for i := range fr.Cells {
		fr.Cells[i].X += 0.3
	}
	b := tr.Assign(Census(fr, CensusCfg{}), 16)
	if a[0].ID == 0 || a[0].ID != b[0].ID {
		t.Fatalf("id %d -> %d", a[0].ID, b[0].ID)
	}
}

func TestConvertCondensation(t *testing.T) {
	prev := Frame{L: 16, Cells: []Cell{
		{X: 8, Y: 8, Z: 8, Em: 0.5, Ee: 1.0},
		{X: 8.4, Y: 8, Z: 8, Em: 0.5, Ee: 0.1},
		{X: 8, Y: 8.4, Z: 8, Em: 0.5, Ee: 0.1},
	}, Links: []Link{{I: 0, J: 1, D: 0.4, A: 1}, {I: 0, J: 2, D: 0.4, A: 1}, {I: 1, J: 2, D: 0.5, A: 1}}}
	cur := prev
	cur.Cells = append([]Cell(nil), prev.Cells...)
	cur.Cells[0].Em = 0.9
	cur.Cells[0].Ee = 0.6 // dEm=+0.4, dEe=-0.4 → I=0.4
	ms := Census(cur, CensusCfg{EmFloor: 0.1})
	cv := MeasureConvert(prev, cur, ms[0], 1, 0, 0, 0, false)
	if math.Abs(cv.I-0.4) > 1e-9 {
		t.Fatalf("I=%g want 0.4", cv.I)
	}
	if cv.O != 0 {
		t.Fatalf("O=%g want 0", cv.O)
	}
}

func TestGlobalJPlaneWave(t *testing.T) {
	// 1D chain, th2 = kx, k>0 ⇒ Jnet_x > 0
	const n, L, k = 8, 16.0, 1.1
	fr := Frame{L: L, Cells: make([]Cell, n), Links: make([]Link, n-1)}
	for i := 0; i < n; i++ {
		x := 2 + float64(i)*1.2
		fr.Cells[i] = Cell{X: x, Y: 8, Z: 8, Em: 1, Th2: k * x}
		if i+1 < n {
			fr.Links[i] = Link{I: i, J: i + 1, D: 1.2, A: 1}
		}
	}
	jx, _, _ := GlobalJ(fr)
	if jx <= 0 {
		t.Fatalf("GJnetX=%g want >0 for k>0", jx)
	}
}

func TestH1EmptyOnPointBook(t *testing.T) {
	rows := make([]Row, 15)
	for i := range rows {
		rows[i] = Row{
			T: float64(i) * 20, N: 50, Em: 100, EmTot: 100,
			Phi: 0.01, I: 0.01, O: 0, Sep: 0.1,
			Jdip: 1, Circ: 20, JdipX: float64(1 - 2*(i%2)),
			CoreEsThru: 0.2, DEmTot: 0,
		}
	}
	rep := GradeH1("point", rows)
	if rep.Inhabited {
		t.Fatalf("point book must be EMPTY\n%s", rep)
	}
	if rep.Best.NPass >= 5 {
		t.Fatalf("expected <5 clauses, got %d\n%s", rep.Best.NPass, rep)
	}
}

func TestH1InhabitedSynthetic(t *testing.T) {
	rows := make([]Row, 12)
	for i := range rows {
		rows[i] = Row{
			T: float64(i) * 20, N: 20, Em: 50, EmTot: 50, DEmTot: 0.01,
			Phi: 2, I: 2.2, O: 2.0, Sep: 1.5,
			Ifore: 2.0, Oaft: 1.8, Iaft: 0.2, Ofore: 0.2,
			Jdip: 8, Circ: 4, JdipX: 8, Jnet: 5,
			CoreEsThru: 3, CoreEsNet: 0.1,
		}
	}
	rep := GradeH1("synth", rows)
	if !rep.Inhabited {
		t.Fatalf("synthetic dipole must be INHABITED\n%s", rep)
	}
	if rep.Best.Span < 100 {
		t.Fatalf("span %g", rep.Best.Span)
	}
}

func TestH1ShortIsEmpty(t *testing.T) {
	rows := []Row{
		{T: 0, N: 10, Em: 10, Phi: 2, I: 2, O: 2, Sep: 2, Jdip: 5, Circ: 1, JdipX: 5, CoreEsThru: 4},
		{T: 40, N: 10, Em: 10, Phi: 2, I: 2, O: 2, Sep: 2, Jdip: 5, Circ: 1, JdipX: 5, CoreEsThru: 4},
	}
	rep := GradeH1("short", rows)
	if rep.Inhabited {
		t.Fatal("short stream must not inhabit")
	}
}

func TestPearson(t *testing.T) {
	xs := []float64{1, 2, 3, 4, 5}
	ys := []float64{2, 4, 6, 8, 10}
	if math.Abs(pearson(xs, ys)-1) > 1e-9 {
		t.Fatalf("pearson=%g", pearson(xs, ys))
	}
	zs := []float64{-1, -2, -3, -4, -5}
	if math.Abs(pearson(xs, zs)+1) > 1e-9 {
		t.Fatalf("pearson=%g", pearson(xs, zs))
	}
}

// --- exact-ledger regression (PLAN.md step 1′; REVIEW_grok_next §F1) ---
//
// The graded I/O on a promo stream must BE the kernel's own door
// ledger. Expected totals come from the kernel RESULT lines of the
// leads runs (runs/leads/p2_space.log):
//   cond=0.524917   evap=0.612923   rad=13.572040  (evap+rad=14.184963)
// If this test fails after a meter edit, the meter regressed to the
// hop-shaped heuristic (or the ledger columns broke).

func openLeads(t *testing.T, name string) []Row {
	t.Helper()
	s, err := Open("../runs/leads/" + name)
	if err != nil {
		t.Skipf("leads stream not present: %v", err)
	}
	defer s.Close()
	rows, err := Analyze(s, Options{Census: CensusCfg{EmFloor: 0.05, MinN: 3}})
	if err != nil {
		t.Fatalf("analyze: %v", err)
	}
	return rows
}

func TestExactLedgerP2Space(t *testing.T) {
	rows := openLeads(t, "p2_space.fcs")
	if len(rows) == 0 {
		t.Fatal("no frames")
	}
	var boxI, boxO, motI, motO, hI, hO float64
	for i, r := range rows {
		if !r.Ledger {
			t.Fatalf("frame %d: ledger missing on a promo stream", i)
		}
		boxI += r.BoxPI
		boxO += r.BoxPO
		if i == 0 {
			continue
		}
		motI += r.I
		motO += r.O
		hI += r.Ih
		hO += r.Oh
		if r.N > 0 && (r.I > r.BoxPI+1e-6 || r.O > r.BoxPO+1e-6) {
			t.Fatalf("frame %d: motif I/O exceeds box ledger (I=%.4g/%.4g O=%.4g/%.4g) — heuristic leaked in",
				i, r.I, r.BoxPI, r.O, r.BoxPO)
		}
	}
	const wantI, wantO = 0.524917, 14.184963
	if math.Abs(boxI-wantI) > 0.02 {
		t.Errorf("Σ box promo_i = %.6g, kernel cond = %.6g", boxI, wantI)
	}
	if math.Abs(boxO-wantO) > 0.30 {
		t.Errorf("Σ box promo_o = %.6g, kernel evap+rad = %.6g", boxO, wantO)
	}
	// The documented failure mode of the old meter: the heuristic
	// inflates I far beyond the door (~40× on this stream). Keep the
	// direction pinned so nobody "fixes" it backwards again.
	if hI < 5*boxI {
		t.Errorf("heuristic ΣI=%.4g not ≫ ledger ΣI=%.4g — inflation direction changed; re-verify against RESULT lines", hI, boxI)
	}
}

func TestHeuristicFallbackP0Ctl(t *testing.T) {
	rows := openLeads(t, "p0_ctl.fcs")
	if len(rows) == 0 {
		t.Fatal("no frames")
	}
	for i, r := range rows {
		if r.Ledger {
			t.Fatalf("frame %d: promo=0 stream must not report a ledger", i)
		}
		if r.BoxPI != 0 || r.BoxPO != 0 {
			t.Fatalf("frame %d: nonzero box promo sums without columns", i)
		}
	}
}

// Dest-only scope: books on the replayed fore/aft picks count; books
// on an unrelated neighbour do not (the dilate-halo contamination).
func TestDestOnlyScope(t *testing.T) {
	// 4 cells on a line: A(0) - M(1) - F(2) - X(3). M is the motif.
	// M's matter J points +x (θ gradient), so fore pick = F, aft = A.
	L := 16.0
	mk := func() Frame {
		fr := Frame{T: 0, L: L, HasPromo: true, Cells: make([]Cell, 4)}
		for i := range fr.Cells {
			fr.Cells[i] = Cell{X: 5 + float64(i), Y: 8, Z: 8, R: 0.8, Em: 1, Es: 1}
		}
		// θ decreasing in +x ⇒ sin(θj−θi) < 0 toward +x… choose the
		// sign so û·Ĵ selects F: θ(x) = +0.5x ⇒ J_MF = 2·sin(+0.5) > 0.
		for i := range fr.Cells {
			fr.Cells[i].Th2 = 0.5 * float64(i)
		}
		fr.Links = []Link{{I: 0, J: 1, A: 1, D: 1}, {I: 1, J: 2, A: 1, D: 1}, {I: 2, J: 3, A: 1, D: 1}}
		return fr
	}
	prev, cur := mk(), mk()
	cur.T = 4
	cur.Cells[2].PromoI = 0.7 // book on the fore pick (F)
	cur.Cells[0].PromoO = 0.3 // book on the aft pick (A)
	cur.Cells[0].PromoI = 9.9 // contamination: intake booked on the AFT cell
	// (dest scope counts only O there; dilate swallows its I too)
	m := Motif{ID: 1, N: 1, Cells: []int{1}, Em: 1, Cx: 6, Cy: 8, Cz: 8}

	cv := MeasureConvert(prev, cur, m, 1, 0, 0, 0, true)
	if math.Abs(cv.I-0.7) > 1e-12 || math.Abs(cv.O-0.3) > 1e-12 {
		t.Fatalf("dest-only I=%.3g O=%.3g, want 0.7/0.3 (X excluded)", cv.I, cv.O)
	}
	// dilate-1 by contrast swallows X's book too
	cd := MeasureConvert(prev, cur, m, 1, 0, 0, 1, false)
	if cd.I < 0.7+9.9-1e-9 {
		t.Fatalf("dilate control should include contamination, I=%.3g", cd.I)
	}
}
