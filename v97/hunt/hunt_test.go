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
	cv := MeasureConvert(prev, cur, ms[0], 1, 0, 0)
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
