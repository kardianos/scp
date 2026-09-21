package calc

import (
	"math"
	"testing"
)

func approx(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6g, want %.6g (±%.2g)", name, got, want, tol)
	}
}

// The selection rule is the DFT of the door profile: closed forms.

func TestVertexDoorIsFlat(t *testing.T) {
	for _, n := range []int{6, 12} {
		d := VertexDoor(n)
		for dm := 0; dm < n; dm++ {
			approx(t, "vertex |M|", DoorAmp(d, dm, 0), 1/float64(n), 1e-12)
		}
	}
}

func TestEdgeSymDoorCosEnvelope(t *testing.T) {
	for _, n := range []int{6, 12} {
		d := EdgeSymDoor(n)
		for dm := 0; dm < n; dm++ {
			want := 2 / float64(n) * math.Abs(math.Cos(math.Pi*float64(dm)/float64(n)))
			approx(t, "edge-sym |M|", DoorAmp(d, dm, 0), want, 1e-12)
		}
	}
}

func TestEdgeDipoleDoorSinEnvelope(t *testing.T) {
	for _, n := range []int{6, 12} {
		d := EdgeDipoleDoor(n)
		for dm := 0; dm < n; dm++ {
			want := 2 / float64(n) * math.Abs(math.Sin(math.Pi*float64(dm)/float64(n)))
			approx(t, "edge-dipole |M|", DoorAmp(d, dm, 0), want, 1e-12)
		}
		// elastic channel closed; max transfer at the zone boundary —
		// the oriented edge door is ANTI-selective in Δm.
		approx(t, "edge-dipole Δm=0", DoorAmp(d, 0, 0), 0, 1e-12)
		if DoorAmp(d, 1, 0) >= DoorAmp(d, n/2, 0) {
			t.Errorf("n=%d: |M(Δm=1)| ≥ |M(Δm=N/2)| — expected rising envelope", n)
		}
	}
}

func TestWoundDoorIsExactPump(t *testing.T) {
	n, k := 12, 2
	d := WoundDoor(n, k)
	for dm := 0; dm < n; dm++ {
		want := 0.0
		if dm == n-k { // M ∝ DFT at Δm = m−m' ⇒ pump lands at Δm ≡ −k (mod N)
			want = 1
		}
		approx(t, "wound |M|", DoorAmp(d, dm, 0), want, 1e-9)
	}
}

func TestCurrentDoorSelectsOnSum(t *testing.T) {
	n := 12
	// |M| = (2/N)|sin(π(m+m')/N)|: zero when m+m' ≡ 0 (mod N),
	// regardless of Δm — sum-selective, not difference-selective.
	approx(t, "J-door m=2,m'=-2", CurrentDoorAmp(n, 2, -2), 0, 1e-12)
	approx(t, "J-door m=2,m'=1",
		CurrentDoorAmp(n, 2, 1), 2/float64(n)*math.Abs(math.Sin(math.Pi*3/float64(n))), 1e-12)
}

// The corrected Φ-bar theorem.

func TestPhiBarPassableAtGZero(t *testing.T) {
	// defaults: k=0.05 p=4 q=1.2 Δt=4, bar 0.015
	x := PhiBarLoadThreshold(0.05, 4, 1.2, 4, 0.015)
	approx(t, "x threshold", x, 0.373, 0.005)
	if x >= 0.62 {
		t.Fatalf("threshold %.3f ≥ xb=0.62 — H1 Φ would need g>0, contradicting the corrected theorem", x)
	}
	// at ambient load the frame ratio clears the bar with margin
	if r := PhiFrameOverEm(0.05, 4, 1.2, 0.62, 4); r < 0.015 {
		t.Fatalf("Φ/Em at xb = %.4g < bar", r)
	}
	// but as a RATE (Δt=1) at ambient it is 0.0208 — and WITHOUT the
	// comp factor it is 0.0119 < 0.015: the review's F2 arithmetic,
	// preserved as the record of what the identification got wrong.
	approx(t, "rate with comp", PhiFrameOverEm(0.05, 4, 1.2, 0.62, 1), 0.0208, 5e-4)
	approx(t, "rate bare", PhiFrameOverEm(0.05, 4, 0, 0.62, 1), 0.0119, 5e-4)
}

func TestBalanceUReproducesThm35(t *testing.T) {
	// qgraph Theorem 3.5 table
	for _, c := range []struct{ g, u float64 }{
		{0, 1.000}, {0.5, 1.153}, {1, 1.380}, {2, 2.107}, {3, 3.036},
	} {
		approx(t, "u*(g)", BalanceU(c.g), c.u, 2e-3)
	}
	// holding 8% above ambient needs g ≈ 0.28 (the H2 threshold)
	approx(t, "g for u=1.08", GForHold(1.08), 0.283, 5e-3)
}
