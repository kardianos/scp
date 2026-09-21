package size

import (
	"math"
	"testing"

	"scp/v97/hunt"
)

func TestJacobiC6Spectrum(t *testing.T) {
	n := 6
	L := make([]float64, n*n)
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		L[i*n+i] += 1
		L[j*n+j] += 1
		L[i*n+j] -= 1
		L[j*n+i] -= 1
	}
	eval, _, ok := jacobiSym(L, n)
	if !ok {
		t.Fatal("jacobi failed")
	}
	sortFloats(eval)
	want := []float64{0, 1, 1, 3, 3, 4}
	for k := range want {
		if math.Abs(eval[k]-want[k]) > 1e-8 {
			t.Fatalf("eval[%d]=%g want %g (%v)", k, eval[k], want[k], eval)
		}
	}
}

func TestIsLineThresholds(t *testing.T) {
	rep := &Report{
		HistN:   map[int]int{6: 80, 7: 20},
		NFrame:  100,
		HistDom: map[int]int{6: 80, 7: 20},
	}
	// top-bin 0.80 but we need N_eff/N too — empty snaps ⇒ mean 0 ⇒ not a line
	line, _ := isLine(rep)
	if line {
		t.Fatal("empty snaps must not count as a compact line")
	}
	rep.Snaps = make([]FrameSnap, 100)
	for i := range rep.Snaps {
		n := 6
		if i >= 80 {
			n = 7
		}
		rep.Snaps[i] = FrameSnap{Motifs: []MotifSnap{{N: n, NEff: float64(n)}}}
	}
	line, why := isLine(rep)
	if !line {
		t.Fatalf("compact 80%% peak should be a line: %s", why)
	}
}

func TestChartMatch6OnCycle(t *testing.T) {
	fr := hunt.Frame{L: 16}
	fr.Cells = make([]hunt.Cell, 6)
	for i := 0; i < 6; i++ {
		fr.Cells[i] = hunt.Cell{Em: 1, Th2: float64(i) * math.Pi / 3}
		j := (i + 1) % 6
		fr.Links = append(fr.Links, hunt.Link{I: i, J: j})
	}
	sp := analyzeChart(fr, []int{0, 1, 2, 3, 4, 5})
	if !sp.OK || !sp.Match6 {
		t.Fatalf("C6 chart %+v", sp)
	}
}
