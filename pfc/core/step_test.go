package core

import (
	"testing"
	"unsafe"
)

func TestStepEnergyStable(t *testing.T) {
	law := DefaultLaw()
	law.L = 8
	w := BuildLattice(law, 512)
	e0 := w.TotalE()
	for k := 0; k < 20; k++ {
		w.Step()
	}
	e1 := w.TotalE()
	rel := (e1 - e0) / e0
	if rel < 0 {
		rel = -rel
	}
	if rel > 1e-9 {
		t.Fatalf("energy drift %g (e0=%g e1=%g)", rel, e0, e1)
	}
}

func TestEdgeColorValid(t *testing.T) {
	w := BuildLattice(DefaultLaw(), 216)
	color, nColor := w.EdgeColor()
	if nColor < 1 {
		t.Fatal("no colors")
	}
	for i := 0; i < len(w.Lo); i++ {
		seen := make(map[int]bool)
		for q := w.Cls[i]; q < w.Cls[i+1]; q++ {
			si := w.Inc[q]
			c := color[si]
			if c < 0 || c >= nColor {
				t.Fatalf("bad color %d", c)
			}
			if seen[c] {
				t.Fatalf("cell %d color clash %d", i, c)
			}
			seen[c] = true
		}
	}
}

func TestSizes(t *testing.T) {
	var lo CellLo
	var hi CellHi
	t.Logf("sizeof CellLo=%d CellHi=%d Slot=%d",
		unsafe.Sizeof(lo), unsafe.Sizeof(hi), unsafe.Sizeof(Slot{}))
	if unsafe.Sizeof(lo) > 200 {
		t.Fatalf("CellLo larger than expected: %d", unsafe.Sizeof(lo))
	}
}
