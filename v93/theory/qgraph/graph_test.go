package qgraph

import (
	"math"
	"testing"
)

func TestLaplacianRowSumsVanish(t *testing.T) {
	for _, g := range []*Graph{NewCycle(6), NewPath(5), NewHypercube(3), NewIcosahedron()} {
		L := g.Laplacian()
		for i := range L {
			s := 0.0
			for j := range L[i] {
				s += L[i][j]
			}
			if math.Abs(s) > 1e-12 {
				t.Fatalf("row %d sum %g (Laplacian must annihilate the constant)", i, s)
			}
		}
	}
}

func TestCycleRankAndForest(t *testing.T) {
	if b := NewCycle(6).CycleRank(); b != 1 {
		t.Fatalf("C6 b1 = %d, want 1", b)
	}
	if b := NewPath(9).CycleRank(); b != 0 {
		t.Fatalf("P9 b1 = %d, want 0", b)
	}
	if !NewPath(9).IsForest() {
		t.Fatal("a path must be a forest")
	}
	if NewCycle(6).IsForest() {
		t.Fatal("a cycle is not a forest")
	}
	// K4: 6 edges, 4 vertices, 1 component -> b1 = 3
	if b := NewComplete(4).CycleRank(); b != 3 {
		t.Fatalf("K4 b1 = %d, want 3", b)
	}
}

// A tree exterior can carry no circulating harmonic at all: the path spectrum
// is non-degenerate, so every eigenspace is one-dimensional, so every mode is
// real up to a global phase and has zero link current. Circulation — and
// therefore charge — needs cycles.
func TestTreeSpectrumIsNonDegenerate(t *testing.T) {
	sh := Shells(Spectrum(NewPath(9)), 1e-9)
	for _, s := range sh {
		if s.Cap != 1 {
			t.Fatalf("path shell %d has capacity %d; a tree must have only singlets", s.K, s.Cap)
		}
	}
	// the ring, by contrast, is built out of counter-rotating doublets
	sh = Shells(Spectrum(NewCycle(9)), 1e-9)
	doublets := 0
	for _, s := range sh {
		if s.Cap == 2 {
			doublets++
		}
	}
	if doublets != 4 {
		t.Fatalf("C9 has %d doublets, want 4", doublets)
	}
}
