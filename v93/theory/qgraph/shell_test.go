package qgraph

import (
	"math"
	"testing"
)

// Shells are degeneracy classes of the Laplacian, so their capacities are
// forced by the graph's symmetry, not chosen. The hypercube Q_d has spectrum
// 2k with multiplicity C(d,k) — an exact, checkable shell ladder.
func TestHypercubeShellCapacitiesAreBinomial(t *testing.T) {
	for d := 2; d <= 5; d++ {
		sh := Shells(Spectrum(NewHypercube(d)), 1e-8)
		if len(sh) != d+1 {
			t.Fatalf("Q_%d: %d shells, want %d", d, len(sh), d+1)
		}
		for k, s := range sh {
			want := binom(d, k)
			if s.Cap != want {
				t.Fatalf("Q_%d shell %d: cap %d, want C(%d,%d)=%d", d, k, s.Cap, d, k, want)
			}
			if math.Abs(s.Lambda-float64(2*k)) > 1e-7 {
				t.Fatalf("Q_%d shell %d: lambda %g, want %d", d, k, s.Lambda, 2*k)
			}
		}
	}
}

// The icosahedron gives 1,3,5,3 — an s/p/d-like ladder with no orbital
// imported anywhere. This is the sharpest demonstration that "shells" here
// are graph symmetry and nothing else.
func TestIcosahedronShellLadder(t *testing.T) {
	sh := Shells(Spectrum(NewIcosahedron()), 1e-6)
	got := make([]int, len(sh))
	for i, s := range sh {
		got[i] = s.Cap
	}
	want := []int{1, 3, 5, 3}
	if len(got) != len(want) {
		t.Fatalf("icosahedron shells %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("icosahedron shells %v, want %v", got, want)
		}
	}
}

// Ring shells are the +-m doublets plus the m=0 (and, for even N, m=N/2)
// singlets: magic numbers 1,3,5,7,... (shells.mac §S4).
func TestRingMagicNumbers(t *testing.T) {
	cases := map[int][]int{
		6:  {1, 3, 5, 6},
		7:  {1, 3, 5, 7},
		10: {1, 3, 5, 7, 9, 10},
	}
	for n, want := range cases {
		sh := Shells(Spectrum(NewCycle(n)), 1e-8)
		w := NewWell(sh, math.Inf(1))
		got := w.MagicNumbers()
		if len(got) != len(want) {
			t.Fatalf("C%d magic %v, want %v", n, got, want)
		}
		for i := range want {
			if got[i] != want[i] {
				t.Fatalf("C%d magic %v, want %v", n, got, want)
			}
		}
	}
}

// Aufbau must fill exactly, must never exceed a capacity, and must reproduce
// the requested count.
func TestAufbauRespectsCapacity(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewIcosahedron()), 1e-6), math.Inf(1))
	for n := 0; n <= 12; n++ {
		occ := w.Aufbau(n)
		tot := 0
		for k, v := range occ {
			if v > w.Shells[k].Cap {
				t.Fatalf("n=%d shell %d over cap: %d > %d", n, k, v, w.Shells[k].Cap)
			}
			if v < 0 {
				t.Fatalf("n=%d negative occupancy", n)
			}
			tot += v
		}
		if tot != n {
			t.Fatalf("aufbau(%d) placed %d quanta", n, tot)
		}
	}
}

// The binding threshold truncates the stack: Z is a property of the well, not
// of what happens to be in it.
func TestWellBindingThresholdSetsZ(t *testing.T) {
	sh := Shells(Spectrum(NewHypercube(3)), 1e-8) // caps 1,3,3,1; omegas sqrt(0,2,4,6)
	full := NewWell(sh, math.Inf(1))
	if full.Z() != 8 {
		t.Fatalf("full Q3 well Z=%d, want 8", full.Z())
	}
	// bind only the two lowest shells
	part := NewWell(sh, math.Sqrt(2)+1e-9)
	if part.NBound != 2 || part.Z() != 4 {
		t.Fatalf("truncated well NBound=%d Z=%d, want 2 and 4", part.NBound, part.Z())
	}
}

func binom(n, k int) int {
	if k < 0 || k > n {
		return 0
	}
	r := 1
	for i := 0; i < k; i++ {
		r = r * (n - i) / (i + 1)
	}
	return r
}
