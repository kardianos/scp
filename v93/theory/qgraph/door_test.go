package qgraph

import (
	"math"
	"testing"
)

// SELECTION RULES, DERIVED. The door profile IS the selection rule: its
// Fourier transform at the winding difference is the matrix element.
//
//	uniform door -> only Delta m = 0   (radiation cannot touch circulation)
//	wound door   -> exactly Delta m = k (organized field pumps circulation)
//	local door   -> every Delta m      (a point contact scrambles winding)
func TestUniformDoorCannotChangeWinding(t *testing.T) {
	n := 12
	got := AllowedTransfers(n, UniformDoor(n), 1e-9)
	if len(got) != 1 || got[0] != 0 {
		t.Fatalf("uniform door allows %v, want only {0}", got)
	}
	mean, variance := CirculationTransfer(n, UniformDoor(n))
	if math.Abs(mean) > 1e-12 || math.Abs(variance) > 1e-12 {
		t.Fatalf("uniform door mean=%g var=%g; it must not move circulation", mean, variance)
	}
}

func TestWoundDoorPumpsExactlyKTurns(t *testing.T) {
	n := 12
	for _, k := range []int{1, 2, 3, 5} {
		got := AllowedTransfers(n, WoundDoor(n, k), 1e-9)
		if len(got) != 1 || got[0] != k%n {
			t.Fatalf("wound door k=%d allows %v, want only {%d}", k, got, k%n)
		}
		mean, variance := CirculationTransfer(n, WoundDoor(n, k))
		wantMean := float64(k)
		if 2*k > n {
			wantMean = float64(k - n)
		}
		if math.Abs(mean-wantMean) > 1e-9 {
			t.Fatalf("wound door k=%d mean transfer %g, want %g", k, mean, wantMean)
		}
		if math.Abs(variance) > 1e-9 {
			t.Fatalf("wound door k=%d variance %g, want 0 (it is a clean pump)", k, variance)
		}
	}
}

// This is the derivation of Face B. A point-contact door is the operator with
// FLAT support over every winding transfer: it moves circulation on every
// single event, with no preferred sign. Mean zero, variance large — the
// signature of dephasing, not of neutrality.
func TestLocalDoorScramblesWinding(t *testing.T) {
	n := 12
	got := AllowedTransfers(n, LocalDoor(n, 0), 1e-9)
	if len(got) != n {
		t.Fatalf("local door allows %d transfers, want all %d", len(got), n)
	}
	mean, variance := CirculationTransfer(n, LocalDoor(n, 0))
	// No preferred sign at all: the mean is exactly zero, while the spread
	// is large. That pair — zero mean, large variance — IS the dephasing
	// signature, and it is why a point-contact door cannot be read as
	// "neutral" just because it transfers no net circulation.
	if math.Abs(mean) > 1e-12 {
		t.Fatalf("local door mean %g, want 0", mean)
	}
	if variance < 1 {
		t.Fatalf("local door variance %g; it must be large (this is the scramble)", variance)
	}
	// contrast, on the same ring and the same normalisation
	_, uvar := CirculationTransfer(n, UniformDoor(n))
	if !(variance > 100*uvar+1) {
		t.Fatalf("local variance %g is not decisively above uniform %g", variance, uvar)
	}
}

// Matrix elements are the same object on a general graph; a uniform door is
// diagonal there too, because the harmonics are orthonormal.
func TestUniformDoorIsDiagonalOnAnyGraph(t *testing.T) {
	modes := Spectrum(NewIcosahedron())
	d := UniformDoor(12)
	for a := range modes {
		for b := range modes {
			m := GraphMatrixElement(modes, a, b, d)
			if a == b && math.Abs(m-1) > 1e-8 {
				t.Fatalf("diagonal element %g, want 1", m)
			}
			if a != b && math.Abs(m) > 1e-8 {
				t.Fatalf("off-diagonal element %g at (%d,%d), want 0", m, a, b)
			}
		}
	}
}
