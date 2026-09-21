package qgraph

import (
	"math"
	"testing"
)

func icoWell() Well { return NewWell(Shells(Spectrum(NewIcosahedron()), 1e-6), math.Inf(1)) }

// Q = 0 exactly at the magic numbers, and nowhere else.
func TestChargeVanishesOnlyAtClosure(t *testing.T) {
	w := icoWell() // caps 1,3,5,3 -> magic 1,4,9,12
	for n := 0; n <= 12; n++ {
		o := NewObject(w, n)
		neutral := o.IsNeutral()
		magic := n == 12 // Q = Z - N = 12 - N, so only full occupancy is neutral
		if neutral != magic {
			t.Fatalf("n=%d neutral=%v want %v (Q=%d)", n, neutral, magic, o.Charge())
		}
	}
}

// The chemically meaningful statement: with the well truncated at a binding
// threshold, Q = 0 exactly when the bound stack is closed, and Q counts the
// deficit otherwise.
func TestChargeIsTheShellDeficit(t *testing.T) {
	sh := Shells(Spectrum(NewIcosahedron()), 1e-6)
	w := NewWell(sh, sh[2].Omega+1e-9) // bind shells 0,1,2 -> Z = 1+3+5 = 9
	if w.Z() != 9 {
		t.Fatalf("Z=%d want 9", w.Z())
	}
	for n := 0; n <= 9; n++ {
		if q := NewObject(w, n).Charge(); q != 9-n {
			t.Fatalf("n=%d Q=%d want %d", n, q, 9-n)
		}
	}
	// one quantum held above closure is the opposite sign
	if q := NewObject(w, 10).Charge(); q != -1 {
		t.Fatalf("over-filled Q=%d want -1", q)
	}
}

// CONSERVATION THEOREM. Emission, absorption and internal promotion all leave
// the quantum count N invariant, so they leave Q invariant — exactly.
func TestPromotionConservesCharge(t *testing.T) {
	w := icoWell()
	o := NewObject(w, 6)
	q0, n0 := o.Charge(), o.N()
	moves := 0
	for _, mv := range [][2]int{{1, 2}, {2, 3}, {3, 1}, {0, 3}, {2, 0}} {
		if _, ok := o.Promote(mv[0], mv[1]); ok {
			moves++
		}
		if o.N() != n0 {
			t.Fatalf("promotion changed N: %d -> %d", n0, o.N())
		}
		if o.Charge() != q0 {
			t.Fatalf("promotion changed Q: %d -> %d", q0, o.Charge())
		}
	}
	if moves == 0 {
		t.Fatal("no promotion succeeded; test is vacuous")
	}
}

// A promotion returns the pitch of the field quantum it exchanges. Emission
// and absorption are the SAME operation with opposite sign — there is no
// second law for intake.
//
// Note the first assertion: an aufbau GROUND state cannot emit at all. Every
// downward move is blocked by cap refusal. That is the law's own statement
// that a cold closed object is dark, derived rather than asserted.
func TestEmissionAndAbsorptionAreOneOperation(t *testing.T) {
	w := icoWell()
	ground := NewObject(w, 5) // occ = [1,3,1,0]
	for a := range ground.Occ {
		for b := 0; b < a; b++ {
			if _, ok := ground.Promote(a, b); ok {
				t.Fatalf("ground state emitted via %d->%d; cap refusal must block it", a, b)
			}
		}
	}

	o := NewObject(w, 5)
	up, ok := o.Promote(2, 3) // climb: absorb
	if !ok {
		t.Fatal("upward promotion blocked")
	}
	if up >= 0 {
		t.Fatalf("a climbing promotion must ABSORB (negative pitch), got %g", up)
	}
	down, ok := o.Promote(3, 2) // fall: emit the same quantum back
	if !ok {
		t.Fatal("downward promotion blocked after excitation")
	}
	if down <= 0 {
		t.Fatalf("a falling promotion must EMIT (positive pitch), got %g", down)
	}
	if math.Abs(down+up) > 1e-12 {
		t.Fatalf("absorbed %g and emitted %g are not time-reverses", up, down)
	}
}

// Charge moves between objects and is conserved in total: ionization of one
// object is capture by another.
func TestTransferConservesTotalCharge(t *testing.T) {
	w := icoWell()
	a, b := NewObject(w, 8), NewObject(w, 3)
	q0 := TotalCharge(a, b)
	for i := 0; i < 4; i++ {
		if !Transfer(a, b) {
			t.Fatalf("transfer %d blocked", i)
		}
		if q := TotalCharge(a, b); q != q0 {
			t.Fatalf("total charge drifted: %d -> %d", q0, q)
		}
	}
	if a.Charge() == 0 || b.Charge() == 0 {
		t.Fatal("after transfer both objects should be charged")
	}
	if a.Charge()+b.Charge() != q0 {
		t.Fatal("charges do not sum back")
	}
}

// A tree can host no charge: its spectrum is non-degenerate, every shell has
// capacity 1, so every occupancy closes a shell and every Q is a pure count
// with no circulating partner. Charge needs cycles.
func TestTreeHostsNoCirculatingCharge(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewPath(8)), 1e-9), math.Inf(1))
	for _, s := range w.Shells {
		if s.Cap != 1 {
			t.Fatalf("tree shell capacity %d", s.Cap)
		}
	}
	// every occupancy leaves the valence shell CLOSED: there is no
	// half-filled multiplet anywhere, hence no residual circulation.
	for n := 1; n <= 8; n++ {
		if w.Incomplete(w.Aufbau(n)) {
			t.Fatalf("n=%d: a tree produced an incomplete shell", n)
		}
	}
}

// SHELL NEUTRALITY (shells.mac §S2, proved identically there; checked
// numerically here): a full +-m doublet carries zero net circulation, and a
// half-filled one carries exactly +-sin(2*pi*m/N).
func TestFullShellIsNeutralHalfShellIsCharged(t *testing.T) {
	for _, n := range []int{5, 6, 7, 8, 9, 12} {
		for _, sh := range RingShells(n) {
			full := make([]int, len(sh))
			for i := range full {
				full[i] = 1
			}
			if c := ShellCirculation(sh, full); math.Abs(c) > 1e-12 {
				t.Fatalf("C%d full shell circulation %g, want 0", n, c)
			}
			if wnd := ShellWinding(sh, full); wnd != 0 && !isSelfConjugate(sh, n) {
				t.Fatalf("C%d full shell winding %d, want 0", n, wnd)
			}
			if len(sh) == 2 { // half-fill the doublet
				half := []int{1, 0}
				if c := ShellCirculation(sh, half); math.Abs(c) < 1e-9 {
					t.Fatalf("C%d half shell circulation %g, want nonzero", n, c)
				}
				if wnd := ShellWinding(sh, half); wnd == 0 {
					t.Fatalf("C%d half shell winding 0, want nonzero", n)
				}
			}
		}
	}
}

// The zone-boundary singlet m = N/2 is a real standing mode: its winding label
// is ambiguous mod N, but its physical circulation is exactly zero, so it is
// neutral like m = 0.
func TestZoneBoundaryModeIsNeutral(t *testing.T) {
	for _, n := range []int{6, 8, 10, 12} {
		modes := RingModes(n)
		found := false
		for _, m := range modes {
			if 2*abs(m.M) == n {
				found = true
				if math.Abs(m.Circ) > 1e-12 {
					t.Fatalf("C%d zone-boundary mode circulation %g, want 0", n, m.Circ)
				}
			}
		}
		if !found {
			t.Fatalf("C%d has no zone-boundary mode", n)
		}
	}
}

func isSelfConjugate(sh []RingMode, n int) bool {
	return len(sh) == 1 && (sh[0].M == 0 || 2*abs(sh[0].M) == n)
}

func abs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}
