package qgraph

import (
	"math"
	"testing"
)

// THE DUALITY IDENTITY. One coupling A drives both directions; at a common
// temperature the residual is zero to machine precision. Emission and
// absorption are one law read forwards and backwards.
func TestDetailedBalanceIsExact(t *testing.T) {
	for _, T := range []float64{0.05, 0.3, 1, 7} {
		for _, mu := range []float64{-1, 0, 0.4, 2} {
			for _, pair := range [][2]float64{{1, 0.5}, {2.4, 0.1}, {0.9, 0.87}} {
				r := DetailedBalance(1.7, pair[0], pair[1], mu, T)
				if math.Abs(r) > 1e-12 {
					t.Fatalf("T=%g mu=%g eps=%v: residual %g", T, mu, pair, r)
				}
			}
		}
	}
}

// Detailed balance must FAIL when the two sectors are at different
// temperatures — otherwise the identity above would be vacuous.
func TestDetailedBalanceBreaksOffEquilibrium(t *testing.T) {
	na := FermiOcc(1.0, 0.4, 0.5)
	nb := FermiOcc(0.3, 0.4, 0.5)
	Nw := Planck(0.7, 2.0) // field hotter than matter
	r := EmitRate(1, na, nb, Nw) - AbsorbRate(1, na, nb, Nw)
	if math.Abs(r) < 1e-6 {
		t.Fatalf("mismatched temperatures gave residual %g; identity is vacuous", r)
	}
	if r > 0 {
		t.Fatalf("a hotter bath must drive net ABSORPTION, got net emission %g", r)
	}
}

// With a random-phase bath (g=0) the fixed point is exactly the bath load:
// matter thermalizes and can NEVER sit above its surroundings. That is the
// box-fill null, derived.
func TestRandomBathFixedPointIsTheBathItself(t *testing.T) {
	for _, xb := range []float64{0.2, 0.62, 1.0, 2.5} {
		x := FixedPointLoad(0.05, 2.5, xb, 4, 0)
		if math.Abs(x-xb) > 1e-6 {
			t.Fatalf("xb=%g gave x*=%g, want x*=xb", xb, x)
		}
	}
}

// With a coherent (phase-matched) bath the fixed point rises above the bath:
// u* solves u^4 - g u^3 - 1 = 0 (detailed_balance.mac §CLAIM 3).
func TestCoherentBathHoldsMatterAboveAmbient(t *testing.T) {
	xb := 0.62
	want := map[float64]float64{0: 1.0, 1: 1.380277568940073, 2: 2.106919340323657, 3: 3.03574411245063}
	for g, u := range want {
		x := FixedPointLoad(0.05, 2.5, xb, 4, g)
		if math.Abs(x/xb-u) > 1e-6 {
			t.Fatalf("g=%g: u*=%g, want %g", g, x/xb, u)
		}
	}
}

// The balanced state RADIATES: throughput is strictly positive and the fixed
// point is restoring. I = O is a live cycle, not a sealed bag.
func TestBalanceIsALiveCycleAndRestoring(t *testing.T) {
	phi := Throughput(0.05, 2.5, 0.62, 4, 0)
	if phi <= 0 {
		t.Fatalf("throughput %g at balance; a balanced object must still emit", phi)
	}
	if s := Restoring(0.05, 2.5, 0.62, 4, 0); s >= 0 {
		t.Fatalf("d(I-O)/dx = %g >= 0; fixed point is not attracting", s)
	}
}

// STEFAN EXPONENT = SPECTRAL DIMENSION + 1. p_rad = 4 is the claim d_s = 3,
// and d_s is a MEASUREMENT of the contact graph, not a modelling choice.
func TestStefanExponentTracksSpectralDimension(t *testing.T) {
	for d := 1; d <= 3; d++ {
		L := 24
		if d == 3 {
			L = 14
		}
		// the counting estimator carries a known +0.25 finite-size bias
		ds := SpectralDimensionFromEigenvalues(TorusSpectrum(d, L), 0.10)
		if ds < float64(d) || ds > float64(d)+0.5 {
			t.Fatalf("d=%d torus: eigenvalue-counting d_s=%.3f out of band", d, ds)
		}
		if p := StefanExponent(float64(d)); math.Abs(p-float64(d+1)) > 1e-12 {
			t.Fatalf("d=%d: p=%g", d, p)
		}
	}
	if p := StefanExponent(3); math.Abs(p-4) > 1e-12 {
		t.Fatalf("d_s=3 must give p_rad=4, got %g", p)
	}
}

// The random-walk estimator is the one that will run on a real free-cell
// contact graph (no diagonalisation, O(|E| * t) per source).
func TestReturnDimensionAgreesOnTori(t *testing.T) {
	for d := 1; d <= 3; d++ {
		L := 60
		if d == 2 {
			L = 30
		}
		if d == 3 {
			L = 14
		}
		g := NewTorus(d, L)
		ds := ReturnDimension(g, 128, []int{0, g.N / 3, 2 * g.N / 3})
		if math.Abs(ds-float64(d)) > 0.05 {
			t.Fatalf("d=%d torus: return-probability d_s=%.3f", d, ds)
		}
	}
}

// The Stefan constant for d_s = 3 is pi^4/15 — the ordinary Stefan-Boltzmann
// number, recovered rather than inserted.
func TestStefanConstantIsPi4Over15(t *testing.T) {
	got := StefanConstant(3)
	want := math.Pow(math.Pi, 4) / 15
	if math.Abs(got-want) > 1e-9 {
		t.Fatalf("C(3)=%g want pi^4/15=%g", got, want)
	}
	if f := StefanFlux(3, 2); math.Abs(f-want*16/4) > 1e-8 {
		t.Fatalf("flux %g", f)
	}
}

// The adopted p_rad = 4 is the SATURATED-load limit. When the object's own
// modes freeze out, the effective exponent softens to 4/(dm+1).
func TestLowLoadRadianceExponentSoftens(t *testing.T) {
	cases := map[float64]float64{1: 2, 2: 4.0 / 3, 3: 1}
	for dm, want := range cases {
		if p := EffectiveRadianceExponent(3, dm); math.Abs(p-want) > 1e-12 {
			t.Fatalf("dm=%g: p_eff=%g want %g", dm, p, want)
		}
	}
	if EffectiveRadianceExponent(3, 3) >= 4 {
		t.Fatal("frozen-stack exponent must be softer than the adopted 4")
	}
}
