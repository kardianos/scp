package qgraph

import (
	"math"
	"testing"
)

// MONODISPERSITY. Cold, the mass distribution collapses onto a magic number —
// a species with one size. Hot, it smears across every occupancy — a
// hierarchy with no line, which is what Face C measured.
func TestColdDistributionCollapsesOntoAMagicNumber(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewCycle(10)), 1e-8), math.Inf(1))
	const a0 = 2 * math.Pi // eps = omega
	mu := 1.7

	pCold := MassDistribution(w, 10, mu, 0.05, a0)
	nCold, probCold := PeakOccupancy(pCold)
	if !w.IsMagic(nCold) {
		t.Fatalf("cold peak at n=%d, which is not magic (magic=%v)", nCold, w.MagicNumbers())
	}
	if probCold < 0.5 {
		t.Fatalf("cold peak probability %.3f; a species must dominate", probCold)
	}

	pHot := MassDistribution(w, 10, mu, 1.5, a0)
	_, probHot := PeakOccupancy(pHot)
	if probHot > 0.30 {
		t.Fatalf("hot peak probability %.3f; a hot stack must be a hierarchy", probHot)
	}
	if probHot >= probCold {
		t.Fatal("heating must broaden the mass distribution")
	}
}

// The sharpness of a species is Delta/T and nothing else. This is the formal
// reason a soft outtake ceiling can never manufacture a line.
func TestSpeciesSharpnessIsTheGapOverTemperature(t *testing.T) {
	prev := 0.0
	for _, r := range []float64{1, 2, 4, 8, 16} {
		tau := ResidenceTime(r, 1)
		if tau <= prev {
			t.Fatalf("residence time not increasing with gap: %g", tau)
		}
		prev = tau
	}
	if math.Abs(ResidenceTime(16, 1)-math.Exp(8)/2) > 1e-6 {
		t.Fatal("residence time must be exp(Delta/2T)/2")
	}
	if s := Sharpness(2, 0.1); math.Abs(s-20) > 1e-12 {
		t.Fatalf("sharpness %g want 20", s)
	}
}

// A larger shell gap must make the species sharper at fixed temperature.
// The hypercube ladder (gaps between sqrt(2k)) provides the comparison.
func TestLargerGapGivesSharperSpecies(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewHypercube(4)), 1e-8), math.Inf(1))
	const a0 = 2 * math.Pi
	g0 := w.Gap(0, a0)
	g2 := w.Gap(2, a0)
	if !(g0 > g2) {
		t.Fatalf("Q4 gaps: g0=%g g2=%g (sqrt spacing must decrease)", g0, g2)
	}
	if ResidenceTime(g0, 0.2) <= ResidenceTime(g2, 0.2) {
		t.Fatal("the wider gap must hold longer")
	}
}

// TEMPERATURE IS AN ATTRIBUTE OF THE BOUND STACK, recoverable from the
// stack's own energy at fixed quantum count. A field with energy but no bound
// stack has no temperature at all — there is nothing to fit.
func TestStackTemperatureRoundTrips(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewIcosahedron()), 1e-6), math.Inf(1))
	const a0 = 2 * math.Pi
	N := 6
	for _, T := range []float64{0.05, 0.2, 0.7, 2.0, 6.0} {
		E := StackEnergy(w, N, T, a0)
		got := StackTemperature(w, E, N, a0)
		if math.IsNaN(got) {
			t.Fatalf("T=%g: inversion failed", T)
		}
		if rel := math.Abs(got-T) / T; rel > 1e-3 {
			t.Fatalf("T=%g recovered as %g (rel %g)", T, got, rel)
		}
	}
}

// Energy at fixed quantum count must be monotone increasing in temperature —
// otherwise "temperature" would not be an attribute at all.
func TestStackEnergyIsMonotoneInTemperature(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewIcosahedron()), 1e-6), math.Inf(1))
	const a0 = 2 * math.Pi
	prev := -1.0
	for _, T := range []float64{0.01, 0.1, 0.5, 2, 10, 100} {
		e := StackEnergy(w, 6, T, a0)
		if e <= prev {
			t.Fatalf("E(T=%g)=%g did not exceed %g", T, e, prev)
		}
		prev = e
	}
}

// The identification x ~ T (which is what makes the law's O ~ x^4 mean
// O ~ T^4) holds ONLY when every mode of the stack is live. Em/T must flatten
// at high T and fall at low T — the freeze-out that softens the exponent.
func TestLoadIsLinearInTOnlyWhenSaturated(t *testing.T) {
	w := NewWell(Shells(Spectrum(NewHypercube(4)), 1e-8), math.Inf(1))
	const a0 = 2 * math.Pi
	Ts := []float64{0.02, 0.1, 0.5, 2, 10, 50}
	r := LoadIsLinearInT(w, 8, a0, Ts)

	if !(r[len(r)-1] < r[len(r)-2]) && !(math.Abs(r[len(r)-1]-r[len(r)-2])/r[len(r)-1] < 0.15) {
		t.Fatalf("Em/T did not settle at high T: %v", r)
	}
	// at the cold end the stack is frozen: Em/T is far larger than at the hot end
	if !(r[0] > 5*r[len(r)-1]) {
		t.Fatalf("no freeze-out signature in Em/T: %v", r)
	}
}
