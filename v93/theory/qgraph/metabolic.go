package qgraph

import "math"

// ---------------------------------------------------------------------------
// The metabolic layer: temperature, monodispersity, and why matter emits
// without evaporating.
//
// The coarse fixed point x* = xb is CONTINUOUS — it is thermalization, not a
// species. Discreteness enters through the shell ladder only: the neutral
// totals are the magic numbers, and a closed configuration survives for
// ~exp(Delta/2T). Monodispersity is therefore a statement about the SHELL GAP
// against the temperature, never about switching radiance off.
// ---------------------------------------------------------------------------

// MassWeight is the grand-canonical weight of an occupancy n:
// P(n) ~ exp(-(E(n) - mu*n)/T).
func MassWeight(w Well, n int, mu, T, a0 float64) float64 {
	if T <= 0 {
		return 0
	}
	e := w.Energy(w.Aufbau(n), a0)
	return math.Exp(-(e - mu*float64(n)) / T)
}

// MassDistribution returns the normalised P(n) for n = 0..nMax.
//
// Cold: it collapses onto a magic number — a SPECIES with one size.
// Hot:  it smears across every n — a HIERARCHY with no line, which is exactly
// what Face C measured (broad, long-tailed hoard spectra).
func MassDistribution(w Well, nMax int, mu, T, a0 float64) []float64 {
	p := make([]float64, nMax+1)
	z := 0.0
	for n := 0; n <= nMax; n++ {
		p[n] = MassWeight(w, n, mu, T, a0)
		z += p[n]
	}
	if z > 0 {
		for n := range p {
			p[n] /= z
		}
	}
	return p
}

// PeakOccupancy returns the most probable n and its probability.
func PeakOccupancy(p []float64) (int, float64) {
	best, bp := 0, -1.0
	for n, v := range p {
		if v > bp {
			best, bp = n, v
		}
	}
	return best, bp
}

// ResidenceTime is the mean dwell at a closed configuration, ~exp(Delta/2T)/2.
//
// This is the ONLY source of sharpness in the law. A soft outtake ceiling can
// erode a hierarchy from above; it can never manufacture a line.
func ResidenceTime(gap, T float64) float64 {
	if T <= 0 {
		return math.Inf(1)
	}
	return math.Exp(gap/(2*T)) / 2
}

// Sharpness is Delta/T, the dimensionless monodispersity number.
// Delta/T >> 1 is a species; Delta/T <~ 1 is a hoard.
func Sharpness(gap, T float64) float64 {
	if T <= 0 {
		return math.Inf(1)
	}
	return gap / T
}

// ---------------------------------------------------------------------------
// Temperature.
//
// T is an attribute of the BOUND stack, not of the field: it is the parameter
// of the maximum-entropy occupancy that reproduces the object's own energy. A
// vacuum with mode energy but no bound stack has no T at all — which is why
// "hot bath" arms box-fill instead of warming a particle.
// ---------------------------------------------------------------------------

// StackEnergy is the equilibrium energy of a bound stack holding N quanta at
// temperature T: the max-entropy occupancy of the object's OWN harmonics.
func StackEnergy(w Well, N int, T, a0 float64) float64 {
	mu := chemicalPotential(w, N, T, a0)
	e := 0.0
	for k := range w.Shells {
		eps := a0 * w.Shells[k].Omega / (2 * math.Pi)
		e += float64(w.Shells[k].Cap) * FermiOcc(eps, mu, T) * eps
	}
	return e
}

// StackTemperature fits T to a bound stack carrying energy E with occupancy N,
// by solving sum_k cap_k * FermiOcc(eps_k, mu, T) = N and matching the energy.
// It returns NaN when the stack cannot host (E,N).
func StackTemperature(w Well, E float64, N int, a0 float64) float64 {
	if N <= 0 || len(w.Shells) == 0 {
		return math.NaN()
	}
	energyAt := func(T float64) float64 { return StackEnergy(w, N, T, a0) }
	lo, hi := 1e-9, 1.0
	for energyAt(hi) < E && hi < 1e9 {
		hi *= 2
	}
	if energyAt(hi) < E {
		return math.NaN()
	}
	for i := 0; i < 200; i++ {
		mid := 0.5 * (lo + hi)
		if energyAt(mid) < E {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// chemicalPotential solves sum_k cap_k * FermiOcc(eps_k, mu, T) = N.
func chemicalPotential(w Well, N int, T, a0 float64) float64 {
	count := func(mu float64) float64 {
		c := 0.0
		for k := range w.Shells {
			eps := a0 * w.Shells[k].Omega / (2 * math.Pi)
			c += float64(w.Shells[k].Cap) * FermiOcc(eps, mu, T)
		}
		return c
	}
	lo, hi := -1e3, 1e3
	for i := 0; i < 200; i++ {
		mid := 0.5 * (lo + hi)
		if count(mid) < float64(N) {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// LoadIsLinearInT reports the regime check behind "O ~ x^4 means O ~ T^4":
// it returns Em(T)/T over a temperature sweep. A flat result is the
// Dulong-Petit (saturated) regime where the identification x ~ T is legitimate;
// a rising result means modes are still freezing out and the adopted quartic
// radiance over-estimates the object's coolness.
func LoadIsLinearInT(w Well, N int, a0 float64, Ts []float64) []float64 {
	out := make([]float64, len(Ts))
	for i, T := range Ts {
		out[i] = StackEnergy(w, N, T, a0) / T
	}
	return out
}
