package qgraph

import "math"

// ---------------------------------------------------------------------------
// The dual I/O law.
//
// ONE coupling A per ordered mode pair drives BOTH directions:
//
//	O(a->b) = A * n_a * (Nu - n_b) * (1 + N_w)     emission   (spont + stim)
//	I(b->a) = A * n_b * (Nu - n_a) * N_w           absorption
//
// with w = eps_a - eps_b and N_w the field occupancy at that pitch. The (Nu-n)
// factors are cap refusal; the (1+N_w) vs N_w asymmetry is the whole content
// of "spontaneous". detailed_balance.mac proves the residual is identically
// zero at a common temperature — emission and absorption are one law read
// forwards and backwards, not two channels that happen to be tuned.
// ---------------------------------------------------------------------------

// Planck is the field-mode occupancy at pitch omega and temperature T.
func Planck(omega, T float64) float64 {
	if T <= 0 || omega <= 0 {
		return 0
	}
	return 1 / math.Expm1(omega/T)
}

// FermiOcc is the matter-mode occupancy under cap refusal (Nu = 1).
func FermiOcc(eps, mu, T float64) float64 {
	if T <= 0 {
		if eps < mu {
			return 1
		}
		return 0
	}
	return 1 / (math.Exp((eps-mu)/T) + 1)
}

// EmitRate is O(a->b) with the coupling A factored in.
func EmitRate(A, na, nb, Nw float64) float64 { return A * na * (Nu - nb) * (1 + Nw) }

// AbsorbRate is I(b->a) with the SAME coupling A.
func AbsorbRate(A, na, nb, Nw float64) float64 { return A * nb * (Nu - na) * Nw }

// DetailedBalance returns O(a->b) - I(b->a) for occupancies drawn from a
// single temperature T and chemical potential mu. It is zero to machine
// precision — this is the identity the whole law rests on.
func DetailedBalance(A, epsA, epsB, mu, T float64) float64 {
	na := FermiOcc(epsA, mu, T)
	nb := FermiOcc(epsB, mu, T)
	Nw := Planck(epsA-epsB, T)
	return EmitRate(A, na, nb, Nw) - AbsorbRate(A, na, nb, Nw)
}

// ---------------------------------------------------------------------------
// Coarse-grained form — the version that maps onto laws_V3a.
//
//	O = k * cap * x^p                       (adopted radiance, p = p_rad)
//	I = k * cap * (xb^p + g * x^(p-1) * xb) (its dual, SAME k and p)
//
// The first intake term is the spontaneous dual of radiance; the second is the
// stimulated cross term, nonzero only when the incoming field is phase-matched
// to what the object already carries. g = 0 is a random-phase bath.
// ---------------------------------------------------------------------------

// Outtake is the adopted radiance channel O = k*cap*x^p.
func Outtake(k, cap, x, p float64) float64 {
	if cap <= 0 || x <= 0 {
		return 0
	}
	return k * cap * math.Pow(x, p)
}

// Intake is the dual of Outtake: same k, same p, plus the coherent cross term.
func Intake(k, cap, xb, x, p, g float64) float64 {
	if cap <= 0 || xb <= 0 {
		return 0
	}
	i := math.Pow(xb, p)
	if g != 0 && x > 0 {
		i += g * math.Pow(x, p-1) * xb
	}
	return k * cap * i
}

// Balance is dM/dt = I - O at load x against a bath at load xb.
func Balance(k, cap, xb, x, p, g float64) float64 {
	return Intake(k, cap, xb, x, p, g) - Outtake(k, cap, x, p)
}

// FixedPointLoad solves I = O for x by bisection on the load ratio.
//
// With g = 0 the root is exactly x = xb: matter merely thermalizes to its bath
// and can never sit above it — which is precisely the measured "box-fill"
// null (MECHANISM.md S08/S12/S19: I ~ O at filled capacity is not a particle).
// With g > 0 the root moves above xb: a coherent bath holds matter up.
func FixedPointLoad(k, cap, xb, p, g float64) float64 {
	f := func(x float64) float64 { return Balance(k, cap, xb, x, p, g) }
	lo, hi := 1e-12, math.Max(xb, 1)*4
	for f(hi) > 0 && hi < 1e9 {
		hi *= 2
	}
	for i := 0; i < 200; i++ {
		mid := 0.5 * (lo + hi)
		if f(mid) > 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// Throughput is the flux at balance: Phi = O(x*) = I(x*).
// It is strictly positive whenever xb > 0, which is the formal statement that
// a balanced object RADIATES. I = O is a live cycle, never a sealed bag.
func Throughput(k, cap, xb, p, g float64) float64 {
	return Outtake(k, cap, FixedPointLoad(k, cap, xb, p, g), p)
}

// Restoring reports d(I-O)/dx at the fixed point; negative means the mass
// fixed point is attracting.
func Restoring(k, cap, xb, p, g float64) float64 {
	x := FixedPointLoad(k, cap, xb, p, g)
	h := 1e-6 * math.Max(x, 1)
	return (Balance(k, cap, xb, x+h, p, g) - Balance(k, cap, xb, x-h, p, g)) / (2 * h)
}

// ---------------------------------------------------------------------------
// Stefan exponent from spectral dimension.
//
//	u(T) = C(d_s) * T^(d_s+1),  C(d_s) = Gamma(d_s+1) * zeta(d_s+1)
//	j(T) = (c/4) * u(T)         => p_rad = d_s + 1
//
// so p_rad = 4 is not an axiom: it is the assertion that the FIELD sector has
// three transport dimensions. On a graph, d_s is measurable (SpectralDimension).
// ---------------------------------------------------------------------------

// StefanExponent returns p = d_s + 1.
func StefanExponent(ds float64) float64 { return ds + 1 }

// StefanConstant returns Gamma(ds+1)*zeta(ds+1); for ds = 3 it is pi^4/15.
func StefanConstant(ds float64) float64 { return math.Gamma(ds+1) * zeta(ds+1) }

// StefanFlux is the emitted flux j ~ C(ds) * T^(ds+1) (units absorbed into c).
func StefanFlux(ds, T float64) float64 {
	if T <= 0 {
		return 0
	}
	return StefanConstant(ds) * math.Pow(T, ds+1) / 4
}

// EffectiveRadianceExponent is the SMALL-LOAD radiance exponent.
//
// O ~ T^4 always (field side, d_s = 3). The law writes O ~ x^4, which requires
// x ~ T — true only when every internal mode of the object is live
// (Dulong-Petit). When the stack freezes out, Em ~ T^(dm+1), so
//
//	T ~ x^(1/(dm+1))  =>  O ~ x^(4/(dm+1)).
//
// The adopted p_rad = 4 is therefore the SATURATED limit only; small hoards
// must radiate with a much softer exponent. This is a falsifiable deviation
// from the adopted table (SIM_PREDICTIONS P3).
func EffectiveRadianceExponent(dsField, dmMatter float64) float64 {
	return (dsField + 1) / (dmMatter + 1)
}

// zeta is the Riemann zeta function for real s > 1, by Euler-Maclaurin.
func zeta(s float64) float64 {
	if s <= 1 {
		return math.Inf(1)
	}
	const nTerm = 40
	sum := 0.0
	for n := 1; n < nTerm; n++ {
		sum += math.Pow(float64(n), -s)
	}
	n := float64(nTerm)
	sum += math.Pow(n, 1-s) / (s - 1)
	sum += 0.5 * math.Pow(n, -s)
	sum += s * math.Pow(n, -s-1) / 12
	sum -= s * (s + 1) * (s + 2) * math.Pow(n, -s-3) / 720
	return sum
}
