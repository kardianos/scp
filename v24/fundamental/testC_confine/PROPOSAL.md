# Test C: Confining Potential — Does Linear V(P) Stop Leaking?

## Thesis

Replace the saturating potential V = (μ/2)P²/(1+κP²) with a CONFINING
potential V = -σ|P| (linear in |P|, like the QCD string tension). The
linear potential has no saturation — the cost of separating the fields
grows without bound.

## Setup

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a + σ·sign(P)·∂|P|/∂φ_a

The force from V = -σ|P|:
    ∂V/∂φ_a = -σ·sign(P)·∂P/∂φ_a

where ∂P/∂φ₁ = φ₂φ₃, etc.

Note: |P| has a cusp at P=0, making the force discontinuous.
Regularize: V = -σ·√(P² + ε²) with small ε = 10⁻⁶.

    ∂V/∂φ_a = -σ·P·(∂P/∂φ_a)/√(P² + ε²)

## What to Test

1. Initialize Gaussian oscillon (A=0.8, σ_init=3.0) with the linear potential
2. Scan string tension σ = {1.0, 5.0, 10.0, 20.0, 50.0}
3. Evolve for t=20000 (long time to test true stability)
4. Measure: E(t), dE/dt, fc, ω, peak amplitude
5. KEY: does dE/dt → 0? Is the oscillon PERFECTLY stable?
6. Compare with the standard saturating potential at the same parameters

Also test: V = -σ|P| + (κ_conf/2)P⁴ (confining + quartic stabilization)
to prevent the amplitude from growing without bound.

## Physics Rationale

The saturating potential V = (μ/2)P²/(1+κP²) reaches maximum |V_max| = |μ|/(2κ)
at large P. Beyond this, increasing P gives no more binding → the fields can
"leak" by reducing P without energy cost.

The linear potential V = -σ|P| has NO maximum — every bit of P gives more
binding. The fields CANNOT reduce P without paying energy → confinement.

## Reference Code

- v21/src/triad1d.c

## Output

- `src/confine.c`, `data/`, `RESULTS.md`

## Parameters

m=1.0, σ scan: {1.0, 5.0, 10.0, 20.0, 50.0}
ε = 1e-6 (regularization)
Nx=4000, xmax=100, t=20000

Compile: `gcc -O3 -Wall -o confine src/confine.c -lm`
