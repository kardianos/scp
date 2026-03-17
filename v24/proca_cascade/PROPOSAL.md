# V24-P3: Cascade — Proca Sources a Massless Scalar

## Thesis

The Proca mediator from V24-ME has a Yukawa tail ~exp(-m_A r)/r around each
oscillon. Add a massless scalar θ sourced by the STATIC component of the
antisymmetric mode energy:

    □θ = g · ⟨(φ₁-φ₂)² + (φ₂-φ₃)² + (φ₃-φ₁)²⟩_static

The antisymmetric energy density is localized (Yukawa) but nonzero over a
region of size ~1/m_A. The Poisson Green's function extends θ BEYOND the
Proca range → longer-range interaction.

## Setup

Four field types: φ₁, φ₂, φ₃ (with pairwise coupling λ) + θ (massless).

EOM:
    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b+φ_c) - ∂V_triple/∂φ_a
    ∂²θ/∂t² = ∂²θ/∂x² - g · S(x,t)

The source S should be the TIME-AVERAGED antisymmetric energy, not the
oscillating part (to avoid the traveling-wave problem from V24-MF).

Implementation: compute S as a running average of the antisymmetric
energy density over one breathing period T = 2π/ω:

    S(x,t) = (1/T) ∫_{t-T}^{t} [(φ₁-φ₂)² + (φ₂-φ₃)² + (φ₃-φ₁)²] dt'

This is the DC component of the antisymmetric energy — nonzero wherever the
oscillon has antisymmetric field content (which it does when pairwise coupling
is present, since the equilibrated fields are not exactly equal).

Actually, for the SYMMETRIC oscillon (φ₁=φ₂=φ₃): the antisymmetric
combination (φ₁-φ₂) = 0 everywhere. There's no antisymmetric content.

So the source S = 0 for a single symmetric oscillon. The cascade only works
for ASYMMETRIC oscillons (UUD, 180°, etc.) or for PERTURBATIONS of the
symmetric state.

**Alternative source**: Use the TOTAL energy density ρ instead:
    □θ = -g · ⟨ρ⟩_static

This is nonzero for any oscillon and gives the V6/gravity-like coupling.
The Proca isn't needed as an intermediate step — θ directly couples to ρ.

**Better alternative**: Source θ from the PAIRWISE energy specifically:
    S(x,t) = ⟨φ₁φ₂ + φ₂φ₃ + φ₃φ₁⟩_time

For the symmetric oscillon: ⟨φ₁φ₂⟩ = ⟨f²cos²ωt⟩ = f²/2. So S = 3f²/2.
This IS nonzero, localized at the oscillon, and static. It sources θ:
    □θ = -g · (3f²/2)  → θ ~ g·Q/(4πr) with Q = ∫(3f²/2)dx

## Method

1. Implement the 4-field system (φ₁,φ₂,φ₃,θ) with pairwise coupling λ
2. Source θ from the time-averaged pairwise energy:
   Simple approach: use instantaneous S = φ₁φ₂+φ₂φ₃+φ₃φ₁ (oscillates at 2ω)
   Better: accumulate running average over one period
3. Run single oscillon at λ=0.99, g=0.1. Measure θ profile.
4. Does θ develop a LINEAR (1D) or 1/r (3D) tail?
5. Two oscillons: measure θ-mediated force at D = 50, 100

## Reference Code

- v24/maxwell_e: pairwise coupling code
- v24/maxwell_f: θ field implementation

## Output

- `src/cascade1d.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ=0.99
g scan: {0.01, 0.1, 1.0}
Nx=8000, xmax=200, tfinal=10000

Compile: `gcc -O3 -Wall -o cascade1d src/cascade1d.c -lm`
