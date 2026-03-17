# V24-PT: Pairwise + Triple Product — 3ω < m Threshold

## Thesis

Combine pairwise coupling V_pw = λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) with the existing
triple product V_triple = (μ/2)P²/(1+κP²). The pairwise coupling lowers
the antisymmetric base mass to m²_A = m² - λ. The triple product provides
nonlinear saturation.

At λ ≈ 0.88m²: m²_A = 0.12m², ω_A ≈ 0.31m, 3ω_A ≈ 0.93m < m.
The 120° oscillon would radiate NOTHING (3rd harmonic below gap).

## EOM

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b + φ_c) - μP·dP/dφ_a/(1+κP²)²

## Method

1. Fine-scan λ near the threshold: {0.80, 0.82, 0.84, 0.86, 0.87, 0.88,
   0.89, 0.90, 0.92, 0.95} (units of m²)
2. Initialize 120° oscillon at each λ (same as V24-PW)
3. Evolve t=20000 (long run to test true stability)
4. Measure: ω, 3ω, whether 3ω < m, dE/dt, lifetime
5. The KEY diagnostic: does dE/dt → 0 when 3ω drops below m?

## Prediction

At the 3ω = m threshold (λ ≈ 0.88):
- Below threshold (λ < 0.88): 3ω > m, oscillon radiates, finite lifetime
- Above threshold (λ > 0.88): 3ω < m, NO radiation, infinite lifetime
- The transition should be SHARP (radiation is exponentially suppressed
  once 3ω < m)

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Nx=4000, xmax=100, tfinal=20000

## Output

- `src/pt120.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -Wall -o pt120 src/pt120.c -lm`
