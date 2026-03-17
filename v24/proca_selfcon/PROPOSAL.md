# V24-P4: Self-Consistent λ from Field Amplitude

## Thesis

Make the pairwise coupling position-dependent:

    λ_eff(x) = λ₀ + α·(φ₁² + φ₂² + φ₃²)

Near the oscillon (φ large): λ_eff large → m_A small → long range.
In vacuum (φ = 0): λ_eff = λ₀ → m_A = √(m²-λ₀) → shorter range.

This creates a "low-mass corridor" around each oscillon where the Proca
mediator propagates freely. Two nearby oscillons have overlapping corridors,
enhancing their interaction.

## EOM

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ_eff(x)·(φ_b+φ_c)
                 - α·2φ_a·(φ_aφ_b + φ_aφ_c + additional)
                 - ∂V_triple/∂φ_a

The derivative of λ_eff with respect to φ_a creates additional force terms:

    ∂/∂φ_a [λ_eff·(φ₁φ₂+φ₂φ₃+φ₃φ₁)]
    = λ_eff·(φ_b+φ_c) + 2αφ_a·(φ₁φ₂+φ₂φ₃+φ₃φ₁)

The second term couples the pairwise energy to the individual field amplitude.

## Method

1. Fix λ₀ = 0.5 (moderate base coupling, m_A,vac = 0.71, range_vac = 1.4)
2. Scan α ∈ {0.0, 0.5, 1.0, 2.0, 5.0, 10.0}
3. For each α: evolve the oscillon for t=10000
4. Measure: effective λ at center (λ₀ + α·3f²), effective m_A at center
5. Add an antisymmetric perturbation at the oscillon edge, measure how far
   it propagates before decaying. This gives the LOCAL range.
6. Compare the range near the oscillon vs in vacuum.

## Key Question

Does the oscillon create a "halo" of low m_A that extends beyond its
core? If the halo extends to r_halo > 1/m_A,vac, the effective range is
r_halo, not 1/m_A,vac.

At the oscillon center with f ≈ 0.5: λ_eff = 0.5 + α·3·0.25 = 0.5 + 0.75α.
For α = 1: λ_eff(0) = 1.25 > m² → locally tachyonic! But stabilized by the
oscillon's nonlinear terms.

For α = 0.5: λ_eff(0) = 0.875, m_A(0) = 0.35 (range 2.8 locally)
At r where f ≈ 0.1: λ_eff = 0.53, m_A = 0.69 (range 1.5)
At r where f ≈ 0: λ_eff = 0.5, m_A = 0.71 (range 1.4)

The range variation is modest at α=0.5. Need larger α for dramatic effect.

## Reference Code

- v24/maxwell_e: `/home/d/code/scp/v24/maxwell_e/src/maxwell_e.c`
- v21 1D: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/selfcon1d.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ₀=0.5
α scan: {0.0, 0.5, 1.0, 2.0, 5.0, 10.0}
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o selfcon1d src/selfcon1d.c -lm`
