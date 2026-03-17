# V24-ME: Proca Field from Mass-Split Sector

## Background

The pairwise coupling V_pw = λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) splits the mass spectrum:
m²_anti = m² - λ (antisymmetric, lighter) and m²_sym = m² + 2λ (symmetric,
heavier). At λ → m², the antisymmetric mode becomes nearly massless.

If the lightest mode is identified as the "photon" (Proca field), we get a
massive vector-like mediator that approaches Maxwell in the m → 0 limit.

V24-PW showed that λ near m² causes condensation. This test finds the
MAXIMUM λ that avoids condensation, giving the LIGHTEST possible Proca
"photon" mass.

## Setup

The Lagrangian already has the pairwise coupling from V24-PW:

    L = Σ_a [½(∂φ_a)² - ½m²φ_a²] - λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) - V_triple(P)

The antisymmetric modes A₁ = (φ₁-φ₂)/√2, A₂ = (φ₁+φ₂-2φ₃)/√6 have:
    m²_A = m² - λ

These modes propagate with dispersion ω² = k² + m²_A. For m²_A > 0, the
range is 1/m_A (Yukawa). For m²_A → 0: range → ∞ (Coulomb-like).

## What to Compute

### Phase 1: Maximum Stable λ

1. Scan λ from 0 to m² in steps of 0.05
2. At each λ: initialize the symmetric 0° oscillon, evolve t=5000
3. Monitor: does the oscillon survive? Does the vacuum condense?
4. Find λ_max: largest λ where the oscillon is stable (fc > 0.9)
5. Record m_A = √(m² - λ_max): the lightest Proca mass

### Phase 2: Antisymmetric Mode Propagation

6. At λ_max: excite an antisymmetric perturbation on the oscillon
7. Measure how far the perturbation propagates before decaying
8. The propagation distance = the Proca range = 1/m_A
9. Compare with the prediction 1/√(m² - λ_max)

### Phase 3: Two-Oscillon Interaction via Proca

10. Two oscillons at separation D, with λ = λ_max
11. The antisymmetric mode mediates a Yukawa interaction with range 1/m_A
12. Measure the force F(D). Compare with F = F₀·exp(-m_A·D)
13. Is the range LONGER than without pairwise coupling (λ=0)?

## Reference Code

- v24/120pw: `/home/d/code/scp/v24/120pw/src/pw120.c` (pairwise coupling)
- v21 1D: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/maxwell_e.c` — solver with pairwise coupling + perturbation analysis
- `data/` — TSV output
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
λ scan: {0.0, 0.1, 0.2, ..., 0.9, 0.95}
Nx=4000, xmax=100, tfinal=5000 (Phase 1), 10000 (Phase 3)

Compile: `gcc -O3 -Wall -o maxwell_e src/maxwell_e.c -lm`
