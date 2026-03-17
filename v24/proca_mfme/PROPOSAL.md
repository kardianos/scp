# V24-P5: Combined MF+ME — Goldstone Backreaction + Pairwise Coupling

## Thesis

V24-MF found: Goldstone θ backreaction at g=1.0 STRENGTHENS the oscillon
(E increases 131%, ω drops from 0.87 to 0.76, deeper below gap).

V24-ME found: pairwise coupling λ creates a tunable Proca mediator, and
the oscillon tolerates any λ < m².

Combine both: the θ backreaction pushes ω deeper below the gap, which means
the oscillon can tolerate EVEN LARGER λ (further above the stability
threshold). This gives even smaller m_A → longer range.

The positive feedback loop:
1. θ strengthens oscillon (deeper ω below gap)
2. Stronger oscillon tolerates larger λ
3. Larger λ → smaller m_A → longer Proca range
4. Longer range → more θ radiation → more backreaction → goto 1

## Setup

Five field arrays: φ₁, φ₂, φ₃, θ, plus the pairwise coupling λ.

EOM:
    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b+φ_c) - ∂V_triple/∂φ_a
                 + 2g·(∂θ/∂t)·P·(∂P/∂φ_a)/(1+κP²)²   [θ backreaction]
    ∂²θ/∂t² = ∂²θ/∂x² - g·∂(P²)/∂t                    [θ source]

## Method

### Phase 1: Baseline — MF alone then ME alone

1. MF alone (g=1.0, λ=0): measure ω_MF, E_MF
2. ME alone (g=0, λ=0.99): measure ω_ME, E_ME, m_A

### Phase 2: Combined MF+ME

3. g=1.0, λ=0.99: the strongest combination
4. Measure: ω_combined, E_combined, effective m_A
5. Does ω drop further below gap?
6. Does the oscillon tolerate λ > 0.99 (which failed without θ)?

### Phase 3: Push the boundary

7. Fix g=1.0. Scan λ from 0.99 to 1.05 (past the tachyonic boundary!)
8. At each λ: does the oscillon survive? What is the effective m_A?
9. If the θ backreaction stabilizes the antisymmetric sector beyond λ=m²,
   we get the condensed-phase Goldstone from V24-P2 with the OSCILLON
   surviving inside it.

### Phase 4: Two-oscillon

10. At the optimal (g, λ): two oscillons at D=50, 80, 120
11. Measure force vs D. Compare with ME-only at the same λ.
12. Is the range EXTENDED by the θ backreaction?

## Reference Code

- v24/maxwell_e: pairwise coupling
- v24/maxwell_f: θ field with backreaction

## Output

- `src/mfme1d.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
Phase 1: (g=1.0,λ=0), (g=0,λ=0.99)
Phase 2: (g=1.0,λ=0.99)
Phase 3: g=1.0, λ scan {0.99, 0.995, 1.0, 1.01, 1.05}
Nx=8000, xmax=200, tfinal=10000

Compile: `gcc -O3 -Wall -o mfme1d src/mfme1d.c -lm`
