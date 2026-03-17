# V24-S3: Multi-Scale Forces — Two Pairwise Couplings

## Thesis

In real physics, the strong force (range ~1 fm) and gravity (range ∞)
coexist. Can the three-field system support TWO different-range forces
simultaneously?

The pairwise coupling V_pw = λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) has the same coupling
strength for all three pairs. But what if different pairs have different
coupling strengths?

### Asymmetric pairwise coupling

    V_pw = λ₁₂·φ₁φ₂ + λ₂₃·φ₂φ₃ + λ₃₁·φ₃φ₁

With λ₁₂ ≠ λ₂₃ ≠ λ₃₁, the mass matrix becomes:

    M² = | m²     λ₁₂    λ₃₁  |
         | λ₁₂    m²     λ₂₃  |
         | λ₃₁    λ₂₃    m²   |

This has THREE distinct eigenvalues (not just two as in the symmetric case).
Each eigenvalue corresponds to a different mediator mass → different range.

Example: λ₁₂ = 0.9 (strong, range 3.2), λ₂₃ = λ₃₁ = 0.99 (weak, range 10).
The mass matrix eigenvalues:
- Mode 1 (mostly symmetric): m² + λ₁₂ + λ₂₃ + λ₃₁ ≈ m² + 2.88
- Mode 2: depends on the specific structure
- Mode 3: depends on the specific structure

The two antisymmetric modes now have DIFFERENT masses → two different ranges.

### Alternative: Nested coupling

Keep the symmetric pairwise coupling but add a SECOND coupling at a
different order:

    V = λ_strong·(φ₁φ₂+φ₂φ₃+φ₃φ₁) + λ_weak·(φ₁φ₂+φ₂φ₃+φ₃φ₁)²

The quadratic term provides a second range scale.

## Method

### Phase 1: Asymmetric pairwise

1. Set λ₁₂ = 0.5 (strong, range 1.4), λ₂₃ = λ₃₁ = 0.99 (weak, range 10)
2. Equilibrate the oscillon with this asymmetric coupling
3. Does the oscillon survive? What is ω?
4. Compute the three eigenvalues of the mass matrix
5. Measure the propagation range for each eigenmode

### Phase 2: Force measurement

6. Two oscillons at D = 5, 10, 15, 20, 30, 40
7. Measure force F(D) at each separation
8. Fit to TWO-Yukawa form: F = F_s·exp(-D/λ_s) + F_w·exp(-D/λ_w)
9. Extract the two ranges λ_s (strong) and λ_w (weak)
10. Are there TWO distinct force scales?

### Phase 3: Scan coupling ratios

11. Fix λ_strong + λ_weak = 1.0 (total coupling constant)
12. Scan ratio: (0.9/0.1), (0.7/0.3), (0.5/0.5), (0.3/0.7), (0.1/0.9)
13. At each ratio: measure the two interaction ranges
14. Can we get both a ~1 range (nuclear) and a ~30 range (long-range)?

## Reference Code

- v24/maxwell_e: pairwise coupling
- v24/proca_force (V24-S1): force measurement

## Output

- `src/proca_multi.c` — asymmetric pairwise solver
- `data/` — force measurements, eigenvalue analysis
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
Phase 1: λ₁₂=0.5, λ₂₃=λ₃₁=0.99
Nx=8000, xmax=200, t_equil=10000, t_run=2000

Compile: `gcc -O3 -Wall -o proca_multi src/proca_multi.c -lm`
