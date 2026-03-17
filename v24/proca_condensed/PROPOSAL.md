# V24-P2: Condensed Phase Goldstone — Cross the Tachyonic Boundary

## Thesis

At λ > m²: m²_A = m² - λ < 0 → antisymmetric modes are tachyonic in vacuum.
They condense to a nonzero VEV. The condensate breaks SO(2) symmetry in the
antisymmetric subspace → Goldstone mode with m = 0 exactly → true 1/r.

The oscillon survives because m²_S = m² + 2λ > 3m² (very stable).

## Setup

Same EOM as V24-ME but with λ > m²:

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b+φ_c) - μP·dP/dφ_a/(1+κP²)²

At λ > 1.0 (with m=1): the vacuum φ=0 is unstable for the antisymmetric
modes. The system should condense to a new vacuum where the antisymmetric
mode has a nonzero expectation value.

### New Vacuum

The pairwise potential for uniform fields φ_a = const:
V_pw = λ(φ₁φ₂+φ₂φ₃+φ₃φ₁)In the antisymmetric direction (e.g., φ₁ = v, φ₂ = -v/2, φ₃ = -v/2 with
φ₁+φ₂+φ₃ = 0), the total potential is:
V_total = ½m²(v² + v²/4 + v²/4) - λ(v²/4) + triple product...
       = ½m²·(3v²/2) - λv²/4 + ...

The new minimum is at v² = (λ - m²) × (some positive factor). This is
the antisymmetric condensate.

Around this new vacuum: fluctuations have one massive (radial) mode and
one massless (Goldstone) mode.

## Method

### Phase 1: Vacuum Structure

1. Initialize uniform fields in the antisymmetric configuration:
   φ₁ = v, φ₂ = -v/2, φ₃ = -v/2 (with small v)
2. Let the system evolve to find the true vacuum at λ > m²
3. Measure the condensate VEV v_eq

### Phase 2: Fluctuation Spectrum

4. Around the condensed vacuum: add small perturbations
5. Measure the oscillation frequencies of the two antisymmetric modes
6. One should be massive (Higgs-like), one should be massless (Goldstone)
7. Verify: ω_Goldstone → 0 as k → 0

### Phase 3: Oscillon in Condensed Vacuum

8. Initialize an oscillon (symmetric Gaussian) ON TOP of the condensate
9. Does the oscillon survive in the new vacuum?
10. Does the Goldstone mode create a 1/r (linear in 1D) tail around it?

### Phase 4: Two-Oscillon Interaction

11. Two oscillons at D=50, 100, 150 in the condensed vacuum
12. Measure force vs D. Is it 1/D² (from Goldstone) or Yukawa?

## Scan

λ scan: {1.01, 1.05, 1.1, 1.2, 1.5, 2.0} (all above m² = 1.0)

## Reference Code

- v24/maxwell_e: `/home/d/code/scp/v24/maxwell_e/src/maxwell_e.c`

## Output

- `src/condensed1d.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
λ scan: {1.01, 1.05, 1.1, 1.2, 1.5, 2.0}
Nx=8000, xmax=200, tfinal=10000

Compile: `gcc -O3 -Wall -o condensed1d src/condensed1d.c -lm`
