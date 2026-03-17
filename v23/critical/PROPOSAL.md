# V23-C: Critical Gravity — Correlation Length Near Gap Edge

## Thesis

The v21 oscillon exists at ω ≈ 0.95m (5% gap margin in 3D). This near-
criticality is not accidental — the oscillon CAN ONLY EXIST near the gap
edge. As ω → m, the correlation length ξ ~ 1/√(m²-ω²) diverges.

This investigation measures the actual correlation length of vacuum
perturbations around the oscillon as a function of the gap margin (m-ω)/m.
If ξ diverges, the oscillon's influence extends to macroscopic distances,
enabling long-range gravity-like interactions even through a massive
(Yukawa) vacuum.

The key question: does the INTERACTION between two oscillons transition
from Yukawa (e^{-mr}/r) to power-law (1/r^n) as the gap margin shrinks?

## Mathematical Setup

### Correlation Length

For a massive field with mass m in vacuum: ξ_vac = 1/m.

Near an oscillon with breathing frequency ω < m, the field has an
exponential tail:

    φ(r) ~ A·e^{-κr}/r,    κ = √(m² - ω²)

The "effective" correlation length is:

    ξ = 1/κ = 1/√(m² - ω²)

At the v21 parameters (ω ≈ 0.95m):

    ξ = 1/√(1 - 0.9025) = 1/√0.0975 = 3.2/m = 3.2 code lengths

As ω → m: ξ → ∞.

### Two-Oscillon Interaction

Two oscillons at separation D interact through tail overlap. The
interaction energy scales as:

    E_int(D) ~ exp(-2κD) / D    (product of two tails)

This is Yukawa with range 1/(2κ) = ξ/2.

For large ξ (near criticality): E_int ~ exp(-D/1.6·ξ) / D, which
approaches 1/D (power-law) as ξ → ∞.

### What Controls the Gap Margin

The gap margin (m-ω)/m depends on the nonlinear frequency shift, which
is set by the oscillon amplitude A and the coupling parameters:

    ω² ≈ m² - C|μ|A⁴/(1+κA⁶)²

Larger |μ| or larger A → smaller ω → larger margin → shorter ξ.
Smaller |μ| or smaller A → ω closer to m → longer ξ.

To scan ξ: vary μ (or equivalently m) while keeping the oscillon alive.

## What to Compute

### Phase 1: Correlation Length vs Gap Margin (1D)

1. Start from the v21 1D triad code (triad1d.c).
2. Scan the coupling μ from -20 (standard) toward -5 (weaker coupling).
3. For each μ: let the oscillon equilibrate (t > 2000), then measure:
   a. Breathing frequency ω (from FFT of central field value)
   b. Tail decay rate κ (fit e^{-κr}/r to the profile at large r)
   c. Gap margin δ = (m-ω)/m
   d. Correlation length ξ = 1/κ
4. Plot ξ vs δ. Expect ξ ~ 1/√(2mδ) = 1/√(2m(m-ω)).
5. Find the minimum |μ| for which the oscillon still exists (critical μ).

### Phase 2: Two-Oscillon Interaction vs Gap Margin (1D)

6. For each μ from Phase 1: initialize two oscillons at separation D.
7. Measure the force (acceleration of separation) vs D.
8. Fit to F(D) = F₀·exp(-D/λ)/D^n to extract range λ and exponent n.
9. Plot λ vs ξ. If λ ∝ ξ, confirms that interaction range is controlled
   by the gap margin.
10. Check: does n approach 0 (pure exponential) or 1 (Yukawa) or 2 (1/r²)?

### Phase 3: Power-Law Transition (1D)

11. At the most critical μ (longest ξ): measure the interaction at many
    separations D = 5, 10, 15, 20, 25, 30, 40, 50.
12. Log-log plot of |F| vs D. Check for a power-law regime at D < ξ.
13. If power-law exists: determine exponent. Is it consistent with
    1/D² (gravitational)?

### Phase 4: Verification in 3D (if Phase 1-3 positive)

14. Repeat the most promising μ value in 3D (triad3d.c).
15. Measure the 3D correlation length and interaction range.
16. Compare with 1D results (dimensionality affects exponents).

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (primary base code)
- v21 3D solver: `/home/d/code/scp/v21/src/triad3d.c` (for Phase 4)
- v22 two-oscillon: `/home/d/code/scp/v22/src/two_oscillon.c`
  (two-body interaction measurement — adapt separation tracking)
- v12 normal modes: `/home/d/code/scp/v2/src/normal_modes.c`
  (eigenvalue solver for ξ measurement)

## Output Structure

- `src/critical1d.c` — 1D oscillon with parameter scanning
- `src/critical_two.c` — two-oscillon interaction measurement
- `data/` — output data files (profiles, spectra, force measurements per μ)
- `RESULTS.md` — results and analysis

## Parameters

Base: κ=20, m=1.0 (fixed)
Scan: μ ∈ [-20, -5] in steps of 1.0 (Phase 1)
      μ ∈ [-8, -5] in steps of 0.5 (fine scan near critical, Phase 3)
Grid: Nx=4000, xmax=100, tfinal=10000 (need large domain for long ξ)

Note: for small |μ|, the oscillon will be weakly bound and may need longer
equilibration time. Use tfinal=10000 to ensure convergence.

Compile: `gcc -O3 -Wall -o critical1d src/critical1d.c -lm`
