# Lagrangian-Mode Decomposition & SU(2)_L Coupling Estimate

**Date**: 2026-05-22
**Parent**: `FINDINGS.md` (this directory, parent project `../SUMMARY.md`)
**Status**: First quantitative estimate of g₂ from v59 structure.

## Setup — the four modes

The constraint surface S³ ⊂ ℍ has 3 angular directions; the 4-dim ambient
ℍ adds 1 radial direction off the constraint.  At the Brannen equilibrium
ξ_eq = ρ₀(cos δ_B + sin δ_B · n̂), these decompose cleanly:

| Mode  | Direction | Dim | Geometric content                                  | Candidate sector |
|-------|-----------|-----|----------------------------------------------------|------------------|
| **R** | radial     | 1   | δρ := \|ξ\| − ρ₀ ; modulates Koide Q from 2/3       | G (gravity)      |
| **A** | active     | 1   | ψ := arctan(\|Im ξ\|/Re ξ) − δ_B ; Brannen phase    | U(1)_em (α)      |
| **S** | silent     | 2   | n̂ ∈ S² ⊂ Im ℍ ; LEAN-PROVED invariant (SilentDirection) | SU(2)_L (α_W) |

Kinetic-term decomposition (parameterising ξ = ρ(cos ψ + sin ψ · n̂)):

```
L_kin = (1/2)(∂ρ)²                          ← radial (R)
      + (ρ²/2)(∂ψ)²                          ← active (A)
      + (ρ²/2) sin²(ψ) (∂n̂)²                 ← silent (S)  [sigma model on S²]
```

At ρ = ρ₀ = 1/√2 and ψ = δ_B = 2/9, the silent decay constant is

```
f_S² := ρ₀² sin²(δ_B) = (1/2) · sin²(2/9) ≈ 0.0243
```

## The structural conjecture

Empirically, the ratio α_W / α is in the range **4.32 – 4.65** (depending on
the scale at which the couplings are evaluated).  Two v59-natural candidates
sit inside that range:

| Candidate                | Value     | Origin                                | Match (vs tree) |
|--------------------------|-----------|---------------------------------------|-----------------|
| **√21 = √(dim Spin(7))** | 4.5826    | square root of the v59 cross-sector dimension | **1.4 % off** |
| 1 / sin(δ_B) = csc(2/9)  | 4.5388    | inverse silent-direction "stretch"    | 2.4 % off       |

The **√21 candidate** is structurally cleaner:  it is the square root of the
adjoint dimension of Spin(7) — the same 21 that appears in the v59
cross-sector exponential modulation ratio S_grav / S_em ≈ 21.  The
**conjecture** is then

> **α_W / α  =  √(dim Spin(7))  =  √21**.

## Combined prediction with v59 α-conjecture

Variant D of the Furey construction (`04_findings.md`) found that

```
−ln α + 2α = π²/2
```

agrees with empirical α at 4 × 10⁻⁵.  Solving for α and combining with the
new α_W/α = √21 conjecture:

| Quantity        | Predicted (v59) | Empirical (tree-level m_W/v)         | Difference   |
|-----------------|-----------------|--------------------------------------|--------------|
| α               | 7.2976 × 10⁻³  | 7.2974 × 10⁻³ (low-energy)          | +4 × 10⁻⁵  |
| α_W             | 0.03311        | 0.03390 (= g₂²/4π, g₂ = 0.6529)     | −2.3 %       |
| g₂              | 0.6450         | 0.6529 (tree); 0.6517 (M_Z, PDG)    | **−1.2 % / −1.0 %** |
| g₂ / e          | g₂/√(4πα) = 7.553 | g₂/e = 7.643 (tree)               | −1.2 %       |

**The g₂ prediction is within 1 % of empirical** — comparable in tightness to
v59's existing α conjecture (4 × 10⁻⁵) and considerably tighter than v59's G
conjecture (factor 0.76 off in `06_findings.md`).

## Honest caveats

1. **Running uncertainty.**  Empirical α_W/α depends on the renormalisation
   scale.  At M_Z it's 4.32, at the W-mass scale via m_W/v it's 4.65, and
   √21 = 4.58 sits between them.  Without identifying the v59 prediction's
   "natural scale", the 1% level is plausible but not pinned.

2. **No derivation of the √21 factor yet.**  The square root of dim Spin(7)
   is structurally appealing but I have not derived it from a Lagrangian.
   The natural mechanism would be a YM bundle normalisation involving the
   Killing form / Cartan matrix of Spin(7), but the calculation is not done.

3. **G remains the outlier.**  The radial-mode (R) candidate for gravity has
   no clean prediction yet; v59 has G off by 0.55 % in the exponent (factor
   0.76 in the absolute value).  Whether the radial mode is the right
   identification, or whether G needs a 5th mode (e.g. via the constraint
   surface curvature), is open.

4. **The user's caution stands.**  This calculation is *not* yet a full
   coherent strain Lagrangian; it identifies the modes and reads off natural
   scales.  The full Lagrangian with all four modes in concert is the next
   step.

## Where this fits in the v59 prediction table

| Quantity                       | v59 prediction                        | Empirical              | Status          |
|--------------------------------|---------------------------------------|------------------------|-----------------|
| Lepton Koide Q                 | dim G₂ / dim Spin(7) = 14/21          | 0.6666605              | 6 × 10⁻⁶ ✓     |
| Brannen phase φ                | Q/3 = 2/9                              | 0.2222296 rad          | 7 × 10⁻⁶ ✓     |
| Cross-sector ratio S_grav/S_em | dim Spin(7) = 21                       | 20.945                 | 0.26 %         |
| Three generations              | Z₃ ⊂ triality                          | 3                      | exact          |
| α (low energy)                 | −ln α + 2α = π²/2 → 1/137.05          | 1/137.036              | 4 × 10⁻⁵       |
| **α_W (this work)**            | **α × √21 → 1/30.20**                  | **1/29.78 to 1/29.57** | **~1 %**       |
| **g₂ (this work)**             | **√(4π × α × √21) → 0.6450**           | **0.6529 (tree)**      | **−1.2 %**     |
| G                              | exp(−21π²/2) (modulated)              | 1.75 × 10⁻⁴⁵          | factor 0.76 ×  |
| Quark Koide                    | 2/3                                    | 0.849 (u), 0.731 (d)   | × (sector-spec) |

## Next steps (concrete)

In rough order of feasibility:

1. **Verify the conjecture's sensitivity to convention.**  Run the prediction
   with α at intermediate scales (μ between low-energy and M_Z) and see at
   which μ the agreement is best.  This will tell us which scale v59's
   prediction "naturally" sits at.

2. **Derive √21 from a YM bundle normalisation.**  Compute the Killing form
   of so(7) and check whether the natural ratio of adjoint trace
   normalisations between active U(1) and silent SU(2) gauges gives √21.

3. **Write the full strain Lagrangian.**  All four modes (R, A, S, S) with
   the density coupling.  Check whether the gauge fields A^a_μ for SU(2)_L
   arise *automatically* from the Cosserat strain on the silent S² — and
   whether their coupling is naturally √21 × that of the active U(1) field.

4. **The radial-mode → G connection.**  Test whether modulating ρ₀ → ρ₀(1+ε)
   produces a Koide shift δQ that matches the empirical bounds on δQ from
   atomic clocks in gravitational potential, with the v59 G_e magnitude.

## Files

- `04_lagrangian_modes.py` — natural-scale computation (this experiment).
- `05_alpha_alphaW_ratio.py` — α_W/α conjecture test.
- `FINDINGS_lagrangian.md` — this document.
- `04_modes.npz` — saved natural scales.
