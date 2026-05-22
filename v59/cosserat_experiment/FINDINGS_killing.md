# Killing-Form Calculation — Refutation of √21, Structural √5 Instead

**Date**: 2026-05-22
**Parent**: `FINDINGS_lagrangian.md` (this directory)
**Status**: NEGATIVE for √21; clean structural number √5 instead.

## Summary

The conjecture from `FINDINGS_lagrangian.md` (g_W/g_parent = √21) is **NOT
supported** by a direct Lie algebra calculation.  The natural Killing-form
embedding index for `so(3) ⊂ so(7)` is

```
X = (7 − 2) / (3 − 2)  =  5
```

so the natural gauge-coupling ratio is

```
g_W / g_parent  =  √X  =  √5  ≈  2.236
```

not √21 ≈ 4.58.

The √21 numerology in `05_alpha_alphaW_ratio.py` matched empirical α_W/α at
mixed scales (low-energy α and tree-level g₂) but was off by 5.7 % at the
consistent M_Z scale.  The clean Lie-algebra calculation gives √5.

## What was computed

`06_killing_form.py` builds so(7) explicitly as 7×7 antisymmetric matrices,
computes the Killing form as a 21×21 symmetric matrix, identifies the
so(3) ⊂ so(7) subalgebra acting on the imaginary-direction subspace
⟨i, j, k⟩ of ℍ, and reads off the embedding index by comparing the
intrinsic so(3) Killing form to the restriction of the so(7) Killing form.

Results:
- Killing form B_{so(7)}(E_a, E_a) = −10  for each basis vector E_a.
- Restriction to so(3) sub-basis: diag(−10, −10, −10).
- Intrinsic so(3) Killing form: diag(−2, −2, −2).
- **Embedding index = (−10)/(−2) = 5** — exactly the predicted (N−2)/(n−2) = 5/1.
- Vector trace form ratio = 1 (the vector reps are naturally inclusive).

## Empirical comparison

| Identification              | Predicted        | Empirical          | Gap   |
|-----------------------------|------------------|--------------------|-------|
| g_W / e  (low-E α)         | √5 ≈ 2.236       | 2.152              | +3.9% |
| g_W / e  (M_Z α)           | √5 ≈ 2.236       | 2.087              | +7.1% |
| g_W / g_Y  (M_Z, hypercharge)| √5 ≈ 2.236      | 1.820              | +23%  |

The closest match is to **g_W/e at low energy** (√5 vs 2.152, 3.9 % gap).
This is meaningful but not at the precision level of the v59 α conjecture
(4 × 10⁻⁵).

The Spin(7) → SO(3) × SO(4) symmetry breaking pattern is NOT how the SM
electroweak structure works (SM has U(1)_Y × SU(2)_L → U(1)_em via Higgs).
So the "g_parent" of our identification is ambiguous — it could be e, g_Y,
or something more abstract.

## Why √21 fails

The Killing form of so(7) is B_{so(7)}(X, Y) = (7−2) Tr_vec(X Y) = 5 Tr_vec.
The dim of so(7) (= 21) does NOT enter the natural coupling-ratio formula
for SO(3) ⊂ SO(7) embedding.  The 21 is the dim of the adjoint rep of
Spin(7), but the embedding ratio uses the trace-form normalization
difference (N−2)/(n−2) = 5/1, not the adjoint dim.

The earlier α_W/α = √21 conjecture was a **numerological coincidence**:
- At mixed scales (low-E α, tree g₂):  empirical ratio 4.65; √21 = 4.58; gap 1.4 %.
- At consistent M_Z scale: empirical ratio 4.32; √21 = 4.58; gap 5.7 %.

Neither match has Lie-algebra justification.  The √5 result (g_W/g_parent
= 2.24) is the structurally correct number; it has its own ~7 % gap with
empirical g_W/e, but it comes from a derivation rather than a fit.

## Implications

1. **The √21 prediction is withdrawn.**  v59's "α_W/α = √21" claim from
   `FINDINGS_lagrangian.md` does NOT survive the Killing-form check.
   Update `SUMMARY.md` and the memory entry accordingly.

2. **The structural prediction is √5.**  The Killing-form ratio
   (7−2)/(3−2) = 5 gives g_W/g_parent = √5 with no fitting.  This is
   within 7 % of empirical g_W/e at M_Z — in the ballpark, not exact.

3. **What's needed for precision.**  A precise g_W requires:
   - Identifying the correct "parent" coupling g_parent (probably not e
     directly, since the SM mixing is U(1)_Y × SU(2)_L not SO(7) → SO(3) × SO(4)).
   - A Lagrangian that produces the YM kinetic terms with the right
     normalisation — *which* gauge group structure (so(7), Spin(7)/G_2 quotient,
     or something else) is the v59 parent of SU(2)_L is not yet determined.

4. **The Cosserat-strain Lagrangian still needs to be written**, but in
   concert with the cross-sector structure, as the user warned.  The
   Killing-form calculation tells us the *natural* coupling ratio under
   one specific embedding hypothesis; the right hypothesis remains open.

## Honest scoring

| Conjecture                  | Source                                | Match (best case) | Status |
|-----------------------------|---------------------------------------|-------------------|---|
| Koide Q = 14/21             | dim G_2 / dim Spin(7)                  | 6 × 10⁻⁶           | ✓ |
| Brannen phase 2/9           | Q/3                                    | 7 × 10⁻⁶           | ✓ |
| −ln α + 2α = π²/2          | instanton on Cl(3,1)                   | 4 × 10⁻⁵           | ✓ |
| α_W/α = √21                | numerology (mixed scales)              | 1.4 % (mixed) / 5.7 % (M_Z) | ✗ refuted |
| g_W = √5 × g_parent        | Killing form so(3) ⊂ so(7)             | 3.9 % – 7.1 %      | ⚠ ballpark only |
| G via 21·π²/2              | cross-sector extrapolation             | factor 0.76        | ⚠ form OK, value off |

## Lessons

- **Numerological fits at the few-percent level can mislead.**  √21 looked
  promising but had no Lie-algebra support.
- **Lie-algebra calculations give clean structural numbers** but those
  numbers aren't always close to empirical (here √5 vs 2.08, 7 % gap).
- **The honest v59 status on gauge couplings**: α at 10⁻⁵, α_W at the
  ~7 % level via √5, g_Y not yet addressed structurally.

## Files

- `06_killing_form.py` — Killing-form calculation (this experiment).
- `06_killing.npz` — saved data.
- `FINDINGS_killing.md` — this document.
