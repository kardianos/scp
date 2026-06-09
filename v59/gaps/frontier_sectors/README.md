# Frontier Sectors (G10 / G11 / G12) — Background

**Cluster**: the untouched Standard-Model sectors of v59 — quark mixing (CKM),
neutrinos + PMNS, and the strong sector (α_s, θ_QCD).
**Date**: 2026-05-25 · **Status**: mostly OPEN reconnaissance.

> **Status tags** (same convention as `UNIFIED_THEORY.md`):
> **[thm]** machine-checked Lean identity; **[emp≈X]** empirical match at precision X;
> **[conj]** structural ansatz, not derived; **[free]** genuine input, no v59 reading.

This cluster is, by design, the part of v59 that is **least reduced**. The
algebraic skeleton (Furey `Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆`, the L⊕F grade bisection, the
exceptional chain `G₂⊂Spin(7)⊂Spin(8)`) pins the lepton + electroweak + Higgs
block down to ~2 inputs (`a_ℓ`, `α`), but the quark-flavour / CKM / neutrino /
strong sectors add **~10–12 free parameters** (`UNIFIED_THEORY.md` §8). This
folder maps what v59 currently *does* and *does not* say about three of them.

---

## G10 — CKM quark mixing  (3 angles + 1 phase)

| quantity | v59 status |
|---|---|
| `sin²θ_C` (Cabibbo) | `≈ 7·α(0) = dim(Im 𝕆)·α` **[conj, emp≈0.9% on sin², 0.45% on sinθ]** |
| `V_cb / sin²θ_C` (Wolfenstein A) | `≈ 7/9 = t²_u = cos²θ_W` **[conj, emp≈6% — weak]** |
| `V_ub`, `V_td` | **[free]** — only Wolfenstein-derived from empirical (ρ̄,η̄) |
| `δ_CP` | **[free]** — no v59 candidate |

**The one real signal** is the Cabibbo conjecture
`sin²θ_C = (dim Im 𝕆)·α = 7·α(0)` (`FINDINGS_ckm_and_selection.md`, `10_cabibbo_test.py`).
The `7` is a deeply recurrent v59 integer: `dim Im 𝕆 = dim S⁷ = dim Λ⁶ℝ⁷` = the
top grade of the lepton ambient `L = Λ²⊕Λ⁶`, and the same `7` that appears in
`cos²θ_W = 7/9`, `m_Z/m_W = 3/√7`, `m_H²/v² = 7/27`. The match is 0.45% on
`sinθ_C` — softer than the lepton sector's `10⁻⁵`, and it uses `α(0)` (α(M_Z)
misses by 4%). The rest of CKM is essentially free.

**Mechanism attempts** (`09_ckm_and_selection.py`): Brannen circulant kernels with
ℂ-valued ξ share the DFT eigenbasis ⇒ `V_CKM` collapses to a *permutation*
(empirically wrong). ℍ-valued ξ in sector-specific complex slices does give
mixing, but ~0.4–0.6 (far larger than the observed ≤0.22) for hand-chosen slices.
A *principled* slice choice that reproduces the near-diagonal CKM is **not found**.

## G11 — Neutrinos + PMNS  (≥7 parameters)

**Completely unaddressed [free].** Three open pieces:
- **PMNS mixing** is *large* (θ₁₂≈34°, θ₂₃≈48°, θ₁₃≈8.5°) — near **tri-bimaximal**
  (`sin²θ ≈ 1/3, 1/2, 0`), the opposite of the small CKM angles. This is a
  *symmetry* pattern; v59 has **no mechanism** for large near-degenerate mixing.
- **The ν-Koide ambient puzzle**: the charged-fermion pattern is `t²_N = 1−14/D_N`
  for `D ∈ {28,35,63}` (`Q ∈ {2/3, 11/15, 23/27}`). The neutrino Koide fits **none**
  of these ambients — a known v59 puzzle (`UNIFIED_THEORY.md` §9-D G11). Neutrinos
  appear to live *outside* the Brannen-circulant charged-fermion structure.
- **Absolute scale / ordering / Majorana-vs-Dirac**: no candidates. A seesaw frame
  is order-plausible but needs an intermediate `M_R ~ 10³–10⁴ GeV` with no v59 origin.

## G12 — Strong sector  (α_s, θ_QCD)

| quantity | v59 status |
|---|---|
| color group `SU(3)` | **[thm]** — constructed in `ColorSU3.lean` from `Cl(6)≅ℂ⊗𝕆` Witt decomposition (8 generators, full A₂ closure, kills lepton singlet, commutes with `J_c`) |
| `α_s(M_Z)` | **[free]** — the gauge GROUP exists but its coupling does not |
| `θ_QCD ≈ 0` | **[free]**, with a STRUCTURAL LEAD (color complex structure `J_c`) |

The color gauge **group** is genuinely built (theorem-grade). What is missing is
its **coupling** `α_s` and the **strong-CP angle** `θ_QCD`. A naive analog of the
EW gauge ansatz `g² = h∨·√α` with the color dual Coxeter `h∨(A₂)=3` gives
`g_s² = 3√α ≈ 0.26`, but the measured `g_s² = 4π·α_s ≈ 1.48` — **off by >80%**
(see `pmns_strong_test.py`, `FrontierSectors.lean:strong_coupling_analog_fails`).
The `√α` form fundamentally cannot produce a *strong* coupling. `θ_QCD≈0` may be
forced by the color complex structure `J_c` reality (see `ALTERNATIVES.md`).

---

## Honest framing

This cluster is the **most prone to spurious numerology**: many free parameters,
each easy to hit with a small-integer×α or simple-ratio form. The accompanying
scripts (`ckm_overfitting_test.py`, `pmns_strong_test.py`) deliberately quantify
*how easy* — and `FINDINGS.md` reports the null results as prominently as the one
live lead. The single clean identity (`sin²θ_C = 7α`) is formalized in
`FrontierSectors.lean`; everything else is `[free]` or a flagged `[conj-lead]`.

## Files
- `README.md` — this file.
- `ALTERNATIVES.md` — candidate structural origins, each with a test + falsifier.
- `FINDINGS.md` — per-sector verdict + honest follow-up assessment.
- `ckm_overfitting_test.py` — CKM conjecture tests + overfitting bound.
- `pmns_strong_test.py` — PMNS / α_s / θ_QCD reconnaissance (mostly null).
- `FrontierSectors.lean` — the clean `sin²θ_C=7α` identity + α_s falsifier
  + ν-ambient-puzzle arithmetic. *Written, NOT built this run.*
