# FINDINGS — Frontier Sectors (G10 / G11 / G12)

**Date**: 2026-05-25 · **Cluster**: CKM, neutrinos+PMNS, strong sector.
**Bottom line up front**: this cluster is, as expected, **mostly genuinely
open**. There is **one** soft live lead worth follow-up (G12 θ_QCD via the color
complex structure `J_c`), **one** previously-claimed signal that survives but
only as a coincidence-grade match (G10 `sin²θ_C=7α`), and a clean **null result**
that should be recorded against the theory (the EW `g²=h∨√α` ansatz does **not**
extend to `α_s`). Everything else is free.

---

## Per-sector verdict

### G10 — CKM  →  ONE soft signal, rest FREE
- **`sin²θ_C = 7α(0)` [conj, emp≈0.9%]**: survives. The `7=dim Im𝕆` is genuinely
  recurrent in v59, and the overfitting scan (`ckm_overfitting_test.py`) shows
  that *among `×α` forms* it is not generic (ranks #1, only 2 forms within 2%).
  **Caveat**: it requires the IR `α(0)` (α(M_Z) misses by 4%), which is
  physically unmotivated for a quark-mixing angle, and the match (0.45% on
  `sinθ_C`) is two orders of magnitude *looser* than the lepton Koide/phase
  (10⁻⁵). **Honest grade: a structurally-flavoured coincidence, not a prediction**
  — until a mechanism explains *why the IR photon coupling sets quark mixing*.
- **`V_cb/V_us² = A ≈ 7/9` (6%) or `4/5=28/35` (3%)**: the `28/35` (L/F grade
  overlap) reading is *closer* than the documented `7/9` and is v59-structural —
  a marginal lead for A2 (sector-overlap), but at 3% it is weak.
- **`V_ub`, `V_td`, `δ_CP`**: **[free]**. The overfitting test found `δ_CP`
  matched by *two* unrelated forms (`π/√7` at 1%, `π·7/18` at 1.8%) — the
  signature of overfitting, not signal. No v59 reading.

### G11 — Neutrinos + PMNS  →  GENUINELY OPEN
- **PMNS angles** are near **tri-bimaximal** (`sin²θ ≈ 1/3, 1/2, 0`) — an `S₃/A₄`
  *symmetry* pattern, not a v59 one. v59 *has* an `S₃` (sedenion generation
  symmetry), so a democratic-mixing reading (B2) is *qualitatively* plausible,
  but it predicts `θ₁₃=0` whereas `θ₁₃≈8.5°≠0`, and no v59 term supplies the
  correction. **Open.**
- **ν-Koide ambient puzzle CONFIRMED**: no v59 ambient `{28,35,63}` (giving
  `Q∈{2/3,11/15,23/27}`, all `>2/3`) reproduces the neutrino Koide (which trends
  to `~1/3` in the normal-ordering massless-lightest limit, *below* the charged
  minimum). Formalized in `FrontierSectors.lean:neutrino_below_charged_minimum`.
  Best hypothesis (B1): neutrinos are a **non-Brannen, symmetric-grade (Majorana)
  sector** `Λ⁰⊕Λ⁴`, not the skew/circulant charged-fermion family — testable.
- **Absolute scale / ordering / Majorana-vs-Dirac**: no candidates. Seesaw needs
  `M_R ~ 5×10³ GeV`; `v·dim(L) = 246·28 ≈ 6.9 TeV` is coincidentally in range
  (B3) but rests on `m_D~m_e` with no v59 justification. **Open.**

### G12 — Strong sector  →  one LIVE LEAD (θ_QCD), `α_s` FREE
- **`θ_QCD ≈ 0` — the best lead in the cluster.** The color sector carries a
  *theorem-grade* complex structure `J_c=γ₀γ₅` (`ColorSU3.lean`). The θ-term is
  CP-odd; the same reality/skewness forcing that pins the lepton `J`
  (`LeptonRealityForcing`) plausibly forbids the CP-odd projection of `⟨F∧F⟩`.
  This would **explain an observed near-zero from existing structure** rather
  than fit a number — the rare case where the structural lead is genuinely
  attractive. **Deserves real follow-up** (C1): compute the `J_c`-grading of the
  color `⟨F∧F⟩` bilinear.
- **`α_s(M_Z)` — clean NULL / FALSIFIER.** The EW ansatz `g²=h∨·√α` with the
  color dual Coxeter `h∨(A₂)=3` gives `g_s²=3√α≈0.256`, vs measured
  `g_s²=4π·α_s≈1.48` — **>80% off** (machine-checkable:
  `FrontierSectors.lean:strong_coupling_analog_fails` proves `3√α<1<1.48`). The
  `√α` form cannot produce a strong coupling at *any* `h∨`. GUT-style
  unification (C3) also fails (`g_s²/g_W²=3/5≠1`). **`α_s` is genuinely free**;
  a v59 reading would need to generate the *dimensionful* `Λ_QCD` via running —
  a different kind of object than v59's ratio-predictions.

---

## Which deserves real follow-up vs which is just free?

| item | grade | follow-up? |
|---|---|---|
| **G12 θ_QCD via `J_c` reality** | structural lead | **YES** — explains observed ≈0 from a theorem-grade structure; concrete test (grade of `⟨F∧F⟩`) |
| G11 ν = symmetric-grade (Majorana) sector | hypothesis | maybe — would *explain* the ambient puzzle; needs the symmetric-bilinear eigenvalue computation |
| G10 `sin²θ_C=7α` | coincidence-grade | only IF a mechanism ties quark mixing to IR `α(0)`; otherwise leave as-is |
| G11 PMNS democratic `S₃` | qualitative | low priority — over-predicts `θ₁₃=0`, no correction term |
| G10 `V_cb`,`V_ub`,`V_td`,`δ_CP` | **free** | no |
| **G12 `α_s` via `h∨√α`** | **FALSIFIED** | no — record the null result |
| G11 absolute ν scale / `M_R` | **free** | no |

## Honest assessment

Three of v59's biggest gaps remain gaps. The headline "≈2 inputs" of the
lepton+EW block does **not** extend here: the full SM still carries ~10–12 free
parameters, and this cluster contains most of them. The most valuable output of
this reconnaissance is **negative**: the `g²=h∨√α` gauge ansatz — one of the
load-bearing v59 conjectures (G2) — **does not generalize to the strong sector**,
and `δ_CP`/PMNS show the easy-to-overfit failure mode explicitly. The one place
where v59 *structure* (not numerology) could genuinely contribute is **strong CP
(θ_QCD=0) via the color complex structure `J_c`** — that is where to look next.

## Artifacts
- `ckm_overfitting_test.py` — runs clean; CKM tests + overfitting scan.
- `pmns_strong_test.py` — runs clean; PMNS / α_s / θ_QCD reconnaissance.
- `FrontierSectors.lean` — `sin²θ_C=7α` identity, ν-ambient arithmetic, and the
  `α_s` falsifier. **Written, NOT built this run** (shared lean project; no
  `lake build` per task rules). Imports `Mathlib`; targets the project's
  Mathlib `v4.29.0`. Open items flagged in comments; no `sorry` used (all stated
  theorems are arithmetic/inequality facts verified numerically in Python).
