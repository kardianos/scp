# 00 — SUMMARY of the v58↔v59 integration session (2026-05-24)

*Capstone for the `integration_v58/` series (`01`–`20`) and the Lean modules added this session.
The one-line result: **the octonion framework (v59 kinematics) + the v58/OFE dynamics fix the
DISCRETE structure and the MEANING of the couplings, but not their continuous VALUES** — and every
route to derive the values is now closed with a specific, computed reason (the kinematics/dynamics
divide, `17_`).  This indexes the whole arc and the honest status of each claim.*

## Headline status table

| quantity | status | where |
|---|---|---|
| **discrete data** (gauge group, electric charges `Q_em=⅓·color#`, generation count `3`, mass-grade `L`, ratios) | **FIXED by the algebra** (rep-theory/counting) | `09_`, `16_`, `17_` |
| **Koide amplitude `t²=½`** | **structural *meaning* = L-grade complex reversion-norm** `‖(1+u)/2‖`; **value = potential parameter (input)** | `10_`,`12_`,`13_`,`14_`,`20_`; `OctoHalf.lean` |
| **Koide phase `φ=2/9`** | **free residual**, near the chiral/electron-massless edge `π/12`; not symmetry-fixed, not π-rational | `15_`; `ChiralPhaseWindow.lean` |
| **`v_Higgs = 28²·a²`** | `= ‖Y_lepton‖²_Frobenius` (so(8)-democratic bilinear); the `28²` = #components.  Dimensionful "vacuum=bilinear" identification **open** | `03_`; `HiggsVevReframe.lean` |
| **`α` / `g_W²=5√α` / `α(0)`** | **value-conjecture, NOT a derivable law** (the `√α` form is `α(M_Z)` in disguise; `α(0)` rides a fitted `+2α`) | `04_`; `ALPHA_SCOPING.md` |
| **quark Koide, selection rule** | soft, RG-dependent; undriven | (`RIGOR_AUDIT.md`) |

**Genuine inputs (lepton+EW+Higgs):** the scale `a_lepton`, the amplitude `t²` (one number;
`φ` follows via the phase law, `19_`), and `α`.  Not the optimistic "1 input."

## The arc, with back-references

1. **v_Higgs** — `01_`,`02_` (reframings) → **`03_`**: the Cl(7) bridge.  `v_Higgs=dim(L)²a²` is the
   Frobenius² of the so(8)-democratic lepton mass bilinear (the "equipartition over directions"
   reading is *wrong* — gives `√28`; the *bilinear/component-count* reading is right).  Lean:
   `HiggsVevReframe.{frobeniusSq_democratic,bilinear_reading,vector_reading,readings_differ}`.
2. **α** — **`04_`** + `ALPHA_SCOPING.md`: three separate conjectures; `g_W²=5√α` is **not a law**
   (single-point coincidence with the definitional `e=g_W sinθ_W`); `α(0)`'s tight match rides a
   reverse-engineered `+2α`; `G_e` numerology.  α is a value-conjecture.
3. **Koide = ⅔ (the root)** — the equipartition story:
   - `EQUIPARTITION_PRINCIPLE.md`: Koide and `v_Higgs` are *one* mystery (maximal mixing).
   - **`05_`**: v58 energy is **flat** on the vacuum manifold (`⟨M²⟩₂≡0` on `L=skew`) → maximal
     symmetry must be primitive, not energy-selected.
   - **`06_`** (+ agent **`07_`**): `G₂` maximal-mixing gives `t²=(D−dimG₂)/D`; *mechanism retracted*.
   - **`08_`**: the `G₂×S₃` generation map — the bridge **fails** (`S₃` leaves `t²` free).
   - **`10_`**: grade-balance gives `t²=½` (right value, but a reframing).
   - **`11_`**: idempotent / criticality / color-singlet — **none force `½`** (recurring value `⅓`).
4. **The "½" itself (octomath)** — **`12_`** (missing-½ = idempotent normalization), **`13_`**: the
   `½` is the **root-product of the L-grade complex idempotent** `(1+u)/2`, `u²=−1` (vs `0` for the
   F-grade `u²=+1`); fixed by the grade square-sign; `mass∈L` proven.  Lean: `OctoHalf.lean`
   (`half_element_law`, `complex_half(_field)`, `real_half(_field)`, `root_product_complex`).
   **`14_`**: the `½` *is* carried to the amplitude (`t²=ξξ̃=½`, L-grade complex norm), but **not the
   phase** (the `45°` is the amplitude balance, not a phase).
5. **The phase (chiral)** — **`15_`**: M4 chiral window.  Electron-massless point `φ=π/12`; physical
   `φ=2/9<π/12` (light-but-massive).  A *constraint*, not a derivation.  Lean: `ChiralPhaseWindow.lean`.
6. **The audit** — **`09_`**: structure/magnitude audit.  Everything fixed is *discrete*; everything
   free is a *continuous magnitude*.  We used only octospace's *symmetry* (magnitude-blind); the
   *multivector* (product/associator/norm) is the magnitude-carrying content.
7. **color-3 vs generation-3 (Task #1)** — **`16_`**,**`17_`**: **distinct** — independent commuting
   factors of `Aut(𝕊)=G₂×S₃` (color block-diagonal, generation block-mixing); not unified by Spin(8)
   triality (the sedenion `S₃` is the Cayley–Dickson doubling automorphism).  So the fixed charge
   sector can't fix the free couplings.
8. **The dynamics (OFE)** — **`18_`**: lifting v58→Cl(7) adds the grade-4 (F) channel… **`19_`**:
   stacking the non-circular constraints collapses the sector to **one number `t²`** (`φ` follows)…
   **`20_`**: but the constrained OFE is **NEGATIVE** — the lepton is the **color singlet**, confined
   to the `ℍ`-slice where the F-channel is **identically zero**, so the OFE reduces to the `{0,2}=ℍ`
   `XiVacuum`, giving `|ξ|²=½` only as the Mexican-hat parameter (input) and the phase as a Goldstone
   (free).  (Corrects `18_`'s optimism.)

## The robust conclusion (`17_`, confirmed by `20_`)

> **Kinematics vs dynamics.**  The octonion algebra fixes the *discrete* data and grounds the
> *meaning* of the couplings (e.g. `½` = the L-grade complex norm).  Neither the kinematics
> (symmetry — magnitude-blind) nor the dynamics in hand (the OFE — flat on the color-singlet `ℍ`
> slice) fixes the *continuous values* `t²=½`, `φ=2/9`, `α`.  These are inputs/residuals — the same
> standing as the SM's Yukawa couplings (the flavor problem).

Every escape route closed with a specific reason: symmetry is magnitude-blind (`09_`); `G₂`/maximal-
mixing leaves `t²` free (`08_`,`11_`); the half-norm is the grade-balance restated (`14_`); the phase
is non-π-rational, near the chiral edge, free (`15_`); color ⊥ generation (`16_`,`17_`); the OFE
F-channel is inert on the `ℍ`-confined singlet (`20_`).

## What survives, and the floor

- **Theorem-grade (Lean, axiom-clean):** the discrete skeleton (Koide ratio = `dimG₂/dimSpin7`, grade
  structure, `5=h∨`, `sin²θ_W=2/9`); the `½` as the L-grade complex root-product (`OctoHalf`); the
  chiral window (`ChiralPhaseWindow`); the maximal-mixing arithmetic (`MaximalMixingKoide`); the
  v_Higgs reframings (`HiggsVevReframe`).  All in `AxiomCheck` (standard trio).
- **The floor:** the coupling *values* (`t²`, `φ`, `α`, `v/a²`) are **not derivable** from the
  octonion+v58 framework — they are inputs, or await a flavor dynamics outside it.  The reduction
  `19_` (sector = one number + `φ` free-rider) is the cleanest statement of what would need fixing.

## Lean modules added this session (all axiom-clean, imported in `AxiomCheck.lean`)
`HiggsVevReframe` (+vector/bilinear), `MaximalMixingKoide`, `OctoHalf`, `ChiralPhaseWindow`.

## Open
- **Task #1 done** (color-3 ≠ generation-3); its sub-lead (Spin(8) triality) also closed.
- The genuinely-open targets are *outside* this framework: a flavor dynamics that fixes `t²` (hence
  `φ` via the phase law), and the dimensionful `v_Higgs` "vacuum = bilinear" identification.
