# v59 Gap Attack — Consolidated Synthesis

**Date**: 2026-05-25. Six agents worked the 12 live gaps of `UNIFIED_THEORY.md` §9 (G13,
the dynamical Lagrangian, was excluded). Each folder holds `README.md` (background),
`ALTERNATIVES.md` (opened solution space, each avenue with test + falsifier), numerical
scripts, a Lean module (written, not built this run — no concurrent `lake build`), and
`FINDINGS.md`. This file is the cross-cutting summary.

## Verdict table

| gap | folder | verdict | status change |
|---|---|---|---|
| G2 `g_W²=5√α` form | `alpha_couplings` | **DEAD as independent law** — = α(M_Z) "in disguise" (one point where a line meets a √). The "5"=h∨ survives [thm]. | demoted: G2 ≡ G3 |
| G3 `α` value | `alpha_couplings` | **OPEN** — genuine dimensionless input. Tight within a form-family, large look-elsewhere across families. | unchanged; best avenue = RG fixed point |
| G1 `v=28²a²` | `ew_scale_bridge` | **CONJECTURE, sharply localized.** 784=dim End(L) is forced [thm, generic]. New **rank tension** is the real obstruction (below). | sharpened |
| G7 `Q=2/3` phase mag | `lepton_phase` | **Magnitude LIVE** (non-geometric: max-mixing G₂-content + equipartition t²=½). `cos(2/3)` as algebra invariant **NULL**. Radian-insert **OPEN**. | upgraded (live readings) + reframed (skewness) |
| G4 selection rule | `quark_flavour` | d→F **[thm]**, u→L⊕F thm-up-to-glue, **lepton→L reduced to one decidable question** (is the ℂ-unit J∈Λ²?). | sharpened |
| G5 quark Koide | `quark_flavour` | **DEMOTED to soft pattern** — the "0.3%" is a scale-convention artifact; single-scale gap 2–5%, RG spread 16–27× the gap, 12 rationals fit. | downgraded (honest) |
| G6 quark phases/scales | `quark_flavour` | **FREE.** φ_q≠Q_q/3; n·α numerology weak; scale ratios are convention artifacts. | unchanged |
| G8 gravity magnitude | `gravity` | **Exponent 21=dimSpin(7) data-forced** (0.0024%); α^(21/2)=det(√α·I₂₁)=top-form=instanton e^{−S}, S=21·ln(1/α). Mechanism **not derived**; prefactor weak (17/13 beats 21/16). | sharpened (instanton/det reading) |
| G9 scalar→tensor | `gravity` | **FATAL & decisive.** v59 gravity is a Lorentz scalar → h=0 only; the internal so(8)/Λ² index is inert under spacetime rotations → no h=±2. **Fails LIGO** without an induced-metric (spacetime-symmetric) mode. | **promoted to top blocker** |
| G10 CKM | `frontier_sectors` | `sin²θ_C=7α(0)` survives 0.45% but **coincidence-grade**; δ_CP hit by 2 unrelated forms (overfit signature). V_ub,V_td,δ_CP free. | unchanged |
| G11 neutrinos/PMNS | `frontier_sectors` | **GENUINELY OPEN.** ν-ambient puzzle confirmed; best hypothesis = non-Brannen Majorana Λ⁰⊕Λ⁴; PMNS≈tri-bimaximal (S₃) but predicts θ₁₃=0 (wrong). | unchanged |
| G12 strong | `frontier_sectors` | **Key null:** `g²=h∨√α` does NOT generalize (g_s²=3√α≈0.26 vs 1.48). **One live lead:** θ_QCD≈0 may be *forced* by the color complex structure J_c [thm] — explains a zero, not fits a number. | new lead + new null |

## Three cross-cutting results (the most important outputs)

1. **The rank tension links G1 and G8/gravity.** The EW bridge wants a **full-rank** democratic
   `Y∈End(L)` so that `‖Y‖²_F = 784 a²` (all 784 entries ~a). But the physical spectrum is **3
   generations** = a **rank-3** `Y`, whose `‖Y‖²_F = Σm = 9Qa²` — which is exactly the **gravity
   charge** (`gravity_charge_test`), and is `(6/784)v` for leptons. So the *same* mass bilinear
   `Y` is read two incompatible ways: full-rank (EW, 784) vs rank-3 (gravity/spectrum, Σm). The
   "25 unaccounted directions" is the symptom, and no single-step `so(8)→H` breaking gives 3 light
   directions (no 25-dim subalgebra). **This is now the sharpest open problem of the lepton/EW
   block**, and it ties the two most-developed sectors (EW bridge and gravity charge) together.

2. **G9 is the decisive blocker and is sector-fatal.** A scalar source can never produce the
   `h=±2` LIGO polarizations, and the internal bivector index does not carry spacetime helicity.
   So *the gravity magnitude (G8) is moot until G9 is solved* — a beautifully-tuned `α²¹` scalar
   force still fails. The only escape is an emergent/induced-metric recast (Plebański-style, with
   `Ω∈Λ²` as the fundamental 2-form `B` and `h_μν` derived). **Prioritize G9 over G8.**

3. **Two "matches" demoted by honest analysis.** G2 (`g_W²=5√α`) is not an independent law (≡
   α(M_Z)), and G5 (quark Koide 11/15, 23/27) is a scale-convention artifact (2–5% at a physical
   scale, not 0.3%). Both were load-bearing-sounding; both are now correctly weak. Conversely, two
   things *strengthened*: the G7 magnitude (live non-geometric readings) and a new G12 lead
   (θ_QCD=0 possibly forced by J_c — explaining an observed zero, the strongest kind of result).

## Revised gap priority (after the attack)

1. **G9 (tensor mode)** — decisive; gravity is scalar and fails LIGO. Test the induced-metric recast.
2. **G1 rank tension** — the EW/gravity-shared `Y` cannot be both full-rank (784) and rank-3 (Σm);
   resolving this is the real content of R1/R2 and connects to the gravity charge.
3. **G3 (α)** — the one genuine dimensionless input; RG fixed-point is the live avenue (subsumes G2).
4. **G7 radian-insert** — why `Q` re-enters as the phase argument (the magnitude readings are live).
5. **G12 θ_QCD via J_c** — a clean "explain a zero" lead worth a dedicated check.

Quark Koide (G5), CKM (G10), neutrinos (G11) remain soft/open; treat as patterns, not predictions.

## Verification status (2026-05-25, actually built/run — no longer "written, not built")

**Core Lean project** (`furey_construction/lean`, Lean 4.29.0 + Mathlib): `lake build`
**succeeds** (8278 jobs); `AxiomCheck` rebuilt after the `OctoBridge` removal — every headline
theorem uses only `{propext, Classical.choice, Quot.sound}` or no axioms.

**Gap modules** (typechecked via `lake env lean` against the project Mathlib):

| module | errors | `sorry`-backed (flagged-open) | fully proved? |
|---|---|---|---|
| `AlphaCouplingIdentities` | 0 | 1 (the open α derivation) | structural identities yes |
| `EwScaleBridge` | 0 | 1 (R2 value-pin) | dim End(L)=784, ratios yes |
| `FrontierSectors` | 0 | **0** | **fully** |
| `G8G9_Gravity` | 0 | 1 (the α²¹ empirical match) | structural identities yes |
| `LeptonPhaseMagnitude` | 0 | **0** | **fully** (max-mixing, equipartition, skewness) |
| `QuarkKoide` | 0 | 3 (G4 selection, G5 RG-refutation, G6 phases) | closed-forms yes |

The `sorry`s sit **only** on the genuine conjectures (α, R2, the α²¹ match, G4/G5/G6) — the
build confirms these are the *sole* unproved steps; all structural identities are machine-checked.
Four proof bugs the agents could not catch (they never built) were fixed during verification:
`gravity_prefactor` (missing `norm_num`), `endL_dim` (`finrank_matrix` API), `R2_component_scale`
(`OfNat (Fin dimL)` on opaque `dimL` → all-pairs-equal form), and 3× `Fin.val` reductions in the
skewness identity. **All 14 gap Python scripts run clean.** So "backed by Lean and Python" is now
literal, not aspirational — the prerequisite for the version-closing summary is met.
