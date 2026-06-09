# Quark-Flavour Cluster — Opening the Solution Space (G4, G5, G6)

**Date**: 2026-05-25 · Companion to `README.md` (background) and `FINDINGS.md` (verdict).

This document proposes avenues *beyond* the conjectures already tried, for each gap.
Each avenue has an explicit **test** (what computation/observation would support it)
and a **falsifier** (what would kill it). The goal is to move G4/G5/G6 from
numerology toward either a derivation or a clean refutation.

---

## A. G4 — DRIVING the assignment N ↦ (Bit-L, Bit-F)

The state of the art: **d-quark→F is [thm]** (color Cartan), **u-quark→L⊕F is
[thm]-ingredients + glue**, **lepton→L is [open], reduced to "where does J live"**.
So the productive avenues all target the lepton half or the dynamical wrapper.

### A1. The complex-structure (J-locating) avenue — *most promising*

`LeptonGradeForcing.lean` proved the open question reduces to: **does the ℂ-imaginary
unit `J` of the Furey minimal left ideal lie in Λ² (⊂ L) or in F?** The whole lepton
assignment hinges on this single fact.

- **Mechanism**: In Furey's `ℂ⊗𝕆` ideal construction, `J` is the `i` of the `ℂ`
  factor, fixed by the choice of primitive idempotent `½(1 + i e₇)` (or analog). The
  idempotent is built from a *bivector* `e₇`-type object; one must compute which
  Clifford grade the resulting `J = ad(idempotent)` action populates.
- **Test**: build the explicit `Cl(7)_even → M₈(ℂ)` isomorphism (the finite 64-element
  computation flagged in `13_fock_mass_forcing_report.md` §6), construct the Witt
  idempotent, extract J, and decompose J into Λ²/Λ⁴/Λ⁶. If `J ∈ Λ²`, lepton=L is
  forced [would become thm]; the Hermitian lepton mass is then `J·(L-antisym block)`.
- **Falsifier**: if `J ∈ Λ⁴ (F)` or J is grade-mixed across L and F, the assignment
  flips or is not grade-pure — and lepton=L is *refuted* as a forced statement (it
  would be a convention, demoting G4 entirely).
- **Why promising**: it is a *finite, well-posed* computation (no Lagrangian needed),
  and the prior work already isolated it as THE single missing input.

### A2. Chirality / Witt grading avenue

`F_chi = γ₃γ₄γ₅γ₆ = diag(1,1,1,1,−1,−1,−1,−1)` is the global Cl(6) chirality (Weyl
split N≤1 vs N≥2). The Fock-N grading and the L/F grading interleave through chirality.

- **Mechanism**: posit that the *mass term* must be chirality-odd (connects the two
  Weyl halves, as a Dirac mass does). L-generators are off-diagonal (chirality-odd in
  blocks); the color-neutral F diagonal `F_chi` is chirality-*even*. A Dirac-mass
  requirement would then prefer L for the color-neutral (lepton) channel and force
  the color structure into the off-diagonal F pieces for quarks.
- **Test**: classify all 28 L and 35 F generators by chirality parity on the {0,7}
  lepton pair and on the color blocks; check whether "Dirac mass = chirality-odd"
  uniquely selects L for leptons and (L⊕F) for the composite.
- **Falsifier**: if chirality-odd generators exist in F that give a color-neutral
  lepton mass, the requirement does not select L. (Note `LeptonGradeForcing` already
  shows F has a *symmetric* lepton block — need to check its chirality parity.)

### A3. Index / anomaly avenue

Color triplets (d, u) and singlets (lepton) differ by an SU(3)_c index. An
Atiyah-Singer / anomaly-flow argument might assign grades by a topological index of a
G₂- or Spin(7)-bundle restricted to each sector.

- **Mechanism**: the G₂→SU(3) branching `7 → 1⊕3⊕3̄` (in `15_su3_branching.py`) gives
  each Fock-N state an SU(3) weight. The number of zero modes of a Dirac operator
  twisted by the L-bundle vs the F-bundle could differ by the SU(3) index, selecting
  which grade hosts a non-vanishing mass.
- **Test**: compute the SU(3)_c index of the L-grade (Λ²⊕Λ⁶ = 8⊕3⊕3̄⊕… under SU(3))
  vs F-grade (Λ⁴) restricted to each Fock-N weight; check if the index vanishes
  exactly on the "forbidden" (Bit=0) combinations.
- **Falsifier**: if the indices are equal across L and F for a given N, no selection.
- **Caveat**: this is the most speculative; the v59 program's prior anchor-avenue
  attempts (the MEMORY "Avenue C anomaly: π₃(S⁷)=0" for gravity) found anomalies
  *don't* constrain free parameters — temper expectations.

### A4. SU(3)_c-content filtering — the H6 route, pushed to completion

H6 (`14_selection_rule_deep_think.py`) was "most physically motivated" and *evolved*
into the d-quark=F color-Cartan theorem. Push it to the lepton half.

- **Mechanism**: a color singlet needs *no* color Cartan, so by the d-quark theorem's
  logic it is "compatible with diagonal-free L" — but compatibility ≠ forcing. Sharpen
  to: a color singlet must couple to the **SU(3)-singlet** part of its grade. L's
  SU(3)-singlet content vs F's must be compared; whichever provides the unique
  color-neutral *mass* (not just any singlet) wins.
- **Test**: decompose L and F under SU(3)_c (not just G₂); count color-singlet mass
  channels in each; check that L's singlet channel is the one compatible with the
  lepton Brannen kernel's measured |ξ|²=1/2.
- **Falsifier**: F also has a color-singlet mass channel (it does — `*φ` and `F_chi`);
  so this avenue must explain why the lepton uses L's singlet *rather than* F's. This
  is exactly where A1 (J-location) becomes necessary — A4 alone is insufficient.

---

## B. G5 — an RG-INVARIANT reformulation of the quark mass relation

The fatal weakness of G5 is that Koide Q is **not RG-invariant** — it drifts 5–7%
across scales (`quark_koide_rg.py`). Any honest quark mass relation must either be
RG-invariant or be pinned to a *physical* scale (not the PDG mixed convention).

### B1. Fixed-scale anchoring (the conservative fix)

State the prediction at ONE physical scale and accept the resulting precision.

- **Mechanism**: choose μ = M_Z (or μ = m_t) as the canonical scale where v59 quantities
  are defined (consistent with α(M_Z), m_W, m_Z all being M_Z-scale in v59). Then the
  honest prediction is `Q_d(M_Z) ≈ 0.749`, `Q_u(M_Z) ≈ 0.888` — and the structural
  targets 11/15, 23/27 are then **2–4% predictions**, not 0.3%.
- **Test**: ask whether *any* pair of structural rationals from the v59 integer roster
  matches Q_d(M_Z), Q_u(M_Z) at <1%. (From `quark_koide_rg.py`: at M_Z, 11/15 is 2.0%
  low for Q_d; candidates like 3/4=0.75 are 0.2% — so the M_Z data prefers DIFFERENT
  rationals than the mixed-scale convention does.)
- **Falsifier**: if no clean rational matches at the chosen physical scale, G5 is dead
  as a prediction at that scale (this is the likely outcome — and an honest null).

### B2. An RG-invariant combination (the ambitious fix)

Find a function of the masses that IS one-loop RG-invariant and equals a structural
number. Quark mass *ratios* run mildly; certain combinations are protected.

- **Mechanism candidates**:
  - The one-loop QCD running multiplies all quark masses of a given type by a common
    factor `η(μ)` (flavor-universal anomalous dimension at leading order). **Koide Q is
    invariant under a common rescaling** `m_i → η m_i` — so at *leading-log, flavor-
    universal* running, Q *would* be invariant! The 5–7% drift is from (a) flavor
    *non*-universal threshold/electroweak corrections and (b) the mixed-scale inputs.
  - **Test**: recompute Q_d, Q_u with masses run by the *common* QCD factor only
    (strip the flavor-dependent and EW pieces). If the residual Q is scale-stable AND
    near 11/15, 23/27, then G5 is rescued: the relation holds for the QCD-universal
    part, and the drift is a known SM correction. This is a concrete, falsifiable
    computation.
  - **Falsifier**: if Q still drifts >1% under common-factor-only running, then the
    drift is intrinsic (flavor non-universal) and no RG-invariant rescue exists.
- **Why this is the single most valuable test in the cluster**: it directly addresses
  the "is it a prediction?" question. Koide's leptonic Q is famously RG-stable
  *because* leptons run flavor-universally (QED); the quark question is precisely
  whether the flavor-*non*-universal QCD/Yukawa running spoils it. (Literature: the
  charged-lepton Koide relation is one-loop RG-stable; the quark analog is known to be
  only approximately stable — this test quantifies "how approximately" for v59's D_N.)

### B3. Running-Koide that IS invariant (reformulate the observable)

Instead of fixing the masses' scale, define a scale-dependent target.

- **Mechanism**: if the v59 structural ratio is attached to a *fixed* algebraic D_N but
  the data runs, then perhaps the correct statement is `Q_N(μ*) = (1+2(1−14/D_N))/3`
  at the scale μ* where the Brannen kernel's vacuum is defined (the EW/compositeness
  scale, ~v). This ties the quark Koide to the *same* scale as v_Higgs = 28²a_l².
- **Test**: compute Q_d, Q_u at μ = v ≈ 246 GeV (close to m_t scale; from the table
  Q_d(1 TeV)≈0.748, Q_u≈0.896). Check against 11/15, 23/27.
- **Falsifier**: the EW-scale values (≈0.748, ≈0.896) are 2% and 5% off — so even the
  "natural v59 scale" does not give 0.3%. This avenue likely *confirms* G5 is soft.

---

## C. G6 — structural origin for φ_d, φ_u and the scale ratios

### C1. Effective-potential / loop-induced phases (the prior conjecture, sharpened)

`FINDINGS_brannen_dynamics.md` §3 proposed `V_eff(φ) = −A cos3φ + B cos6φ + …` with
critical points `cos(3φ) = A/(4B)`, the sector-specific ratio set by fermion loops.

- **Mechanism**: compute A_X, B_X from one-/two-loop fermion contributions in a
  complete Brannen-Yukawa Lagrangian. Different sector mass spectra → different A/B →
  different φ.
- **Test**: does the *measured* ratio `cos(3φ_X) = A_X/(4B_X)` come out of the loop
  integral with the sector's own masses as input (self-consistent)? Lepton:
  cos(2/3)=0.786; d: cos(1/3)=0.945; u: cos(−3/14)=0.978.
- **Falsifier**: if the loop integral gives A/4B independent of sector (or wildly off),
  the V_eff origin is dead. **Blocker**: requires the executable Lagrangian (G13),
  which does not yet exist — so this is currently untestable, not promising near-term.

### C2. The α-suppression as a *threshold/mixing* effect, not a phase law

The numerology `φ_d ≈ 14α, φ_u ≈ −10α` (dim G₂·α, dim Spin(5)·α) is suggestive but
overfit-prone. Reinterpret: a phase of order α is what you'd expect if the quark
Brannen phase is *radiatively generated* by EW/EM mixing (CKM-like), whereas the
lepton phase (2/9, order 1) is tree-level.

- **Mechanism**: the lepton phase is O(1) (tree, set by the G₂ ratio); the quark phases
  are O(α) (radiative). The integers (14, 10) would then be the multiplicities of the
  loop (number of gauge generators in the loop: dim G₂=14 for down, dim Spin(5)=10 for
  up). This *predicts the SIZE* (O(α)) even if the integer is fuzzy.
- **Test**: is `|φ_quark|/|φ_lepton| ≈ α·(integer)/(2/9)`? Numerically
  φ_d/φ_l = 0.1086/0.2222 = 0.489; if = 14α(M_Z)/(2/9) = 0.492 — match to 0.6%. And
  it predicts the *hierarchy* φ_quark ≪ φ_lepton, which is observed.
- **Falsifier**: if a full EW-loop calculation of the radiatively-induced Brannen phase
  does NOT come out O(α), the radiative reading is wrong. The integer's structural
  identity (14 vs 10) remains [conj] regardless.
- **Honest caveat**: this is a *size* argument (order-α), which is more defensible than
  a *value* argument (exactly 14α). Promote the size claim, demote the integer.

### C3. Scale ratios from the same Frobenius² bridge as v_Higgs

a_u²/a_l² ≈ 72 and a_u²/a_d² ≈ 35. The "35" = dim F = D_d, and 72 appeared in v59 as
`m_top = (1+2√(7/9))²·72·a_l²`. The 72 = 2·36 = 8·9 might tie to dim End or a Frobenius²
count like the 784 = dim End(L) of the v_Higgs bridge.

- **Mechanism**: if a_X² = (Frobenius² norm of the sector's mass bilinear)/(component
  count), the ratio a_u²/a_d² would be the ratio of ambient dims (D_u·something /
  D_d·something). 35 = D_d exactly — so a_u²/a_d² = D_d is the conjecture that the
  up-scale-squared/down-scale-squared = the down ambient dimension.
- **Test**: check a_u²/a_d² against D_d=35 at a *fixed physical scale* (not mixed). From
  `quark_phases_scales.py`: at M_Z it's 49, not 35 — so the "35" is a mixed-scale
  artifact, like everything else in this cluster.
- **Falsifier**: already falsified at single scales — the ratio is not 35 at M_Z or
  2 GeV. **Verdict: dead** as a single-scale relation. Same convention artifact as G5.

---

## Priority ranking (across all three gaps)

| rank | avenue | gap | why | near-term feasible? |
|---|---|---|---|---|
| 1 | **A1 — locate J in the grade decomposition** | G4 | reduces lepton=L to one finite computation; the prior work already isolated it | **YES** (finite 8×8) |
| 2 | **B2 — common-QCD-factor RG-invariance test** | G5 | directly answers "prediction or pattern?"; falsifiable now | **YES** (RG running) |
| 3 | A2/A3 — chirality / index selection | G4 | could upgrade the lepton half if A1 is ambiguous | partial |
| 4 | C2 — radiative O(α) phase *size* (not value) | G6 | defensible size argument; predicts φ_q≪φ_l hierarchy | needs EW loop |
| 5 | B1/B3 — fixed-scale anchoring | G5 | honest fallback; likely a null (2–5%) | YES (already shown) |
| — | C1 (V_eff), C3 (scale ratios=35) | G6 | C1 blocked on Lagrangian; C3 already falsified | no / dead |

**Bottom line**: the two highest-value, *near-term-feasible* moves are
**(A1)** finish the J-location computation to settle the lepton half of G4, and
**(B2)** test whether common-factor QCD running leaves Q invariant — which would
either rescue G5 as a genuine (if loop-corrected) prediction or kill it cleanly.
