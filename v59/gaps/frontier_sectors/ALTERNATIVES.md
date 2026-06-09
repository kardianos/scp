# Opening the Solution Space — Frontier Sectors

Candidate v59-structural origins for the open quantities in G10/G11/G12.
Each candidate carries a **test** (what would confirm it) and a **falsifier**
(what would kill it). This is reconnaissance: most candidates are expected to
fail, and saying so cleanly is the point.

---

## G10 — CKM

### A1. Is `sin²θ_C = 7α` robust, or fitted? *(the one live lead)*
**Claim**: the Cabibbo mixing is set by `dim(Im 𝕆)·α = 7α(0)`.
- **Test (passed at 0.45%)**: `√(7/137.036) = 0.22601` vs `sinθ_C = 0.22500`.
- **Robustness test (`ckm_overfitting_test.py`)**: scan all `n·α`, `(n/m)·α`
  forms with small/v59 integers. **Result**: because `√α` makes `sin²θ_C≈0.05`
  *small*, very few `×α` forms land near it — `7·α` ranks #1 and only 2 forms
  sit within 2%. So among `×α` forms it is *not* generic. **But** this is a weak
  guarantee: widen to all `n/m` ratios (no α) and the target is easily hit.
- **Falsifier**: if `α(M_Z)` (the scale where weak mixing is defined) were the
  right α, the match degrades to 4% — so the conjecture *requires* the IR `α(0)`,
  which is physically odd for a quark-mixing quantity. A derivation must explain
  *why α(0)*. If no mechanism ties the Cabibbo angle to the IR photon coupling,
  treat `7α` as a **suggestive coincidence**, not a prediction.

### A2. `V_cb`, `V_ub`, `V_td` from sector overlaps of Brannen eigenvectors
**Claim**: the full CKM = `U_u† U_d` from the up/down Brannen kernels; the
hierarchy `λ : λ² : λ³` comes from the *grade overlap* of the u-quark ambient
`L⊕F` (dim 63) with the d-quark ambient `F` (dim 35).
- **Test**: Wolfenstein `A = V_cb/V_us² ≈ 0.826`; candidate `7/9 = t²_u = cos²θ_W`
  gives 0.778 (6% off); `4/5 = 28/35` (the L/F overlap ratio!) gives 0.80 (3% off).
  The `4/5 = dim L_λ⁶... / dim F`-type reading is *closer* than `7/9` and is
  itself v59-structural (28/35). Worth a principled-slice computation.
- **Falsifier**: build the ℍ-slice Brannen kernels with the slices *fixed by the
  μ-eigenspace projection* (not by hand) and compute `V_cb`. If the result is not
  `~0.04` (it was ~0.4–0.6 for arbitrary slices), the Brannen-overlap mechanism
  for the sub-Cabibbo elements is dead.

### A3. `δ_CP` from the Brannen phase
**Claim**: the CKM CP phase descends from the *quark* Brannen phases `φ_u, φ_d`
(the analog of the lepton `φ=2/9`).
- **Test**: `δ_CP ≈ 1.20 rad`. Candidates `π/√7 = 1.187` (1.0%), `π·7/18 = 1.222`
  (1.8%) are *suspiciously close* — but with the quark phases `φ_u,φ_d` already
  **[free]** (G6: `φ_q ≠ Q_q/3`), any of these is a *post-hoc fit*.
- **Falsifier**: the overfitting analysis flags BOTH `π/√7` and `π·7/18` within
  2% — two unrelated forms hitting the same target is the signature of
  overfitting, not signal. **Verdict: δ_CP is free.** A real claim would need
  `δ_CP` to be *predicted by* an independently-fixed `φ_u−φ_d`.

---

## G11 — Neutrinos + PMNS

### B1. Why neutrinos differ: a non-Brannen (non-circulant) sector
**Claim**: the ν ambient `D_ν` fits no v59 number (28/35/63) because neutrinos
are **not** in the Brannen-circulant charged-fermion family at all — they live
in a different algebraic slice (e.g. the *symmetric* `F=Λ⁴` reality sector, or a
Majorana bilinear that is order-2 not circulant).
- **Test**: a Majorana mass is a *symmetric* bilinear `ν^T C ν`; in v59 the
  symmetric sector is `Λ⁰⊕Λ⁴` (the `μ=+1` reality eigenspace), structurally
  distinct from the *skew* `L=Λ²⊕Λ⁶` that carries the complex (Dirac) structure.
  **Prediction to test**: the ν mass matrix should be the symmetric-grade
  bilinear, giving a *different* (non-Koide) eigenvalue structure — explaining
  why no `14/D` works.
- **Falsifier**: if the symmetric-grade bilinear still has a circulant/Koide form
  (i.e. `Q_ν` lands on some `14/D`), this distinction is empty.

### B2. Large PMNS angles from near-degeneracy (seesaw analog)
**Claim**: PMNS angles are large because the *light* neutrino masses are
near-degenerate (small splittings), and mixing of near-degenerate states is
generically maximal — the opposite of the strongly-hierarchical quarks.
- **Test**: tri-bimaximal `sin²θ = (1/3, 1/2, 0)` is the symmetric-group `S₃ / A₄`
  mixing pattern. v59 *has* an `S₃` (the sedenion `Aut(𝕊)=G₂×S₃` generation
  symmetry, `UNIFIED_THEORY.md` §2). **Hypothesis**: PMNS ≈ the `S₃`-symmetric
  (democratic) mixing matrix, while CKM is the *broken* `S₃` (the breaking that
  splits the quark masses also makes CKM near-diagonal).
- **Falsifier**: compute the mixing from the *exact* `S₃` democratic mass matrix.
  If it gives tri-bimaximal (`θ₁₃=0` exactly), it *over*predicts (real `θ₁₃≈8.5°`
  is nonzero) — so a *correction* term is needed; if no v59 term supplies the
  right `θ₁₃`, the `S₃` reading is only qualitative.

### B3. Seesaw scale `M_R` from a v59 dimension
**Claim**: the right-handed Majorana scale is `M_R = (v59 integer)·v` or a power
of `α` times `v`, giving `m_ν = m_D²/M_R` in the observed `0.05 eV` range.
- **Test**: observed needs `M_R ~ m_e²/m_ν ~ 5×10³ GeV` (if `m_D ~ m_e`). Is
  `M_R ~ v/√α`, `v·dim(L)`, …? `v·28 ≈ 6.9 TeV` is in range — **testable**.
- **Falsifier**: `m_D ~ m_e` is itself an assumption; if `m_D` runs with sector
  (u-quark-like vs lepton-like), `M_R` shifts by orders of magnitude. Without a
  v59 reason for `m_D`, this is dimensional analysis, not a prediction.

---

## G12 — Strong sector

### C1. `θ_QCD = 0` from the color complex structure `J_c` (best structural lead)
**Claim**: the strong-CP angle vanishes *structurally* because the color sector
carries a genuine **complex structure** `J_c = γ₀γ₅` (`ColorSU3.lean`: `J_c²=−I`,
skew, color-invariant, makes `ℝ⁸=ℂ⁴=1⊕3`). The θ-term `θ Tr(F∧F)` is **CP-odd**;
if the color field strength is built holomorphically w.r.t. `J_c`, the CP-odd
projection is forbidden by the *same reality/skewness forcing* that pins the
lepton complex structure (`LeptonRealityForcing`: "no symmetric matrix is a
complex structure").
- **Test**: explicitly compute the `J_c`-grading of the color bilinear `⟨F∧F⟩`.
  If `Tr(F∧F)` lies entirely in the `J_c`-*odd* (CP-even) eigenspace, `θ` has no
  home → `θ_QCD = 0` forced.
- **Falsifier**: if `⟨F∧F⟩` has a nonzero `J_c`-even (CP-odd) component, the
  mechanism fails and `θ_QCD` is free (the usual strong-CP problem persists).
- **Status**: this is the **most promising single lead in the whole cluster** —
  it would *explain* an observed near-zero, not fit a number, and it rests on an
  already-theorem-grade structure (`J_c`). Deserves real follow-up.

### C2. `α_s` from a color Killing/dual-Coxeter index — **TESTED, FAILS**
**Claim**: by analogy with `g_W² = 5√α = h∨(Spin(7))·√α`, the strong coupling is
`g_s² = h∨(SU(3))·√α = 3√α`.
- **Test (`pmns_strong_test.py`, `FrontierSectors.lean`)**: `3√α(0) ≈ 0.256` vs
  measured `g_s² = 4π·α_s(M_Z) ≈ 1.48`. **Off by >80%** (factor ~5.8).
- **Verdict: FALSIFIED.** The `√α` form gives couplings of order `√α ≈ 0.085` —
  it structurally *cannot* produce a strong (`O(1)`) coupling. Even `h∨=5` gives
  `0.43`, still 3× too small. The EW `g²=h∨√α` ansatz does **not** generalize.
- **What would be needed**: `α_s` is large *because* of asymptotic-freedom
  running, not a tree-level Killing index. A v59 reading would have to predict
  the QCD scale `Λ_QCD` (a *dimensionful*, RG-generated quantity), which is a
  qualitatively different kind of object than the v59 ratio-predictions. **No
  current v59 structure does this.** `α_s` is genuinely free.

### C3. `α_s` as a high-scale boundary condition (GUT-style unification)
**Claim**: at the unification scale all three gauge `g²` meet at `h∨·√α`-type
values, and `α_s(M_Z)` is the RG-*run* result.
- **Test**: do `g_W², g_R², g_s²` evaluated as `h∨·√α` agree at *some* scale?
  `g_W²=5√α`, `g_s²=3√α` differ by a fixed ratio `5/3` at the matching scale —
  compare to SM unification ratios.
- **Falsifier**: the SM `g₃²/g₂²` at `M_GUT` is `~1` (they unify), not `3/5`. So
  `h∨·√α` with `h∨=(5,3)` does **not** reproduce gauge unification either.
  Likely dead, but cheaper to check than C1.
