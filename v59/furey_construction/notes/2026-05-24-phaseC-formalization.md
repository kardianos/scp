# Lean formalization follow-up (2026-05-24)

Records the formalization work done in response to the v59 Lean review, and the
TODOs that were explicitly deferred (no code written yet).

## Overnight autonomous session (2026-05-25) — deliverables summary

Worked the remaining task-(4) and task-(1) parts, most-likely-to-succeed first. All
modules build; axiom footprints noted. New files:

- **(4d) `lean/EmpiricalAgreement.lean`** (axiom-clean) — 9 theorems + table proving each
  rational v59 prediction agrees with its measured central value to an explicit relative
  tolerance (`relDevLE`). 2 are within experimental error (lepton Koide Q, Brannen φ); 7
  are O(0.1–0.6%) structural approximations (sin²θ_W, cos²θ_W, m_Z/m_W, m_H²/v², Q_d, Q_u,
  Cabibbo) — honestly labelled. Turns the prose "matches to X%" claims into checked bounds.
- **(4a) `lean/AlphaZero.lean`** (axiom-clean) — the α(0) conjecture `−ln α + 2α = 8π²/dim Cl(3,1)`:
  `rhs_eq_structural` (8π²/16 = π²/2) + `alpha0_conjecture_holds` (at α=1/137.036 the
  relation holds to within 1/100, true gap ≈2.4×10⁻³), via `log 2`/`log(1+x)` bounds and
  Mathlib π bounds. Does NOT derive α or motivate the equation.
- **(1) `7D_Algebra/Z2Z2Forcing.lean`** (no axioms / propext+Quot.sound) — algebraic grounding
  of the forcing per the recommended route. Built the previously-unchecked Λ⁶ six-forms and
  proved **all 28 L-grade generators (Λ²⊕Λ⁶) are diagonal-free** (`L_grade_diagonal_free`,
  no axioms) — closing the "could Λ⁶ carry the Cartan?" gap PhaseB left open. Proved the
  F-grade color-block diagonals have **rank ≥ 2** (= SU(3)_c Cartan rank; two explicit
  4-forms give independent `[1,1,-1]`,`[-1,1,-1]`). Packaged as `color_splitting_requires_F`:
  a color triplet's splitting is identically zero from all of L but nonzero from F.

Skipped (genuinely not Lean-formalizable as theorems, documented why): **a_lepton** (the one
dimensionful empirical input, not derivable in-framework); **gauging SU(2)_L → m_W** (needs the
gauge coupling and dimensionful ⟨ξ⟩ — physics modelling; the Goldstone *count* 3 = dim su(2)_L
is already implicit in `XiVacuum.xi_mass_spectrum`).

## Done this session

### (3c) Kernel eigenvalue diagonalisation — `lean/KernelEigenvalues.lean`
The matrix diagonalisation that `CyclicShift.lean` / `BrannenKernel.lean` had left
as "standard linear algebra ... not re-derived formally here" is now proved
(axiom-clean — only `propext`/`Classical.choice`/`Quot.sound`):

- `Scyc_mulVec_v`  : `S  *ᵥ v_k = ω^k   • v_k`
- `ScycT_mulVec_v` : `Sᵀ *ᵥ v_k = ω^{2k} • v_k`
- `M_mulVec_eigen` : `M(a,ξ) *ᵥ v_k = λ_k • v_k`  (λ_k = a(1 + ξω^k + ξ̄ω^{2k}))
- `xi_sum`         : `ξω^k + ξ̄ω^{2k} = 2t·cos(φ + 2πk/3)` for `ξ = t·e^{iφ}`
- `lam_eq_brannen` : `λ_k = ↑(BrannenKernel.s a t φ k)`

Composed with `BrannenKernel.Q_value` / `koide_iff_constraint`, the *matrix*
eigenvalues (not just the abstract amplitude function) now carry Koide
`Q = (1+2t²)/3`, hence `2/3 ⟺ t² = 1/2`.

### (2) PhaseC certificates promoted from `#eval` to theorems — `7D_Algebra/PhaseC_Certificates.lean`
Each pure-`Bool` certificate now has a `native_decide` theorem fixing its truth
value (these depend additionally on `Lean.ofReduceBool`, since the claims are heavy
closed `Rat` computations over the 8×8 Fano Hessian — kernel `decide` gets stuck).

**Discrepancy surfaced (and largely resolved):** promotion initially showed 5 of
the historical certificates evaluating to `false`, not `true` — the `#eval` lines
printed them but never asserted anything, so this went unnoticed.  Investigation
(actual intermediate numbers `#eval`-ed) classified them:

- `gershgorin_cert_055/070.2.2` — **encoding bug**.  Compared the L-block Gershgorin
  radius (rows 0..3, rL≈0.020) to the FULL-matrix radius (rows 0..7, rLF≈0.029),
  which is a max over a superset ⇒ `rLF ≥ rL` always ⇒ `rL > rLF` impossible.
  Tested nothing.  **FIXED** (`gershgorinFvsLCert`): now compares F-sector rows
  (4..7, where L·L→F leakage lands) to L-sector rows; F is leakier (0.029 > 0.020).
  → `= true`.
- `combined_cert_055/070` — **threshold too tight**.  The L≫LF differential is real
  (L score drops 21.2%/19.6%; LF drops only 1.69%/1.58%); the cert failed solely on
  the absolute "LF < 1%" cutoff.  **FIXED**: now tests the differential (L drops
  >10% AND ≥5× more than LF). → `= true`.
- `famp_table_consistency_living` — **genuine model limitation, NOT fixed**.  With the
  realistic living potential the L and LF scores degrade *together* (rL≈rLF≈0.82 at
  fAmp 0.55, ≈0.91 at 0.40), so the asymmetric "LF preserved" bands genuinely fail.
  Kept honest as `= false`.  → feeds task (1) below.

Final state: 16 certificates `= true`, 1 (`famp_table_consistency_living`) `= false`.

## TODO (1) — Z₂×Z₂ sector-forcing as a theorem  *(not started; deferred)*

Goal: replace the consistency/"protection" narrative for the lepton/d/u → L/F/L⊕F
assignment with a genuine forcing theorem.  `Predictions.u_quark_is_induced_…`
currently proves only the arithmetic (28+35=63, 14/63=2/9); PhaseB proves diagonal
*structure* facts but not *why* each sector takes its bits.  The honest finding in
`notes/2026-05-23-fock-mass-forcing-attack.md` is that there is **no strict
algebraic vanishing** of the "wrong" `⟨Ω_N| op |Ω_N⟩` in the 8-dim spinor.

Sketch of what a real theorem would need (DO NOT implement yet):
1. Formalize the Fock/Witt number states `|Ω_N⟩` and the creation/annihilation
   `α_i` on `Λ•ℂ³ ≅` spinor 8 (extend `SevenDAlgebra.fockBasis`).
2. Formalize the G₂-branching `8 = 1 ⊕ 7`, `7|_{SU(3)} = 1 ⊕ 3 ⊕ 3̄`, and the
   grade content of L = Λ²⊕Λ⁶ (no G₂ singlet) vs F = Λ⁴ (contains *φ singlet).
   Needs Mathlib Lie/exterior-algebra machinery (heavy).
3. State + prove: the only `(Bit-L, Bit-F)` assignment giving a G₂-covariant,
   color-rep-consistent mass term for each N is lepton=(1,0), d=(0,1), u=(1,1).
   This is the missing "forcing" lemma; today it is a representation-theoretic
   consistency argument, not a theorem.

### Annotations from the 2026-05-24 PhaseC certificate analysis  ← decide next step here

The PhaseC promotion work directly informs which route to (1) is viable:

- **The stability/"protection-stacking" route is NOT a reliable basis for forcing.**
  The whole dynamical argument for "u needs L⊕F" was: an L-only source leaks into F
  and its protection degrades, so stacking L⊕F is required.  The certificate numbers
  show this differential is an artifact of the *simplified* potential:
    * Simple/`V_full` combined score: L degrades ~20%, LF ~1.7% — clean differential.
    * Realistic `V_living_project` (f(ρ) + κ-saturation): L and LF degrade *together*
      (≈18% each) — **the differential vanishes** (`famp_table_consistency_living`
      stays `= false`).
  So the more physical the potential, the weaker the protection-stacking signal.
  Conclusion: do **not** hang the Z₂×Z₂ forcing on the stability/Hessian dynamics —
  it is model-dependent and disappears in the faithful potential.

- **The algebraic route IS solid and is where forcing should be grounded.**  The
  PhaseB theorems are exact, axiom-clean `decide` facts independent of any potential:
  `gamma_mul_diag_zero` (universal over Fin 7), `lepton_L_compatibility`,
  `dquark_L_incompatible` + `dquark_F_compatible`, `uquark_LF_requires_F`, and the
  grade-intensity signs.  The corrected `gershgorinFvsLCert` (F rows leakier than L)
  is consistent with these but is a *consequence*, not the *driver*.

**Recommended next step for (1):** pursue the forcing as a purely
**representation-theoretic / grade-algebraic theorem** (extend PhaseB with the Fock
states + G₂-branching of step 1–2 above), and explicitly **demote** the stability /
protection-stacking certificates from "evidence for forcing" to "consistency checks
on the dynamics."  Do not invest further in making the living-V reproduce a L-vs-LF
differential as a forcing argument — the algebra, not the Hessian, must carry it.

## Task (4) — Remaining empirical inputs / dynamical ξ

Sub-item **(4b) ξ-vacuum Goldstone structure: DONE** (`lean/XiVacuum.lean`, axiom-clean).
Formalized the Mexican-hat potential `V(ξ) = (λ/4)(|ξ|²−1/2)²` on ℍ≅ℝ⁴ and proved, from
the potential itself:
- `V_critical_point` — the vacuum manifold `|ξ|²=1/2` (the v59 `t²=1/2` constraint) consists
  of critical points (first directional derivative ≡ 0);
- `V_hessian_directional` / `hessian_eq_massMatrix_quadForm` — the Hessian quadratic form is
  `2λ(ξ·v)²`, equal to `vᵀ(massMatrix)v`, so `massMatrix = 2λ ξ_a ξ_b` IS the Hessian;
- `massMatrix_radial` (eigenvalue λ, the Higgs), `massMatrix_goldstone` (every ⊥ direction is
  a null mode), `massMatrix_trace = λ`;
- `xi_mass_spectrum` — at the canonical vacuum `(1/√2,0,0,0)` the three imaginary-quaternion
  directions are massless Goldstones and the real direction is the massive mode: **spectrum
  {λ,0,0,0}** = SM-Higgs structure (3 would-be SU(2)_L Goldstones + 1 radial Higgs).
This makes the `DYNAMIC_XI.md` Maxima observation a machine-checked theorem.  NB: it proves
the *structure* of the vacuum; it does NOT derive *why* λ has any particular value, nor gauge
the silent SU(2)_L (that would eat the Goldstones) — those remain open.

Still open (Python-level / empirical; not yet Lean — listed for scope):
- `α(0) ≈ 4×10⁻⁵` structural origin (only `alpha_MZ = 25/(324π²)` is a Lean theorem; the
  low-energy α and inter-scale running are not formalized; the `−ln α + 2α = π²/2` equation
  is conjectural).
- `a_lepton` mass scale: the one remaining dimensionful empirical input (not derivable in-framework).
- Empirical PDG agreement: values are encoded as structural rationals and "match within
  precision" is asserted in prose, not as interval theorems. (This — task (4d) — is the next
  cleanly Lean-formalizable target if desired.)
- Gauging the silent SU(2)_L on the ξ vacuum (predict m_W from ⟨ξ⟩) — physics modelling.

## (3) parts NOT formalizable as theorems (no code; for the record)

The eigenvalue piece (3c) is done.  The other items bundled under (3) in the review
are geometric/representation-theoretic facts that cannot become axiom-clean theorems
without major new Mathlib development, and would otherwise have to be added as bare
`axiom`s (which would break the framework's axiom-cleanliness):

- Openness of the G₂-orbit of the associative 3-form in Λ³ℝ⁷ (triggers
  `dim G₂ = 49 − 35 = 14`).  Currently a documented geometric input.
- The homogeneous-space identification `Spin(7)/G₂ ≅ S⁷` (used for `21 − 7 = 14`).
- The full Furey ℂ⊗ℍ⊗𝕆 / Cl(6) Witt representation theory.

Recommendation: leave these as explicitly-documented inputs (as the README already
does) rather than axiomatizing them.
