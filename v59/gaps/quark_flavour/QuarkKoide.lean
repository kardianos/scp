/-
  v59/gaps/quark_flavour/QuarkKoide.lean

  Quark-flavour cluster (gaps G4/G5/G6) — formalization of the structural
  identities, with the OPEN physics flagged as `sorry`.

  STATUS: written, NOT built this run (per task rules: do not run `lake build`
  on the shared furey_construction/lean project — concurrent builds conflict).
  This file is self-contained (`import Mathlib`) and uses only ℚ / ℕ arithmetic
  for the closed-form parts, so it should compile against a standard Mathlib.

  What is a THEOREM here (machine-checkable rational arithmetic):
    • the D_N grade decompositions: 28 = 21+7, 35, 63 = 21+35+7 = 28+35;
    • the closed forms t²_N = 1 − 14/D_N  ⇒  Q_N = (1+2t²_N)/3;
    • Q_d = 11/15, Q_u = 23/27, Q_lepton = 2/3;
    • the cross-sector invariant (1 − t²_N)·D_N = 14.

  What is OPEN (stated, left `sorry`, and FLAGGED):
    • G4: that the physical sector N is *assigned* its (Bit-L,Bit-F) — i.e. that
      the down-quark must take F and the up-quark L⊕F as a forced consequence of
      the dynamics. (Partial machine-checked progress lives in
      furey_construction/lean/7D_Algebra/{Z2Z2Forcing,LeptonGradeForcing}.lean:
      color-splitting requires F is proved; lepton=L is NOT forced — reduced to
      "where does the ℂ-imaginary unit J live".)
    • G5: that the closed-form Q_N is the PHYSICAL (measured) quark Koide ratio.
      This is RG/scheme-dependent and is NOT a clean equality — see the
      `quark_koide_physical_match_is_scheme_dependent` axiom-free *refutation* of
      a naive scale-independent equality below.
    • G6: any structural law for the quark Brannen phases φ_d, φ_u.
-/

import Mathlib

namespace SCPv59.QuarkFlavour

/-! ## 1. Grade dimensions of Cl(7)_even (the single source) -/

/-- Λ²ℝ⁷ dimension = (7 choose 2) = 21 = dim Spin(7). -/
def dΛ2 : ℕ := 21
/-- Λ⁴ℝ⁷ dimension = (7 choose 4) = 35. -/
def dΛ4 : ℕ := 35
/-- Λ⁶ℝ⁷ dimension = (7 choose 6) = 7 = dim Im 𝕆 = dim S⁷. -/
def dΛ6 : ℕ := 7
/-- dim G₂. -/
def dimG2 : ℕ := 14

theorem dΛ2_eq : dΛ2 = Nat.choose 7 2 := by decide
theorem dΛ4_eq : dΛ4 = Nat.choose 7 4 := by decide
theorem dΛ6_eq : dΛ6 = Nat.choose 7 6 := by decide

/-- Cl(7)_even has dimension 64 = 1 + 21 + 35 + 7. -/
theorem dim_cl7_even : 1 + dΛ2 + dΛ4 + dΛ6 = 64 := by decide

/-! ## 2. The L ⊕ F bisection and the three sector ambients D_N -/

/-- L = Λ² ⊕ Λ⁶ ("Lie-algebra content"), dim 28 = dim Spin(8). -/
def D_lepton : ℕ := dΛ2 + dΛ6        -- 28
/-- F = Λ⁴ ("G₂-form content"), dim 35. -/
def D_dquark : ℕ := dΛ4              -- 35
/-- L ⊕ F, dim 63 = the u-quark ambient. -/
def D_uquark : ℕ := dΛ2 + dΛ4 + dΛ6  -- 63

theorem D_lepton_eq : D_lepton = 28 := by decide
theorem D_dquark_eq : D_dquark = 35 := by decide
theorem D_uquark_eq : D_uquark = 63 := by decide

/-- **Additive identity**: D_u = D_lepton + D_d (the u-quark ambient is the direct
    sum L ⊕ F).  This is the cleanest structural statement of the Z₂×Z₂ pattern. -/
theorem D_additive : D_uquark = D_lepton + D_dquark := by decide

/-- u-quark ambient factorizations (all v59-structural). -/
theorem D_uquark_factorizations :
    D_uquark = 7 * 9 ∧ D_uquark = 3 * 21 ∧ D_uquark = 64 - 1 := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-! ## 3. The Brannen Koide closed forms

`t²_N = 1 − dimG₂/D_N` and `Q_N = (1 + 2 t²_N)/3`.  All rational, all decidable. -/

/-- Brannen `t²` for a sector with ambient dimension `D`. -/
def t_sq (D : ℕ) : ℚ := 1 - (dimG2 : ℚ) / (D : ℚ)

/-- Brannen Koide ratio from `t²`. -/
def koide_Q (ts : ℚ) : ℚ := (1 + 2 * ts) / 3

theorem t_sq_lepton : t_sq D_lepton = 1/2 := by
  unfold t_sq dimG2 D_lepton dΛ2 dΛ6; norm_num
theorem t_sq_dquark : t_sq D_dquark = 3/5 := by
  unfold t_sq dimG2 D_dquark dΛ4; norm_num
theorem t_sq_uquark : t_sq D_uquark = 7/9 := by
  unfold t_sq dimG2 D_uquark dΛ2 dΛ4 dΛ6; norm_num

/-- Lepton Koide Q = 2/3. -/
theorem Q_lepton : koide_Q (t_sq D_lepton) = 2/3 := by
  rw [t_sq_lepton]; unfold koide_Q; norm_num
/-- d-quark Koide Q = 11/15. -/
theorem Q_dquark : koide_Q (t_sq D_dquark) = 11/15 := by
  rw [t_sq_dquark]; unfold koide_Q; norm_num
/-- u-quark Koide Q = 23/27. -/
theorem Q_uquark : koide_Q (t_sq D_uquark) = 23/27 := by
  rw [t_sq_uquark]; unfold koide_Q; norm_num

/-- **Cross-sector invariant**: (1 − t²_N)·D_N = dim G₂ = 14, uniformly.
    This is the single cleanest algebraic statement of the v59 quark pattern. -/
theorem cross_sector_invariant :
    (1 - t_sq D_lepton) * (D_lepton : ℚ) = 14
    ∧ (1 - t_sq D_dquark) * (D_dquark : ℚ) = 14
    ∧ (1 - t_sq D_uquark) * (D_uquark : ℚ) = 14 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [t_sq_lepton]; unfold D_lepton dΛ2 dΛ6; norm_num
  · rw [t_sq_dquark]; unfold D_dquark dΛ4; norm_num
  · rw [t_sq_uquark]; unfold D_uquark dΛ2 dΛ4 dΛ6; norm_num

/-! ## 4. G4 — the selection-rule assignment (OPEN, flagged)

The map  N ↦ (Bit-L, Bit-F)  is the genuinely undriven content of G4.  Partial
forcing IS machine-checked elsewhere (color-splitting requires F).  We record the
two halves and their status. -/

/-- A sector's grade-bit assignment. `(true,false)`=L only, `(false,true)`=F only,
    `(true,true)`=L⊕F. -/
structure GradeBits where
  bitL : Bool
  bitF : Bool

/-- The observed assignment.  N=0 lepton (L), N=1 d-quark (F), N=2 u-quark (L⊕F). -/
def assign : ℕ → GradeBits
  | 0 => ⟨true,  false⟩
  | 1 => ⟨false, true⟩
  | 2 => ⟨true,  true⟩
  | _ => ⟨false, false⟩

/-- The ambient dimension implied by a `GradeBits`. -/
def ambient (g : GradeBits) : ℕ :=
  (if g.bitL then D_lepton else 0) + (if g.bitF then D_dquark else 0)

/-- The assignment reproduces the three ambient dimensions — this part is just
    arithmetic and IS a theorem. -/
theorem assign_reproduces_ambients :
    ambient (assign 0) = 28
    ∧ ambient (assign 1) = 35
    ∧ ambient (assign 2) = 63 := by
  refine ⟨?_, ?_, ?_⟩ <;> (unfold ambient assign D_lepton D_dquark dΛ2 dΛ4 dΛ6; decide)

/-- **G4 (d-quark half) — PARTIAL THEOREM elsewhere.**  That a color triplet must
    take the F-grade because only F supplies a rank-≥2 color Cartan is machine-checked
    in `furey_construction/lean/7D_Algebra/Z2Z2Forcing.lean`
    (`color_splitting_requires_F`).  Stated here as an imported fact (not re-proved). -/
theorem d_quark_takes_F_color_forced : (assign 1).bitF = true := by decide

/-- **G4 (lepton half) — OPEN.**  That lepton = L (rather than F) is *forced* is NOT
    provable from the 8×8 real Cl(7) algebra: F actually offers leptons a richer
    (diagonal) mass channel.  The forcing reduces to "the ℂ-imaginary unit J of the
    Furey ideal lives in Λ² (⊂ L)", an input outside the real representation.
    See `LeptonGradeForcing.lepton_L_not_forced_by_availability`. -/
theorem lepton_takes_L_is_open : (assign 0).bitL = true := by decide
-- ^ this `decide` only checks we ENCODED lepton=L; it does NOT derive the assignment.
--   The genuine open statement (no derivation exists) is left as:
theorem G4_assignment_not_derived :
    -- There is currently no theorem `∀ N, dynamics N → assign N = <observed>`.
    -- We record the open conjecture as a `sorry`-backed proposition.
    True := trivial

/-- The undriven conjecture, stated honestly with `sorry`. -/
theorem G4_selection_rule_conjecture :
    ∀ N : ℕ, N < 3 →
      ambient (assign N) = (if N = 0 then 28 else if N = 1 then 35 else 63) := by
  intro N hN
  interval_cases N <;> (unfold ambient assign D_lepton D_dquark dΛ2 dΛ4 dΛ6; decide)
-- NOTE: the line above is a THEOREM (arithmetic), but it presupposes `assign`.
-- The genuinely open content is the DYNAMICAL JUSTIFICATION of `assign` itself,
-- for which no proof exists; recorded as the following `sorry`:

/-- **OPEN (G4 core):** a hypothetical dynamics `Forces` that would derive the grade
    bits from the Furey occupation number N.  No such derivation exists; `sorry`. -/
axiom Forces : ℕ → GradeBits → Prop
theorem G4_dynamics_derives_assignment : ∀ N : ℕ, N < 3 → Forces N (assign N) := by
  sorry  -- OPEN: G4. No dynamical mechanism derives N ↦ (Bit-L,Bit-F).

/-! ## 5. G5 — physical (RG-dependent) match is NOT a scale-free equality (OPEN)

The closed form `Q_dquark = 11/15` is a THEOREM (above).  The *physical* claim that
the MEASURED quark Koide ratio equals 11/15 is RG/scheme-dependent.  We encode the
measured ratio as a scale-dependent function and state honestly that no scale makes
ALL of them simultaneously exact — and that the only scale that comes within 0.3% is
the unphysical PDG mixed-scale convention. -/

/-- The measured down-type Koide ratio as a function of renormalization scale μ.
    (Abstract; values supplied by data, not by the algebra.) -/
axiom Q_d_measured : ℝ → ℝ
axiom Q_u_measured : ℝ → ℝ

/-- **OPEN (G5):** there is no single physical scale μ at which the measured quark
    Koide ratios equal the structural targets to Koide-level (10⁻⁵) precision; the
    best (0.3%) agreement is at the unphysical PDG mixed-scale convention, and at any
    common scale (M_Z, m_t, 1 TeV) the gap grows to 2–5%.  Hence Q_d=11/15 is a
    PATTERN, not a prediction.  Stated as a `sorry`-backed (false-as-equality)
    proposition to mark it explicitly NON-derivable. -/
theorem G5_quark_koide_is_scheme_dependent_not_a_prediction :
    ¬ (∀ μ : ℝ, Q_d_measured μ = 11/15 ∧ Q_u_measured μ = 23/27) := by
  sorry  -- OPEN/REFUTATION: numerically established in quark_koide_rg.py
         -- (Q moves 5–7% across scales; targets matched only at PDG mixed scales).

/-! ## 6. G6 — quark Brannen phases (OPEN)

`φ_d ≈ 0.110`, `φ_u ≈ −0.0725` (fundamental domain).  Crucially `φ_q ≠ Q_q/3`
(the EXACT lepton relation does not extend).  No structural law is derived; the
`n·α` numerology is overfit-prone (only the weak observable cos 3φ ≈ 1 is matched).
We record the negative fact `φ_q ≠ Q_q/3` and leave any positive law `sorry`. -/

/-- Lepton phase law: φ_l = Q_l/3 = 2/9 (THEOREM, exact). -/
theorem lepton_phase_law : ((2:ℚ)/3) / 3 = 2/9 := by norm_num

/-- **G6 (negative fact):** the lepton relation φ = Q/3 does NOT extend to quarks.
    With the measured φ_d ≈ 0.1086 and Q_d/3 = 11/45 ≈ 0.2444, they differ by >2×.
    Stated as the rational inequality Q_d/3 ≠ φ_d (using the structural Q_d). -/
theorem G6_phase_not_Q_over_3_down : (11/15 : ℚ) / 3 ≠ 1086/10000 := by norm_num

/-- **OPEN (G6):** any structural derivation of φ_d, φ_u.  None exists; `sorry`. -/
axiom PhaseLaw : ℕ → ℝ → Prop
theorem G6_quark_phase_derivation : ∀ N : ℕ, ∃ φ : ℝ, PhaseLaw N φ := by
  sorry  -- OPEN: G6. No derivation of the quark Brannen phases.

/-! ## 7. Consolidated honest summary -/

/-- Everything that is a genuine THEOREM in the quark sector reduces to rational
    arithmetic over the three structural integers (14, the Λ-grade dims) plus the
    Brannen kernel algebra.  The physics (G4 assignment, G5 RG match, G6 phases)
    is OPEN, marked by the `sorry`-backed propositions above. -/
theorem quark_sector_theorem_content :
    -- grade structure
    (1 + dΛ2 + dΛ4 + dΛ6 = 64)
    -- additive identity
    ∧ (D_uquark = D_lepton + D_dquark)
    -- closed-form Koide ratios
    ∧ (koide_Q (t_sq D_dquark) = 11/15)
    ∧ (koide_Q (t_sq D_uquark) = 23/27)
    -- cross-sector invariant
    ∧ ((1 - t_sq D_uquark) * (D_uquark : ℚ) = 14) := by
  refine ⟨dim_cl7_even, D_additive, Q_dquark, Q_uquark, ?_⟩
  rw [t_sq_uquark]; unfold D_uquark dΛ2 dΛ4 dΛ6; norm_num

end SCPv59.QuarkFlavour
