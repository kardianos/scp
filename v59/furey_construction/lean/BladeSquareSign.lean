/-
  v59/furey_construction/lean/BladeSquareSign.lean

  The **complex-structure grade law**, proved universally (representation-independent),
  completing the lepton = L forcing.

  `7D_Algebra/LeptonComplexStructure.lean` proved — by exhaustive `decide` over the
  genuine Cl(7) — that every L-grade blade (Λ²⊕Λ⁶) squares to `−I` (is a complex
  structure) while every F-grade blade (Λ⁴) squares to `+I` (a real structure), so the
  complex structure `J` of the Brannen lepton phase `ξ = e^{iφ}` must be L-grade.

  This module proves the *reason*, for **any** ring, not just the 8×8 Cl(7) realization:
  the product of `k` pairwise-anticommuting elements each squaring to `−1` squares to
  `(−1)^{k(k+1)/2}`.  Hence among even grades the complex structures (square `−1`) are
  exactly grades `k ≡ 2 (mod 4)` — i.e. `Λ²` and `Λ⁶` (= L), never `Λ⁴` (= F).  This is
  the structural, basis-free counterpart of the enumerative Cl(7) check, in the same
  spirit as `CliffordBladeGrade.lean` (which makes the *composite*-half forcing universal).
-/
import Mathlib

namespace SCPv59.BladeSquareSign

variable {A : Type*} [Ring A]

/-! ## Centrality of `(−1)^n` -/

/-- `(−1)^n` is central in any ring. -/
lemma neg_one_pow_central (n : ℕ) (x : A) : Commute ((-1 : A) ^ n) x :=
  ((Commute.one_left x).neg_left).pow_left n

/-! ## Moving a generator past a product of anticommuting generators -/

/-- If `a` anticommutes with every element of `l`, then moving `a` across the whole
    product picks up the sign `(−1)^{|l|}`:  `a · ∏ l = (−1)^{|l|} · (∏ l) · a`. -/
lemma comm_past (a : A) (l : List A) (h : ∀ b ∈ l, a * b = -(b * a)) :
    a * l.prod = (-1 : A) ^ l.length * (l.prod * a) := by
  induction l with
  | nil => simp
  | cons b t ih =>
    have hb : a * b = -(b * a) := h b (by simp)
    have ht : ∀ x ∈ t, a * x = -(x * a) := fun x hx => h x (by simp [hx])
    have hc : Commute ((-1 : A) ^ t.length) b := neg_one_pow_central _ b
    rw [List.prod_cons, ← mul_assoc, hb, neg_mul, mul_assoc, ih ht, List.length_cons, pow_succ]
    rw [show b * ((-1 : A) ^ t.length * (t.prod * a))
          = (-1 : A) ^ t.length * (b * (t.prod * a)) from by
          rw [← mul_assoc, ← hc.eq, mul_assoc]]
    noncomm_ring

/-! ## The triangular exponent -/

/-- The triangular number `k(k+1)/2`, as a binomial coefficient (clean Pascal recurrence). -/
def tri (n : ℕ) : ℕ := (n + 1).choose 2

@[simp] lemma tri_zero : tri 0 = 0 := rfl
lemma tri_two : tri 2 = 3 := rfl
lemma tri_four : tri 4 = 10 := rfl
lemma tri_six : tri 6 = 21 := rfl

/-- `tri (n+1) = (n+1) + tri n`  (Pascal's rule). -/
lemma tri_succ (n : ℕ) : tri (n + 1) = (n + 1) + tri n := by
  unfold tri; rw [Nat.choose_succ_succ (n + 1) 1, Nat.choose_one_right]

/-! ## The universal square-sign law -/

/-- **The blade square-sign law (representation-independent).**  For a list `l` of
    pairwise-anticommuting elements, each squaring to `−1`, the product squares to
    `(−1)^{tri |l|} = (−1)^{|l|(|l|+1)/2}`.

    In a Clifford algebra `l` is the list of distinct generators of a simple blade, so
    this is exactly `B² = (−1)^{k(k+1)/2} I` for a grade-`k` blade. -/
theorem prod_sq (l : List A)
    (hsq : ∀ x ∈ l, x * x = (-1 : A))
    (hanti : l.Pairwise (fun a b => a * b = -(b * a))) :
    l.prod * l.prod = (-1 : A) ^ (tri l.length) := by
  induction l with
  | nil => simp
  | cons b t ih =>
    rw [List.pairwise_cons] at hanti
    obtain ⟨hb_anti, ht_anti⟩ := hanti
    have hsq_t : ∀ x ∈ t, x * x = (-1 : A) := fun x hx => hsq x (by simp [hx])
    have hbb : b * b = (-1 : A) := hsq b (by simp)
    have hc : Commute ((-1 : A) ^ t.length) b := neg_one_pow_central _ b
    have key : t.prod * b = (-1 : A) ^ t.length * (b * t.prod) := by
      have hcp := comm_past b t hb_anti
      have h2 : (-1 : A) ^ t.length * (b * t.prod)
          = (-1 : A) ^ t.length * ((-1 : A) ^ t.length * (t.prod * b)) := by rw [hcp]
      rw [h2, ← mul_assoc, ← pow_add, ← two_mul, pow_mul]; simp
    rw [List.prod_cons, List.length_cons]
    calc (b * t.prod) * (b * t.prod)
        = b * (t.prod * b) * t.prod := by noncomm_ring
      _ = b * ((-1 : A) ^ t.length * (b * t.prod)) * t.prod := by rw [key]
      _ = (-1 : A) ^ t.length * (b * (b * t.prod) * t.prod) := by
            rw [show b * ((-1 : A) ^ t.length * (b * t.prod))
                  = (-1 : A) ^ t.length * (b * (b * t.prod)) from by
                  rw [← mul_assoc, ← hc.eq, mul_assoc]]
            noncomm_ring
      _ = (-1 : A) ^ t.length * ((b * b) * (t.prod * t.prod)) := by noncomm_ring
      _ = (-1 : A) ^ t.length * ((-1 : A) * (-1 : A) ^ (tri t.length)) := by
            rw [hbb, ih hsq_t ht_anti]
      _ = (-1 : A) ^ (t.length + 1 + tri t.length) := by rw [pow_add, pow_add]; noncomm_ring
      _ = (-1 : A) ^ (tri (t.length + 1)) := by rw [tri_succ]

/-! ## Grade corollaries: L = Λ²⊕Λ⁶ are complex structures, F = Λ⁴ is not

The even grades of `Cl(7)_even` are `{0, 2, 4, 6}`.  These three corollaries cover the
nontrivial ones and split them exactly into L (square `−1`) and F (square `+1`). -/

/-- **Λ² blades are complex structures.**  A simple bivector (grade 2) squares to `−1`. -/
theorem prod_sq_grade_two (l : List A) (hlen : l.length = 2)
    (hsq : ∀ x ∈ l, x * x = (-1 : A))
    (hanti : l.Pairwise (fun a b => a * b = -(b * a))) :
    l.prod * l.prod = (-1 : A) := by
  rw [prod_sq l hsq hanti, hlen, tri_two]; exact Odd.neg_one_pow (by decide)

/-- **Λ⁴ blades are real structures.**  A simple 4-form (grade 4) squares to `+1` —
    so it is an involution, *not* a complex structure.  This is why `F = Λ⁴` cannot
    carry the lepton phase's `J`. -/
theorem prod_sq_grade_four (l : List A) (hlen : l.length = 4)
    (hsq : ∀ x ∈ l, x * x = (-1 : A))
    (hanti : l.Pairwise (fun a b => a * b = -(b * a))) :
    l.prod * l.prod = (1 : A) := by
  rw [prod_sq l hsq hanti, hlen, tri_four]; exact Even.neg_one_pow (by decide)

/-- **Λ⁶ blades are complex structures.**  A simple 6-form (grade 6) squares to `−1` —
    the other half of `L = Λ²⊕Λ⁶`. -/
theorem prod_sq_grade_six (l : List A) (hlen : l.length = 6)
    (hsq : ∀ x ∈ l, x * x = (-1 : A))
    (hanti : l.Pairwise (fun a b => a * b = -(b * a))) :
    l.prod * l.prod = (-1 : A) := by
  rw [prod_sq l hsq hanti, hlen, tri_six]; exact Odd.neg_one_pow (by decide)

/-- **The universal complex-structure / grade dichotomy.**  Among the even-grade simple
    blades of `Cl(7)_even`:

      grade 2 (Λ²)  and  grade 6 (Λ⁶)  square to `−1`  — **complex structures** (= L),
      grade 4 (Λ⁴)                     squares to `+1` — a **real structure** (= F).

    Representation-independent: it is the parity of `k(k+1)/2`, not any feature of the
    8×8 octonion matrices.  This is the structural reason the Brannen lepton phase's
    complex structure `J` (`J² = −1`) is forced into `L`, never `F`.

    (Over `ℝ`/`ℤ`-matrices — the actual Cl(7) realization, char 0 — the two types are
    genuinely distinct since `−1 ≠ +1`; cf. `LeptonComplexStructure.id8_ne_negId8`.
    In characteristic 2 they would coincide, which is why distinctness is left to the
    concrete realization rather than asserted here.) -/
theorem even_grade_complex_structure_dichotomy
    (l₂ l₄ l₆ : List A)
    (h₂ : l₂.length = 2) (hsq₂ : ∀ x ∈ l₂, x * x = (-1 : A))
      (ha₂ : l₂.Pairwise (fun a b => a * b = -(b * a)))
    (h₄ : l₄.length = 4) (hsq₄ : ∀ x ∈ l₄, x * x = (-1 : A))
      (ha₄ : l₄.Pairwise (fun a b => a * b = -(b * a)))
    (h₆ : l₆.length = 6) (hsq₆ : ∀ x ∈ l₆, x * x = (-1 : A))
      (ha₆ : l₆.Pairwise (fun a b => a * b = -(b * a))) :
    (l₂.prod * l₂.prod = (-1 : A))      -- Λ² ∈ L : complex structure
    ∧ (l₆.prod * l₆.prod = (-1 : A))    -- Λ⁶ ∈ L : complex structure
    ∧ (l₄.prod * l₄.prod = (1 : A)) :=  -- Λ⁴ = F : real structure
  ⟨prod_sq_grade_two l₂ h₂ hsq₂ ha₂,
   prod_sq_grade_six l₆ h₆ hsq₆ ha₆,
   prod_sq_grade_four l₄ h₄ hsq₄ ha₄⟩

end SCPv59.BladeSquareSign
