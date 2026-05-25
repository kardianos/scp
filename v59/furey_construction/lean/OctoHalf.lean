/-
  v59/furey_construction/lean/OctoHalf.lean

  **The "missing 1/2" as the L-grade complex-structure root-product** (octomath grounding of
  the Koide amplitude `t²=1/2`).

  Pushing connection (3) of `integration_v58/12_missing_half_idempotent.md` with "octomath":
  the `1/2` is NOT an arbitrary ordinary-algebra value — it is the **root-product (constant term
  of the monic quadratic)** of the *half-element* `P=(1+u)/2`, and its value is fixed entirely by
  the **square-sign of `u`** (a grade/octonionic fact, `BladeSquareSign`):

    * `u² = −1`  (a COMPLEX structure — the **L-grade** `Λ²⊕Λ⁶`, where the lepton mass is
      *proven* to live, `LeptonRealityForcing`):  `P² = P − 1/2`, monic poly `X²−X+1/2`,
      **root-product = 1/2**.
    * `u² = +1`  (a REAL structure — the **F-grade** `Λ⁴`):  `P² = P` (genuine idempotent),
      monic poly `X²−X`, **root-product = 0**.

  So "the lepton mass is L-grade (complex)" ⟹ "its half-element carries root-product `1/2`",
  whereas an F-grade (real) mass would carry `0`.  The `1/2` is the **complex-structure
  signature**.  This is signature-independent (holds in any ring; the octonionic/Clifford content
  is *which* elements are complex structures — the grade square-sign law).

  HONEST SCOPE.  This formalises *why the value is `1/2`* (the L-grade root-product), grounding it
  in the proven `mass ∈ L`.  It does NOT by itself transport this internal `1/2` to the *generation*
  Brannen amplitude `t²` — that is the still-open generation↔internal map (`08_…md`).  But it
  removes the "arbitrary `1/2`" objection: the value is forced by the grade.
-/
import Mathlib

namespace SCPv59.OctoHalf

/-! ## The signature-independent half-element law (any ring) -/

variable {R : Type*} [Ring R]

/-- **Half-element law.**  In *any* ring, `(1+u)² = 2(1+u) + (u²−1)`.  Everything about the
    half-element `(1+u)/2` is controlled by `u²`. -/
theorem half_element_law (u : R) :
    (1 + u) * (1 + u) = 2 * (1 + u) + (u * u - 1) := by
  noncomm_ring

/-- **L-grade (complex structure, `u²=−1`):** `(1+u)² = u + u  (= 2(1+u) − 2)`.  The constant
    `−2` (≠ 0) is the root-product signature (`P=(1+u)/2 ⇒ P² = P − 1/2`). -/
theorem complex_half (u : R) (h : u * u = -1) :
    (1 + u) * (1 + u) = u + u := by
  have e : (1 + u) * (1 + u) = 1 + u + u + u * u := by noncomm_ring
  rw [e, h]; abel

/-- **F-grade (real structure, `u²=+1`):** `(1+u)² = 1 + 1 + u + u  (= 2(1+u))` — idempotent
    doubling, root-product 0. -/
theorem real_half (u : R) (h : u * u = 1) :
    (1 + u) * (1 + u) = 1 + 1 + u + u := by
  have e : (1 + u) * (1 + u) = 1 + u + u + u * u := by noncomm_ring
  rw [e, h]; abel

/-! ## The `1/2` explicitly (over any field, char ≠ 2): the root-product of the half-element -/

variable {K : Type*} [Field K] [CharZero K]

/-- **Complex (`u²=−1`): `P=(1+u)/2` satisfies `P² = P − 1/2`** — the monic quadratic
    `X² − X + 1/2`, whose constant term (root-product) is the **`1/2`**. -/
theorem complex_half_field (u : K) (h : u * u = -1) :
    ((1 + u) / 2) * ((1 + u) / 2) = (1 + u) / 2 - 1 / 2 := by
  field_simp
  linear_combination h

/-- **Real (`u²=+1`): `P=(1+u)/2` is idempotent** (`P² = P`) — monic `X² − X`, root-product `0`. -/
theorem real_half_field (u : K) (h : u * u = 1) :
    ((1 + u) / 2) * ((1 + u) / 2) = (1 + u) / 2 := by
  field_simp
  linear_combination h

/-- **The root-product is the complex/real (L/F) signature.**  For `P=(1+u)/2`:
    `u²=−1 ⇒ P²−P+1/2 = 0` (root-product **1/2**); `u²=+1 ⇒ P²−P = 0` (root-product **0**). -/
theorem root_product_complex (u : K) (h : u * u = -1) :
    ((1 + u) / 2) * ((1 + u) / 2) - (1 + u) / 2 + 1 / 2 = 0 := by
  rw [complex_half_field u h]; ring

theorem root_product_real (u : K) (h : u * u = 1) :
    ((1 + u) / 2) * ((1 + u) / 2) - (1 + u) / 2 = 0 := by
  rw [real_half_field u h]; ring

/-- The two root-products are genuinely different: `1/2 ≠ 0`.  So the complex (L-grade) half
    carries `1/2` and the real (F-grade) half carries `0` — the `1/2` is the L-grade signature. -/
theorem root_products_differ : (1 / 2 : ℚ) ≠ 0 := by norm_num

end SCPv59.OctoHalf
