/-
  v76 Approach A — LocalityBudget (Round 1 fragment)

  Minimal monist ledger formalization (no fixed metric, no T→G):
  - MediumBudget: free + bound components
  - Budget identity: free + bound = total
  - Depletion lemma: fixed total + more bound ⇒ less free
  - Mass relation: m * c² = E_bound  (c = free-field locality scale)
  - Zero-mass free medium; lock-step depletes free

  Standalone Lean 4 (no Mathlib). Uses Int + omega.
  Not a full dynamics / inertia proof — ledger half only.

  Build (from this directory): lake build
-/

namespace V76.LocalityBudget

/-! ## Medium budget (ledger of one continuum) -/

/-- Local continuum budget: free vs bound are states of one medium.
    Values are integer ledger units (discrete / scaled continuum). -/
structure MediumBudget where
  free  : Int
  bound : Int
  free_nonneg  : 0 ≤ free
  bound_nonneg : 0 ≤ bound

/-- Total budget at a point / region. -/
def MediumBudget.total (b : MediumBudget) : Int :=
  b.free + b.bound

/-- Strong budget identity holds by definition of `total`. -/
theorem budget_identity (b : MediumBudget) :
    b.free + b.bound = b.total :=
  rfl

/-! ## Depletion: monist content of free + bound = total -/

/-- If two budgets share the same total and one has strictly more bound,
    it has strictly less free. Ledger form of
    "mass-form consumes free fabric." -/
theorem free_depleted_of_bound_gt
    (b₁ b₂ : MediumBudget)
    (hTot : b₁.total = b₂.total)
    (hBound : b₁.bound < b₂.bound) :
    b₂.free < b₁.free := by
  simp only [MediumBudget.total] at hTot
  omega

/-- Non-strict form: more bound at fixed total ⇒ free is ≤. -/
theorem free_le_of_bound_ge
    (b₁ b₂ : MediumBudget)
    (hTot : b₁.total = b₂.total)
    (hBound : b₁.bound ≤ b₂.bound) :
    b₂.free ≤ b₁.free := by
  simp only [MediumBudget.total] at hTot
  omega

/-! ## Rest mass from locality scale c

  Propositional form avoids division:
    IsMass c E_bound m  ↔  m * c * c = E_bound
  i.e. E_bound = m c² with c = free-field locality (Ax3).
-/

/-- `IsMass c E m` means rest mass m for bound energy E at locality c:
    E = m c². -/
def IsMass (c E m : Int) : Prop :=
  m * c * c = E

/-- Recover E = m c² from IsMass. -/
theorem energy_eq_mass_mul_c_sq {c E m : Int} (h : IsMass c E m) :
    E = m * c * c :=
  h.symm

/-- Trivial unit choice c = 1: mass equals bound energy. -/
theorem mass_eq_energy_at_c_one {E m : Int} (h : IsMass 1 E m) :
    E = m := by
  simp only [IsMass] at h
  omega

/-- Uniform free medium: no bound ⇒ zero rest mass works for any c. -/
def freeOnly (ρ : Int) (hρ : 0 ≤ ρ) : MediumBudget where
  free := ρ
  bound := 0
  free_nonneg := hρ
  bound_nonneg := by omega

theorem mass_freeOnly (c ρ : Int) (hρ : 0 ≤ ρ) :
    IsMass c (freeOnly ρ hρ).bound 0 := by
  simp [IsMass, freeOnly]

/-- Fully locked budget (extreme): free = 0, bound = total. -/
def fullyLocked (ρ : Int) (hρ : 0 ≤ ρ) : MediumBudget where
  free := 0
  bound := ρ
  free_nonneg := by omega
  bound_nonneg := hρ

theorem free_zero_of_fullyLocked (ρ : Int) (hρ : 0 ≤ ρ) :
    (fullyLocked ρ hρ).free = 0 :=
  rfl

/-- Bound energy of a fully locked region of ledger ρ satisfies mass relation. -/
theorem fullyLocked_mass (c ρ m : Int) (hρ : 0 ≤ ρ)
    (h : IsMass c ρ m) :
    IsMass c (fullyLocked ρ hρ).bound m := by
  simpa [fullyLocked] using h

/-- Converting free → bound by δ at fixed total: free drops, bound rises. -/
theorem lock_step_depletes
    (free₀ bound₀ δ : Int)
    (hf : 0 ≤ free₀) (hb : 0 ≤ bound₀) (hd : 0 ≤ δ)
    (hfit : δ ≤ free₀) :
    let b₁ : MediumBudget :=
      { free := free₀, bound := bound₀
        free_nonneg := hf, bound_nonneg := hb }
    let b₂ : MediumBudget :=
      { free := free₀ - δ, bound := bound₀ + δ
        free_nonneg := by omega, bound_nonneg := by omega }
    b₁.total = b₂.total ∧ b₂.free ≤ b₁.free ∧ b₁.bound ≤ b₂.bound := by
  intro b₁ b₂
  refine And.intro ?_ (And.intro ?_ ?_)
  · simp only [MediumBudget.total, b₁, b₂]; omega
  · simp only [b₁, b₂]; omega
  · simp only [b₁, b₂]; omega

/-- Lock step with δ > 0 strictly decreases free budget. -/
theorem lock_step_strict
    (free₀ bound₀ δ : Int)
    (hf : 0 ≤ free₀) (hb : 0 ≤ bound₀)
    (hδ : 0 < δ) (hfit : δ ≤ free₀) :
    let b₁ : MediumBudget :=
      { free := free₀, bound := bound₀
        free_nonneg := hf, bound_nonneg := hb }
    let b₂ : MediumBudget :=
      { free := free₀ - δ, bound := bound₀ + δ
        free_nonneg := by omega, bound_nonneg := by omega }
    b₂.free < b₁.free := by
  intro b₁ b₂
  simp only [b₁, b₂]
  omega

/-- Instance: free_depleted follows from a positive lock step. -/
theorem lock_step_depletes_via_lemma
    (free₀ bound₀ δ : Int)
    (hf : 0 ≤ free₀) (hb : 0 ≤ bound₀)
    (hδ : 0 < δ) (hfit : δ ≤ free₀) :
    let b₁ : MediumBudget :=
      { free := free₀, bound := bound₀
        free_nonneg := hf, bound_nonneg := hb }
    let b₂ : MediumBudget :=
      { free := free₀ - δ, bound := bound₀ + δ
        free_nonneg := by omega, bound_nonneg := by omega }
    b₁.total = b₂.total ∧ b₂.free < b₁.free := by
  intro b₁ b₂
  have htot : b₁.total = b₂.total := by
    simp only [MediumBudget.total, b₁, b₂]; omega
  have hbound : b₁.bound < b₂.bound := by
    simp only [b₁, b₂]; omega
  exact ⟨htot, free_depleted_of_bound_gt b₁ b₂ htot hbound⟩

end V76.LocalityBudget
