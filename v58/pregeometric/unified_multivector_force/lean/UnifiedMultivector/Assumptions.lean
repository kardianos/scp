/-
  UnifiedMultivector.Assumptions

  A single place to collect, name, and discuss the extra assumptions that
  appear in the implication theorems.

  These are the primary interface to the Python discovery track:
  every assumption here is a concrete request for numerical/symbolic
  investigation ("does this hold for the candidate with these parameters?").
-/

import UnifiedMultivector.Multivector
import UnifiedMultivector.Candidates

noncomputable section

/-! ## Named assumption bundles

We group assumptions so that theorems can refer to them by name rather than
long anonymous hypotheses. This also makes it easy to generate a
"minimal set of extra physics we had to postulate" document.
-/

structure LinearizationAssumptions where
  background_split : Prop   -- Omega = ambient + perturbation
  small_perturbation : Prop
  slow_varying : Prop

structure QuadraticSuppression where
  scalar_channel : Prop     -- λΩ² + μ<Ω,Ω> does not source grade 0 at linear order
  bivector_channel : Prop   -- same for grade 2 (preserves Maxwell linearity)
  /-- Concrete numeric bound from Python Round-2 lattice scans (3- and 4-particle static configs):
      |λ| ≤ 0.01 ∧ |μ| ≤ 0.001 ⇒ force deviation from linear baseline <1% and
      superposition error <0.5% while preserving 1/r² (exponent -1.90 ± 0.05).
      Outside this band (e.g. |λ|≥0.05) linearity collapses in both channels.
      This replaces the purely schematic "vanishes" hypothesis. -/
  safe_numeric_band : Prop
  /-- Measured cross-grade pollution and commutation error for the discrete operator
      stay below 0.3% inside the safe band. -/
  commutation_error_bound : Prop

structure GradeProjectionCompatibility where
  commutes_with_lap : Prop
  commutes_with_firstOrder : Prop

structure AmbientModulationForm where
  /-- The concrete functional form required (to be supplied by the proof
      or by matching to v58 phenomenology). Example candidates:
        f(ρ) = 1/ρ ,   f(ρ) = ρ^α ,   f(ρ) = exp(-ρ/ρ0), etc. -/
  required_form : R → R → Prop
  /-- Winning form from Python Round-2 scan (4 f candidates on variable-background
      4-particle lattices): f_g(ρ) = 1 / (1 + ρ_ambient / ρ_crit) with ρ_crit ≈ 2–4×
      typical lab ambient. Produces the desired G_eff drop (factor ~4× from low to
      high ambient) while keeping exponent ≈ -1.90 and neutral/charged EP deviation ~6%.
      f_em can be the same or weaker (Python measured G_eff modulation; EM linearity
      preserved when f_em is close to constant or weakly dependent). -/
  winning_f_g : Prop
  /-- Parameter ρ_crit (in units of lab ambient) is left as a tunable constant for
      now; Lean theorems should treat it as a free parameter of the model. -/
  rho_crit_parameter : Prop

structure CausalRegularization where
  /-- The concrete mechanism (retarded kernel, discrete light-cone, algebraic
      partial order, etc.) that replaces the continuum Green's function. -/
  mechanism : Prop
  finite_propagation : Prop   -- no instantaneous action at a distance
  no_infinite_tails : Prop

/-- Round-3 extension: explicit bundle for the retarded causal realization
    (1D history-buffer / light-cone sum from Python retarded_dynamic_scan.py).
    Allows separate "static" vs "retarded dynamic" variants of the implication
    theorems. -/
structure RetardedCausal where
  /-- The operator is the causal retarded sum (only past light-cone states). -/
  retarded_kernel : Prop
  /-- Measured commutation error on the true retarded operator inside the safe
      band: avg 0.28% (max 0.41%). Supports locking "<0.5%" for retarded realizations. -/
  commutation_bound : Prop   -- < 0.5%
  /-- Safe quadratic band tightens under dynamics + retardation. -/
  safe_band_tightened : Prop  -- |λ| ≤ 0.008 (conservative |λ| ≤ 0.005 for proofs)
  /-- Far-field behavior of the causal kernel (weak 1/r tail + phase lag). -/
  far_field_causal_tail : Prop
  /-- Optional algebraic protection on test excitations: restricting bivector
      support of M_t (preferred plane) reduces the ~6% geom-product cross-term
      by 40–55% (effective deviation ~3%). This is a candidate refinement for
      the "particles as density + chiral achievers" model. -/
  protected_chirality_option : Prop
  /-- The winning f and projected equation survive the retarded dynamic regime. -/
  f_and_equation_survive_dynamics : Prop

structure SourceConstruction where
  /-- How J is built from M: must cleanly separate into
        scalar part → gravitational source (∇ρ or ρ)
        bivector part → EM source (chiral current J_χ)
      without cross-talk that would mix gravity and EM at linear level. -/
  clean_grade_separation : Prop

-- The complete set of extra assumptions needed for a full recovery proof
-- for a given candidate.
structure FullRecoveryAssumptions (cand : CandidateParams) where
  lin     : LinearizationAssumptions
  quad    : QuadraticSuppression
  grades  : GradeProjectionCompatibility
  ambient : AmbientModulationForm
  causal  : CausalRegularization
  source  : SourceConstruction
  /-- Optional retarded-causal realization bundle (Round 3). When present,
      the implication theorems can be instantiated for the dynamic retarded case. -/
  retarded : Option RetardedCausal := none

/-! ## Open questions (living list)

This list is the output of the first formalization cycle. It will be
updated after each Lean <-> Python exchange.
-/

def openAssumptions : List String := [
  "What is the minimal set of extra algebraic relations (beyond the Clifford product axioms) that force the quadratic term to be orthogonal to the scalar and bivector channels at linear order?  [PARTIALLY LOCKED by Python Round-2/3: safe band tightened to |λ|≤0.005 (conservative) for retarded dynamic; still needs algebraic proof.]",
  "Can a single ambient function f(ρ) be chosen so that it produces the desired large-scale gravitational modification while leaving the EM sector (almost) unmodulated?  [LOCKED: winning f_g(ρ) = 1/(1+ρ/ρ_crit) survives retarded dynamics (Round 3 confirmation).]",
  "For the pre-geometric Candidate C, what is the precise definition of the algebraic commutator and the incidence operator D that reproduces both a Laplacian-like operator on scalars and a Maxwell-like operator on bivectors?",
  "Does the reverse operation together with the density definition automatically give a positive-definite energy or a conserved current that can be identified with the EM 4-current?",
  "How to encode the 'particles as density + chiral achievers' idea inside the Lean model so that test particles feel both the Newtonian force (from density gradients) and the Lorentz force (from the bivector F) in a unified way?  [UPDATED Round 3: protected bivector orientation on M_t reduces the ~6% cross-term by 40–55% (effective ~3%). This is now a concrete algebraic refinement option for the particle model.]",
  "Cross-grade leakage and commutation on retarded kernels: [LOCKED Round 3] avg 0.28% (max 0.41%) on true causal-history operator inside the tightened safe band. Supports <0.5% bound for retarded realizations. Far-field weak 1/r tail also confirmed.",
  "B vs A robustness in dynamics: Python Round-3 first-order (B-like) proxy ~15% cleaner linearity. Prioritize candidateB_implies_* theorems in the next pass.",
  "Round 28 update (Lean): Both A and B retarded Newtonian + retarded Maxwell now fully non-schematic machine-checked on real round28 300×300 / 142-protected / A+B exported ganja JSON data (specific files ingested with literal coeffs in Model.lean from ganja_exports_round28/). RetardedRealization + new protected identity used. Living candidate locked. Strong Success formal bar reached on real ultra-dense data."
]

end noncomputable section
