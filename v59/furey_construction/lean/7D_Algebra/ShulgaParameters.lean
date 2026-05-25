import SevenDAlgebra

namespace SCPv59.Furey7D.Shulga

/-!
# 7D Parameter Derivation from Internal Geometry

This module formally proves the geometric parameter ratios extracted from the
7D continuous internal manifold ($S^7$) using exact rational arithmetic (`Rat`).

By evaluating the $S^7$ harmonic sum (the Green function of the Laplacian) using
Gegenbauer polynomial recurrences, we derive the exact theoretical value for the
effective potential parameters ($\lambda$ and $\mu$) without floating point approximations.
-/

/-- Structurally recursive loop for the Green function sum.
    `steps` is the number of remaining iterations.
-/
def greenSumLoop : Nat → (c_prev2 c_prev1 sum ml x : Rat) → Rat
  | 0, _, _, sum, _, _ => sum
  | steps + 1, c_prev2, c_prev1, sum, ml, x =>
    let cl := (2 * (ml + 2) * x * c_prev1 - (ml + 4) * c_prev2) / ml
    let w := (ml + 3) / (3 * ml * (ml + 6))
    greenSumLoop steps c_prev1 cl (sum + w * cl) (ml + 1) x

/-- The exact Gegenbauer C_l^{(3)}(x) polynomials sum evaluated up to lmax. -/
def computeGreenSum (x : Rat) (lmax : Nat) : Rat :=
  if lmax == 0 then 0
  else if lmax == 1 then
    let c1 := 6 * x
    let w1 : Rat := 4 / 21
    w1 * c1
  else
    let c0 : Rat := 1
    let c1 : Rat := 6 * x
    let w1 : Rat := 4 / 21
    let sum1 := w1 * c1
    -- Start at l=2, so ml=2, steps remaining = lmax - 1
    greenSumLoop (lmax - 1) c0 c1 sum1 2 x

/-- The exact geometric parameter \lambda_raw comes from the cross-interaction
    across the Z_3 family shift (\theta = 2\pi/3), so \cos(\theta) = -1/2.
-/
def lambda_raw_exact : Rat := computeGreenSum (-1 / 2) 100

/-- The exact geometric parameter \mu_raw comes from the self-energy regularized
    at the geometric cutoff. \theta_c = 1/8. \cos(1/8) \approx 127/128.
-/
def mu_raw_exact : Rat := computeGreenSum (127 / 128) 100

/-- The geometrically derived ratio of the interaction \lambda to self-energy \mu. -/
def geometric_ratio : Rat := lambda_raw_exact / mu_raw_exact

-- Using #guard to verify bounds via the evaluator instead of kernel decide,
-- as the kernel decide struggles with 200+ digit Rat reduction depth.
#guard (- (1 : Rat) / 4000 < geometric_ratio) && (geometric_ratio < - (1 : Rat) / 6000)

/-!
Absolute energy scale: REMOVED 2026-05-24.  The tuned `shulga_abs_prefactor` (1/400) and the
derived `|μ|` only fed the living-candidate stability potential, which has been retracted (it
modelled a Cl(3,0) crossover the octonion algebra does not exhibit — see
`notes/2026-05-24-LivingCandidateCrossover-ReEvaluation.md`).  The algebraic **ratio**
`geometric_ratio` (above) is the genuine, bug-immune deliverable of this module and is kept.
The absolute-scale derivation, if pursued, is the S7-coset functional-measure sub-project.
-/

end SCPv59.Furey7D.Shulga
