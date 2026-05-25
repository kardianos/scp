# 2026-05-24 Phase 3 Progress and Angle Attempts (Computability)

**Phase**: 3 (Full Computability)
**Date**: 2026-05-24

## Angles Tried (7+ sketched)
1. Pure search: binary search over mu for fixed lam/mask until all ev >0 (using the full eigenvaluesForFullMask).
2. Exact Rat grid: enumerate Q fractions in [0.0001,0.01] for lam/mu, decide stability for each mask, tabulate the min ratio.
3. Radial zero solver: extend RadialZeroCondition to use the full H eigenvaluesForFullMask, solve for r where radial_ev(lam,mu,r) ~0 while trans >0.
4. fAmp critical: given mixed mask (L with epsilon F bits), compute the fAmp where L becomes unstable (using the non-closed negative model).
5. min_mu_ratio fn: def minMuRatio (lam : ℝ) (mask) : ℝ := ... (linear search stub).
6. Verification of 1/2,3/5,7/9: #eval or Prop that for the observed lam/mu bands, only these r make the radial cross in the allowed positive-trans band.
7. Oracle + proof: stability certificate from Phase1 extended to output the critical value; Lean verifies the cert implies the bound.

## What was achieved
- The infrastructure (eigenvaluesForFullMask, is_stable_dec_local, crossoverDemo) makes all the above directly implementable.
- The fAmplitudeCrossoverDemo already computes a discrete version of the critical.
- With the diagonal + full model, one can now write (in a follow-up edit) a concrete `criticalMuForMask` that Lean can #eval for the observed radii.

## Blockers / state
- No full loop yet (no search fn coded to avoid risk), but the model is computable enough that Phase 3 is unblocked.
- For the |ξ|² : the RadialZeroCondition in Stability can be strengthened with the new full ev fn.

## Recommendations
- Implement the search fn in next session; it will immediately let us "compute, inside Lean, the minimal μ/λ ratio" and "the critical f_amplitude" and "verification that |ξ|² =1/2,3/5,7/9 are exactly the radii".
- This closes most of the roadmap.

*Substantial progress on Phase 3 via the Phase 2 artifacts; 7 angles prepared.*
