# Lean Formalization of the Stability Bounds — Current State and Roadmap (2026-05-24)

## Goal
Take the concrete mathematical bounds that emerged from the Maxima 8D sensitivity analysis (Hessians, eigenvalue conditions, protection projectors, radial vs transverse modes, the discrete |ξ|² values, the f_amplitude crossover, etc.) and turn them into statements that Lean can check or compute.

This is the "pave the way for their computation" part of the question.

## What Now Exists

- `lean/DensityForcing.lean` — the original high-level `Prop` statements (S³ radius as unique trade-off, LF selection forced by budgets, 21 as internal cost).
- `lean/StabilityBounds.lean` (new, 2026-05-24) — a more concrete module that directly encodes the Hessian-derived conditions:
  - 8-component masks as projectors
  - Simplified but explicit `EffectiveHessian` with the anisotropic form discovered in Maxima:
      - scalar mode: 2(μ + λ)
      - transverse modes (inside the projector): 2(μ − λ)
  - `eigenvalues_for_mask` and `is_stable` — the computable core
  - Example stability bound for concrete (λ, μ)
  - `RadialZeroCondition` + axiom that only the three observed radii work
  - Clear roadmap comments for the next steps

## Realistic Path Forward (Phased)

### Phase 1 — Axiomatic + Numerical Checks (current)
- Keep the key algebraic identities as `axiom` (the exact 8×8 Hessian coming from the octonion table, the precise form of the living-candidate quadratic terms, the force-separation/commutation bound).
- Make the eigenvalue check for *concrete numerical* (λ, μ, mask) tuples decidable inside Lean (using `Float` matrices + a tactic, or `Rat` + exact arithmetic).
- This already lets us *compute* the critical μ/λ ratio for a given mask, or verify that a particular f_amplitude crossover from the Python sweeps satisfies the stability condition.

### Phase 2 — Internalize More Algebra
- Replace the simplified diagonal `EffectiveHessian` with the real matrix whose entries come from the Fano multiplication table (the same table used in `octonion_sensitivity_analysis.mac`).
- Prove (or axiomatize + compute) the 7-fold degeneracy of the (μ − λ) eigenspace for quaternion-triple subalgebras.
- Turn the statement "only the observed |ξ|² = 1/2, 3/5, 7/9 make the radial eigenvalue cross zero while transverse modes stay positive" into a theorem with explicit bounds on λ and μ.

### Phase 3 — Full Computability
- Once we have a computable way to evaluate the Hessian for any mask and any (λ, μ), we can write a function that, given a mask, returns the minimal μ/λ ratio that makes the well stable.
- We can then prove (or at least machine-check for dense grids of values) that the observed discrete structures are exactly the ones that survive inside the physically allowed band (small |λ| ≤ 0.005, |μ| ≤ 0.001, commutation error < 0.5 %, etc.).
- Link this directly to the Python `f_amplitude` sweep results (the point ~0.4 where pure-L protection starts losing well depth on mixed sources becomes a theorem: "for f_amplitude > 0.4 the L-only Hessian has a negative eigenvalue on the F-contamination directions").

## Immediate Next Actions (suggested order)

1. Make `eigenvalues_for_mask` and `is_stable` actually computable for `Float` (or `Rat`) inputs so `lake build` succeeds and we can evaluate concrete cases inside Lean.
2. Add a small `example` or `#eval` that reproduces one of the concrete eigenvalue spectra the Maxima agent produced.
3. Write a short tactic or `norm_num` instance that can prove "all eigenvalues > 0" for a given numerical mask + (λ, μ).
4. Import the exact 8×8 matrix from the Maxima output (or hard-code the table) and repeat the above with the real (non-diagonal) Hessian.
5. State the f_amplitude crossover as a theorem once we have a model of "mixed source = L background + ε F contamination".

## Relation to Existing Work

- Directly continues the `DensityForcing.lean` skeleton.
- Is the Lean counterpart of `maxima/octonion_sensitivity_analysis.mac` (the file the Maxima agent heavily developed).
- The numerical side (Python `improved_density_protection_scan.py` + the `f_amplitude` sweep) supplies the concrete numbers and crossover points that the Lean statements should eventually certify.
- The high-level physical forcing relationships in `HYPOTHESES.md` and `CONSTRAINTS_AND_TROUGHS.md` are the informal version of what we are trying to make rigorous here.

This is the path from "we have observed these bounds symbolically and numerically" to "Lean can state them precisely and, for concrete parameters, compute and verify the critical values."

---

*Written after the two background agents completed their work on the Maxima and Python sides (May 2026).*