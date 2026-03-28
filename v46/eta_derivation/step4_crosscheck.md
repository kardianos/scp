# Step 4: Cross-Check and Analysis

## The Virial Approach Failed (as applied)

The virial residual from global diagnostics gives R_sol = -252,254 — a
huge number dominated by background subtraction errors. The background
energy (~136k per component) dwarfs the soliton contribution (~90 for
E_pot). Subtracting two 136k numbers to get a ~200 residual produces
massive errors.

**The virial approach requires per-voxel background subtraction on the
SFA grid data, not global diagnostic subtraction.** The global diagnostics
are not precise enough for this calculation.

## The Force Balance Approach is More Reliable

The force balance estimate (Step 1, section 1g) doesn't require background
subtraction — it operates at the braid core where the background is
irrelevant:

    η₁ ≈ (1/η₀ - η₀) / A_core³

This gives:

| A_core | η₁ | η_eff at |P|=0.1 |
|--------|-----|---------|
| 0.20 | 187 | 19.2 |
| 0.25 | 96 | 10.1 |
| 0.30 | 56 | 6.1 |
| 0.35 | 35 | 4.0 |

The V34/V43 data shows core phi_rms ≈ 0.09-0.32 depending on the
measurement (V43 proton formation: core depleted to 0.33× initial).
The effective A_core is in the range 0.2-0.3.

## The Conjecture η₁ = |μ|/η₀ = 82.7

This falls within the force balance range (corresponds to A_core ≈ 0.22).
However, the virial cross-check shows it's NOT derivable from global
energy balance — the required ⟨|P|⟩ = 57 is unphysical (|P| is bounded
at ~0.15).

**The conjecture is numerologically suggestive but NOT confirmed by
the virial analysis.** It may be approximately correct but the exact
relationship η₁ = |μ|/η₀ is not established.

## What IS Established

1. **η₁ is in the range 35-200** from the force balance argument.
   This is a robust estimate that doesn't depend on background
   subtraction or global diagnostics.

2. **The force balance derivation IS a first-principles result**:
   at the braid core, the curl force must balance the gradient
   restoring force. This gives η₁ in terms of η₀ and A_core.

3. **A_core is itself determined by the other parameters** (m², μ, κ).
   The braid profile is the solution of the nonlinear ODE:
   ∇²φ - m²φ - V'(P)∂P/∂φ = 0 (static equilibrium).
   So A_core = f(m², μ, κ), making η₁ = g(m², μ, κ, η₀).

4. **η₁ is NOT a free parameter** — it's determined by the requirement
   that the braid is in force equilibrium with the modified coupling.
   The numerical value depends on the precise braid profile, which we
   can compute from simulation or from the radial ODE.

## Conclusion

η₁ is constrained to 35-200 by force balance. The exact value requires
either:
- (a) Solving the self-consistent radial ODE with η(P) coupling
- (b) Per-voxel virial computation from SFA data (not global diagnostics)

The conjecture η₁ = |μ|/η₀ ≈ 83 is in the correct range and has an
appealing physical interpretation, but is not rigorously derived.

**For simulation testing, η₁ = 80 is a well-motivated starting point.**
A sweep η₁ = {40, 80, 120, 200} would bracket the force balance range
and identify the binding threshold.
