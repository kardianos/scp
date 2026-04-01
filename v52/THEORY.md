# V52 Theoretical Results: Chirality-Dependent Hardening

## Summary

Algebraic analysis (Maxima) and formal proofs (Lean) show that the
V50/C4 equations produce a chirality-dependent interaction barrier
between composite baryons. The proton-neutron (UUD+UDD) pair experiences
19-69% less hardening repulsion than the proton-proton (UUD+UUD) pair,
while feeling identical binding attraction. This emerges purely from the
chirality geometry — no new terms or parameters.

## The Three Results

### Result 1: Triple Product P is Chirality-Independent

The binding potential V(P) = (μ/2)P²/(1+κP²) depends on P = φ₀φ₁φ₂.
At the midpoint between two composites approaching along x:

  P(UUD+UUD) = P(UUD+UDD) at all separations D

**Ratio = 1.000 exactly**, verified numerically at D = 1..50.

**Mechanism**: At the midpoint (y=z=0), chirality enters the field as
cos(χ_b × k × coord_b). For braids 1,2 (y,z axes), coord=0, so
cos(χ×k×0) = cos(0) = 1 regardless of chirality. For braid 0 (x-axis),
both UUD and UDD have χ₀ = +1 (both "Up" along x).

**Lean proof**: `P_independent_of_chirality` in `v52/UUDvsUDD.lean`

**Implication**: Both pair types feel the SAME attractive binding force.
The difference must come from the repulsive channel.

### Result 2: Hardening Energy is 19-69% Lower for UUD+UDD

The hardening term β|∇×φ|²|θ|² acts as a repulsive barrier. At the
midpoint between two composites:

| D | |curl(φ)|²(UU) | |curl(φ)|²(UD) | Ratio UD/UU | Reduction |
|---|---------------|---------------|-------------|-----------|
| 50 | 0.437 | 0.352 | 0.80 | -20% |
| 25 | 0.442 | 0.356 | 0.81 | -19% |
| 20 | 0.182 | 0.058 | 0.32 | -68% |
| 15 | 0.268 | 0.166 | 0.62 | -38% |
| 12 | 0.184 | 0.057 | 0.31 | -69% |
| 10 | 0.285 | 0.185 | 0.65 | -35% |
| 8  | 0.431 | 0.345 | 0.80 | -20% |
| 6  | 0.207 | 0.095 | 0.46 | -54% |
| 4  | 0.187 | 0.058 | 0.31 | -69% |
| 2  | 0.240 | 0.134 | 0.56 | -44% |

The ratio oscillates with D (wavelength ≈ 2π/k ≈ 4.2 code units) but
is ALWAYS less than 1.0. The UD pair never has more hardening than UU.

**Lean proof**: `hardening_reduction` and `UD_hardening_leq_UU` in
`v52/UUDvsUDD.lean`

### Result 3: The Mechanism — Chirality Cancellation in curl_z

The curl of the total field has three components. The x-component
(from braids 1,2) is position-independent:

  curl_x(single UUD) = -0.1543
  curl_x(single UDD) = -0.0819

  curl_x(UU pair) = -0.3085  (2 × UUD, constructive)
  curl_x(UD pair) = -0.2361  (UUD + UDD, partially cancel)

The critical difference is in **curl_z**, which receives a contribution
from ∂φ_x/∂y (braid 1, y-axis). Braid 1 has **opposite chirality** in
UUD vs UDD:

  UUD braid 1: χ₁ = +1 → ∂φ_x/∂y = -Ak sin(δ₀ + Δ₁) = -C
  UDD braid 1: χ₁ = -1 → ∂φ_x/∂y = +Ak sin(δ₀ + Δ₁) = +C

For the UU pair: sum = (-C) + (-C) = -2C (constructive)
For the UD pair: sum = (-C) + (+C) = 0   (DESTRUCTIVE CANCELLATION)

**Lean proof**: `curl_z_cancellation_in_UD` and
`curl_z_constructive_in_UU` in `v52/UUDvsUDD.lean`

This single cancellation removes the entire braid-1 contribution to
curl_z² from the hardening energy of the UD pair. The reduction is
4C² = 4A²k²sin²(δ₀+Δ₁), which accounts for most of the 19-69% gap.

## Result 4: Breathing Invariance During Collision

The proton breathing mode (P_int oscillation, period ≈ 2.2 time units)
is UNPERTURBED by the collision. V51 zoom data shows phi_max cycling
between 0.31-0.86 with identical amplitude before, during, and after
the proton-proton first contact at t≈25.

**Mechanism**: The carrier phase cancellation Σ_b cos(δ + Δ_b) = 0 for
Δ = {0, 2π/3, 4π/3} is a property of the carrier phases alone. When two
protons overlap, each proton's contribution independently cancels by
linearity: α×0 + β×0 = 0.

**Lean proof**: `overlap_preserves_cancellation` in
`v52/CoreInteraction.lean`

The breathing is an internal mode protected by the carrier phase identity.
External overlap modulates P, curl, and M at the midpoint, but does not
couple into the internal oscillation because the cancellation holds at
all amplitudes.

## Physical Prediction

Same attraction (P identical) + weaker repulsion (|curl|² lower) →
**UUD+UDD can approach more closely and bind more strongly than UUD+UUD.**

The hardening shell that prevents merger is 19-69% thinner for the
proton-neutron pair. This predicts:

1. **Deuterium (UUD+UDD) binds** — the reduced barrier allows the
   depletion attraction to capture the pair into a bound state.

2. **Diproton (UUD+UUD) may not bind** — the full hardening barrier
   may be sufficient to prevent capture, or the bound state is weaker
   and more easily disrupted. (V51 shows UU binding, but with weaker
   binding energy than expected for a real nucleus.)

3. **The binding energy ratio** should satisfy:
   |E_bind(UD)| > |E_bind(UU)|
   because the UD pair reaches smaller equilibrium separation.

4. **The equilibrium separation** D_eq(UD) < D_eq(UU) because the
   hardening barrier is shifted inward for the UD pair.

## Source Files

- `v52/core_interaction.mac` — Two-proton overlap algebra (P, curl, M)
- `v52/core_interaction_output.txt` — Raw Maxima output
- `v52/uud_udd_interaction.mac` — UU vs UD comparison
- `v52/uud_udd_interaction_output.txt` — Raw Maxima output
- `v52/CoreInteraction.lean` — Breathing invariance proofs
- `v52/UUDvsUDD.lean` — Hardening reduction proofs
- `v52/uud_vs_udd.mac` — Initial chirality analysis
