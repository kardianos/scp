# T11: Field Depletion Gravity

## The Hypothesis
The braid doesn't impact the field — it REDUCES it. The braid consumes
background field energy to sustain itself. The missing field (depletion zone)
propagates at c·dt, creating a gradient. Other braids fall into this gradient.
That IS gravity.

## Why This Is Different from V6/V7
- V6: Added a separate density field ρ with □ρ = -source. Got 1/r but
  coupling was 10⁴⁰× too strong (nuclear, not gravitational).
- V7: Disformal self-trapping. Same problem.
- THIS: The depletion IS the field itself. No separate mediator.
  The braid is MADE OF the field, and what it consumes creates the well.

## The Key Insight: Mass = Field Consumption
The braid's mass is not a parameter — it's how much field the braid consumes
per unit time. The depletion rate IS the mass. This makes mass emergent.

## Mathematical Formulation

Background field density: ρ(x) = ½ Σ_a [(∂_t φ_a)² + (∂_i φ_a)² + m²φ_a²]
(total energy density of the field)

Without braids: ρ = ρ₀ (uniform vacuum energy)
With braid: ρ < ρ₀ near braid (field consumed to form the braid structure)

Depletion: δρ(x) = ρ₀ - ρ(x) > 0 near braid

The braid fields propagate with effective speed that depends on local ρ:
    c_eff²(x) = c² × ρ(x)/ρ₀

Where ρ < ρ₀ (near braid): c_eff < c → waves slow down → geodesics curve
toward braid → GRAVITY.

The depletion propagates outward at c → gravitational waves.

## Implementation

### Phase 1: Measure the natural depletion profile
Run the bimodal braid and measure ρ(r) = energy density as function of
distance from braid center. This is already available from existing runs
but needs to be extracted specifically.

Question: is ρ < ρ₀ outside the braid? Or does the braid ADD energy
everywhere (E > 0 total)? If the latter, we need a different mechanism.

### Phase 2: Add background field and measure depletion
Initialize the domain with uniform ρ₀ background PLUS the braid.
Run dynamics. Does the braid deplete the background?
How fast does the depletion spread?

Concretely:
- φ_a(x) = braid(x) + A_bg × cos(random phases)
  where A_bg sets the background energy density
- Run dynamics, measure ρ(r, t)
- Does ρ decrease near the braid over time?

### Phase 3: Two braids in depleted field
If Phase 2 shows depletion:
- Two braids at D=30, each depleting the field
- The combined depletion profile: ρ(x) = ρ₀ - δρ₁(x) - δρ₂(x)
- Each braid moves along -∇ρ (toward the other braid's well)
- Measure: do they attract? Is F ∝ 1/r²?

### Phase 4: Speed-of-light coupling
Make the braid EOM depend on local ρ:
    ∂²φ_a/∂t² = (ρ/ρ₀) × ∇²φ_a - m²φ_a - ∂V/∂φ_a

This couples the depletion back into the braid dynamics.
When ρ < ρ₀: waves slow down → braid curves toward depleted region.

### Phase 5: Gravitational wave extraction
Accelerate a braid (give it a kick), measure the time-varying δρ at far field.
Is it quadrupolar (spin-2)? Does it propagate at c?

## The Hidden Dimension Variant
If the field lives in 4+1 dimensions (3 space + 1 hidden + 1 time):
- The braid wraps around the hidden dimension
- Depletion in 3+1 is a projection → naturally 1/r in 3D
- The hidden dimension sets the scale: G ∝ 1/(ρ₀ × L_hidden)

## What Would Make This Work vs V6
V6's problem: G_eff/G_N ≈ 10⁴⁰ (too strong). This was because:
- Source was ½|ω|² (field gradient energy) with nuclear amplitude
- ρ₀ was O(1) in code units

For correct coupling: need ρ₀ ≫ (braid energy). If ρ₀ ~ 10²⁰ × E_braid,
then G_eff ~ E_braid² / ρ₀² ~ 10⁻⁴⁰ (correct ratio to Newton).

This means: the background field must be VASTLY more energetic than the
braid. The braid is a tiny ripple on an enormous sea. The depletion is
a tiny fractional change in ρ, but it's enough to curve geodesics.

## Grid & Runtime
- Phase 1: Reanalyze existing data (~0 min)
- Phase 2: N=128, L=20, T=500 with background (~10 min)
- Phase 3: N=192, L=40, T=500 with two braids (~30 min)
- Phase 4: Modified solver (~20 min)
- Phase 5: N=192, L=60, T=1000 (~45 min)
