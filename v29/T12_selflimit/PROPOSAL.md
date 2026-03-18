# T12: Self-Limiting Sink — Emergent Mass from Field Consumption

## The Problem

T11 showed the braid depletes the background field (-6.6% at r≈19).
But the current model has no mechanism to LIMIT this consumption.
The braid will keep eating until the background is gone.

In reality, the consumption must self-limit:
- The braid moves (V>0 always), releasing field behind and consuming ahead
- Steady state: release + radiation = consumption
- The steady-state consumption RATE is the emergent mass
- Heavier braids consume FASTER (higher rate), not more volume

We need a modified Lagrangian where the sink is self-limiting and the
equilibrium consumption rate IS the gravitational mass.

## Status: NOTHING CONFIRMED YET

The depletion was observed (T11) but:
- Not confirmed as 1/r
- Not confirmed as self-limiting
- Not confirmed as producing attraction between braids
- Not confirmed as propagating at c in steady state

This proposal defines the mechanisms to test and the equations to try.

## Step 0: Mechanism Screening

Before committing to a specific model, test 8 candidate mechanisms for
self-limiting consumption. Each is a different modification to the EOM.
Run each for T=200 at N=96, measure whether the depletion stabilizes.

### Mechanism 1: Density-Dependent Consumption Rate
The sink weakens as the local field density drops. Like a pump that
weakens as the well runs dry.

    Sink rate: S(x) = -α × ρ(x)/ρ₀ × B(x)

where B(x) is the braid intensity (e.g., |P|² or energy density in core),
and ρ(x) is local background density. When ρ→0, S→0: no more consumption.

**Implementation**: After each Verlet step, multiply velocities in the
core by a factor (1 - α × dt × ρ_local/ρ₀). This extracts kinetic
energy from the background proportional to the local density.

**Self-limiting because**: consumption → ρ drops → consumption drops.

### Mechanism 2: Back-Pressure from Gradient
The depletion gradient creates a restoring force that opposes further
consumption. Deeper well = stronger back-pressure.

    F_back = +β × ∇ρ

This force pushes field BACK toward the depleted zone. At equilibrium,
the gradient-driven inflow balances the consumption.

**Implementation**: Add β × ∂ρ/∂x_i to the acceleration of each field.
Compute ρ(x) as the local energy density, take its gradient.

**Self-limiting because**: deeper depletion → steeper gradient → more
back-flow → limits further depletion.

### Mechanism 3: Speed-Dependent Consumption (Most Physical?)
The braid propagates through the field. Its consumption rate depends on
how fast it moves through the background. In depleted regions, the braid
slows down → consumes less → equilibrium.

    c_eff²(x) = c² × ρ(x)/ρ₀

EOM: ∂²φ_a/∂t² = c_eff²(x) × ∇²φ_a - m²φ_a - ∂V/∂φ_a

The braid in a depleted zone propagates slower → sweeps through less
background per unit time → lower consumption rate → equilibrium.

**Implementation**: Replace the Laplacian coefficient (currently 1.0)
with ρ(x)/ρ₀ where ρ is smoothed over a few grid cells.

**Self-limiting because**: consumption → depletion → slower propagation
→ less consumption. The feedback loop has a fixed point.

**Mass = c_eff at the braid location × geometric factor.**
Heavier braids: deeper depletion → lower c_eff → different fixed point.

### Mechanism 4: Saturating Potential
Modify V(P) to include a depletion cost:

    V_total = V(P) + λ_depl × (ρ₀ - ρ)²

This penalizes departures from the preferred density ρ₀. The braid
can deplete up to a maximum set by the balance V(P) vs λ_depl.

**Implementation**: Add a term to the force that pushes ρ toward ρ₀.
This requires computing ρ(x) and its functional derivative w.r.t. φ_a.

**Self-limiting because**: depletion is energetically expensive. The
braid depletes until the marginal cost exceeds the marginal benefit.

### Mechanism 5: Conservation (Redistribution, Not Consumption)
Total field energy is strictly conserved. The braid doesn't consume
field — it REDISTRIBUTES it. Energy concentrates at the core, creating
a deficit in the surroundings. Total integral = 0.

    ∫ δρ(x) d³x = 0    (exactly, at all times)

No modification to EOM needed — this is a DIAGNOSTIC constraint.
If the current model already conserves energy (it does in periodic BC),
then the depletion IS redistribution.

**Self-limiting because**: zero-sum by construction. The depletion
is exactly balanced by the concentration at the core.

**BUT**: Does redistribution produce 1/r? In electrostatics, a point
charge redistributes the field: excess at the charge, deficit at 1/r.
The same should happen here.

### Mechanism 6: Nonlinear Wave Speed
The field's own wave speed depends on local amplitude:

    c² = c₀² + γ × Σ_a φ_a²

In depleted regions (low |φ|): c lower → waves slow down → field
accumulates → depletion fills back in. Natural regulation.

**Implementation**: Replace Laplacian coefficient with
(1 + γ × Σφ²/Σφ²_vac). Requires specifying γ.

**Self-limiting because**: depletion → lower c → field accumulates
→ depletion heals.

### Mechanism 7: Two-Component Field
Split into "structural" (S, forms braid) and "background" (B, mediates gravity):

    L = ½(∂S)² + ½(∂B)² - V_S(S) - V_B(B) - g × S² × B

The coupling g×S²×B converts background B into structural S at the braid.
The reverse process (S→B) happens through radiation.

**Implementation**: Add 3 more field arrays (B_a). Double the grid memory.
The coupling S²B is like a Yukawa interaction.

**Self-limiting because**: as B depletes near the braid, conversion
rate drops. Braid radiation replenishes B far from the core.

### Mechanism 8: Preferred-Density Condensate (Mexican Hat for ρ)
The background is a condensate with preferred density ρ₀:

    V_bg(ρ) = λ(ρ² - ρ₀²)²    where ρ = √(Σ φ_a²)

This is a Mexican hat potential for the field MAGNITUDE. Depletion
(ρ < ρ₀) and excess (ρ > ρ₀) both cost energy. The braid sits at ρ
slightly below ρ₀ in its surroundings.

**Implementation**: Replace current V(P) with V(P) + V_bg(ρ).
This adds a radial restoring force toward ρ = ρ₀.

**Self-limiting because**: depletion below ρ₀ costs energy quadratically.
Maximum depletion set by V_bg vs V(P) balance.

**Connection to Higgs**: This IS the Higgs mechanism for the background
field. The braid lives in a Higgs condensate and depletes it locally.

## Ranking and Selection

| # | Mechanism | Simplicity | Physicality | Self-limit | 1/r? |
|---|-----------|-----------|-------------|------------|------|
| 5 | Conservation | trivial | high | by construction | maybe |
| 3 | Speed-dependent | moderate | **highest** | feedback loop | likely |
| 1 | ρ-dependent sink | simple | moderate | direct | likely |
| 8 | Condensate (Mexican hat) | moderate | high | energetic | likely |
| 4 | Saturating potential | moderate | moderate | energetic | maybe |
| 6 | Nonlinear wave speed | simple | moderate | feedback | maybe |
| 2 | Back-pressure | simple | low | gradient | maybe |
| 7 | Two-component | complex | moderate | conversion | likely |

### Recommended Order

**First**: Mechanism 5 (conservation). Check if the EXISTING model
already conserves total energy in periodic BC. If so, the depletion IS
redistribution, and we just need to verify 1/r profile. Zero code changes.

**Second**: Mechanism 3 (speed-dependent). Most physical: the braid
slows in depleted regions, creating a natural feedback. Connects motion
to mass. One line change in EOM (Laplacian coefficient → ρ/ρ₀).

**Third**: Mechanism 8 (condensate). Adds preferred density via Mexican
hat. Connects to Higgs physics. Provides both self-limiting AND a
background that can support 1/r depletion.

**Fourth**: Mechanism 1 (ρ-dependent sink). Simplest explicit sink.
Good for calibrating the consumption rate.

## Step 1-5: After Mechanism Selection

Once a self-limiting mechanism is identified:

1. **Verify 1/r**: Large domain (L=100), measure δρ(r) profile at late time
2. **Measure rate**: Track dρ/dt at fixed r, confirm it stabilizes
3. **Moving braid**: Asymmetric wake (leading consumption, trailing release)
4. **Two braids**: Attraction via depletion gradient
5. **Gravitational waves**: Accelerated braid → time-varying depletion → ripples

## The Target Equation

If Mechanism 3 works, the complete model would be:

    L = ½(ρ/ρ₀)(∂_t φ_a)² - ½(ρ/ρ₀)(∂_i φ_a)² - ½m²φ_a² - V(P)

    ρ(x) = smoothed energy density = ⟨½(∂φ)² + ½m²φ²⟩_local

    EOM: ∂²φ_a/∂t² = (ρ/ρ₀)∇²φ_a - m²φ_a - ∂V/∂φ_a
                       + (1/2ρ₀)(∂_i ρ)(∂_i φ_a)  [gradient coupling]

Mass emerges as: M_grav = (consumption rate at equilibrium) ∝ ∫ δρ d³x

If Mechanism 8 works, add to the above:

    V_total = V(P) + λ(ρ² - ρ₀²)²

where ρ = √(Σ φ_a²) and ρ₀ is the condensate density.

## Grid & Runtime
- Step 0 (screening): N=96, L=20, T=200, 8 mechanisms × ~3 min = ~24 min
- Steps 1-5: N=128-192, L=30-100, T=500-1000, ~30 min each
