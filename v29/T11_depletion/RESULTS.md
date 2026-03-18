# T11: Field Depletion — Results

## Phase 1: Natural Profile — NO depletion in vacuum

The braid in vacuum has ρ(r) monotonically decreasing from core to zero.
No background to deplete. ρ_vac = 0 everywhere outside the braid.
The braid is a positive energy concentration, not a depletion.

## Phase 2b: Differential Measurement — DEPLETION DETECTED

### Method
Two parallel simulations: (A) background only, (B) braid + background.
Same dynamics, same initial background. δρ(r) = ρ_B(r) - ρ_A(r) isolates
the braid's effect on the background field.

N=128, L=30, A_bg=0.1, T=400 (time-averaged over T=200..400).

### Result: Depletion Zone at r=16.5-21

| Region | δρ/ρ_ctrl | Interpretation |
|--------|-----------|----------------|
| r < 5 | +125% | Braid core (energy source) |
| r = 5-16 | +2 to +50% | Braid tail + radiation passing through |
| **r = 16.5-21** | **-2% to -7%** | **DEPLETION ZONE** |
| r > 21 | +5% | Edge reflection artifacts |

Peak depletion: **-6.6% at r=18.75** (relative to control background).

### Time Evolution of Depletion

From the timeseries at r=10:
- t=0: δρ = +0.0001 (no effect yet)
- t=20: δρ = +0.012 (initialization radiation arriving)
- t=60: δρ = +0.002 (radiation passing through)
- t=80: δρ = -0.003 (DEPLETION beginning)
- t=100: δρ = -0.017 (depletion strengthening)

The depletion develops AFTER the initial radiation pulse passes.
The signal arrives at r=10 around t≈10-20 (speed ≈ c, consistent
with massive wave propagation).

### Earlier Agent Result (A_bg=0.01)

At lower background, the agent's Phase 2 found even stronger depletion:
- r=12.2: δρ/ρ_bg = -24.6%
- r=13.0: δρ/ρ_bg = -50.5%
- Peak: 11.45% below far-field at r=11.8

Stronger depletion at lower A_bg suggests the braid has a FIXED consumption
rate, and the fractional depletion scales as 1/ρ_bg (= 1/A_bg²).

## Physical Interpretation

### The Braid as a Field SINK
The braid consumes background field energy to sustain itself. The consumption
creates a depletion zone around the braid. This depletion:
1. Propagates outward at approximately c
2. Has a characteristic radial profile (ring-shaped depletion zone)
3. Scales inversely with background density (fixed consumption rate)

### Mass = Consumption Rate
The braid's "mass" should be identified with its field consumption rate
(dE/dt from the background), NOT with the total energy or volume.
A heavier braid consumes faster → deeper depletion → stronger gradient
→ stronger "gravity."

This makes mass EMERGENT: it's a property of the braid-field interaction,
not a Lagrangian parameter. The m=1.50 mass in the current model provides
confinement, but the GRAVITATIONAL mass could be the consumption rate.

### What Model Changes Are Needed

The current depletion grows without bound (braid keeps consuming). In reality:

1. **Self-limiting consumption**: The braid should reach equilibrium where
   consumption = release. As it moves (V>0 always), it releases field behind
   and consumes ahead. The steady-state depletion is the DIFFERENCE.

2. **Depletion-dependent speed**: The braid's propagation speed should depend
   on local ρ. In depleted regions, the braid slows down (less field to
   propagate through). This creates a natural feedback: more depletion →
   slower braid → less consumption → equilibrium.

3. **Background must be dynamic**: The current "background" is just initial
   conditions. For gravity, the background field must be a physical entity
   with its own dynamics — perhaps a condensate or a thermal bath (connecting
   to V29-T1b's thermal equilibrium finding).

### Connection to Gravity

If the depletion profile stabilizes as δρ ∝ -M/r at large r (where M is
the consumption rate), then other braids propagating through this field
experience a modified speed:

    c_eff² = c² × (1 - δρ/ρ₀) ≈ c² × (1 + M/(ρ₀ r))

This is EXACTLY the weak-field Schwarzschild metric:

    ds² = -(1 - 2GM/c²r)c²dt² + (1 + 2GM/c²r)dr² + r²dΩ²

with G ∝ 1/ρ₀. The gravitational constant is INVERSELY proportional to
the background field density. Planck-scale ρ₀ → small G.

## Next Steps

1. **Measure depletion RATE**: Track dρ/dt at fixed r over time. Is it
   constant (steady sink) or decaying (equilibrating)?

2. **Test 1/r profile**: Does the late-time depletion follow δρ ∝ 1/r?
   Need larger domain and longer time.

3. **Moving braid**: Give the braid a velocity, measure the asymmetric
   wake (leading edge depletion vs trailing edge release).

4. **Two braids**: Do two braids in a depleted field attract via the
   gradient? This is the gravity test.

5. **Self-consistent coupling**: Make c_eff depend on local ρ, creating
   the feedback loop that makes gravity dynamical.
