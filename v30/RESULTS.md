# V30 Results: Rotating FRW Expansion

## Gate Test: Does Rotation Create Particles? — YES

| n_wind | n_particles | max_rho | E_total | max_rho/avg_rho |
|--------|------------:|--------:|--------:|----------------:|
| 0 (control) | 0 | 0.059 | 1.5e5 | ~1× (uniform) |
| 1 | 409 | 0.325 | 1.3e5 | ~6× |
| 2 | 432 | 0.322 | 1.4e5 | ~6× |
| 3 | 548 | 0.334 | 1.4e5 | ~6× |
| 5 | 500 | 0.336 | 1.4e5 | ~6× |

Rotation spontaneously creates 400-550 localized energy peaks during
cosmological expansion. n_wind=0 creates ZERO — confirming that rotation
(angular momentum) is the essential ingredient.

### CORRECTION: v1 particles were boundary artifacts (atan2 + periodic BC)

All 400-550 "particles" in v1 were at r_perp=96-137 (domain boundary).
The atan2 phase winding is incompatible with periodic BC → artificial
singularities at the boundary.

### v2: Smooth velocity rotation → NO effect

Rotation applied as Lie derivative (δv = Ω×(-y∂φ/∂x + x∂φ/∂y)). But
the initial field has no xy structure (z-only wave), so the Lie derivative
is exactly ZERO. Rotation has no effect. nw=0 and nw=3 are identical.

### v3: Gaussian blob + atan2 vortex lines → shell fragmentation

| n_wind | Particles (late) | Location | Mechanism |
|--------|-----------------|----------|-----------|
| 0 (control) | 175 | r≈45-55 (blob edge) | Shell fragmentation |
| 3 (rotation) | 90 | r≈51-60 (blob edge) | Same + vortex disruption |

Particles form at the EXPANDING FRONT of the Gaussian blob, not in the
bulk. This is modulational instability (Rayleigh-Taylor), not topological
defect formation. The control (nw=0) produces MORE particles than rotation
(nw=3) because the vortex lines disrupt the shell.

**Pin 7 remains OPEN.** The expansion + rotation approach does not produce
braids. The particles are shell fragments, not vortex-line braids.
They form at the expanding front, have no winding, and dissipate slowly.

## V30 c-field: Density-Dependent Speed of Light

### The Physics
c_eff²(x) = 1/(1 + (ρ/ρ_ref)^0.5)
- Dense field → SLOWER c (Schwarzschild analog)
- Depleted field → c approaches c₀=1
- Waves curve TOWARD depleted regions → gravitational lensing

### The Approach
- Large domain: N=512, L=5000 (dx≈20)
- Concentrated Gaussian blob (A₀=1.0, R=300) in empty space
- Pinned boundaries (φ=0 at edges)
- Quartic potential V=(λ/4)(φ²-v²)² replaces triple product (stable at coarse dx)
- Massless (m=0) — mass emerges from c(ρ), not as Lagrangian parameter
- Small random asymmetry (1%) to seed instabilities

### Status: RUNNING (stable, slow)
- t=0: E=2.6e4, R_rms=379, fc=0.87 (concentrated blob)
- t=50: E=1.5e8, R_rms=1468, fc=0.04 (expanded, energy released from quartic)
  - max_rho increased 100× (peak concentration forming)
  - 57K peaks detected (grid-scale noise, not resolved structures)
- Running to T=10000, diagnostics every T=50
- ~40 min per diagnostic at N=512 with smoothing pass
- Will take several hours to complete

### Issues Encountered
1. **v1 (FRW + atan2 rotation)**: Boundary artifacts from atan2 + periodic BC
2. **v2 (smooth velocity rotation)**: Zero effect (z-only field has no xy gradient)
3. **v3 (Gaussian blob + atan2)**: Shell fragmentation, not braids
4. **c-field v1 (c∝ρ^α)**: Blowup — positive feedback (dense→fast→denser)
5. **c-field v2 (c∝1/(1+ρ^α))**: Blowup — dt too large for mass oscillation
6. **c-field v3 (quartic, massless, dt capped)**: STABLE, running

### Key Design Decisions
- Triple product replaced by quartic: the triple product requires dx<1 resolution
  and is unstable at dx=20. The quartic V=(λ/4)(φ²-v²)² is stable at any dx.
- Mass set to zero: m=1.5 creates evanescent modes at dx>>1/m that grow
  exponentially. Mass should emerge from c(ρ), not be a parameter.
- c SLOWER in dense regions: this creates attraction (geodesic lensing toward
  mass). The opposite (c faster in dense) creates runaway positive feedback.
- Smoothing pass on ρ: removes grid-scale noise before computing c_eff.
  Accounts for ~50% of compute cost. Could optimize with separable filter.

### Result: MONOTONIC DISPERSAL — No Structure

| t | R_rms | fc | n_peaks | max_rho |
|---|-------|------|---------|---------|
| 0 | 379 | 0.865 | 8483 | 1.6e-4 |
| 100 | 1840 | 0.017 | 158K | 1.9e-2 |
| 250 | 2841 | 0.004 | 381K | 1.9e-2 |
| 350 | 3365 | 0.002 | **1** | 1.8e-2 |
| 400 | 3603 | 0.002 | **0** | 1.8e-2 |

The blob expands at constant speed until it fills the domain. By T=350,
the field is completely uniform (zero peaks). No structure formation.

**Root cause**: c SLOWER in dense regions → the core FREEZES while the
edges fly off at c=1. This is NOT gravity — it's dispersion. The dense
core doesn't attract anything; it just stops evolving. The edges peel
away freely.

**The fundamental insight**: For gravity, the SOURCE should propagate
normally while OTHER objects (test particles) curve toward it. Making
c(ρ) affect everything equally (including the source) creates freezing,
not attraction. This is the equivalence principle problem — the source
shouldn't feel its own gravity the same way test particles do.

## V30 Lessons Learned

### What DOESN'T work for c(ρ) gravity:

1. **c faster in dense regions** (c∝ρ^α): positive feedback → blowup
2. **c slower in dense regions** (c∝1/(1+ρ^α)): core freezes, edges disperse
3. **Uniform c(ρ) on all fields**: source and test particles both affected
   → freezing or blowup, never attraction

### What MIGHT work (untested):

1. **M7 two-component + c(ρ_B)**: The BACKGROUND field B determines c_eff.
   The braid (S) propagates at c_eff set by B. The braid depletes B locally
   → c_eff drops near braid → other braids (also propagating on B) slow
   down and curve toward the depletion. The braid itself is NOT affected
   by its own depletion because S and B are separate fields.

2. **Self-interaction excluded**: c_eff for field a depends on ρ of fields
   b,c (not itself). Each field propagates on a metric set by the others.

3. **Retarded c(ρ)**: c_eff depends on the TIME-AVERAGED ρ (not instantaneous).
   This prevents the oscillating braid from freezing at its peak density.
   The time-averaged ρ is smooth and slowly varying → stable c_eff.

## Pin Status

V30 is PINNED. The c-field concept (c depends on field density) is
physically sound but the implementation needs the M7 two-component
framework to separate source from test particle. This connects back
to V29's T12 M7 depletion results (which already showed 1/r^1.2
power-law depletion with M7).

The natural next step: M7 + c(ρ_B) at large scale. The braid depletes
B, B sets c_eff, other braids propagate on c_eff → attraction.


## Setup
- N=256, L=100, dx=0.78
- FRW expansion: H₀=0.05, t_inf=50 → a_final=12.2
- Phase 2 (formation): T=300 with H=0 (static spacetime)
- Periodic BC in all directions
- Initial: uniform dense rotating field A₀=3.0 filling entire box
- Triple product: μ=-41.3, κ=50, m=1.5
