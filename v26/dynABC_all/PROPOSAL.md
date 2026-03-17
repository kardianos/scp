# V26-DynABC: Full Dynamic Braid (Propagating + Rotating + Massless)

## Thesis

The complete dynamic braid: propagating along z at c (massless), rotating
around z with angular momentum, all energy from dynamics. No mass parameter.
This is the "particle as process" realized as a helical wave that both
translates and rotates through a braided topology.

## Initialization (3D, N=128, L=20)

    phase_a = θ + kz + 2πa/3
    φ_a(x,0) = A(r_⊥) · cos(phase_a)
    v_a(x,0) = (k + Ω) · A(r_⊥) · sin(phase_a)

where k=2π/L (axial), Ω=0.1 (rotation), m=0 (massless).
v_g = c along z, angular velocity Ω around z.

## Constraints

- V > 0: both axial (c) and angular (Ω) velocities nonzero
- T > 0: from combined oscillation k+Ω
- δV > 0: both components maintained

## Why This Is the Ultimate Test

If this configuration survives:
- The braid is dynamical (never static, always moving)
- The mass is emergent (no m² parameter, energy from dynamics)
- The angular momentum is spin (from rotation, quantized by topology)
- The spatial structure is aspherical (l=2 from braid geometry)
- The metric is self-consistent (strain = geometry)
- Causality is automatic (propagation at c)
- Three forces emerge: gravity (symmetric strain), EM (antisymmetric
  torsion), strong (triple product)

This would be the COMPLETE theory: matter, forces, and spacetime from
three dynamically braided displacement fields.

## Lagrangian

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - (μ/2)P²/(1+κP²)

The simplest possible Lagrangian. No mass, no extra couplings.
Everything from the braided initial condition + continuous dynamics.

## What to Measure (in order of importance)

1. SURVIVAL: does the dynamical braid persist at t=500?
   (fc > 0.3, |P| > 0.01, energy bounded)

2. NON-BREATHING: is ρ(center,t) constant/non-oscillating?
   (DFT peak at ω < 0.1, not at ω≈1.4 like the oscillon)

3. QUADRUPOLE: l=2 fraction on shell at R=8?
   (target: > 10%, ideally > 40% like V26 Phase 4 baseline)

4. PROPAGATION: does the pattern move at v_g?
   (track phase of φ₁ maximum along z vs time)

5. ROTATION: is angular momentum conserved?
   (L_z(t) = const)

6. Compare with STATIC braid (V26 Phase 4): same Lagrangian,
   different initial conditions. Does dynamics help?

## Parameters

μ=-20, κ=20, m=0.0, A₀=0.8, R_tube=3.0
k=2π/20=0.314, Ω=0.1
N=128, L=20, periodic BC in z, absorbing in x,y
t=500

## Reference

v26/src/v26.c, v26/src/v26_phase4.c

## Output

`src/dynABC.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynABC src/dynABC.c -lm`
