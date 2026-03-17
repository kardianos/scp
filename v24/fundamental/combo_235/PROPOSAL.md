# Combo 2+3+5: Confined Lattice with Inertial Deformation

## Thesis

Combine the three strongest positive results:
- Test B (Inertia): oscillon deforms under acceleration, a=F/M verified
- Test C (Confine): confining potential at σ=1 produces deeper binding (ω=0.69)
- Test E (Lattice): pairwise coupling λ=0.5 stabilizes the chain

Build a chain of CONFINED oscillons with pairwise coupling. Perturb one
oscillon and measure whether the response propagates as a SHEAR mode.

Combo 2+5 found the antisymmetric mode propagates and is amplified in the
standard lattice. Does the confining potential change this? The confining
potential has ω=0.69 (lower than standard ω=0.87), so the antisymmetric
optical branch gap may be different.

## Setup

The confining+quartic potential from Test C (σ=1, κ_c=0.01):

    V = -σ√(P² + ε²) + (κ_c/2)P⁴

Combined with pairwise coupling λ=0.5.

Total force on φ_a:
    acc = ∇²φ_a - m²φ_a + σ·P·(∂P/∂φ_a)/√(P²+ε²) - 2κ_c·P³·(∂P/∂φ_a)
          - λ(φ_b + φ_c)

## Method

### Phase 1: Equilibrate single confined+pairwise oscillon

1. Initialize Gaussian (A=0.8, σ_init=3.0) with confining potential + pairwise
2. Evolve t=10000 for equilibration
3. Save profile. Record ω, E, peak amplitude
4. Compare with standard potential (same λ=0.5, no confinement)

### Phase 2: Build chain

5. Place 8 equilibrated profiles in periodic chain at d=16
6. Add small random displacements δ∈[-0.5,0.5]
7. Evolve t=10000
8. Track positions, compute phonon spectrum
9. Compare stability with Test E (standard potential, same λ)

### Phase 3: Antisymmetric perturbation (from Combo 2+5)

10. At equilibrium: apply antisymmetric kick at oscillon #4
    (φ₁ += ε, φ₂ -= ε with ε=0.01)
11. Evolve t=10000
12. Track: does the antisymmetric mode propagate?
13. Is it still amplified (as in Combo 2+5)?
14. What is the optical branch frequency? Lower than 0.51 (standard)?

### Phase 4: Inertial test on confined chain

15. Apply a gentle gradient force (from Test B method) to the CHAIN
16. Does the chain respond as a rigid body or does it deform?
17. Measure the deformation pattern — is it tensorial?

## Key Questions

1. Does the confining potential make the lattice MORE stable (deeper well)?
2. Does the lower ω (0.69 vs 0.87) change the antisymmetric branch gap?
3. Is the parametric amplification (from Combo 2+5) preserved or changed?
4. Does the confined chain resist deformation differently than the standard?

## Reference Code

- v24/fundamental/testC_confine/src/confine.c (confining potential)
- v24/fundamental/testE_lattice/src/lattice.c (lattice with pairwise)
- v24/fundamental/combo_25/src/combo25.c (antisymmetric perturbation)
- v24/fundamental/testB_inertia/src/inertia.c (momentum injection)

## Output

- `src/combo235.c`, `data/`, `RESULTS.md`

## Parameters

σ=1.0, κ_c=0.01, λ=0.5, m=1.0, ε=1e-6
N_osc=8, d=16, periodic BC
Nx = 8*320 = 2560, t_equil=10000, t_chain=10000

Compile: `gcc -O3 -Wall -o combo235 src/combo235.c -lm`
