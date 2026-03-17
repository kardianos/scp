# Locality Option 2: Antisymmetric Mode AS Causal Gravity

## Thesis

Don't add Φ at all. The near-gapless antisymmetric mode from the confining
potential (Combo 2+3+5, ω≈0.001) already propagates causally through the
lattice. Use it AS the gravitational degree of freedom.

The chain: oscillon deformation (from acceleration or perturbation) sources
the antisymmetric mode → mode propagates at finite speed → other oscillons
feel a force from the arriving antisymmetric field.

No new fields, no new parameters, no Poisson. Causality is automatic.

## The Key Test

Does the antisymmetric mode actually mediate a FORCE between oscillons?

Combo 2+5 showed it PROPAGATES and is AMPLIFIED. But does it create an
ATTRACTIVE force? Or repulsive? Or no net force?

## Method

### Phase 1: Build confined lattice

1. Use confining potential V=-σ√(P²+ε²)+(κ_c/2)P⁴ with σ=1, κ_c=0.01
2. Add pairwise coupling λ=0.5
3. Build 8-oscillon periodic chain at d=16
4. Equilibrate t=5000

### Phase 2: Measure antisymmetric force

5. Apply an antisymmetric perturbation at oscillon #1 ONLY
   (φ₁ += ε, φ₂ -= ε, ε=0.01)
6. Track: does oscillon #5 (opposite side of ring) feel a FORCE?
7. Specifically: does the CENTER of oscillon #5 shift?
8. If yes: the antisymmetric mode creates a force at distance d·4 = 64
   (4 spacings away, much longer than direct tail overlap)
9. Measure the time delay: does the force arrive at t = 64/c_group?

### Phase 3: Force sign and distance dependence

10. Measure the displacement of oscillons #2,3,4,5,6,7,8 after the
    perturbation at #1
11. Plot displacement vs distance (number of spacings from #1)
12. Is the force attractive (displacement toward #1) or repulsive?
13. Does it decay with distance?

### Phase 4: Two isolated oscillons (no lattice)

14. Two confined oscillons at D=40 (no lattice, just two objects)
15. Give one an antisymmetric kick
16. Track: does the other respond? After what delay?
17. This isolates the antisymmetric-mode force from lattice phonon effects

## No Free Parameters

This test has NO tunable parameters — it uses the existing confining
potential, pairwise coupling, and field dynamics. The "gravitational
constant" is determined by the coupling between symmetric (oscillon
shape) and antisymmetric (mediator) modes.

## Reference Code

- v24/fundamental/combo_235/src/combo235.c (confined lattice + antisym)
- v24/fundamental/combo_25/src/combo25.c (antisymmetric perturbation)

## Output

- `src/locality_antisym.c`, `data/`, `RESULTS.md`

## Parameters

σ=1.0, κ_c=0.01, λ=0.5, m=1.0, ε_pert=0.01
N_osc=8, d=16, periodic BC
Nx=2560, t_equil=5000, t_test=10000
Two-oscillon: Nx=8000, xmax=200, D=40

Compile: `gcc -O3 -Wall -o locality_antisym src/locality_antisym.c -lm`
