# Combo ALL 7: The Complete Theory

## Thesis

All seven ideas combined into one system:

1. EOS: the confined condensate has a definite P(ρ)
2. Inertia: deformation of the condensate is tensorial
3. Confinement: linear potential prevents leaking
4. Gauge: local breathing frequency creates a gauge structure
5. Condensate: oscillon lattice is the vacuum
6. Time crystal: phase-locked breathing defines time
7. Self-reference: fields evolve on the metric they create

## Method

This is the FINAL test, to be run only if the individual components work.

1. Use confining potential (from Test C)
2. Build stable condensate lattice (from Test E) with pairwise coupling
3. Add the self-consistent metric (from Test F: Φ from Poisson)
4. Add the gauge field Ω for local breathing frequency (from Test D)
5. Phase-lock the oscillons (automatic from triple product)
6. Compute the EOS (from Test A)
7. Apply a perturbation and measure the FULL response:
   - Compression phonon (scalar, spin-0): speed c_s
   - Shear phonon (tensor, spin-2): speed c_shear
   - Gauge wave (vector, spin-1): speed c_gauge
   - Gravitational response (from self-consistent Φ change)

## What This Tests

Does the COMBINED system naturally produce three different force carriers:
- Spin-0 (scalar Proca) from the antisymmetric pairwise mode
- Spin-1 (vector gauge) from the local breathing gauge
- Spin-2 (tensor shear) from the condensate deformation

And does the self-consistency condition fix the coupling constants?

## Key Output

The MASS SPECTRUM of the combined system:
- How many distinct modes?
- What are their spins?
- What are their masses (ranges)?
- Are the coupling constants determined or free?

This would be the complete low-energy effective theory of the three-field
oscillon condensate.

## Reference Code

- ALL prior tests (A through G, combos 25, 1257, 235)

## Output

- `src/all7.c`, `data/`, `RESULTS.md`

## Parameters

From all prior tests: best parameters for each component
Full complexity: 3 scalar fields + Φ + Ω = 5 field arrays
Lattice: N_osc=8 periodic, t=10000

Compile: `gcc -O3 -Wall -o all7 src/all7.c -lm`
