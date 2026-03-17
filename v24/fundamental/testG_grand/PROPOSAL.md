# Test G: Grand Test — Confined Condensate Deformation

## Thesis

Combine ideas 2 (inertia), 3 (confinement), and 5 (condensate) into one
test. Build a chain of oscillons with CONFINING potential. Apply a
perturbation and measure whether the lattice response is spin-2 (tensor).

This is the most ambitious test — it requires confinement to work (Test C),
the lattice to be stable (Test E), and the deformation to be measurable
(Test B). Run this AFTER the individual tests, using their results.

## Method

1. Use confining potential V = -σ√(P²+ε²) with σ from Test C
2. Build 8-oscillon periodic chain at appropriate spacing
3. Add pairwise coupling λ from Test E
4. Equilibrate the chain
5. Apply a LOCALIZED perturbation: push oscillon #4 by +δx
6. Track how the perturbation propagates through the chain
7. Decompose into compression (longitudinal) and shear (transverse) modes
8. In 1D: only compression exists. But the ASYMMETRY of the deformation
   (how oscillon shapes change) reveals the tensor structure.

For full tensor analysis: need 2D lattice.

## 2D Version (if 1D positive)

9. Build a 4×4 oscillon lattice in 2D
10. Perturb one oscillon
11. Measure the angular pattern of the perturbation wave
12. Decompose into monopole (compression), dipole, quadrupole (shear)
13. The quadrupole fraction = spin-2 content

## Reference Code

- Tests B, C, E results (required before running)
- v23/phonon/src/chain1d.c
- v24/vortex/src/vortex2d.c (2D grid)

## Output

- `src/grand.c`, `data/`, `RESULTS.md`

## Parameters

Use best parameters from Tests B, C, E.
1D: N_osc=8, periodic, t=10000
2D: N_osc=16 (4×4), Nx=Ny=512, t=5000

Compile: `gcc -O3 -Wall -o grand src/grand.c -lm`
