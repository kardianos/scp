# Combo 2+5: Inertia + Condensate — Lattice Deformation = Spin-2

## Thesis

The HIGHEST-rated combination. The condensate (oscillon lattice) has
collective modes (phonons). When the lattice is DEFORMED (perturbed),
the deformation propagates as elastic waves. The SHEAR component of
these elastic waves is spin-2 (symmetric traceless tensor).

This is the elastic interpretation (post-v22 Proposal 2) applied to the
oscillon LATTICE, not to the vacuum. The lattice provides the material
that has shear rigidity — something the empty vacuum lacks.

## What Makes This Different from V23-D

V23-D built a lattice WITHOUT pairwise coupling → too weakly bound → melted.

This test uses:
- Pairwise coupling λ to stabilize the lattice
- Confining potential (if Test C succeeds) to prevent leaking
- The deformation (Test B) applied to the LATTICE, not a single oscillon

## Method

1. Build a stable oscillon lattice (from Test E, with λ and/or confinement)
2. Apply a TRANSVERSE perturbation: push oscillon #4 perpendicular to
   the chain direction (in 2D)

   In 1D: apply a perturbation that breaks the φ₁=φ₂=φ₃ symmetry
   at oscillon #4 (e.g., boost φ₁ component only). This creates an
   ANTISYMMETRIC wave that propagates along the chain.

3. Track the perturbation propagation speed and spatial pattern
4. The antisymmetric perturbation in the chain = the "shear" mode
5. Measure: speed, amplitude decay, polarization

## Key Question

Does the antisymmetric lattice perturbation propagate as a COHERENT
wave (phonon), or does it scatter/damp?

If coherent: the lattice supports a shear-like mode = candidate spin-2.
If damped: the antisymmetric content is expelled (as in V24-DG at the
single-oscillon level).

## Reference Code

- Test E (lattice) results
- Test B (deformation) results

## Output

- `src/combo25.c`, `data/`, `RESULTS.md`

## Parameters

From Test E: best stable lattice parameters
Perturbation: δφ₁ = +ε at oscillon #4, δφ₂ = -ε (antisymmetric)
Nx = N_osc * 320, t=10000

Compile: `gcc -O3 -Wall -o combo25 src/combo25.c -lm`
