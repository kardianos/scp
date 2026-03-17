# Test B: Inertia — Oscillon Deformation Under Acceleration

## Thesis

Apply a constant external force to the oscillon and measure its deformation.
If the deformation is tensorial (quadrupolar), the oscillon's inertia has
spin-2 character — connecting inertia to gravity via the equivalence principle.

## Method

Add an external potential V_ext = -F·x to the equation of motion. This
applies a constant force F in the x-direction, accelerating the oscillon.

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - ∂V_triple/∂φ_a + F

The oscillon accelerates at a = F/M (where M = E/c²). As it accelerates,
it Lorentz-contracts on the leading edge and stretches on the trailing edge.

## What to Measure

1. Equilibrate the oscillon for t=5000 (no force)
2. Turn on force F at t=5000
3. Track the oscillon center x_c(t). Verify a = F/M.
4. At each timestep, measure the oscillon's SHAPE:
   - Width σ_+(t) on the leading side (x > x_c)
   - Width σ_-(t) on the trailing side (x < x_c)
   - Asymmetry: Δσ = σ_+ - σ_-
   - Energy density profile: ρ(x-x_c) in the co-moving frame
5. The deformation: δρ(x) = ρ_accelerated(x) - ρ_rest(x)
6. Multipole decomposition of δρ: monopole (size change), dipole (shift),
   QUADRUPOLE (spin-2 deformation)
7. Scan F = {0.001, 0.01, 0.1} to check linear response

In 1D: only monopole and "even/odd" decomposition is meaningful (no angular
structure). The deformation IS measurable as leading/trailing asymmetry.

In 2D (if time permits): the deformation has angular structure, and the
quadrupole component can be directly measured.

## Key Prediction

The deformation δρ should be:
- Linear in F (at small F): δρ ∝ F
- Quadrupolar in shape: compressed leading + stretched trailing
- The deformation tensor = the strain = the metric perturbation

## Reference Code

- v21/src/triad1d.c

## Output

- `src/inertia.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, Nx=4000, xmax=100
t_equil=5000, t_accel=3000
F scan: {0.001, 0.01, 0.1}

Compile: `gcc -O3 -Wall -o inertia src/inertia.c -lm`
