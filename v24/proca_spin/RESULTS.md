# V24-S6 Results: The Spin Question

## Summary

The Proca mediator (antisymmetric mode A = (phi_1 - phi_2)/sqrt(2)) radiates
as a **spin-0 scalar** with pairwise coupling alone. Adding the cross-gradient
term induces a weak l=2 quadrupolar component (~11%) but does NOT produce
dominant spin-2 radiation. The fundamental field is scalar; no coupling
repackages it into a tensor mediator.

## Test 1: 1D Baseline (Forward/Backward Symmetry)

Parameters: N=2000, xmax=60, tfinal=1000, lambda=0.99, mu=-20, kappa=20, m=1.

| Quantity | Value |
|----------|-------|
| <\|A\|^2> at x=+40 | 8.514e-06 |
| <\|A\|^2> at x=-40 | 8.483e-06 |
| Forward/backward ratio | **1.004** |
| Classification | **Isotropic (spin-0)** |

The antisymmetric mode radiates equally in both directions, confirming
scalar (spin-0) character in 1D.

## Test 2: 2D Radiation Pattern (Pairwise Only)

Parameters: Nx=Ny=256, L=20, lambda=0.99, eta=0, tfinal=1000, R_meas=15.

Angular multipole decomposition of <|A|^2> on circle at r=15:

| l | |C_l| | |C_l|/C_0 | Interpretation |
|---|-------|-----------|----------------|
| 0 | 3.47e-07 | 1.0000 | monopole (spin-0) |
| 1 | 3.02e-20 | 0.0000 | dipole (spin-1) |
| 2 | 4.69e-22 | 0.0000 | quadrupole (spin-2) |
| 3 | 6.74e-23 | 0.0000 | octupole |
| 4 | 8.24e-10 | 0.0024 | l=4 (lattice artifact) |

Max/min anisotropy ratio: **1.05** (essentially isotropic).

**Classification: ISOTROPIC (spin-0 scalar)**

The pairwise coupling alone does NOT produce any angular structure. The
antisymmetric Proca mediator is a pure scalar field regardless of the
"elastic interpretation" where field index = spatial direction. The
field-space antisymmetric combination (phi_1 - phi_2) is a DIFFERENT object
from the spatial-strain tensor (d_x phi_1 - d_y phi_2).

## Test 3: 2D with Cross-Gradient (Pairwise + Cross-Gradient)

Parameters: Nx=Ny=256, L=20, lambda=0.99, eta=0.5, tfinal=1000, R_meas=15.

| l | |C_l| | |C_l|/C_0 | Interpretation |
|---|-------|-----------|----------------|
| 0 | 3.60e-07 | 1.0000 | monopole (spin-0) |
| 1 | 7.84e-21 | 0.0000 | dipole (spin-1) |
| 2 | 3.92e-08 | **0.109** | quadrupole (spin-2) |
| 3 | 1.15e-20 | 0.0000 | octupole |
| 4 | 1.83e-08 | **0.051** | l=4 |

Max/min anisotropy ratio: **1.72** (max at 138 deg, min at 68 deg).

**Classification: MIXED — predominantly spin-0 with ~11% quadrupolar admixture**

The cross-gradient term eta*(d_a phi_b)(d_b phi_a) couples field indices to
spatial indices, introducing a genuine l=2 component. However:

1. The quadrupole is a **perturbation** (11%), not the dominant channel
2. The l=4 component (5%) is comparable — this is lattice-scale coupling, not
   clean spin-2
3. The odd multipoles (l=1, l=3) remain zero — consistent with the parity
   symmetry of the initial perturbation

## Physical Interpretation

### Why pairwise coupling gives spin-0

The pairwise coupling V_pw = (lambda/2)(phi_1 - phi_2)^2 is a potential
in field space. It depends on field VALUES, not field DERIVATIVES. Therefore
it cannot couple to spatial directions. The antisymmetric mode A satisfies a
massive Klein-Gordon equation with no directional preference:

    (d_t^2 - nabla^2 + m^2 + 2*lambda) A = nonlinear(S, A)

This is a scalar equation — no tensor structure.

### Why cross-gradient introduces (weak) anisotropy

The cross-gradient eta*(d_x phi_2)(d_y phi_1) explicitly couples spatial
direction x to field 2 and spatial direction y to field 1. In the (S,A) basis:

    E_cg = (eta/4) [(d_x S)(d_y S) - (d_x A)(d_y A) + (d_x S)(d_y A) - (d_x A)(d_y S)]

The (d_x A)(d_y A) term acts as an anisotropic kinetic energy for A. This
breaks rotational invariance of the A equation, producing the observed l=2
component. But the effect is perturbative: the dominant isotropic kinetic
energy (1/2)|grad A|^2 always wins.

### Why spin-2 is impossible from scalar fields

A fundamental result from representation theory: the spin of a mediator is
determined by its transformation under spatial rotations. A scalar field
phi(x) transforms trivially (spin-0). No coupling to other scalars can change
this — couplings preserve the representation.

To get spin-2 radiation, the mediating field must BE a rank-2 tensor h_{ij}(x)
that transforms as a symmetric traceless tensor under rotations. The three
scalar fields can form composite rank-2 objects like (d_i phi_a)(d_j phi_b),
but these are NOT independent propagating degrees of freedom — they are
determined by the scalar field equations.

## Conclusion

The Proca mediator is **spin-0**. No coupling within the three-scalar-field
framework produces spin-2 radiation. This is a representation-theory
obstruction, not a parameter-tuning problem.

For gravity (spin-2), one needs either:
1. A fundamental rank-2 tensor field (general relativity)
2. A mechanism that makes the COMPOSITE operator (d_i phi_a)(d_j phi_b)
   propagate independently (emergent graviton — requires new physics)

The cross-gradient coupling produces weak (11%) quadrupolar structure, but this
is anisotropic scattering, not spin-2 graviton radiation. The distinction: spin-2
requires the l=2 pattern to be the DOMINANT and UNIVERSAL radiation channel,
not a perturbative admixture.
