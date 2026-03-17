# V24-S6: The Spin Question — Is the Proca Mediator Scalar or Tensor?

## Thesis

The antisymmetric Proca mediator is a scalar (spin-0) in field space. For
gravity, we need spin-2 (tensor). For electromagnetism, we need spin-1
(vector). What IS the spin of the Proca mediator in the elastic
interpretation where field index = spatial direction?

## Background

### The Elastic Interpretation (Post-V22 Proposal 2)

If φ_a(x) is the displacement of the medium in direction a, then:
- The strain tensor ε_{ij} = ½(∂_iφ_j + ∂_jφ_i) is a rank-2 tensor
- The antisymmetric part ω_{ij} = ½(∂_iφ_j - ∂_jφ_i) is the rotation
- The symmetric traceless part is spin-2

The antisymmetric MODE of the three fields (A₁ = (φ₁-φ₂)/√2) is a
scalar in field space but its STRAIN has tensor character.

### The V23-A Hessian Result

V23-A showed the triple product HARDENS the antisymmetric sector
(λ_anti > 0 when μ < 0). This killed the shear-softening proposal.

But the pairwise coupling changes this: it SOFTENS the antisymmetric
sector (m²_A = m² - λ < m²). The pairwise coupling does what the triple
product couldn't — it makes the antisymmetric (shear) modes lighter.

### Does the Pairwise Coupling Change the Hessian Analysis?

With BOTH triple product and pairwise coupling, the effective antisymmetric
mass is:
    m²_A,eff = m² - λ + |μ|f⁴/(1+κf⁶)²

The triple product contribution is POSITIVE (+0.74 at f=0.5 from V23-A).
The pairwise contribution is NEGATIVE (-λ).
At λ = 0.99: m²_A,eff = 1.0 - 0.99 + 0.74 = 0.75 (still positive).

So the antisymmetric mode has effective mass 0.87 inside the oscillon
and 0.1 in vacuum. The mode is LIGHTER in vacuum than in the oscillon
— the opposite of what's needed for trapping.

## What to Compute

### Test 1: Spin Classification in 1D

In 1D, the antisymmetric mode is a scalar (one component). There's no
spin-2 in 1D. This test is purely a consistency check.

1. Compute the Proca mediator's radiation pattern in 1D
2. It should be isotropic (scalar) — forward and backward equally

### Test 2: Spin Classification in 2D

In 2D (x,y) with fields φ₁ (x-displacement) and φ₂ (y-displacement):
- The symmetric strain ε₊ = ∂_xφ₁ + ∂_yφ₂ (compression, scalar)
- The antisymmetric strain ε₋ = ∂_xφ₁ - ∂_yφ₂ (shear, spin-2)
- The rotation ω = ∂_xφ₂ - ∂_yφ₁ (vorticity, spin-0)

The Proca mediator in the elastic interpretation:
A₁ = (φ₁-φ₂)/√2. Its strain: ∂_x A₁ = (∂_xφ₁-∂_xφ₂)/√2.

This is NOT the same as the shear strain ε₋. The field-space
antisymmetric mode is a DIFFERENT object from the spatial-strain
antisymmetric mode.

**To get spin-2**: need the Proca mediator to couple to the TRACELESS
SYMMETRIC strain, not to the field-space antisymmetric combination.

3. Set up a 2D solver with pairwise coupling
4. Initialize a 2D oscillon with a small antisymmetric perturbation
5. Measure the radiation pattern of the perturbation
6. Is it isotropic (scalar), dipolar (vector), or quadrupolar (tensor)?

### Test 3: Cross-Gradient + Pairwise Combined

The cross-gradient term η(∂_iφ_j)(∂_jφ_i) from V24-B couples field
indices to spatial indices. Combined with the pairwise coupling, it might
route the field-space antisymmetric mode into the spatial-strain tensor
channel.

7. Add BOTH pairwise coupling λ AND cross-gradient η to the 2D solver
8. Measure the radiation polarization
9. Does the combination produce spin-2 radiation?

## The Deeper Question

Is there ANY coupling that makes the Proca mediator spin-2?

The answer from representation theory: a scalar field cannot be made
spin-2 by coupling alone. The spin is determined by the field's
transformation under spatial rotations. A scalar transforms trivially
(spin-0). To get spin-2, the mediating field must BE a rank-2 tensor.

The three scalar fields CAN form rank-2 tensors through PRODUCTS of
DERIVATIVES: (∂_iφ_a)(∂_jφ_b). The traceless symmetric part of this
is spin-2. But this is a COMPOSITE operator, not a fundamental field.

The question: can the Proca mediator's equation of motion be rewritten
as an equation for the composite spin-2 tensor? If so, the effective
mediator is spin-2 even though the fundamental fields are scalar.

## Reference Code

- v23/hessian: Hessian analysis (verify with pairwise)
- v24/vortex: 2D solver
- v24/crossgrad: cross-gradient coupling

## Output

- `src/proca_spin.c` — 2D solver with pairwise + spin analysis
- `data/` — radiation patterns, polarization measurements
- `RESULTS.md`

## Parameters

2D grid: Nx=Ny=256, L=20, μ=-20, κ=20, m=1.0
λ=0.99, η=0.5 (cross-gradient, for Test 3)
t_run=2000

Compile: `gcc -O3 -Wall -o proca_spin src/proca_spin.c -lm`
