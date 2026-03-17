# V24-MC: Topological Current Coupling (Baryon Current as EM Source)

## Background

Three real scalar fields in 3D have a natural topological current:

    j^μ = ε^{μνρσ} ε_{abc} ∂_ν φ_a ∂_ρ φ_b ∂_σ φ_c

This is the Skyrme baryon current — automatically conserved (∂_μj^μ = 0).
Coupling this current to a Maxwell field gives EM from topology.

## The Problem: This Requires 3+1D

The topological current j^μ involves three spatial derivatives of three fields
(ε^{μνρσ} has 4 indices, ε_{abc} has 3). In 1D, all spatial derivatives are
∂_x, so ε_{ijk}∂_x∂_x∂_x = 0 (antisymmetric tensor with repeated indices).

**This test MUST be done in 2D or 3D.**

## 2D Reduction

In 2+1D, the topological current reduces to:

    j^0 = ε_{abc} (∂_x φ_a)(∂_y φ_b) φ_c    (charge density)
    j^x = -ε_{abc} (∂_t φ_a)(∂_y φ_b) φ_c
    j^y = +ε_{abc} (∂_t φ_a)(∂_x φ_b) φ_c

Wait — in 2+1D with only x,y spatial, the 4D Levi-Civita contracts to:
    j^μ = ε^{μνρ} ε_{abc} ∂_ν φ_a ∂_ρ φ_b φ_c  ??? No, this isn't right.

Actually, the baryon current in d spatial dimensions needs d field derivatives.
In 3D: j^0 = ε_{ijk}ε_{abc}∂_iφ_a∂_jφ_b∂_kφ_c (3 derivatives, 3 fields).
In 2D: we need j^0 = ε_{ij}ε_{ab}∂_iφ_a∂_jφ_b (2 derivatives, 2 fields).
This uses only TWO fields, not three.

For THREE fields in 2D: we can use a different topological quantity:
    ρ_top = ε_{ij} φ_c ∂_i φ_a ∂_j φ_b    (with ε_{abc} understood)

This is a scalar density (not a current), and it IS nonzero for a Skyrmion-
like configuration in 2D. The baby Skyrmion.

## What to Compute

### Phase 1: 2D Oscillon + Topological Charge

1. Create a 2D solver (adapt from v21 triad3d.c using a 2D slice)
2. Initialize a 2D symmetric oscillon (φ₁=φ₂=φ₃ = A·exp(-r²/2σ²))
3. Compute ρ_top = ε_{ij}ε_{abc}φ_c∂_iφ_a∂_jφ_b
4. Does the oscillon carry nonzero topological charge Q_top = ∫ρ_top dxdy?

### Phase 2: Couple to Maxwell

5. If Q_top ≠ 0: add A_μ with □A_μ = e·j_μ
6. The EM field should develop a log(r) potential (2D Coulomb)
7. Measure the EM energy and field structure

### Phase 3: Two-Oscillon EM Interaction

8. Two oscillons with topological charges → Coulomb interaction

**Expected**: The standard breathing oscillon has ρ_top = 0 (no topology).
A Skyrmion-like initial condition would have ρ_top ≠ 0 but is a DIFFERENT
object from the oscillon. This test checks whether the oscillon acquires
any topological charge through its dynamics.

## Reference Code

- v21 3D solver: `/home/d/code/scp/v21/src/triad3d.c` (adapt to 2D)
- v24 vortex: `/home/d/code/scp/v24/vortex/src/vortex2d.c` (2D grid code)

## Output

- `src/maxwell_c.c` — 2D solver with topological current
- `data/` — field snapshots, topological charge vs time
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
Grid: Nx=Ny=256, L=20 (dx=0.156)
tfinal=2000

Compile: `gcc -O3 -Wall -o maxwell_c src/maxwell_c.c -lm`
