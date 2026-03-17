# V24-MD: Elastic Dislocation as Gauge Field (Emergent EM)

## Background

In the elastic interpretation (field index = spatial direction), dislocations
in the medium are described by a gauge theory. The dislocation density tensor
acts as the field strength of a gauge field. If the oscillon has non-trivial
strain topology, it naturally carries a "charge" (Burgers vector).

This is the deepest integration: EM is not added but EMERGES from the
existing three-field system. No new degrees of freedom.

## Mathematical Setup

### Strain and Dislocation

The displacement field u_a(x) = φ_a(x) gives strain:
    ε_{ij} = ½(∂_i u_j + ∂_j u_i) = ½(∂_i φ_j + ∂_j φ_i)

The dislocation density (incompatibility) tensor:
    α_{ij} = ε_{ikl} ∂_k ε_{lj}

In a perfect crystal: α = 0. A dislocation has α ≠ 0 at its core.

The Burgers vector (total dislocation content):
    b_j = ∮ ε_{ij} dl_i = ∫∫ α_{ij} dS_i

This is the analog of electric charge: b is quantized (lattice vector) and
conserved (dislocations can't end in the bulk).

### The Elastic "Photon"

The gauge field is the "stress potential" or "incompatible strain":
    ε_{ij}^{inc} = part of the strain that comes from dislocations

The "Maxwell equations" for dislocations:
    ε_{ikl} ∂_k α_{lj} = 0    (no magnetic monopoles for dislocations)
    ∂_i α_{ij} = s_j            (source = dislocation line density)

These are EXACTLY the structure of Maxwell's equations in the linearized
theory of elasticity.

## What to Compute

### Phase 1: Strain Analysis of the Oscillon

This is a 1D diagnostic (no new dynamics needed).

1. Load the equilibrated oscillon profile from v21
2. Compute the strain ε_{ij} = ½(∂_i φ_j + ∂_j φ_i)
   In 1D with three fields: ε_{1x} = ∂_x φ₁, ε_{2x} = ∂_x φ₂, ε_{3x} = ∂_x φ₃
3. Compute the dislocation density:
   α = ε_{ikl} ∂_k ε_{lj}  — but in 1D this is trivially zero (only one
   spatial direction, ε_{ikl} needs 3 indices)

**In 1D: no dislocations.** The elastic EM only exists in 2D+ where the
strain field can have topological defects.

### Phase 2: 2D Strain Topology

4. Create a 2D oscillon (adapt from v24/vortex/src/vortex2d.c)
5. Compute the 2D strain field ε_{ij} (2×2 matrix at each point)
6. Compute the incompatibility: κ = ε_{ij,kl}·ε_{kl,ij} - ... (2D Riemann)
7. Does the breathing oscillon have nonzero strain incompatibility?

### Phase 3: Dislocation Content

8. If incompatibility ≠ 0: compute the Burgers vector
9. Is it quantized? Is it conserved in time?
10. Compute the "EM field" from the incompatible strain

**Expected**: The breathing oscillon likely has ZERO dislocation content
(it's a smooth, spherically symmetric deformation). Only TOPOLOGICAL
defects (dislocations, disclinations) carry Burgers vectors. The oscillon
is not topological.

However: the VORTEX configuration (from V24-A) WOULD have dislocation
content. If a vortex could be stabilized (which V24-A showed it cannot
with m²φ² mass term), it would carry an elastic "charge."

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (profile loading)
- v24 vortex: `/home/d/code/scp/v24/vortex/src/vortex2d.c` (2D grid)

## Output

- `src/maxwell_d.c` — strain analysis code (mostly diagnostic, not dynamics)
- `data/strain_profile.tsv` — strain components vs x
- `data/incompatibility.tsv` — dislocation density (2D)
- `RESULTS.md`

## Parameters

1D analysis: use v21 equilibrated profile
2D analysis: Nx=Ny=256, L=20, μ=-20, κ=20, m=1.0

Compile: `gcc -O3 -Wall -o maxwell_d src/maxwell_d.c -lm`
