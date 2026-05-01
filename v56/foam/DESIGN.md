# V55 Stage 1: Voronoi Foam Cosserat Kernel

## Core idea

Replace the uniform N³ voxel grid with an unstructured **Voronoi cell complex**:

- Cells live at positions {x_c} sampled with Poisson-disk on the 3D box [-L, L]³
- Each cell has volume V_c (set by Voronoi geometry, not chosen)
- Each cell has a list of neighbors {n} with shared face areas A_{cn} and
  centroid distances d_{cn}
- The 6 Cosserat fields (φ₀, φ₁, φ₂, θ₀, θ₁, θ₂) are stored at cell
  centroids — one float per field per cell

**Density-from-geometry**: there is no separate density field. Local density
is defined as ρ_local = 1/V_c (cells with smaller volume = denser region).
For a uniform Poisson-disk sample, V_c ≈ V_total / N_cells (approximately
constant), giving uniform density. To introduce density gradients we
non-uniformly sample cell positions — the geometry IS the density.

## Why this might help

The voxel kernel has an axis-aligned grid bias: the 6-point Laplacian
stencil is anisotropic (6 nearest neighbors only) and the carrier wave
cos(kz) locks to a coordinate axis. The hedgehog seed in v52 produced
cubic structure instead of spherical because of this.

A Poisson-disk Voronoi mesh has isotropic statistics (each cell has ~14
neighbors distributed around the sphere). Stencils are inherently
non-axis-aligned, so a spherical structure stays spherical.

## Pipeline

All-C/C++. Voro++ vendored under `voro_src/` (BSD, from chr1shr/voro on
GitHub) — built into a static library `libvoro++.a` linked at compile
time. No external dependencies beyond a C++ compiler and zstd.

```
voro_src/libvoro++.a    # vendored 3D Voronoi library (BSD)
                        #
gen_foam_mesh.cpp       # C++: Poisson-disk + voro++ container_periodic_3d
                        #   one-shot, writes foam_mesh.bin
                        #
foam_sim.c              # C: reads mesh, runs 6-field Cosserat dynamics
                        #   per-cell field arrays, finite-volume operators
                        #   output: foam.sfa (per-cell field time series)
                        #
foam_to_voxel.c         # C: resamples per-cell fields to NxNxN voxel grid
                        #   output: voxel SFA for volview comparison
```

Build:
```
make -C voro_src libvoro++.a
g++ -O3 -o gen_foam_mesh gen_foam_mesh.cpp -Ivoro_src -Lvoro_src -lvoro++
gcc -O3 -fopenmp -o foam_sim foam_sim.c -lzstd -lm
gcc -O3 -fopenmp -o foam_to_voxel foam_to_voxel.c -lzstd -lm
```

## Mesh binary format (foam_mesh.bin)

Little-endian, all fields packed:

```
Header (32 bytes):
  uint32  magic        = 'FOAM'
  uint32  version      = 1
  double  L            box half-extent
  uint32  N_cells      number of cells (interior to [-L,L]³)
  uint32  N_faces      total face count (each face counted once)
  uint32  reserved[2]

Cell records (N_cells × 32 bytes):
  double  pos[3]       centroid (x, y, z)
  double  volume       V_c

Face records (N_faces × 40 bytes):
  uint32  cell_a       cell index
  uint32  cell_b       neighbor cell index
  double  area         shared face area A_{ab}
  double  delta[3]     b.pos - a.pos (shortest periodic image, baked in)

Cell→face index (N_cells+1 × 4 bytes):
  uint32  face_offsets[N_cells+1]  CSR-style offsets into face list

Face index list (sum of degrees × 4 bytes):
  uint32  face_indices[]            face IDs adjacent to each cell
```

The CSR offsets let `foam_sim.c` iterate `for face in cell.faces` cheaply.
`delta[3]` already includes periodic wrap so the kernel doesn't need to
know the box exists.

## Discrete operators (cell-centred finite volume)

For a scalar field f(c) at cell centroid:

**Gradient** (volume-averaged):
```
grad_f(c) = (1/V_c) Σ_{n ∈ neighbors} A_{cn} × (f(n) - f(c))/2 × n̂_{cn}
```
where n̂_{cn} = δ_{cn}/|δ_{cn}|.

**Laplacian** (divergence of gradient, two-point flux):
```
lap_f(c) = (1/V_c) Σ_{n ∈ neighbors} A_{cn} × (f(n) - f(c)) / |δ_{cn}|
```

**Curl** of vector field F:
```
curl_F(c)_a = ε_{abc} ∂_b F_c
            = (1/V_c) Σ_n A_{cn} (F_c(n) - F_c(c))/2 × ε_{abc} n̂_{cn,b}
```

These are all consistent finite-volume discretizations on the
Delaunay-Voronoi dual mesh. They're isotropic to leading order and reduce
to standard central differences on a regular grid.

## Time stepping

Velocity Verlet, identical to voxel kernel. Stable timestep:
```
dt < dt_factor × min(d_{cn})  (CFL on smallest face)
```
For uniform Poisson-disk, min(d) ≈ avg(d), so dt is similar to voxel case.

## Stage 1a: Eulerian (cells fixed)

This first prototype keeps cells stationary. Field values evolve, mesh
doesn't. This isolates "does Voronoi geometry alone change dynamics?"
from "do moving cells help?"

## Stage 1b (future): Lagrangian (cells flow)

Cells move with the field — small cells migrate to high-|φ| regions,
large cells to vacuum. This is where "density-as-geometry" becomes
fully active: the mesh adaptively concentrates resolution.

## Validation target

Run the same single-braid seed (v34 parameters) on a ~50k-cell foam mesh
in a [-20, 20]³ box (matches voxel N=128). Compare:

- Energy conservation drift (target < 1% over T=300, vs voxel -0.6%)
- E_pot oscillation range
- θ_rms saturation level
- Visual structure after voxel resampling

If foam ≈ voxel: discretization isn't the v54 villain.
If foam survives where voxel doesn't: new physics regimes accessible.
If foam dies where voxel survives: foam adds error not symmetry.

## Cell count vs voxel resolution

| Voxel grid | Voxel count | Equivalent foam cells |
|------------|-------------|----------------------|
| 64³        | 262 144     | ~100k                |
| 96³        | 884 736     | ~300k                |
| 128³       | 2 097 152   | ~700k                |

Foam needs fewer cells than voxels for similar resolution because each
cell averages over a sphere (volume ~ V_voxel × √3 for same minimum
distance). Start with 50k for fast iteration.
