# Multi-Resolution Nested Grid Simulation

## Overview

This document describes the design for running multiple simulation domains at
different resolutions within a single `scp_sim` process, with interpolated
boundary coupling between domains. The grid hierarchy is fixed at initialization
(no dynamic adaptive refinement).

The goal: simulate a large coarse domain containing the full environment, with
one or more fine domains zoomed in on regions of interest (braid cores, interaction
zones). The fine domains see the coarse domain's long-range effects through
interpolated boundary conditions. The coarse domain sees the fine domain's
resolved structure through restriction (volume averaging).

---

## Grid Hierarchy

### Terminology

- **Domain**: A single rectangular grid with its own (N, L, center) parameters.
- **Root domain**: The coarsest domain, covering the largest spatial extent. Parent of all others.
- **Child domain**: A finer domain spatially contained within its parent.
- **Ghost zone**: Extra cells around a domain's boundary, filled by interpolation from the parent.
- **Prolongation**: Interpolating coarse data to fill fine ghost zones (coarse → fine).
- **Restriction**: Volume-averaging fine data back onto the overlapping coarse cells (fine → coarse).
- **Subcycling**: The fine domain takes multiple timesteps per coarse timestep (because dt ∝ dx).

### Nesting Rules

1. Every non-root domain has exactly one parent (the next coarser domain that contains it).
2. A domain's spatial extent must be fully contained within its parent's extent (no partial overlap).
3. Multiple children can share the same parent (siblings at the same level).
4. Siblings at the same level do NOT directly exchange boundary data — they communicate through their common parent.
5. Nesting depth is unlimited in principle but practically 2-3 levels.

### Example: Two Braids at Different Resolutions

```
Root (N=64, L=50):  covers [-50,50]^3, dx=1.587
  ├─ Child A (N=128, L=10, center=(15,0,0)):  covers [5,25]×[-10,10]^2, dx=0.157
  └─ Child B (N=128, L=10, center=(-15,0,0)): covers [-25,-5]×[-10,10]^2, dx=0.157
```

The root domain resolves long-range depletion and wave propagation between the braids.
Each child domain resolves its braid's internal structure at 10× higher resolution.

---

## Boundary Coupling

### Prolongation (Coarse → Fine Ghost Zones)

Each fine domain has `ghost` layers of cells outside its physical boundary (typically
ghost=2 for the 7-point Laplacian stencil). These cells are filled by trilinear
interpolation from the parent domain's field values.

For a fine ghost cell at world position (x, y, z):

1. Map to parent grid indices: `i_p = (x + L_parent) / dx_parent`
2. Find the 8 surrounding parent cells (floor/ceil in each dimension)
3. Trilinear interpolate all 12 field arrays (phi, theta, velocities)

This is done for all 6 faces of the fine domain (or fewer if a face coincides
with the root domain's boundary, where absorbing BC applies instead).

**Temporal interpolation**: During subcycling, the parent's field values change.
To avoid stale boundary data, we store the parent's state at both the start
and end of the coarse timestep and linearly interpolate in time:

```
ghost_data(t_fine) = (1 - α) × parent(t_coarse) + α × parent(t_coarse + dt_coarse)
where α = (t_fine - t_coarse) / dt_coarse
```

This requires double-buffering the parent's boundary region (the subset of
parent cells that overlap the fine domain's ghost zone).

### Restriction (Fine → Coarse)

After the fine domain completes its subcycling steps, its interior data is
volume-averaged back onto the parent grid. For each parent cell that overlaps
the fine domain's interior (excluding ghost zones):

1. Identify all fine cells contained within the parent cell
2. Average their values (arithmetic mean for field values and velocities)
3. Overwrite the parent cell with the averaged value

This ensures the parent "sees" the fine domain's resolved structure. Without
restriction, the parent would evolve independently and drift out of sync.

**Restriction is applied to all 12 arrays** (6 fields + 6 velocities). The
accelerations are NOT restricted — they are recomputed from the restricted
fields on the next parent force evaluation.

### Boundary Coupling Order (Berger-Oliger)

For each coarse timestep:

```
1. STORE parent state at t_n (for temporal interpolation)
2. Parent: half-kick (velocity += 0.5 * dt_coarse * acceleration)
3. Parent: drift (position += dt_coarse * velocity)
4. PROLONGATE: interpolate parent t_n and t_{n+1} boundary → fine ghost zones
5. Fine: subcycle ratio steps:
   for k = 1 to substeps:
     a. Update ghost zones (temporal interpolation at t_n + k*dt_fine)
     b. Fine: half-kick
     c. Fine: drift
     d. Fine: compute forces (uses ghost zone data for boundary stencils)
     e. Fine: half-kick
     f. Fine: absorbing BC (applied at fine domain edges NOT adjacent to parent)
6. RESTRICT: volume-average fine interior → overwrite parent overlap cells
7. Parent: compute forces (now sees restricted fine data)
8. Parent: half-kick
9. Parent: absorbing BC
```

The subcycle count is: `substeps = round(dx_parent / dx_fine)` (integer).

### What Crosses the Boundary

- **Long-range depletion**: The braid in child A depletes ρ in the root domain.
  This depletion propagates across the root to child B's ghost zone, where
  child B sees it as a modified background. This is how gravity works across
  domain boundaries.

- **Theta waves**: The braid sources theta waves (EM analog) that propagate
  outward. These cross into the root domain at the fine/coarse boundary,
  propagate through the root, and enter child B's ghost zone. Resolution
  drops at the boundary (fine → coarse → fine), so high-frequency theta
  content is lost in transit.

- **Reflected waves**: Waves hitting the fine/coarse boundary from inside
  partially reflect due to the resolution change. This is the main artifact
  of the method. Mitigation: place fine domain boundaries in low-activity
  regions (at least 5 dx_fine from the braid surface).

---

## SFA Frame Interleaving

### Frame Layout

Each simulation step produces frames from multiple domains. These are written
to a single SFA file, interleaved, with each frame tagged by its domain.

For a 3-domain hierarchy (root + 2 children) with the root snapping every
5.0 time units and children every 1.0:

```
Frame 0: root,    t=0.0
Frame 1: child_A, t=0.0
Frame 2: child_B, t=0.0
Frame 3: child_A, t=1.0
Frame 4: child_B, t=1.0
Frame 5: child_A, t=2.0
Frame 6: child_B, t=2.0
...
Frame 11: root,    t=5.0
Frame 12: child_A, t=5.0
Frame 13: child_B, t=5.0
...
```

Children may snapshot more frequently than the root (they resolve faster dynamics).

### GDEF Chunk (Grid Definitions)

A new SFA chunk type `GDEF` is written between CDEF and KVMD. It defines all
grids in the file:

```
Chunk header: "GDEF" + size (standard 12-byte SFA chunk header)

Payload:
  version:  u8  (1)
  n_grids:  u16
  reserved: u8[5]
  --- 8-byte sub-header ---

  Per grid (repeated n_grids times):
    grid_id:   u16
    parent_id: u16  (0xFFFF = root, no parent)
    Nx, Ny, Nz: u32 × 3
    Lx, Ly, Lz: f64 × 3  (half-domain size)
    cx, cy, cz: f64 × 3  (center in world coordinates)
    ghost:      u16       (ghost zone width)
    kvmd_set:   u16       (which KVMD set has this grid's parameters)
    reserved:   u8[4]
    --- 64 bytes per grid ---
```

Old SFA readers skip the GDEF chunk (they jump directly to JTOP via the
stored offset, same as KVMD). For single-grid files, no GDEF is written.

### Frame Grid Tagging

The existing SFA frame index entry has a `reserved` u32 field:

```c
typedef struct {
    double time;
    uint64_t offset;
    uint64_t compressed_size;
    uint32_t checksum;
    uint32_t reserved;  // ← repurposed as grid_id
} SFA_L2Entry;
```

For multi-grid files, `reserved` stores the `grid_id` of the domain that
produced this frame. For single-grid files, it remains 0 (backward compatible).

### Variable Frame Size

Different grids have different N, so their frames have different byte sizes.
The SFA reader determines frame_bytes per frame by looking up the grid_id
in the GDEF and computing `N_total × column_bytes`. The header's `frame_bytes`
field is set to the root grid's frame size (largest).

The reader flow:
1. Read frame index entry → get grid_id from reserved field
2. Look up grid_id in GDEF → get (Nx, Ny, Nz)
3. Compute frame_bytes = Nx×Ny×Nz × sum(column_dtype_sizes)
4. Decompress frame data using this size

### KVMD Per Domain

Each domain can have its own KVMD set. The GDEF entry's `kvmd_set` field
points to the KVMD set containing that domain's physics parameters.

For most simulations, all domains share the same physics (same m, eta, mu, kappa)
and only differ in grid parameters (N, L, center). In this case, one KVMD set
suffices (all GDEF entries point to set 0).

For multi-physics simulations (different coupling modes in different regions),
each domain points to a different KVMD set.

---

## Multi-Resolution Seeding

### Seed File Structure

A multi-resolution seed is a single SFA file containing one frame per domain,
with GDEF and KVMD metadata describing the hierarchy.

```
seed.sfa:
  SFAH: version, flags
  CDEF: 12 columns (phi/theta/vel)
  GDEF: 3 grids (root + 2 children)
  KVMD set 0: shared physics (m, eta, mu, kappa, ...)
  KVMD set 1: root grid params (N=64, L=50)
  KVMD set 2: child A params (N=128, L=10, cx=15)
  KVMD set 3: child B params (N=128, L=10, cx=-15)
  JTOP/JMPF: index
  Frame 0: root grid data (64^3 × 12 columns)
  Frame 1: child A data (128^3 × 12 columns)
  Frame 2: child B data (128^3 × 12 columns)
```

### Seed Generator Design

A multi-resolution seed generator takes a hierarchy description and produces
the seed SFA:

```bash
gen_multi_seed \
  --physics "m=1.5 eta=0.5 mu=-41.345 kappa=50" \
  --grid "id=0 N=64 L=50" \
  --grid "id=1 N=128 L=10 cx=15 parent=0 init=braid A=0.8" \
  --grid "id=2 N=128 L=10 cx=-15 parent=0 init=braid A=0.8" \
  -o two_braids.sfa
```

Each `--grid` argument specifies:
- Grid geometry: id, N, L, center, parent
- Initialization mode: braid, oscillon, or empty (background only)
- Init-specific parameters: A, R_tube, ellip, delta, sigma

The generator:
1. Creates each grid's field arrays independently (no coupling at init)
2. For child grids: also initializes the background field (A_bg oscillation) so the
   child is consistent with the root at its spatial location
3. Writes all frames + GDEF + KVMD to a single SFA

### Seeding the Overlap Region

When a child domain overlaps the root, both contain field data for the same
spatial region. At initialization, they should be CONSISTENT — the child's
data should match the root's data (at the root's resolution) plus additional
fine-scale detail.

Procedure:
1. Generate root domain first (full extent, including braid if applicable)
2. For each child domain:
   a. Sample the root domain's field values at the child's grid points (trilinear interpolation)
   b. Add the child's own braid/structure on top of the interpolated root data
   c. This ensures the child starts consistent with the root

If the child has a braid and the root doesn't (braid exists only at the fine scale),
the root is initialized with just the background, and the child is initialized with
background + braid. After the first restriction step, the root will see a coarse
representation of the child's braid.

---

## Config Format for Multi-Resolution

```ini
# Shared physics
m          = 1.5
m_theta    = 0.0
eta        = 0.5
mu         = -41.345
kappa      = 50.0

# Shared simulation params
T          = 200.0
dt_factor  = 0.025
damp_width = 3.0
damp_rate  = 0.01

# Grid hierarchy
n_grids    = 3

grid.0.N      = 64
grid.0.L      = 50.0
grid.0.center = 0,0,0
grid.0.snap_dt = 5.0

grid.1.N      = 128
grid.1.L      = 10.0
grid.1.center = 15,0,0
grid.1.parent = 0
grid.1.ghost  = 2
grid.1.snap_dt = 1.0

grid.2.N      = 128
grid.2.L      = 10.0
grid.2.center = -15,0,0
grid.2.parent = 0
grid.2.ghost  = 2
grid.2.snap_dt = 1.0

# Initialization
init       = sfa
init_sfa   = two_braids_seed.sfa

# Output
output     = multi_res.sfa
precision  = f32
diag_dt    = 2.0
diag_file  = multi_res_diag.tsv
```

---

## Viewer / Analysis Integration

### Per-Frame Mode

The viewer reads the SFA's GDEF to identify grids. For each frame, it checks
the grid_id and renders at that grid's native resolution and position. The user
can filter by grid_id to see only one domain.

### Virtual Frame Mode

The viewer composites multiple domains into a single virtual volume at a
user-specified resolution and viewport:

1. User specifies: target resolution (e.g., N_virtual=256), viewport center and extent
2. For each virtual voxel at world position (x, y, z):
   a. Find the finest grid that contains this point
   b. Sample that grid's field value (trilinear interpolation if needed)
3. Render the composited volume

Priority: finer grids override coarser grids in overlap regions. This gives
a seamless view where fine detail appears where available and coarse data
fills the rest.

### Analysis

The diagnostics TSV already contains per-timestep energy data. For multi-grid,
each diagnostic line is tagged with the grid_id. Analysis tools can:
- Track energy per domain separately
- Compute total energy across all domains (summing non-overlapping contributions)
- Monitor inter-domain energy transfer (energy leaving one domain's boundary
  should appear in the parent)

---

## Known Limitations

### Interpolation Artifacts

- **Reflections at fine/coarse boundary**: Fast waves hitting the resolution
  boundary partially reflect. Amplitude depends on wavelength relative to
  dx_coarse. Waves with λ < 4×dx_coarse reflect significantly.
  Mitigation: place boundaries far from sources (buffer of 5-10 fine cells).

- **Temporal interpolation lag**: Ghost zone data is linearly interpolated in
  time between coarse steps. If the coarse field changes nonlinearly within
  one coarse step (e.g., during a rapid event), the interpolation is inaccurate.
  Mitigation: keep the subcycle ratio moderate (≤10).

### Conservation

- **Energy is not exactly conserved at boundaries**: The interpolation/restriction
  cycle introduces O(h_coarse²) energy error per crossing. For typical
  simulations this is <1% and comparable to the absorbing BC energy drain.

- **The restriction step is lossy**: Fine structure below the coarse resolution
  is lost when restricted. This is by design — the coarse domain cannot
  represent fine detail. But it means the coarse domain's representation
  of the fine region is always a smoothed approximation.

### Static Grid Hierarchy

- The grid positions and sizes are fixed at initialization. If a braid
  drifts out of its fine domain, it will be represented only at the coarse
  resolution. The simulation doesn't detect this or adapt.

- For long-running simulations where braids move significantly, either:
  (a) Make the fine domains large enough to contain the expected trajectory
  (b) Periodically checkpoint and restart with repositioned grids

### No Sibling-to-Sibling Coupling

- Two child domains at the same level don't exchange boundary data directly.
  They communicate through the parent. If two fine domains are adjacent with
  no gap, waves crossing between them go: fine A → restrict to parent → prolongate
  to fine B. This adds one level of interpolation compared to direct exchange.

- For closely spaced braids, use a single fine domain covering both rather
  than two adjacent fine domains.

---

## Implementation Sequence

1. GDEF chunk in sfa.h (write + read)
2. Frame grid_id tagging (use reserved field in L2Entry)
3. Variable frame_bytes per grid_id in reader
4. Multi-grid config parser
5. Domain struct with parent/child pointers and overlap indices
6. Ghost zone allocation (N + 2×ghost per axis)
7. Prolongation kernel (trilinear, with temporal interpolation buffer)
8. Restriction kernel (volume average)
9. Berger-Oliger subcycling orchestrator
10. Multi-grid SFA writer (interleaved frames with grid_id tags)
11. Multi-resolution seed generator (`sfa/seed/gen_multi_seed.c`)
12. Viewer: per-frame grid-aware rendering
13. Viewer: virtual frame compositor
