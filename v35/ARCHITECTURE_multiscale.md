# Multi-Scale Atom Simulation Architecture

## Overview

Two nested 3D regions coupled at a spherical boundary, with optional
GPU acceleration. Simulates the full atom: braid (nucleus) at the
center, electron shell (θ standing wave) at the Bohr radius.

---

## Grid Architecture

### Core Region (r < R_match ≈ 20)

**Purpose**: Resolve the braid's helical structure, V(P) binding,
curl coupling, θ source generation.

**Grid**: Cartesian, clipped to a sphere of radius R_match.
- N = 128-256 (256 with GPU)
- dx ≈ 0.3-0.5 (resolves the helix, ~10 points across core)
- Spherical mask: points with r > R_match are ghost zone
- BC at sphere: interpolated from wedge inner boundary
- Total points: ~1.1M (N=128) to ~8.7M (N=256)
- Memory: 300 MB (N=128) to 2.4 GB (N=256)

**Fields**: 6 components (φ_x, φ_y, φ_z, θ_x, θ_y, θ_z) stored in
Cartesian. Velocities and accelerations: another 12 arrays. Total: 18
arrays per grid point.

**Equation**: Full Cosserat (Eq. 10):
    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P) + η×curl(θ)_a
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a       + η×curl(φ)_a

### Wedge Region (R_match < r < R_max ≈ 500,000)

**Purpose**: Resolve the electron orbital (θ standing wave) with
angular structure for p, d, f orbitals.

**Grid**: Spherical coordinates (r, θ_polar, φ_azimuth), pie-slice.
- N_r = 1000 (logarithmic spacing: dr_min ≈ 0.5 at r=20, dr_max ≈ 5000 at r=500,000)
- N_θ = 40 (angular resolution for p/d orbitals)
- N_φ = 40 (azimuthal resolution for d/f orbitals)
- Angular extent: 30° × 30° wedge (1/12 of sphere, mirror BC on faces)
- Or: full hemisphere (180° × 360°) for complete orbital structure
- Total points: 1.6M (wedge) or 40M (full hemisphere at N_θ=200)
- Memory: 230 MB (wedge) to 5.8 GB (full hemisphere)

**Fields**: Same 6 components, stored in Cartesian (x,y,z) at each
grid point. The spherical (r,θ,φ) layout is for INDEXING only — the
physics stays Cartesian. Each point knows its (x,y,z) position;
Laplacian/curl computed from Cartesian distances between neighbors.

**Equation**: Same Cosserat, but with ℏ_eff quantum confinement term:
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a + η×curl(φ)_a
                 - (ℏ²/m_eff) × quantum_pressure_a

The quantum pressure term prevents classical collapse:
    quantum_pressure_a = -∇²(√ρ_θ)/√ρ_θ × θ_a  (Bohm potential)

Or equivalently: evolve the Schrödinger equation for θ instead of
the classical wave equation at large r. The Schrödinger equation
naturally includes ℏ and produces standing waves (orbitals) instead
of propagating waves.

### Boundary Conditions

**Core outer boundary** (sphere at r = R_match):
- θ values interpolated from wedge inner boundary
- φ values: the braid is entirely inside, so φ ≈ background at R_match
- Updated each coupling step

**Wedge inner boundary** (r = R_match):
- θ source: curl(φ) evaluated at R_match from core simulation
- Decompose into spherical harmonics Y_lm at R_match
- Each (l,m) mode propagates independently in the wedge
- This cleanly separates the s, p, d, f components

**Wedge outer boundary** (r = R_max):
- Free/absorbing: outgoing waves exit without reflection
- Or: zero (ψ → 0 at infinity for bound states)

**Wedge angular boundaries** (mirror):
- For wedge mode: f(θ_edge) = f(reflected θ)
- Vector components: flip perpendicular, keep parallel
- Effectively simulates the full sphere with 1/12 the points

---

## Coupling Protocol

Each global timestep:

```
1. Advance CORE by Δt_core (many fine steps)
   - Full 3D Cosserat, N=128-256
   - Fine dt from CFL condition

2. EXTRACT at R_match:
   - Sample θ(R_match, θ_polar, φ_azimuth) from core
   - Decompose into spherical harmonics Y_lm
   - Compute curl(φ) source at R_match

3. Advance WEDGE by Δt_wedge (coarser steps OK at large r)
   - Schrödinger-like equation for θ with ℏ_eff
   - Source term from core's curl(φ) at R_match
   - Logarithmic r-spacing handles the scale naturally

4. INJECT back to core:
   - Wedge provides θ at R_match → core outer ghost zone
   - Only needed if back-reaction is significant

5. RECORD to SFA:
   - FGRP with two grids (core + wedge)
   - Core as Cartesian FRMD
   - Wedge as spherical FRMD (new grid type in SFA)
```

---

## The Quantum Term

The classical wave equation (∂²θ/∂t² = ∇²θ + ...) produces waves
that propagate and disperse. The Schrödinger equation produces
STANDING WAVES (orbitals) that are stable.

For the wedge, we switch from the wave equation to Schrödinger:

**Classical (core, r < R_match)**:
    ∂²θ_a/∂t² = ∇²θ_a + η×curl(φ)_a

**Quantum (wedge, r > R_match)**:
    iℏ ∂ψ_a/∂t = -ℏ²/(2m_eff) ∇²ψ_a + V_eff(r) ψ_a

Where ψ_a is the quantum wavefunction for the θ sector, and
V_eff(r) = -1.27/(r/5)^1.189 is the potential from the braid.

This transition is physically motivated: at small r (core), the
field dynamics are classical and nonlinear. At large r (wedge),
the θ field is weak (linear regime) and the quantum nature of ℏ
becomes relevant. The transition happens at R_match where both
descriptions agree (the classical θ wave matches the quantum ψ).

---

## GPU Architecture (V100 32GB)

### Memory Layout

| Component | Size | Location |
|-----------|------|----------|
| Core fields (N=256, 18 arrays) | 2.4 GB | GPU |
| Wedge fields (1.6M × 18) | 230 MB | GPU |
| Spherical harmonic buffers | 50 MB | GPU |
| SFA write buffer | 100 MB | CPU (pinned) |
| Total GPU | ~3 GB | Fits easily in 32 GB |

### CUDA Kernels

1. **core_forces**: One thread per core grid point. Computes Laplacian,
   V'(P), curl coupling. ~100 lines of CUDA. Embarrassingly parallel.

2. **wedge_forces**: One thread per wedge grid point. Computes Laplacian
   in Cartesian (using actual distances between spherical-grid neighbors),
   Schrödinger evolution. ~80 lines.

3. **boundary_extract**: Extract θ at R_match, compute Y_lm decomposition.
   Small kernel, runs on ~N_θ × N_φ = 1600 points.

4. **boundary_inject**: Interpolate wedge values to core ghost zone.
   Small kernel.

5. **verlet_step**: Half-kick, drift, recompute, half-kick. Standard
   template for both core and wedge.

### Performance Estimate (V100)

| Operation | Time per step | Bottleneck |
|-----------|--------------|------------|
| Core forces (N=256) | ~2 ms | Memory bandwidth |
| Wedge forces (1.6M pts) | ~0.5 ms | Compute |
| Boundary coupling | ~0.1 ms | Latency |
| Verlet integration | ~1 ms | Memory bandwidth |
| Total per step | ~4 ms | |
| Steps for T=100 | ~2000 | |
| Wall time for T=100 | **~8 seconds** | |

A V100 would make this interactive. A full orbital formation
simulation (T=1000) would take ~80 seconds.

### With 2× V100

- GPU 0: core simulation
- GPU 1: wedge simulation
- NVLink: boundary exchange (~200 MB/s, negligible latency)
- ~2× throughput for independent runs (parameter sweeps)
- Marginal benefit for a single coupled run (boundary exchange is fast)

One V100 is sufficient. Two helps for parallel parameter exploration.

---

## SFA Format Extension

New grid type for the wedge:

```
GRID chunk for spherical wedge:
  grid_type: SPHERICAL_WEDGE (new enum value)
  Nr, Ntheta, Nphi (instead of Nx, Ny, Nz)
  r_min, r_max (instead of Lx, Ly, Lz)
  theta_min, theta_max, phi_min, phi_max
  r_spacing: LOGARITHMIC (new flag)
  center: (cx, cy, cz) — center of the sphere
```

The FRMD chunk for this grid stores data in (r, θ, φ) order:
flat_index = ir × (Nθ × Nφ) + itheta × Nφ + iphi

The viewer composites both grids: inner Cartesian volume + outer
spherical shell rendering.

---

## Implementation Phases

### Phase 1: Spherical core boundary (CPU, no GPU)
- Modify v33_cosserat to clip the grid to a sphere
- Implement spherical ghost zone interpolation
- Test braid stability with spherical BC
- Estimate: 2-3 days

### Phase 2: Wedge grid with Schrödinger evolution (CPU)
- Implement logarithmic spherical grid
- Cartesian physics on spherical mesh (finite differences)
- Schrödinger time stepper for θ in the wedge
- Test: eigenstate formation in V_eff(r)
- Estimate: 3-5 days

### Phase 3: Core-wedge coupling (CPU)
- Boundary extraction: θ at R_match
- Spherical harmonic decomposition
- Injection back to core ghost zone
- Test: coupled evolution, orbital formation
- Estimate: 2-3 days

### Phase 4: CUDA port (GPU)
- Port core_forces kernel
- Port wedge_forces kernel
- Async boundary exchange
- Test: reproduce CPU results, measure speedup
- Estimate: 3-5 days

### Phase 5: Visualization
- SFA spherical grid support
- Volume viewer: nested scales (core + wedge)
- Electron shell rendering at Bohr radius
- Estimate: 2-3 days

### Phase 6: Physics campaigns
- Formation of s, p, d orbitals
- Orbital transitions (emission of θ waves = photons)
- Two-atom interaction
- Estimate: ongoing

**Total to first orbital visualization: ~2-3 weeks**
**With V100: physics campaigns run in real time**

---

## Hardware Recommendation

| Option | Cost | Capability |
|--------|------|-----------|
| Current (8-core CPU) | Free | Phases 1-3 (slow), Phase 6 limited |
| V100 32GB × 1 | ~$1-2/hr | All phases, interactive |
| V100 32GB × 2 | ~$2-4/hr | Parallel campaigns, marginal single-run benefit |
| L4 24GB × 1 | ~$0.50-1/hr | All phases except N=512 core |

**Recommendation: 1× V100 32GB** when ready for Phase 4+.
Phases 1-3 are engineering (CPU is fine for development).
Phase 4 onward benefits enormously from GPU.
