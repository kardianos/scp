# V35 Planning: Electrons, Quantization, and Multi-Scale Simulation

**Status**: Planning. Builds on V34's confirmed results (gravity, EM-like
θ force, charge-dependent winding interaction).

---

## What V34 Established

1. **Gravity**: φ-depletion, power-law 1/r^1.2, F = -C×∇ρ (R²=0.9998)
2. **EM carrier**: massless θ field, curl-coupled, dt-converged
3. **Charge**: winding number W=±1, θ_φ reverses with W
4. **Force**: same-winding attracts 27% more, opposite 57% less
5. **Mechanism**: wave-mediated (QFT-like), not static Biot-Savart

## What V35 Needs to Address

### Problem 1: The Electron

The θ shell peaks at r≈8 (2× the braid radius). The real electron
orbits at 50,000× the proton radius. The missing ingredient is ℏ —
a minimum action that creates a confinement cost preventing collapse.

### Problem 2: Quantization (ℏ)

The theory is purely classical. ℏ must emerge from the substrate's
discreteness: the field quantum ε (minimum field step per cell) or
the grid spacing dx. Without ℏ: no atoms, no chemistry, no electrons.

### Problem 3: Scale Bridging

Can't simulate braid (r≈5) and electron (r≈250,000) on a uniform grid.
Need multi-scale techniques: AMR patches, effective potentials, or
1D radial reduction.

---

## Architecture: Multi-Scale AMR Simulation

### Patch-Based Adaptive Mesh Refinement

Multiple independent grid patches at different resolutions and locations,
coupled at boundaries:

```
Patch 0 (braid):    N=128, L=5,     center=(0,0,0)      — fine, resolves helix
Patch 1 (electron): N=128, L=5,     center=(0,0,50000)   — fine, resolves orbital
Patch 2 (global):   N=256, L=100000                       — coarse, connects patches
```

Each patch evolves the same equation (Eq. 10). Coupling:
- Fine → coarse: RESTRICT (average fine cells onto overlapping coarse cells)
- Coarse → fine: PROLONGATE (interpolate coarse boundary into fine ghost zones)
- Each timestep: evolve fine patches, restrict, evolve coarse, prolongate

### SFA Extension for Multi-Resolution

New chunk types for the archive format:

```
GRID — Grid definition (name, Nx/Ny/Nz, Lx/Ly/Lz, center cx/cy/cz)
FRMD — Frame data (references a GRID by index)
```

Each FRMD has a grid_index field linking it to a GRID definition.
Multiple FRMDs per timestep (one per patch). The viewer shows all
patches simultaneously, color-coded by resolution level.

Frame sequence in the SFA file:
```
[GRID 0: "braid",    N=128, L=5,     c=(0,0,0)]
[GRID 1: "electron", N=128, L=5,     c=(0,0,r_orbit)]
[GRID 2: "global",   N=256, L=100000, c=(0,0,0)]
[JMPT / JMPF]
[FRMD: grid=0, t=0]  [FRMD: grid=1, t=0]  [FRMD: grid=2, t=0]
[FRMD: grid=0, t=dt] [FRMD: grid=1, t=dt] [FRMD: grid=2, t=dt]
...
```

### Bidirectional Coupling

The braid's φ-depletion and θ-radiation propagate outward (fine→coarse).
The global field's boundary conditions flow inward (coarse→fine).
The electron, if it exists as a bound perturbation at large r, influences
the braid through the coarse grid.

Standard AMR time-stepping (Berger-Oliger):
1. Advance coarse grid by Δt_coarse
2. Advance fine grids by multiple Δt_fine (Δt_fine = Δt_coarse / refinement_ratio)
3. Restrict fine → coarse
4. Reflux: correct conservation at fine/coarse boundaries

---

## Phase 1: Extract Effective Potential (from V34 data)

Use existing V34 radial profiles to construct V_eff(r):
- V_grav(r) from the δρ(r) profile (power-law depletion)
- V_em(r) from the θ force measurement (27% enhancement)
- Combined: V_eff(r) = V_grav(r) + V_em(r)

This is pure analysis on existing data. No new simulation needed.

## Phase 2: 1D Radial Eigenvalue Problem

Solve the radial Schrödinger equation with V_eff(r) and ℏ_sim:

    -ℏ²/(2m) d²ψ/dr² + [V_eff(r) + ℏ²l(l+1)/(2mr²)] ψ = E ψ

On a 1D logarithmic grid with N_r = 100,000 points.
Test different ℏ values to find the Bohr radius as function of ℏ.
If bound states exist: we have electron energy levels.

## Phase 3: Field Quantization (ε)

Add the field quantum to the 3D Cosserat simulation:
    After each Verlet step: φ_a = round(φ_a / ε) × ε
    (and similarly for θ_a)

Test at ε = 0.001, 0.01, 0.1.
Does the quantization create a minimum orbital radius?
Does it stabilize previously unstable configurations?

## Phase 4: θ Self-Interaction (Option A)

Add a nonlinear self-coupling to the θ equation:
    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a + η×curl(φ)_a - ∂V_θ/∂θ_a

    V_θ = (μ_θ/2)(θ₀θ₁θ₂)² / (1 + κ_θ(θ₀θ₁θ₂)²)

Can θ form its own solitons (θ-braids)? If yes: these could be electrons.

## Phase 5: CUDA GPU Port

Port the Cosserat force computation to CUDA for 50-100× speedup:
- Each grid point = one GPU thread
- Memory: 6 fields × N³ × 8 bytes per array
- Target: N=256 on 12 GB GPU (RTX 4070), N=512 on 24 GB (RTX 4090)

The force computation is embarrassingly parallel — no data dependencies
between grid points (periodic BC handles neighbors). Perfect for GPU.

## Phase 6: Full AMR Multi-Scale

Combine everything: AMR patches + GPU acceleration + field quantization.
Simulate a braid + electron orbital at realistic scale separation.

---

## GPU Hardware Requirements

| Task | Grid | Memory | Min GPU | Cost |
|------|------|--------|---------|------|
| Current sims | N=128 | 0.3 GB | Any (4+ GB) | $150 |
| High-res braid | N=256 | 2.4 GB | RTX 3060 (12 GB) | $300 |
| Multi-patch AMR | 3×N=128 | 0.9 GB | Any (4+ GB) | $150 |
| Large field | N=512 | 19 GB | RTX 4090 (24 GB) | $1600 |
| Extreme | N=1024 | 154 GB | A100 (80 GB) | $10000+ |

The sweet spot is RTX 4070 (12 GB, ~$500) or RTX 4090 (24 GB, ~$1600).
Multi-patch AMR makes the GPU choice less critical — multiple small
patches fit on any GPU.

---

## Priority Order

1. Extract V_eff(r) from V34 data (analysis only, hours)
2. Solve 1D radial eigenvalue problem (new code, hours)
3. Frequency analysis of θ field (Option B, from existing SFA data)
4. θ self-interaction experiment (Option A, code mod + sim, days)
5. Field quantization ε (Phase 3, code mod + sim, days)
6. CUDA port (Phase 5, significant engineering, week)
7. Full AMR (Phase 6, major engineering, weeks)
