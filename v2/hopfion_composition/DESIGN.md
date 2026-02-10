# Hopfion Composition Simulator — Design Document

## Goal

Build a 3D field simulator that can:
1. Initialize composed hopfion configurations (multiple linked/nested topological objects)
2. Evolve them dynamically on a **spherical** domain (not Cartesian box)
3. Track both Skyrmion charge (B) and Hopf invariant (H)
4. Include L₆ (sextic/BPS) term for topology-sensitive dynamics
5. Build toward simulating simple atoms from fundamental hopfions

## Why Spherical Grid

The existing `scatter.c` uses a Cartesian grid [-L,L]³ with N³ points. Problems:
- Wastes ~48% of memory on corner regions (|x| > L in sphere vs cube)
- Boundary effects are anisotropic (closer along faces, farther along diagonals)
- Absorbing BCs harder to implement cleanly on a cube

A spherical grid computes on a ball of radius R:
- Natural for isolated soliton systems (field → vacuum at r → ∞)
- Uniform boundary distance in all directions
- Absorbing/transparent BCs are a single radial condition
- Memory scales as (4/3)πR³ vs (2R)³ = 8R³ → saves ~48%

## Grid Structure: Yin-Yang Spherical

We use a **Cartesian grid masked to a sphere** (simplest approach that works):
- Underlying uniform grid with spacing h
- Only allocate and compute cells where r = |x| ≤ R
- Index mapping: sparse (list of active cells) or dense (full cube, skip r > R)
- For N=256, L=8: full cube = 16.7M cells, sphere = 8.7M cells (52%)

**Why not spherical coordinates (r,θ,φ)?**
- Coordinate singularity at r=0 and poles
- Non-uniform angular resolution (fine at center, coarse at edge)
- Tensor Laplacian in spherical coords is messy for 8-component field
- Hopfion topology has no preferred axis — spherical coords impose artificial symmetry

**The masked-Cartesian approach** keeps the simple finite-difference stencils from
`field.c` while cutting memory and compute by ~48%.

## Field Representation

Same 8-component Cl⁺(3,0,1) multivector as existing code:
```c
typedef struct {
    double s, f1, f2, f3;   /* bulk quaternion q */
    double j1, j2, j3, p;   /* degenerate weight sector */
} Multivector;
```

## Lagrangian Terms

### Existing (from field.c):
- **L₂**: ½|∂q|² (quadratic gradient)
- **L₄**: (1/4e²)|[R_μ, R_ν]|² (Skyrme, 4th-order)
- **V**: (λ/4)(|q|²-ρ₀²)² (Mexican hat)

### New — L₆ (Sextic/BPS term):
```
L₆ = -λ₆ × (B^μ B_μ)
```
where B^μ is the topological (baryon) current:
```
B^0 = (1/24π²) ε^{ijk} Tr(R_i R_j R_k)     [= baryon density]
B^i = -(1/24π²) ε^{0ijk} Tr(R_0 R_j R_k)    [= baryon current]
```

In terms of the quaternion field, the spatial baryon density is:
```
B⁰ = -(1/2π²) ε^{ijk} (q̃ ∂_i q)(q̃ ∂_j q)(q̃ ∂_k q) / |q|⁶
```

The L₆ energy density is proportional to (B⁰)². This is the most
topology-sensitive local term possible — it literally squares the
winding density.

**Why L₆ matters**: In the BPS limit (L₆ + V only), E = E_FB exactly.
Near BPS (small L₂ + L₄ perturbation), binding appears perturbatively.
This dramatically changes multi-soliton structure.

### New — Pion mass:
```
V_π = m_π² ρ₀² (1 - q·s/|q|)
```
Already implemented in scatter.c, needs porting.

## Hopf Charge Computation

The Hopf invariant H for a map φ: S³ → S² is:
```
H = (1/4π²) ∫ A ∧ F
```
where F = φ*(ω_S²) is the pullback of the area form on S², and dA = F.

**CP¹ projection**: From the bulk quaternion q = s + f·σ, extract the unit
vector n̂ = (f₁,f₂,f₃)/|f| (defined where |f| ≠ 0). This is the CP¹ map.

The area form pullback is:
```
F_ij = n̂ · (∂_i n̂ × ∂_j n̂)
```

The vector potential A satisfying dA = F is found numerically by solving
∇²A = ∇×(n̂ × ∇n̂) or via a direct formula involving the Hopf map coordinates.

**Practical approach**: Use the Whitehead integral formula:
```
H = (1/4π²) ∫∫ [J(x) × J(y)] · (x-y) / |x-y|³  d³x d³y
```
where J = (1/4π) n̂ · (∂n̂/∂x_i × ∂n̂/∂x_j) ε^{ijk} is the topological current.

This is an O(N⁶) computation on the full grid — expensive but feasible for
diagnostic snapshots (not every timestep).

**Fast approximation**: Compute F_ij on the grid, solve ∇²A_i = ε_{ijk}∂_j F_{kl}
via FFT, then H = ∫ A·(∇×A) d³x / (4π²). This is O(N³ log N).

## Hopfion Initialization

### Method 1: Hopf map (axially symmetric, H=1)
The standard Hopf map S³ → S² sends (z₁, z₂) ∈ C² with |z₁|²+|z₂|²=1 to
n̂ = (2 Re(z₁z̄₂), 2 Im(z₁z̄₂), |z₁|²-|z₂|²).

For a hopfion centered at origin with size parameter a:
```
z₁ = (a² - r² + 2iaz) / (a² + r²)
z₂ = 2a(x + iy) / (a² + r²)
```
where r² = x² + y² + z².

### Method 2: Rational map hopfions (H = nm)
Following [19] Balakrishnan: (n,m) parameterization gives Hopf index H = nm.
Preimages form torus knots for general (n,m).

### Method 3: Composed hopfions (linked)
Place two or more hopfions at different positions. Their preimage curves
will link, contributing to the total Hopf charge. This is the key new
capability for studying composition.

### Method 4: Hopfion-Skyrmion hybrid
Initialize a Skyrmion (hedgehog, B=1) with a hopfion ring threaded
around it, as observed experimentally in FeGe [01].

## Boundary Conditions

### Absorbing (sponge) layer
At radius R_inner < r < R:
```
∂q/∂t → (1 - γ(r)) ∂q/∂t
```
where γ(r) smoothly increases from 0 at R_inner to 1 at R.

This damps outgoing radiation without reflection. The sponge width
should be ~2-3 wavelengths of the radiation being absorbed.

### Vacuum clamping at r = R
Beyond the sponge: q(r ≥ R) = ρ₀(1, 0, 0, 0) (vacuum).

## Time Evolution

Leapfrog (Störmer-Verlet), same as scatter.c:
```
v(t+dt/2) = v(t-dt/2) + dt × c² × F(q(t))
q(t+dt)   = q(t)      + dt × v(t+dt/2)
```

CFL condition: dt < h/c for L₂, dt < h²/(some constant) for L₄.
With L₆: need to check — may require dt < h³ scaling.

## Diagnostics (per snapshot)

| Quantity | Method | Cost |
|----------|--------|------|
| B (Skyrmion charge) | Volume integral of B⁰ | O(N³) |
| H (Hopf charge) | FFT-based A∧F integral | O(N³ log N) |
| Energy (E₂,E₄,E₆,V) | Volume integral | O(N³) |
| Soliton positions | Weighted centroids | O(N³) |
| Field snapshots | Slice dumps (z=0 plane) | O(N²) |

## Implementation Plan

### Phase 1: Core infrastructure (`src/spherical_grid.h`, `src/spherical_grid.c`)
- Grid allocation with spherical mask
- Cell enumeration (active cell list for sparse iteration)
- Neighbor lookup (with boundary handling)
- Field I/O (save/load snapshots)

### Phase 2: Physics (`src/hopfion_field.c`)
- Energy functional: L₂ + L₄ + L₆ + V + V_π
- Force computation: -δE/δq for each term
- Gradient verification (finite-difference check)

### Phase 3: Topology (`src/topology.c`)
- Skyrmion charge B (existing, port from field.c)
- CP¹ projection q → n̂
- Hopf charge H via FFT method
- Linking number computation for composed states

### Phase 4: Initialization (`src/hopfion_init.c`)
- Single hopfion (Hopf map, (n,m) parameterization)
- Composed hopfions (product/superposition at offsets)
- Skyrmion-hopfion hybrid states
- Profile loading from 1D solvers

### Phase 5: Time evolution (`src/hopfion_evolve.c`)
- Leapfrog integrator
- Sponge absorbing boundary
- Sigma-model projection (optional)
- Diagnostics and output

### Phase 6: Analysis and atom building
- Equilibrium finding (gradient flow)
- Stability analysis (small perturbation evolution)
- Multi-hopfion binding energy
- Build toward H₁ atom analog (proton + electron as composed hopfions)

## File Structure

```
hopfion_composition/
├── DESIGN.md              ← this file
├── Makefile
├── literature/
│   └── INDEX.md           ← paper catalog
├── src/
│   ├── clifford.h         ← symlink or copy from hopfion_search
│   ├── spherical_grid.h   ← grid types and inline functions
│   ├── spherical_grid.c   ← grid allocation, I/O
│   ├── hopfion_field.h    ← energy, force, params
│   ├── hopfion_field.c    ← L₂+L₄+L₆+V computation
│   ├── topology.h         ← B, H charge computation
│   ├── topology.c         ← topological invariants
│   ├── hopfion_init.h     ← initialization routines
│   ├── hopfion_init.c     ← hopfion/skyrmion init
│   ├── hopfion_evolve.c   ← main time-evolution driver
│   └── hopfion_test.c     ← verification tests
├── data/
│   └── profiles/          ← 1D radial profiles
└── results/
    └── (simulation output)
```

## Parameters

| Parameter | Symbol | Typical Value | Role |
|-----------|--------|---------------|------|
| Vacuum density | ρ₀ | 1.0 | Sets energy scale |
| Skyrme coupling | e | 1.0–4.0 | Soliton size |
| Bulk coupling | λ | 0–∞ | σ-model limit |
| Sextic coupling | λ₆ | TBD | BPS limit control |
| Pion mass | m_π | 0–0.5 | Asymptotic decay |
| Degenerate mass | μ | 0–5.0 | Weight sector |
| Grid radius | R | 8–15 | Domain size |
| Grid spacing | h | 0.05–0.1 | Resolution |
| Sponge width | Δ_sponge | 2.0 | Absorbing layer |
| Speed of light | c | 1.0 | Time scale |

## Memory Estimate

For R=10, h=0.08 (N_eff ≈ 250):
- Full cube: 250³ = 15.6M cells
- Sphere: ~8.2M active cells
- Per cell: 8 doubles (field) + 8 doubles (velocity) + 8 doubles (force) = 192 bytes
- Total: 8.2M × 192 = 1.57 GB
- With L₆: need baryon current storage: +4 doubles/cell = +0.26 GB
- **Total: ~1.8 GB** (fits in RAM)
