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
With L₆: dt ~ h³ (the sextic term involves 3 spatial derivatives cubed).
For h=0.1: dt ~ 0.001 (vs 0.01 for L₄ alone). Manageable for gradient flow;
prohibitive for long-time dynamics. Mitigation: arrested Newton flow
(evolve, freeze when E increases) or skip L₆ force every M timesteps.

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
- Energy functional: L₂ + L₄ + L₆ + V + V_π (COMPLETE)
- Force computation: -δE/δq for each term (COMPLETE)
  - L₂: consistent 9-point Laplacian {1,-16,64,16,-130,16,64,-16,1}/(144h²)
  - L₄: analytical 3-pass Skyrme force (ported from field.c)
  - L₆: numerical finite differences of (B⁰)² energy
  - V, V_π, V_D: analytical forces
- Effective metric computation (BLV acoustic metric)
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
- With Skyrme force: SkyrmePre (7 doubles) + pi (12 doubles) per cell = ~312 MB temporary
- **Total: ~2.1 GB** (fits in RAM)

## Effective Metric and Emergent Gravity

### BLV Acoustic Metric

The Babichev-Langlois-Vernieri (BLV) effective metric for small fluctuations
around a soliton background φ₀ is:

```
g_eff^{μν}(x) = -∂²L/∂(∂_μφ^A)∂(∂_νφ^B)
```

For L₂ + L₄ Skyrme theory, this gives:
- **L₂ contribution**: g^{μν} = η^{μν} (flat Minkowski, always)
- **L₄ contribution**: depends on background currents A_d = q̃∂_dq

### Sturm-Liouville Coefficients

For a hedgehog soliton f(r), the effective metric reduces to:
```
P(r) = 2r² + 4c₄ sin²f(r)   (radial stiffness)
m(r) = r² + 2c₄ sin²f(r)    (temporal inertia)
```
where c₄ = 2ρ₀²/e².

### P/m = 2 Identity (Sigma Model)

**Critical result**: For the sigma model (|q| = ρ₀ = const), we have
c₄ = 2ρ₀²/e² and therefore:

```
P(r)/m(r) = (2r² + 4c₄ sin²f) / (r² + 2c₄ sin²f) = 2
```

This ratio is **identically 2** at every radius. This means the effective
line element is conformally flat — there is NO Schwarzschild-like time
dilation in the sigma model. The soliton does not generate emergent gravity.

### What Breaks P/m = 2

1. **Finite-λ (varying ρ(r))**: When the Mexican hat potential has finite λ,
   ρ(r) < ρ₀ near the soliton core. The effective c₄(r) = 2ρ(r)²/e² now
   varies with r, and P(r)/m(r) ≠ 2.

2. **L₆ term**: The sextic term adds additional metric contributions that
   are not proportional to the L₄ ones.

3. **Dynamic sector**: If the degenerate (J,P) sector becomes dynamical
   (with explicit kinetic term L₂_D), its coupling to the bulk could
   break the identity.

### Gravity Implementation Approach

To test "gravity as field density reduction":
1. Compute equilibrium at finite-λ (ρ(r) varies)
2. Compute effective metric g_00(r), g_rr(r)
3. Compare to Schwarzschild: g_00 = -(1-r_s/r), g_rr = (1-r_s/r)⁻¹
4. Extract effective r_s and compare to soliton mass

This requires the finite-λ solver (already working in hopfion_search) to
be ported, or loading a finite-λ profile from file.

## Consistent Laplacian

The energy E₂ = (1/2) Σ |D_d q|² uses 4th-order central differences D_d.
The exact discrete gradient of E₂ is Σ_d D_d†(D_d q), where D_d† is the
transpose of the derivative operator. This gives a 9-point stencil:

```
{1, -16, 64, 16, -130, 16, 64, -16, 1} / (144 h²)
```

This is NOT the same as the standard 5-point Laplacian {-1, 16, -30, 16, -1}/(12h²).
Using the wrong Laplacian creates an energy-force inconsistency that prevents
gradient flow from converging to the true minimum. The boundary skip must be
increased from ±2 to ±4 grid points (4.5h from sphere edge).
