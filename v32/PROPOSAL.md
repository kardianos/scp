# V32: SPH Field Theory — Particle-Based Simulation

## The Insight

Fixed grids can't represent metric contraction. When the field concentrates,
the COORDINATES should contract with it. On a fixed grid, we fake this with
c(ρ), which creates artifacts (freezing, energy pumping, sign problems).

An SPH (Smoothed Particle Hydrodynamics) approach: the simulation points
ARE the field. They move, cluster, spread. Where the field concentrates,
points cluster → shorter inter-point spacing → faster local dynamics.
This IS metric contraction. No c(ρ) formula needed.

## The Equation (unchanged)

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a

    V(P) = (μ/2)P²/(1+κP²),  P = φ₀φ₁φ₂

Same physics as V28. The difference is HOW we compute ∇²φ:
- Grid: finite differences with fixed dx
- SPH: kernel sum over neighbors with VARIABLE spacing

## SPH Formulation

### Each particle i carries:
    x_i = (x, y, z)    — position (MOVES)
    φ_i = (φ₀, φ₁, φ₂) — field values
    v_i = (v₀, v₁, v₂)  — field time derivatives
    m_i                   — SPH mass (fixed, sets resolution)
    h_i                   — smoothing length (adapts)
    ρ_sph_i               — SPH density (number density of particles)

### SPH kernel
    W(r, h) = (1/π^{3/2} h³) × exp(-r²/h²)   [Gaussian]
    or cubic spline (faster, compact support)

### Density
    ρ_sph(x_i) = Σ_j m_j × W(|x_i - x_j|, h_i)

### Smoothing length (adapts to local density)
    h_i = η × (m_i / ρ_sph_i)^{1/3}
    η ≈ 1.2 (target ~40-60 neighbors in 3D)

### Laplacian (second derivative)
    ∇²φ_a(x_i) = Σ_j (m_j/ρ_sph_j) × (φ_a_j - φ_a_i) × F(r_ij, h_i)
    F(r, h) = 2 × (∂²W/∂r² + 2/r × ∂W/∂r)

    Or the Brookshaw (1985) formulation:
    ∇²φ(x_i) = 2 × Σ_j (m_j/ρ_j) × (φ_i - φ_j)/(|r_ij|² + 0.01h²) × r_ij·∇W(r_ij, h)

### Time stepping
    Velocity Verlet (same as grid version, symplectic)
    dt set by: min_i(h_i / c_signal_i) × CFL
    c_signal ~ max eigenvalue of local wave equation ≈ 1 + |v_transport|

## Particle Transport (what moves the particles)

### Option 2: Energy-Flux Transport (primary)

Particles flow toward energy concentrations:

    dx_i/dt = v_transport_i = -β × ∇ρ_energy(x_i) / ρ_energy(x_i)

    ρ_energy = Σ_a [½v_a² + ½|∇φ_a|² + ½m²φ_a²] at each particle

This is a relaxation toward energy concentrations. β controls the transport
speed. Dense field → particles cluster → inter-particle spacing shrinks →
METRIC CONTRACTION.

With β=0: particles fixed (non-uniform grid, no advantage).
With β>0: particles migrate toward energy, creating adaptive resolution.

The transport velocity is SEPARATE from the field velocity. The field
oscillates (v_a); the particles drift (v_transport). Two different dynamics.

### Option 3a: Mesh Quality Minimization — Equal Energy per Particle

Move particles so each carries approximately equal energy:

    dx_i/dt = -γ × ∇(E_i - E_target)

    E_i = (volume_i) × ρ_energy_i ≈ (h_i³) × ρ_energy_i

This equalizes energy per particle → more particles where energy is dense.

### Option 3b: Mesh Quality Minimization — Centroidal Voronoi

Move each particle toward the centroid of its Voronoi cell, weighted by
energy density:

    dx_i/dt = -γ × (x_i - x_centroid_i)

This produces a smoothly graded mesh that tracks energy concentrations.
Classic approach in computational geometry.

### Option 3c: Mesh Quality Minimization — Spring Network

Connect each particle to its ~20 nearest neighbors with springs:

    dx_i/dt = Σ_j k × (|r_ij| - r_target) × r_hat_ij

    r_target = h_ideal × (ρ₀/ρ_energy)^{1/3}

Shorter target spacing where energy is high → particles cluster.
This is the simplest mesh-quality approach.

## Implementation Plan

### Phase 1: Basic SPH Wave Equation (Option 2)
1. Initialize 200K particles (uniform + extra near braid center)
2. Each particle carries φ_a, v_a, position
3. Neighbor search: spatial hash (O(N) with hash table)
4. SPH Laplacian for the wave equation
5. Velocity Verlet time stepping
6. Energy-flux transport (Option 2) for particle motion
7. Single braid: does it survive in SPH?

### Phase 2: Two Braids
If the braid survives, place two braids and measure separation.
The particles should cluster around each braid (metric contraction).
Other braid's particles entering the cluster → geodesic deflection → attraction?

### Phase 3: Alternative Transports (Options 3a-c)
Test the three mesh-quality approaches. Compare to Option 2.
Which gives the most stable braid? Which gives the strongest gravity signal?

## Grid vs SPH Comparison

| Aspect | Grid (V28-V31) | SPH (V32) |
|--------|---------------|-----------|
| Resolution | Fixed (dx=0.3) | Adaptive (h~0.1 near braid, h~1 far) |
| Memory | N³ × 9 × 8 bytes = 150MB at N=128 | N_part × 20 × 8 bytes = 32MB at 200K |
| Metric | Fixed coordinates, c(ρ) proxy | Moving coordinates, metric from spacing |
| Force cost | O(N³) per step | O(N_part × N_neighbor) per step |
| Neighbor | Trivial (grid) | Hash table or tree (O(N log N)) |
| c(ρ) | Needed (source of problems) | NOT needed (geometric) |

## Parameters

    N_particles = 200,000 (start), up to 500K
    Particle mass m_sph = (L/N_part^{1/3})³ × ρ₀
    Smoothing: η = 1.2 (gives ~50 neighbors)
    Transport: β = 0.01 (gentle drift toward energy)
    CFL: 0.3
    Field: m=1.50, μ=-41.3, κ=50 (same as V28)
    Domain: L=20, periodic in z, free in x,y (or fully periodic)

## Output

Each snapshot saves: particle positions + field values + SPH density.
Format: binary (N_part, then N_part × {x,y,z,φ₀,φ₁,φ₂,v₀,v₁,v₂,ρ_sph,h}).
Loadable in Python for scatter-plot visualization.

## Risks

- SPH Laplacian is noisy (kernel approximation errors)
- The triple-product V(P) needs all three fields aligned → sensitive to noise
- Particle transport might disrupt the braid (particles drift out of the braid)
- The braid's helical structure needs aligned particles to represent
- Neighbor search is the bottleneck for large N_part

## Expected Runtime

200K particles, ~50 neighbors each:
- Force: 200K × 50 × 30 FLOP ≈ 300M → 0.006s at 50GFLOP/s
- Neighbors: 200K hash lookups → 0.002s
- Transport: 200K → 0.001s
- Total: ~0.01s per step (10× faster than N=128 grid)
- T=300 at dt≈0.1: 3000 steps → 30s total

Very fast. Can iterate quickly on the approach.
