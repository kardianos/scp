# V24-A: Vortex Ring (Circulation-Stabilized Bound State)

## Thesis

A vortex ring is a topological object that EXISTS only when moving. At v=0
it collapses; at v>0 it has definite size, energy, and angular momentum.
The three displacement fields φ₁, φ₂, φ₃ form a circulating pattern where
each field drives the next: φ₁ → φ₂ → φ₃ → φ₁.

This requires the elastic interpretation: field index a = spatial direction a.
A vortex ring is a divergence-free (∇·φ = 0) circulating structure.

Since vortex rings require 3D, this test uses the v21 3D infrastructure.
However, 3D is expensive, so we first test the 2D cross-section (a vortex
PAIR in 2D) to validate the concept.

## Mathematical Setup

### 2D Vortex Pair (Testable Proxy)

In 2D (x,y), with two displacement fields φ₁ (x-displacement) and
φ₂ (y-displacement), a vortex has:

    φ₁ = -Γ·(y-y₀)/r²,  φ₂ = +Γ·(x-x₀)/r²

where r² = (x-x₀)² + (y-y₀)², Γ is the circulation strength.

This is divergence-free: ∂φ₁/∂x + ∂φ₂/∂y = 0.

For the three-field system: set φ₃ = 0 (or use it as the propagation
direction in the full 3D ring).

### The Problem

The vortex profile above has 1/r singularity at the core and 1/r decay
at infinity — it has infinite energy. A physical vortex needs:
1. A regularized core (finite φ at r=0)
2. A mass gap to contain the field (exponential decay at infinity)

With the mass term m²φ²: the vortex profile must satisfy:
    ∇²φ_a - m²φ_a - ∂V/∂φ_a = 0

For a single massive scalar vortex: there's NO finite-energy solution
(the mass term prevents the 1/r tail from existing).

UNLESS the coupling V(P) provides a mechanism for vortex stabilization.
With three fields and the triple product: the vortex circulation creates
P ≠ 0 in the core, which provides the binding.

### What to Test

The key question is whether the three-field system with triple-product
coupling supports VORTEX-like solutions (divergence-free, circulating).
This is fundamentally different from the oscillon (breathing, compressive).

## What to Compute

### Phase 1: 2D Vortex Initialization

1. Use a 2D grid (Nx × Ny) with two active fields φ₁, φ₂ (φ₃ = 0 or small).
2. Initialize a vortex-antivortex pair:
   - Vortex at (x₀, y₀) with circulation +Γ
   - Antivortex at (-x₀, -y₀) with circulation -Γ
   Using the regularized Lamb-Oseen profile:
     φ₁ = -Γ·(y-y₀)·(1-exp(-r²/a²))/r²
     φ₂ = +Γ·(x-x₀)·(1-exp(-r²/a²))/r²
   where a is the core radius.

3. Set φ₃ = small random noise (seed for potential triple-product effects).

4. Evolve the 2D system (equations same as 3D but on a 2D grid):
   ∂²φ_a/∂t² = ∂²φ_a/∂x² + ∂²φ_a/∂y² - m²φ_a - ∂V/∂φ_a

5. Measure: does the vortex persist? Does it propagate? Does φ₃ grow?

### Phase 2: 2D Vortex with Third Field

6. If the vortex survives Phase 1: initialize with φ₃ ≠ 0 in the vortex core.
7. The triple product P = φ₁φ₂φ₃ is nonzero in the core → binding.
8. Measure whether the three-field coupling HELPS or HURTS the vortex.

### Phase 3: Propagation Velocity

9. A vortex-antivortex pair propagates at velocity v ∝ Γ/(2πd) where
   d is the pair separation. Measure the propagation velocity.
10. Does the pair speed up, slow down, or maintain velocity?

## Reference Code

- v21 3D solver: `/home/d/code/scp/v21/src/triad3d.c` (adapt for 2D)
- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (simpler base)

## Output

- `src/vortex2d.c` — 2D three-field solver with vortex initialization
- `data/vortex_ts.tsv` — time series
- `data/vortex_profile_t{T}.tsv` — field snapshots
- `RESULTS.md` — analysis

## Parameters

μ=-20, κ=20, m=1.0
Grid: Nx=Ny=512, L=40 (dx=0.156)
Vortex: Γ=1.0, a=2.0 (core radius), pair separation d=10
φ₃_noise = 0.01 amplitude
t_run = 2000

Compile: `gcc -O3 -Wall -o vortex2d src/vortex2d.c -lm`

Note: 2D grid 512² = 262k points × 3 fields × 3 arrays = ~19 MB. Fast.
