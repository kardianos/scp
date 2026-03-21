# Structure Contrast Analysis: What Works vs What Fails

## Side-by-Side Comparison

| Metric (t=0)        | braid3(z) | trefoil_R | figure8  |
|---------------------|-----------|-----------|----------|
| E_pot               | **-74.3** | -45.5     | -51.2    |
| P_integrated        | **142**   | 76        | 80       |
| R_rms               | **9.3**   | 22.6      | 22.5     |
| R_core (50% E_pot)  | **3.0**   | 5.0       | 4.4      |
| Phase coherence     | **0.135** | 0.006     | 0.006    |
| Concentration       | 0.52      | **0.65**  | 0.58     |
| Grad/mass ratio     | 0.083     | 0.106     | **0.136**|
| Box size L          | 15        | 25        | 25       |
| Boundary            | periodic  | absorbing | absorbing|

| Metric (t=300)      | braid3(z) | trefoil_R | figure8  |
|---------------------|-----------|-----------|----------|
| E_pot               | **-35.4** | -0.003    | -0.003   |
| E_pot retention     | **47.6%** | 0.007%    | 0.005%   |
| P_integrated        | **40.9**  | 1.6       | 1.5      |
| Phase coherence     | **0.193** | 0.479     | 0.483    |
| Survived?           | **YES**   | NO        | NO       |

## Why Braid3 Works — Five Critical Features

### 1. Phase coherence at initialization (0.135 vs 0.006)
The braid3 fields are initialized with specific phase offsets δ = {0, 3.0, 4.4}
along the z-axis. At every point, the three fields have a FIXED angular
relationship: φ₁/φ₀ ≈ cos(3.0)/cos(0) = -0.99. This means V(P) = V(φ₀φ₁φ₂)
is consistently negative throughout the tube.

The knots have phase coherence of 0.006 — essentially random. The background
field (A_bg=0.1) fills the entire box with uncorrelated phases, diluting the
knot tube's phase structure. The knot tube is a small perturbation on an
incoherent background.

**Lesson: The phase relationship between fields must be maintained everywhere
the structure exists, not just at discrete crossing points.**

### 2. Compact energy distribution (R_rms=9.3 vs 22.6)
Despite using a smaller box (L=15 vs L=25), braid3 has R_rms=9.3 — the energy
is concentrated within ~60% of the box radius. The knots have R_rms=22.6 —
energy fills 90% of the box, mostly from the background.

The braid3 amplitude (A=0.8) dominates the background (A_bg=0.1) in the core,
giving a signal-to-noise ratio of ~8:1. The knot tube also uses A=0.8 but the
tube volume is smaller relative to the box, so the background contributes more.

**Lesson: The structure must dominate its surroundings. The signal (bound
fields) must be much larger than the noise (background).**

### 3. Traveling wave maintains phase lock
The braid3 initial velocity is ω·A·sin(kz + δ_a) — a traveling wave propagating
in +z. This means the phase relationship between fields is DYNAMICALLY
maintained: the wave propagates, and V(P) keeps the three fields locked
together as they move.

The knots start stationary (vel=0). The V(P) forces create oscillations, but
there's no propagation to maintain coherence. The standing wave breathes and
radiates, with each breathing cycle losing phase coherence.

**Lesson: A propagating mode can maintain phase coherence dynamically.
A stationary mode must rely entirely on V(P) restoring forces, which may
not be strong enough.**

### 4. Periodic BC preserves the traveling wave
The braid3 uses periodic boundary conditions. The traveling wave exits one side
and re-enters from the other, creating a continuous loop. The wave NEVER hits
a boundary — it just wraps around. Energy is perfectly conserved (drift 0.37%).

The knots use absorbing BC, which drains 78% of the total energy. Each breathing
cycle radiates waves that get absorbed, steadily weakening the structure.

**Lesson: For compact objects, we need a structure that doesn't radiate at all —
a true stationary soliton. Or: the structure must be its own boundary (closed,
self-contained).**

### 5. Smooth radial profile (no sharp features)
The braid3 has a smooth Gaussian radial profile: rho drops from 9.9 at center
to 0.015 at the background level, with no sharp edges or cusps.

The knot curves have sharp bends (high curvature at crossings). These create
high gradient energy that wants to smooth out. The initial condition is far
from any equilibrium state.

**Lesson: Sharp spatial features radiate. The initial condition must be
smooth and close to equilibrium.**

## What a Self-Sustaining Compact Structure Needs

From this analysis, a compact structure that survives must have:

1. **Phase-locked fields everywhere** (not just at crossings)
2. **Dominant amplitude** over background in the core region
3. **Smooth spatial profile** (no sharp bends or cusps)
4. **Low radiation** (close to an equilibrium state)
5. **Self-contained** (doesn't rely on periodic wrapping)

The braid3 achieves 1-4 but not 5 (it wraps around periodically).
The knots attempt 5 but fail on 1-4.

## Implication for Backprop Approach

Instead of designing a shape (trefoil, figure-eight, etc.) and hoping it works,
we should let gradient descent find the optimal field configuration:

- **Weights**: φ_a(x,y,z,t=0) at every grid point (~100k parameters at N=32)
- **Forward pass**: N timesteps of the PDE
- **Loss**: L = -|E_pot(T)|/|E_pot(0)| + λ₁·R_rms² + λ₂·(1-phase_coherence)
  - Maximize binding retention
  - Penalize spatial spreading
  - Reward phase coherence
- **Backward pass**: adjoint method (same PDE reversed)
- **Optimizer**: Adam

This is tractable at N=32 (32768 grid points per field, ~100k total parameters).
The PDE is differentiable, and the adjoint method gives exact gradients.

The key advantage: backprop doesn't assume any particular geometry. It can
discover structures we wouldn't think to try — potentially novel topologies
that naturally satisfy all five requirements.

## Files

| File | Contents |
|------|----------|
| `data/braid3_ref.sfa` | braid3(z) reference run (f64, 102 frames) |
| `data/struct_braid3.tsv` | Structure analysis time series |
| `data/struct_trefoil_R.tsv` | Trefoil R structure analysis |
| `data/struct_figure8.tsv` | Figure-eight structure analysis |
| `src/sfa_structure.c` | Structure analysis tool |
