# Locality Option 4: Retarded Integration (Exact Causality)

## Thesis

Replace instantaneous Poisson with the RETARDED Green's function:

    Φ(x,t) = α ∫ G(x-x') · ρ(x', t - |x-x'|/c) dx'

This uses the source ρ at the retarded time t_ret = t - |x-x'|/c.
Exactly causal, exactly correct, no free parameters.

## Implementation Challenge

Requires storing the full history of ρ(x,t) at every grid point for
the time window t - x_max/c ≤ t' ≤ t. With x_max=200 and c=1: need
200 time units of history. At dt=0.025 and Nx=8000: 8000 snapshots ×
8000 points = 512M doubles (4 GB). Too much.

**Optimization**: store ρ only at a COARSE time resolution. With
Δt_store = 1.0 (40× coarser): 200 snapshots × 8000 points = 12.8M
doubles (100 MB). Interpolate linearly between stored snapshots.

**Alternative**: for 1D, the retarded Green's function is:
    G_ret(x,t) = (c/2)·θ(ct - |x|)

So Φ(x,t) = (αc/2) ∫_{-ct}^{ct} ρ(x', t - |x-x'|/c) dx'

In 1D, the retarded integral is over a finite interval [x-ct, x+ct].
We only need ρ history within the light cone.

## No Free Parameters

The retarded solution has zero tunable parameters — it's the exact
physical answer. The speed c and coupling α are the only inputs.

## Method

### Phase 1: Single oscillon retarded Φ

1. Equilibrate oscillon t=5000 (no gravity)
2. Start storing ρ history at t=5000
3. At each timestep t > 5000: compute Φ(x,t) from retarded integral
4. Apply backreaction: m²_eff = m²(1+2Φ), c²_eff = 1+4Φ
5. Evolve t=5000 with retarded gravity
6. Measure: Φ(0) profile, oscillation, fc

### Phase 2: Causality verification

7. At t=7500: boost the oscillon (give velocity v=0.01)
8. Track Φ at x=50: it should respond at t = 7500 + 50/c = 7550
9. Verify the delay is exactly 50 time units (c=1)

### Phase 3: Compare with other options

10. Compare Φ profile with EMA-wave (Option 1), telegraph (Option 3),
    and Poisson (Test F)
11. The retarded solution is the REFERENCE — others should approximate it

## Implementation Details

Store ρ in a circular buffer: ρ_history[t_index % N_history][x_index]
where N_history = (int)(2*xmax/c / Δt_store) + 1.

At each timestep, compute Φ(x) by integrating over the past light cone:
    Φ(x) = (α/2) Σ_{x'=x-c(t-t_start)}^{x+c(t-t_start)} ρ_hist(x', t-|x-x'|/c) · dx'

Use linear interpolation in time for ρ between stored snapshots.

## Reference Code

- v24/fundamental/testF_selfref/src/selfref.c (gravity framework)

## Output

- `src/locality_retarded.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, α=-1e-4, c=1.0
Δt_store = 0.5 (history resolution)
N_history = 800 (400 time units of history)
Nx=4000, xmax=100, t_equil=5000, t_test=5000

Compile: `gcc -O3 -Wall -o locality_retarded src/locality_retarded.c -lm`
