# Soliton Search v2: Design Document

## Diagnosis: Why Everything Has Failed

### The core problem in one sentence
Random seeds have zero binding (E_pot = -0.000), and gradient descent
cannot bootstrap binding from nothing because V(P) has a flat landscape
at P=0.

### Evidence
- `results_n32_partial`: 500 epochs, 32 seeds, best E_pot(T) = -0.000024.
  `alive_seeds = 0` for all 500 epochs. The optimizer shrinks R_rms
  (compactness loss) but never develops binding.
- `results_n48`: 1247 epochs, best E_pot(T) = -0.000000. Same failure.
- The loss is dominated by the compactness term; the binding term
  contributes nothing because V(P) is near zero everywhere.

### Why the braid3 seed works (and random doesn't)
The braid3 seed starts with E_pot = -70 to -125. The optimizer has
a signal to work with: it can reduce radiation by 3-4% tweaks to the
tail. Random seeds start at E_pot ~ 0. The V(P) potential is:
```
V(P) = (mu/2) P^2 / (1 + kappa P^2)
```
where P = phi_0 * phi_1 * phi_2. For small fields (A ~ 0.3-0.5 after
smoothing), P ~ 0.003 and V ~ -0.0002. The gradient is proportional
to P, so at P~0 the gradient is also ~0. **There is no gradient signal
to create binding from nothing.**

### The traveling wave requirement
The braid3 survives because it is a traveling wave. Standing waves
breathe and radiate. Every compact structure tested (knots, rings,
spherical harmonics, random spheres) starts stationary and radiates
to death within ~30 time units. The truncated helix survives because
its breathing mode is a feature, not a bug — but only with periodic
BC. With absorbing BC, breathing = radiation = death.

### What backprop actually optimized
When seeded with braid3, backprop made tiny changes (3-4% in core,
12-27% at boundary) that reduced radiation at the braid-vacuum
interface. It did NOT discover a new structure. It polished an
existing one.


## What Must Be Different This Time

### Principle 1: Seed with binding, not with noise
Every successful structure has initial E_pot in the range -45 to -170.
We must start there. Seeds must have P = phi_0 * phi_1 * phi_2
significantly nonzero in a coherent region.

### Principle 2: The search space is PARAMETERS, not voxels
Optimizing 3*32^3 = 98k raw voxels is the wrong abstraction. The
structures that work are described by ~10-20 parameters: amplitude,
wavelength, phase offsets, envelope shape, tube radius, propagation
direction. Searching in parameter space is 5000x more efficient.

### Principle 3: The compact soliton may not exist
The braid3 only survives with periodic BC (wrapping traveling wave).
With absorbing BC, nothing has survived beyond t=200. This might be
a fundamental property of the PDE — the wave equation with m^2>0
might simply not support compact solitons. The search must be
designed to produce a definitive answer: either find one or prove
the search was thorough enough to conclude there are none.

### Principle 4: Evolutionary selection at low res, gradient polish at high res
The Game of Life analogy is correct for the exploration phase. But
once a candidate is found, gradient refinement is the right tool to
optimize it. Hybrid = correct approach.


## The Strategy: Three-Phase Hierarchical Search

### PHASE 1: Parametric Evolutionary Search (N=16, ~2 hrs)

**Goal**: Find parameter combinations that produce binding at low resolution.

**Key difference from field_evolve.py**: Search in a ~30-dimensional
parameter space, not in a 3*16^3 = 12288 voxel space. Every individual
is a parametric field generator, not raw field values.

#### Parameter genome (30 floats)

```python
genome = {
    # Core structure (12 params)
    'type':        int,    # 0=braid, 1=torus, 2=crossed, 3=spherical, 4=shell
    'A':           float,  # amplitude [0.3, 1.5]
    'k':           float,  # wavenumber [0.1, 2.0]
    'delta_1':     float,  # phase offset field 1 [0, 2*pi]
    'delta_2':     float,  # phase offset field 2 [0, 2*pi]
    'R_tube':      float,  # tube/envelope radius [0.5, 5.0]
    'sigma_z':     float,  # z-envelope width [1.0, 8.0]
    'n_osc':       int,    # oscillations per wavelength [1, 4]
    'R_helix':     float,  # helix radius [0.3, 3.0]
    'n_strands':   int,    # 1 or 3
    'R_major':     float,  # torus major radius [2.0, 6.0]
    'chirality':   int,    # +1 or -1

    # Envelope (6 params)
    'env_type':    int,    # 0=gaussian, 1=super-gaussian, 2=tanh, 3=cosine
    'env_sigma':   float,  # envelope width [1.0, 6.0]
    'env_power':   float,  # super-gaussian power [2, 8]
    'env_offset':  float,  # center offset from origin [0, 3.0]
    'env_asym':    float,  # z-asymmetry [-0.3, 0.3]
    'env_smooth':  float,  # transition width [0.5, 3.0]

    # Velocity (6 params)
    'v_type':      int,    # 0=traveling, 1=standing, 2=breathing, 3=rotating
    'v_amp':       float,  # velocity amplitude [0, 1.5]
    'v_dir':       float3, # propagation direction (3 components, normalized)

    # Composition (6 params)  -- for crossed/multi-structure types
    'n_axes':      int,    # 1, 2, or 3
    'cross_A':     float,  # relative amplitude of secondary axes [0, 1]
    'cross_phase': float,  # phase between axes [0, 2*pi]
    'A_bg':        float,  # background amplitude [0, 0.3]
    'bg_corr':     float,  # background correlation with structure [-1, 1]
    'mix_mode':    int,    # 0=additive, 1=multiplicative, 2=max-envelope
}
```

#### Generator function

```python
def genome_to_field(genome, X, Y, Z, R, N, L):
    """Convert 30-parameter genome to (phi[3,N,N,N], vel[3,N,N,N])."""

    if genome['type'] == 0:  # braid
        # Helical braid along z-axis (or rotated direction)
        # with optional z-truncation envelope
        ...
    elif genome['type'] == 1:  # torus
        # Single torus tube with oscillation along toroidal angle
        # Three fields share one tube with phase offsets
        ...
    elif genome['type'] == 2:  # crossed
        # 1-3 braids along different axes, superimposed
        ...
    elif genome['type'] == 3:  # spherical
        # Radial oscillation with Y_lm angular structure
        ...
    elif genome['type'] == 4:  # shell
        # Spherical shell with tangential oscillation
        ...

    # Apply envelope
    env = make_envelope(genome, R, Z)
    phi = phi_core * env[None, :, :, :]

    # Apply velocity
    vel = make_velocity(genome, phi, X, Y, Z)

    return phi, vel
```

#### Evolution strategy: CMA-ES + islands

CMA-ES (Covariance Matrix Adaptation) is the right optimizer here
because:
- 30 continuous dimensions with unknown landscape
- Evaluation is cheap (~0.05s per individual at N=16)
- We want MULTIPLE diverse solutions, not just one
- CMA-ES naturally adapts its search distribution

**Island model**: Run 5 independent CMA-ES populations of size 64,
one per structure type (braid, torus, crossed, spherical, shell).
Every 20 generations, migrate the top 3 from each island to all
others. This prevents premature convergence while allowing
specialization.

```
5 islands x 64 population x 100 generations = 32,000 evaluations
At 0.05s each (N=16, T=20, ~200 steps): 1600s = 27 minutes

+ longer T_eval for top candidates: another ~30 minutes
Total Phase 1: ~1 hour
```

#### Fitness function (multi-objective)

```python
def fitness(genome, grid_params):
    phi, vel = genome_to_field(genome, ...)

    # Quick sanity check: does this config have binding AT ALL?
    epot_0 = compute_epot(phi, dx)
    if epot_0 > -0.5:
        return -1e6  # dead on arrival, don't waste time evolving

    # Evolve for T_eval with absorbing BC
    phi_T, vel_T = evolve(phi, vel, ...)

    # Primary: absolute binding at end (not relative!)
    epot_T = compute_epot(phi_T, dx)

    # Secondary: binding retention
    retention = epot_T / (epot_0 - 1e-10)

    # Tertiary: compactness
    rrms = compute_rrms(phi_T, ...)

    # Quaternary: coherence (is it one blob or fragments?)
    # Simple proxy: what fraction of |P| is in the central region?
    P_T = phi_T[0] * phi_T[1] * phi_T[2]
    P_center = sum(|P| where R < L/3)
    P_total = sum(|P|) + 1e-30
    concentration = P_center / P_total

    # Weighted fitness (higher = better)
    fit = (
        -5.0 * epot_T          # maximize binding (epot is negative)
        + 2.0 * retention       # maximize retention (0 to 1)
        - 0.5 * (rrms / L)**2  # penalize spreading
        + 1.0 * concentration   # reward localization
    )
    return fit
```

**Critical: early termination.** If epot_0 > -0.5 (no binding at
initialization), skip the evolution entirely and return minimum
fitness. This happens for ~70% of random parameter combos and saves
enormous time. At N=16 with T=20, each eval without early-term takes
~0.05s; with it, ~0.001s.

#### Increasing T_eval (curriculum within evolution)

Generations 0-30:   T_eval = 10  (200 steps at N=16)
Generations 30-60:  T_eval = 20  (400 steps)
Generations 60-100: T_eval = 40  (800 steps)

This is crucial: early generations select for "has binding at all"
(fast, cheap). Later generations select for "binding survives"
(slower, more expensive). Don't waste 800-step evals on configs
that don't even bind at t=0.


### PHASE 2: Validation and Upscaling (N=32/48, ~2 hrs)

**Goal**: Take the top 20 candidates from Phase 1, validate at higher
resolution, and refine with gradient descent.

#### Step 2A: Upscale and re-evaluate

Take each Phase 1 winner (parameterized as a genome). Regenerate
at N=32 (8x more voxels). Evolve for T=50 with absorbing BC.

Expected: many Phase 1 winners will fail at N=32 because:
- Low-res artifacts masked instability
- The actual wavelength structure doesn't survive discretization
  changes
- Breathing modes are resolution-dependent

This is fine — it's a filter. Keep the ones that survive.

```
20 candidates x T=50 at N=32 ≈ 20 x 1.5s = 30 seconds
```

#### Step 2B: Gradient refinement of survivors

For each surviving candidate (expect 3-8 of 20), use backprop to
optimize the RAW FIELD VALUES at N=32. But seed with the parametric
config, not random noise. This gives the optimizer a strong starting
signal (E_pot << 0).

```python
for candidate in phase1_survivors:
    phi, vel = genome_to_field(candidate.genome, N=32, L=15)

    # Freeze core, optimize envelope
    # (Based on backprop insight: core is already good, tail matters)
    core_mask = (R < 3.0)
    phi_fixed = phi * core_mask
    phi_learn = phi * (1 - core_mask)  # only optimize outer region

    # Or: optimize the full field but start from the parametric IC
    # which already has E_pot ~ -50 to -150
    params = {'phi': phi, 'vel': vel}
    for epoch in range(200):
        loss, grads = loss_and_grad(params, T_opt=30, absorbing=True)
        params = adam_step(params, grads)
```

Key differences from the failed field_search.py:
1. Seeded with binding (E_pot ~ -50), not random (E_pot ~ 0)
2. T_opt = 30-50 (not 5-10 from curriculum warm-up)
3. Optional: freeze core, optimize only envelope/tail

```
8 candidates x 200 epochs x 0.15s/epoch = 240s ≈ 4 minutes
```

#### Step 2C: Head-to-head tournament at N=48, T=100

Take the top 5 gradient-refined candidates. Upscale to N=48.
Evolve for T=100 with absorbing BC. Rank by binding retention.

```
5 candidates x T=100 at N=48 ≈ 5 x 15s = 75 seconds
```


### PHASE 3: Long-horizon validation (N=64/128, ~2 hrs)

**Goal**: Determine whether ANY candidate is a true soliton (survives
T=500) or just a slow-decaying transient.

#### Step 3A: N=64, T=300, absorbing BC

The definitive test. If a structure retains >30% of E_pot at T=300
with absorbing BC at N=64, it is a genuine compact soliton.

```
3 candidates x T=300 at N=64 ≈ 3 x 10 minutes = 30 minutes
```

#### Step 3B: C simulation at N=128, T=500

For any that pass 3A, hand off to the C code (v37_crossed_v2.c or
new v37_candidate.c) for full production run. Generate SFA for
visualization and fragmentation analysis.

**If no candidate passes 3A**, this is a definitive null result: the
Cosserat PDE with these parameters does not support compact solitons
with absorbing BC. Document and move on.


## Critical Implementation Details

### 1. Phase correlation in random parametric seeds

The phase offsets delta_1, delta_2 control the sign and magnitude of
P = phi_0 * phi_1 * phi_2. Not all phases give binding:

```
cos(0) * cos(d1) * cos(d2) = 1 * cos(d1) * cos(d2)
```

For V(P) < 0 (binding), we need P^2 large, which requires
|cos(d1) * cos(d2)| large. The known-good values are:
  d1 = 3.0005 (cos = -0.990)
  d2 = 4.4325 (cos = -0.307)
  product = +0.304

Random (d1, d2) gives E[|cos(d1)*cos(d2)|] = (2/pi)^2 ≈ 0.41.
So random phases are actually OK on average, but many will be near
zero (e.g., d near pi/2). The CMA-ES will quickly find the good
regions.

**Important**: The initial CMA-ES distribution for phases should be
seeded near the known-good values, with sigma ~ 1.0 to allow
exploration. Don't start uniform [0, 2pi].

### 2. Velocity initialization matters enormously

Standing wave (v=0) always breathes and radiates.
Traveling wave survives with periodic BC but walks off with absorbing.

**New idea: ROTATING field**. Instead of translating, have the fields
rotate in the internal (phi_0, phi_1, phi_2) space:

```python
# Rotating initial velocity: phi_a_dot = omega * epsilon_{abc} * phi_b * some_function
vel[0] = omega * (phi[1] * cos(angle) - phi[2] * sin(angle))
vel[1] = omega * (phi[2] * cos(angle) - phi[0] * sin(angle))
vel[2] = omega * (phi[0] * cos(angle) - phi[1] * sin(angle))
```

This gives angular momentum in field space without spatial
translation. If V(P) locks the phases, the rotation could provide
gyroscopic stability (the backprop result showed 6% z-asymmetry
helps — this is a hint toward rotation).

### 3. Multi-resolution consistency

N=16 has dx = 2*15/15 = 2.0. The braid wavelength is
lambda = 2*pi/k ~ 2*pi/0.21 ~ 30. At N=16 (dx=2.0), this is
15 grid points per wavelength — adequate for the fundamental.

But the binding potential V(P) depends on P = prod(phi), which is
a product of three oscillating fields. If each oscillates at k,
then P oscillates at ~3k, requiring 3x the resolution. At N=16,
this is only 5 points per P-oscillation — borderline.

**Resolution test**: Before committing to Phase 1 at N=16, verify
that the known-good braid3 seed produces qualitatively similar
binding at N=16 vs N=32 vs N=64. If N=16 is too coarse for P,
use N=20 or N=24 instead.

### 4. Absorbing BC calibration

The absorbing BC has two parameters:
- R_damp: where damping starts
- damp_rate: damping strength per step

At N=16 with L=15, R_damp = L - 4 = 11 leaves only a 4-unit buffer.
The braid structure has R_core ~ 3 and waves propagate at c ~ 1.5,
so the first radiation reaches the boundary at t ~ (11-3)/1.5 = 5.3.

**Calibration protocol**: Use the braid3 seed. Measure energy loss
rate vs (R_damp, damp_rate). Find the setting that:
- Absorbs >90% of outgoing radiation (no PBC wrap-around artifacts)
- Drains <5% of bound-state energy per crossing time

This is a 2D parameter scan: ~50 combinations, 0.05s each = 2.5s.

### 5. VRAM budget

V100 has 16GB VRAM.

Single field state: 3 * N^3 * 4 bytes (f32) each for phi and vel,
plus theta and theta_vel.
- N=16: 3 * 4096 * 4 * 4 = 196 KB
- N=32: 3 * 32768 * 4 * 4 = 1.5 MB
- N=48: 3 * 110592 * 4 * 4 = 5.3 MB
- N=64: 3 * 262144 * 4 * 4 = 12 MB

For JAX scan (no checkpointing), memory is O(n_steps * state_size):
- N=16, 200 steps: 200 * 196KB = 38 MB -- trivial
- N=32, 500 steps: 500 * 1.5MB = 750 MB
- N=48, 500 steps: 500 * 5.3MB = 2.6 GB
- N=64, 500 steps: 500 * 12MB = 6 GB

With gradient checkpointing (jax.checkpoint every 50 steps):
- N=64, 500 steps: 10 segments * 12MB + 50 * 12MB = 720 MB

**Verdict**: All phases fit in 16GB with checkpointing. N=64 with
backprop needs checkpointing. N=16/32 with CMA-ES (forward only,
no backprop) is negligible.

For CMA-ES at N=16 with pop=64: 64 * 196KB = 12.5 MB for the
entire population. Trivially fits. Can potentially vmap over
the population for massive parallelism.


## Pseudocode: Main Search Loop

```python
# =====================================================
# PHASE 1: Parametric CMA-ES (N=16, ~1 hr)
# =====================================================

N, L = 16, 15.0
grid = make_grid(N, L)

# Step 0: Verify N=16 is adequate resolution
braid_16 = make_braid3(N=16, L=15)
braid_32 = make_braid3(N=32, L=15)
epot_16 = compute_epot(braid_16.phi, grid.dx)
epot_32 = compute_epot(braid_32.phi, grid32.dx)
assert abs(epot_16/epot_32 - 1) < 0.3, f"N=16 too coarse: {epot_16} vs {epot_32}"

# Step 1: Calibrate absorbing BC
best_damp = calibrate_absorbing_bc(braid_16, grid)

# Step 2: Initialize 5 island populations
islands = {
    'braid':     CMA_ES(dim=30, pop=64, sigma=0.5, center=braid_prior),
    'torus':     CMA_ES(dim=30, pop=64, sigma=0.5, center=torus_prior),
    'crossed':   CMA_ES(dim=30, pop=64, sigma=0.5, center=crossed_prior),
    'spherical': CMA_ES(dim=30, pop=64, sigma=0.5, center=spherical_prior),
    'shell':     CMA_ES(dim=30, pop=64, sigma=0.5, center=shell_prior),
}

# Prior centers: use known-good parameters where available
braid_prior.delta = [0, 3.0005, 4.4325]
braid_prior.A = 0.8
braid_prior.k = pi/15
# Others: reasonable defaults from physical intuition

# Step 3: Evolution loop
T_evals = [10, 10, 20, 20, 40]  # curriculum by generation bracket
all_top = []

for gen in range(100):
    T_eval = T_evals[min(gen // 20, 4)]

    for name, island in islands.items():
        genomes = island.ask()  # get 64 candidates
        fitnesses = []

        for g in genomes:
            phi, vel = genome_to_field(g, grid)
            epot_0 = compute_epot(phi, grid.dx)

            if epot_0 > -0.5:
                fitnesses.append(-1e6)  # early termination
                continue

            phi_T, _ = evolve(phi, vel, T_eval, grid, best_damp)
            fit = compute_fitness(phi_T, grid, L)
            fitnesses.append(fit)

        island.tell(genomes, fitnesses)

    # Migration every 20 generations
    if gen % 20 == 0 and gen > 0:
        for name, island in islands.items():
            for other_name, other in islands.items():
                if name != other_name:
                    other.inject(island.top_k(3))

    # Track best per island
    for name, island in islands.items():
        log(gen, name, island.best_fitness, island.best_genome)

    # Save global top 20
    all_candidates = [(name, island.best_genome, island.best_fitness)
                      for name, island in islands.items()]
    all_candidates.sort(key=lambda x: -x[2])
    all_top = all_candidates[:20]


# =====================================================
# PHASE 2: Upscale + Gradient Refinement (N=32, ~30 min)
# =====================================================

N = 32
grid32 = make_grid(N, L)
damp32 = calibrate_absorbing_bc(make_braid3(N=32, L=15), grid32)

# 2A: Upscale and filter
survivors = []
for name, genome, phase1_fit in all_top:
    phi, vel = genome_to_field(genome, grid32)
    phi_T, _ = evolve(phi, vel, T=50, grid32, damp32)
    epot_T = compute_epot(phi_T, grid32.dx)
    retention = epot_T / compute_epot(phi, grid32.dx)
    if retention > 0.1:  # kept >10% of binding
        survivors.append((name, genome, phi, vel, epot_T))

print(f"Phase 2A: {len(survivors)}/{len(all_top)} survived upscale")

# 2B: Gradient refinement
refined = []
for name, genome, phi, vel, _ in survivors:
    params = {'phi': phi, 'vel': vel}
    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)

    for epoch in range(200):
        loss, aux, grads = loss_and_grad(params, T_opt=30)
        grads = clip_grads(grads)
        updates, opt_state = optimizer.update(grads, opt_state)
        params = optax.apply_updates(params, updates)

    phi_T, _ = evolve(params['phi'], params['vel'], T=100, grid32, damp32)
    epot_T = compute_epot(phi_T, grid32.dx)
    refined.append((name, params, epot_T))

refined.sort(key=lambda x: x[2])  # most negative E_pot first

# 2C: Head-to-head at N=48
N = 48
grid48 = make_grid(N, L)
top5 = refined[:5]
results_48 = []
for name, params, _ in top5:
    phi48, vel48 = upscale(params['phi'], params['vel'], 32, 48)
    phi_T, _ = evolve(phi48, vel48, T=100, grid48)
    epot_T = compute_epot(phi_T, grid48.dx)
    results_48.append((name, phi48, vel48, epot_T))


# =====================================================
# PHASE 3: Definitive validation (N=64, ~2 hrs)
# =====================================================

top3 = sorted(results_48, key=lambda x: x[3])[:3]

for name, phi, vel, _ in top3:
    phi64, vel64 = upscale(phi, vel, 48, 64)
    # Run for T=300 with diagnostics every dt=5
    for t in range(0, 300, 5):
        phi64, vel64 = evolve(phi64, vel64, T=5, grid64)
        epot = compute_epot(phi64, grid64.dx)
        rrms = compute_rrms(phi64, grid64)
        log(name, t, epot, rrms)

        if epot > -0.1:
            print(f"  {name}: DEAD at t={t}")
            break
    else:
        print(f"  {name}: SURVIVED T=300 with E_pot={epot:.1f}")
        # WINNER — save for N=128 C simulation
```


## Implementation: What to Build

### File 1: `parametric_seeds.py` (~300 lines)

The genome-to-field generator. Contains:
- `genome_to_field(genome, X, Y, Z, R, N, L) -> (phi, vel)`
- Five structure-type generators (braid, torus, crossed, spherical, shell)
- Envelope generators
- Velocity initializers (traveling, standing, breathing, rotating)
- `genome_bounds()` returning min/max for each parameter

Test: verify each structure type produces |E_pot| > 1 for
hand-chosen parameters.

### File 2: `cmaes_search.py` (~400 lines)

The Phase 1 search engine. Contains:
- CMA-ES implementation (or use `cmaes` pip package)
- Island model with migration
- Fitness evaluation (with early termination)
- Curriculum T_eval schedule
- Logging and checkpointing

Dependencies: jax, cmaes (pip), numpy

### File 3: `refine_candidates.py` (~200 lines)

Phase 2 gradient refinement. Reuses the PDE code from field_search.py
but seeds from parametric configs. Contains:
- Upscaling (trilinear interpolation)
- Backprop optimization loop
- Head-to-head tournament

### File 4: `validate_long.py` (~150 lines)

Phase 3 long-horizon validation. Contains:
- N=64 evolution with periodic diagnostics
- Convergence test (is E_pot still decaying or stabilized?)
- Final pass/fail determination
- Export to C-readable format


## Time and Resource Budget

| Phase | Resolution | Wall time | GPU memory | Evaluations |
|-------|-----------|-----------|------------|-------------|
| 1     | N=16      | ~1.0 hr   | < 1 GB     | ~32,000     |
| 2A    | N=32      | ~5 min    | < 2 GB     | 20          |
| 2B    | N=32      | ~10 min   | < 4 GB     | 8 x 200     |
| 2C    | N=48      | ~5 min    | < 4 GB     | 5           |
| 3     | N=64      | ~1.5 hr   | < 8 GB     | 3           |
| Setup | all       | ~30 min   | —          | calibration |
| **Total** |       | **~3.5 hr** | **< 8 GB** |            |

Leaves ~2.5 hours of slack in the 6-hour budget for reruns,
parameter adjustments, and the N=128 C simulation if we find a winner.


## What Success Looks Like

**Strong positive**: A structure at N=64 that retains >30% of E_pot
at T=300 with absorbing BC. This would be a genuine compact soliton.

**Weak positive**: A structure that decays but with a clear separation
between fast radiation (first 50 time units) and slow secular decay.
If the decay rate is < 1% per time unit after the transient, this is
a quasi-soliton that might survive with parameter tuning.

**Definitive negative**: All 32,000+ Phase 1 candidates, including
the best parameter regions found by CMA-ES, fail to retain >10%
of binding at N=32, T=50. Combined with the known failures of
hand-designed structures (knots, rings, crossed braids with absorbing
BC), this constitutes strong evidence that the PDE does not support
compact solitons at these parameters.


## Differences from Previous Attempts

| Previous approach | Problem | This design |
|-------------------|---------|-------------|
| Random voxels + backprop | E_pot starts at 0, no gradient signal | Parametric seeds with guaranteed binding |
| Single structure type | Limited to braids | 5 structure families, CMA-ES explores |
| Gradient only | Gets stuck in local optima | CMA-ES global + gradient local |
| No early termination | Wastes time on dead configs | Skip if E_pot > -0.5 |
| Fixed T_eval | Either too short (misses decay) or too expensive | Curriculum: 10 -> 20 -> 40 |
| N=32 only | Too expensive for population search | N=16 for exploration, upscale winners |
| Single seed or few seeds | Limited coverage | 320 per generation, 32,000 total |
| field_evolve.py (correlated random) | Still starts at E_pot ~ 0 | Parametric structure guarantees P >> 0 |
