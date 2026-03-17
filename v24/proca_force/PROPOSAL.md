# V24-S1: Clean Force Measurement with Pairwise Coupling

## Thesis

V24-ME and P1 showed the Proca mediator exists and the oscillon survives at
extreme λ. But the actual FORCE between two oscillons through this channel
has never been cleanly measured — radiation artifacts from Gaussian
initialization overwhelmed the signal in every prior attempt (V23-C, V23-F,
V24-P1 Phase 2).

This test uses the V23-F methodology (pre-equilibrated profiles) applied
to the pairwise-coupled system. This is the most important missing
measurement in the entire V24 program.

## Method

### Step 1: Equilibrate with pairwise coupling

1. Initialize symmetric oscillon (A=0.8, σ=3.0) with pairwise coupling λ
2. Evolve for t_equil=10000 on Nx=4000, xmax=80
3. Save the FULL state (φ_a, v_a) at a breathing peak
4. Record: ω, E, peak amplitude, profile shape
5. Do this for λ = {0.0, 0.5, 0.9, 0.99, 0.999}

### Step 2: Two-oscillon force measurement

For each λ from Step 1, at separations D = {15, 20, 30, 40, 60, 80, 100}:

1. Place two equilibrated profiles at ±D/2
2. Evolve for t_run=2000 on Nx=16000, xmax=400
3. Track separation via energy-weighted centroids (split at x=0)
4. Measure INITIAL acceleration: fit sep(t) to a+bt+ct² over t∈[0,500]
   (short window before significant displacement)
5. Force F = M·c (where M is oscillon mass from equilibration)

### Step 3: Extract interaction range

For each λ: fit F(D) to Yukawa form F = F₀·exp(-D/λ_range)
Compare λ_range with predicted 1/m_A = 1/√(m²-λ)

### Step 4: Verify scaling

Does F₀ depend on λ? Does λ_range track 1/m_A as predicted?
At λ=0.999: predict range = 32. Is F detectable at D=60?

## Key Diagnostics

The CRITICAL measurement: at λ=0.99 (range 10), is the force at D=30
significantly larger than at λ=0 (range ~1)? This tests whether the
Proca mediator actually enhances the interaction range.

## Reference Code

- v23/critical2: `/home/d/code/scp/v23/critical2/src/critical2.c`
  (equilibrated profile + interaction measurement)
- v24/maxwell_e: `/home/d/code/scp/v24/maxwell_e/src/maxwell_e.c`
  (pairwise coupling implementation)

## Output

- `src/proca_force.c` — equilibration + force measurement
- `data/profile_lam{L}.dat` — equilibrated profiles
- `data/force_lam{L}_D{D}_ts.tsv` — separation vs time
- `data/force_summary.tsv` — F(D) table for each λ
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
λ scan: {0.0, 0.5, 0.9, 0.99, 0.999}
D scan: {15, 20, 30, 40, 60, 80, 100}
Equil: Nx=4000, xmax=80, t_equil=10000
Interact: Nx=16000, xmax=400, t_run=2000

Compile: `gcc -O3 -Wall -o proca_force src/proca_force.c -lm`
