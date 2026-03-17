# V24-DG: Rotating Triad (Three-Frequency + Gyroscopic Stabilization)

## Thesis

The v21 oscillon has all three fields in phase: φ₁=φ₂=φ₃=f(r)cos(ωt).
This is maximally aligned but also maximally degenerate — all radiation
channels are active.

Split the three frequencies: ω₁=ω₀+δ, ω₂=ω₀, ω₃=ω₀-δ. The triad
ROTATES in internal field space at rate δ. The triple product P acquires
beat frequencies at δ and 2δ instead of ω₀. If δ < m, these beats are
trapped below the mass gap → dramatically reduced radiation.

At δ=0: standard oscillon (degenerates, decays at dE/dt ~ -10⁻⁵)
At δ>0: rotating triad with continuous internal energy cycling

The three fields MUST all participate — removing any one destroys the
rotation. This is the self-alignment process operating continuously.

## Mathematical Setup

### Initial Condition

    φ₁(x, 0) = A·exp(-x²/2σ²) · cos(0)     = A·g(x)
    φ₂(x, 0) = A·exp(-x²/2σ²) · cos(0)     = A·g(x)
    φ₃(x, 0) = A·exp(-x²/2σ²) · cos(0)     = A·g(x)

    v₁(x, 0) = -ω₁ · A·g(x) · sin(0) = 0   ... wait, need phase offset

Better: initialize with phase offsets that produce the rotation:

    φ₁(x, 0) = A·g(x)
    φ₂(x, 0) = A·g(x) · cos(2π/3·0)  ... no, we need FREQUENCY difference

The frequency splitting emerges from giving the three fields different
INITIAL VELOCITIES:

    φ₁(x, 0) = A·g(x),     v₁(x, 0) = -(ω₀+δ)·A·g(x)·0 = 0
    φ₂(x, 0) = A·g(x),     v₂(x, 0) = 0
    φ₃(x, 0) = A·g(x),     v₃(x, 0) = 0

This won't work — identical initial conditions give identical evolution.

**Correct approach**: Two-step initialization.
1. Equilibrate the standard oscillon (t=5000, all fields equal)
2. At equilibrium, PERTURB the phases by shifting the velocities:
   v₁ → v₁ + δ·φ₁  (advance phase by δ·dt at each point)
   v₃ → v₃ - δ·φ₃  (retard phase by δ·dt)
   v₂ unchanged (reference)

This creates a small frequency splitting that may grow or decay. If the
rotating state is an ATTRACTOR, the splitting grows to a finite value.
If it's UNSTABLE, it collapses back to the degenerate state.

### Triple Product Analysis

With φ_a = f·cos(ω_a·t + θ_a):

P = f³ cos(ω₁t+θ₁) cos(ω₂t+θ₂) cos(ω₃t+θ₃)

Using product-to-sum formulas, P contains frequencies:
  ω₁+ω₂+ω₃ = 3ω₀  (high, above gap)
  ω₁+ω₂-ω₃ = ω₀+2δ  (near ω₀)
  ω₁-ω₂+ω₃ = ω₀  (at ω₀)
  -ω₁+ω₂+ω₃ = ω₀-2δ  (near ω₀)

So P still has components at ω₀ ± 2δ. The radiation is NOT eliminated
by frequency splitting alone — it's shifted by ±2δ.

For radiation suppression: need ω₀ + 2δ < m and ω₀ - 2δ > 0, i.e.,
δ < (m - ω₀)/2. With ω₀=0.87, m=1.0: δ < 0.065.

The splitting must be SMALL (δ < 0.065) to avoid opening new radiation
channels. The internal rotation rate is slow.

## What to Compute

### Phase 1: Equilibrate + Perturb

1. Equilibrate standard oscillon (μ=-20, κ=20, m=1.0) for t=5000
2. Save profile. Record ω₀, E₀, A₀.
3. Apply phase perturbation: v₁ += δ·φ₁, v₃ -= δ·φ₃
4. Scan δ ∈ {0.0, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20}

### Phase 2: Evolve and measure

5. For each δ: evolve for t=10000
6. Measure:
   a. Energy E(t) and energy loss rate dE/dt
   b. Per-field frequencies ω₁(t), ω₂(t), ω₃(t) from DFT of each field
   c. Frequency splitting Δω(t) = ω₁ - ω₃
   d. Does the splitting GROW, DECAY, or remain constant?
   e. Peak amplitudes per field
   f. Oscillon lifetime
7. Compare dE/dt at each δ with the δ=0 baseline

### Key Question

Does the rotating triad (δ>0) have a LOWER radiation rate than the
symmetric state (δ=0)?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`
- v24 boost code: `/home/d/code/scp/v24/src/boost1d.c` (profile saving)

## Output

- `src/rotating1d.c` — code
- `data/rotating_d{delta}_ts.tsv` — time series per δ
- `data/rotating_summary.tsv` — summary table
- `RESULTS.md` — analysis

## Parameters

μ=-20, κ=20, m=1.0
Nx=4000, xmax=100, t_equil=5000, t_run=10000
δ scan: {0.0, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20}

Compile: `gcc -O3 -Wall -o rotating1d src/rotating1d.c -lm`
