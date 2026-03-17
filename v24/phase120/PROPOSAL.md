# V24-B: 120° Phase-Separated Oscillon

## Thesis

The symmetric oscillon (φ₁=φ₂=φ₃, all in phase) has P = f³cos³(ωt),
which radiates at ω, 2ω, 3ω through the nonlinear coupling.

The 120° phase-separated oscillon:
    φ₁ = f·cos(ωt)
    φ₂ = f·cos(ωt + 2π/3)
    φ₃ = f·cos(ωt + 4π/3)

has P = (f³/4)·cos(3ωt). The fundamental and second harmonic VANISH in the
triple product — only the third harmonic remains.

If 3ω < m (i.e., ω < m/3 ≈ 0.33), the oscillon radiates NOTHING through
the coupling term. This would be a perfectly stable bound state.

The challenge: ω < 0.33 requires much deeper binding than v21 (ω ≈ 0.87).
But even at v21 parameters, the 120° state eliminates two radiation channels,
potentially extending the lifetime dramatically.

## Mathematical Setup

### Triple Product Identity

cos(θ)·cos(θ+2π/3)·cos(θ+4π/3) = (1/4)·cos(3θ)

Proof: use cos(A+B) = cosAcosB - sinAsinB to expand, terms cancel pairwise.

So P_120 = (f³/4)cos(3ωt) vs P_sym = f³cos³(ωt) = f³[(3/4)cos(ωt) + (1/4)cos(3ωt)]

The 120° state has P_max = f³/4 (vs f³ for symmetric). Binding energy is
WEAKER (P² is 16× smaller). The oscillon must compensate with larger
amplitude or stronger coupling.

### Stability of the 120° Configuration

The symmetric state (0° phase) MAXIMIZES P = f³, so it minimizes V(P²)
when μ < 0 (attractive). The 120° state has P = f³/4, which gives
HIGHER potential energy (weaker binding).

The 120° state is a SADDLE POINT: stable against perturbations that
maintain the 120° pattern, but unstable against collapse to the 0° state.

The question: how fast does the 120° state collapse to 0°? If the collapse
timescale is LONGER than the radiation timescale of the 0° state, then the
120° state is effectively more stable.

## What to Compute

### Phase 1: Initialize 120° State

1. Start with the equilibrated symmetric oscillon profile f_eq(x)
2. At a moment when the symmetric state is at phase θ=0 (maximum):
   - Set φ₁ = f_eq · cos(0) = f_eq
   - Set φ₂ = f_eq · cos(2π/3) = -f_eq/2
   - Set φ₃ = f_eq · cos(4π/3) = -f_eq/2
   - Set v₁ = 0
   - Set v₂ = -ω·f_eq·sin(2π/3) = -ω·f_eq·√3/2
   - Set v₃ = -ω·f_eq·sin(4π/3) = +ω·f_eq·√3/2

   where ω is the breathing frequency from equilibration.

### Phase 2: Evolve and Compare

3. Evolve the 120° state for t=10000
4. Also evolve the 0° (symmetric) control for t=10000
5. Measure for both:
   a. Energy E(t) and loss rate dE/dt
   b. Phase differences Δθ₁₂(t), Δθ₂₃(t) — do the phases stay at 120°?
   c. Collapse time: when do phases snap to 0° (if ever)?
   d. Per-field amplitudes A₁(t), A₂(t), A₃(t)
   e. DFT of P(t) — what frequencies appear?

### Phase 3: Strong Binding Regime

6. Try parameters where ω < m/3:
   - Increase |μ| to deepen the binding: try μ=-60, -100, -200
   - Or decrease m: try m=0.5 with μ=-20
   - Find ω at each parameter set
   - If ω < m/3: the 120° state should be PERFECTLY stable (no radiation)

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`
- v23 hessian: `/home/d/code/scp/v23/hessian/src/hessian.c` (profile saving)

## Output

- `src/phase120_1d.c` — code
- `data/phase120_ts.tsv` — time series for 120° state
- `data/phase120_control_ts.tsv` — time series for 0° control
- `data/phase120_summary.tsv` — comparison table
- `RESULTS.md` — analysis

## Parameters

Base: μ=-20, κ=20, m=1.0
Strong binding: μ=-60, κ=20, m=1.0 AND μ=-20, κ=20, m=0.5
Nx=4000, xmax=100, t_equil=5000, t_run=10000

Compile: `gcc -O3 -Wall -o phase120_1d src/phase120_1d.c -lm`
