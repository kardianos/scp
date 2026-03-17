# V24-E: Dissipative Soliton (Gain-Loss Balance)

## Thesis

The v21 oscillon is a CONSERVATIVE system — energy is either in the
oscillon or radiated to infinity. There is no energy input. The oscillon
slowly dies as it radiates harmonics.

In nature, many stable structures are DISSIPATIVE — they exist as a
balance between energy input and dissipation. Laser pulses, nerve impulses,
flame fronts, convection cells. These structures REQUIRE energy flow to
exist and CANNOT exist in equilibrium (T=0, v=0).

Add a nonlinear gain term that feeds energy to ALIGNED field configurations:

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a - γ·∂φ_a/∂t + G·S_a

where:
  γ = dissipation coefficient (removes energy uniformly)
  G = gain coefficient
  S_a = gain source term (feeds energy to aligned states)

The gain source S_a should reward configurations where P = φ₁φ₂φ₃ is large:

  S_a = P² · φ_a  (gain proportional to alignment squared, times the field)

This means: where all three fields are large and aligned, energy flows IN.
Where fields are small or misaligned, energy flows OUT (from the γ term).
The oscillon exists at the balance point.

## Mathematical Setup

### Modified EOM

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a - γ·(∂φ_a/∂t) + G·P²·φ_a

The gain term G·P²·φ_a acts as a NEGATIVE mass squared when P is large:
  m²_eff = m² - G·P² = m² - G·f⁶·cos⁶(ωt)

Time-averaged: m²_eff = m² - G·⟨f⁶cos⁶⟩ = m² - G·(5/16)·f⁶

For the oscillon core (f≈0.5): ⟨P²⟩ ≈ (5/16)·0.5⁶ ≈ 0.005.
Gain contribution: G·0.005. Needs G ~ 200 to match m² = 1.0.

### Energy Balance

Power input from gain: P_in = G · ∫ P² φ_a (∂φ_a/∂t) dx  (summed over a)
Power lost to damping: P_out = γ · ∫ (∂φ_a/∂t)² dx

At steady state: P_in = P_out. This determines the oscillon amplitude as
a function of G/γ.

### Limit Cycle Behavior

The dissipative oscillon is a LIMIT CYCLE, not a fixed point. It breathes
at a frequency and amplitude determined by the gain-loss balance, not by
initial conditions. Any perturbation that changes the amplitude is corrected:
- Amplitude too high → dissipation exceeds gain → amplitude decreases
- Amplitude too low → gain exceeds dissipation → amplitude increases

This is SELF-CORRECTING — the defining property of a dissipative soliton.

## What to Compute

### Phase 1: Find the Gain-Loss Balance

1. Start with the equilibrated v21 oscillon (μ=-20, κ=20, m=1.0).
2. Turn on damping γ and gain G simultaneously.
3. Scan γ ∈ {0.001, 0.01, 0.1} and for each γ, scan G to find the
   balance point where dE/dt ≈ 0.

4. The balance condition: G_balance = γ · ⟨(∂φ/∂t)²⟩ / ⟨P² φ (∂φ/∂t)⟩

### Phase 2: Stability of the Dissipative Oscillon

5. At the balanced (γ, G): perturb the amplitude by ±20%.
6. Does the oscillon return to the original amplitude (limit cycle)?
7. How fast? (Measure the relaxation rate.)

### Phase 3: From Nothing

8. Start from NOISE (no oscillon, just random small fields).
9. With gain turned on: does an oscillon spontaneously FORM from noise?
10. This would demonstrate that the dissipative oscillon is a true attractor
    — it forms from any initial condition, not just carefully prepared ones.

### Phase 4: Comparison with Conservative

11. Compare the dissipative oscillon (γ>0, G>0) with the conservative
    one (γ=0, G=0). Which lives longer? Which is more robust against
    perturbation?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/dissipative1d.c` — code
- `data/dissipative_g{gamma}_G{gain}_ts.tsv` — time series
- `data/dissipative_summary.tsv` — balance points
- `data/dissipative_fromNoise_ts.tsv` — spontaneous formation test
- `RESULTS.md` — analysis

## Parameters

μ=-20, κ=20, m=1.0
γ scan: {0.001, 0.01, 0.1}
G scan: adaptive to find balance
Nx=4000, xmax=100
t_equil=5000 (conservative), t_dissipative=10000

Compile: `gcc -O3 -Wall -o dissipative1d src/dissipative1d.c -lm`
