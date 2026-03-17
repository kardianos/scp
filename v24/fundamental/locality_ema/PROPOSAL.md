# Locality Option 1: Time-Averaged Wave Equation (EMA Source)

## Thesis

Replace the Poisson equation (instantaneous) with a wave equation sourced
by the time-averaged energy density (causal, static limit matches Poisson):

    □Φ = ∂²Φ/∂t² - c²∂²Φ/∂x² = α · ⟨ρ⟩

where ⟨ρ⟩ is an exponential moving average (EMA) of ρ with time constant τ.

The EMA filters out the oscillating 2ω component of ρ, leaving only the
static (DC) part as the source. Changes to ⟨ρ⟩ (from oscillon motion)
propagate at c — causally.

## Free Parameter: τ (EMA timescale)

τ must be tuned: too small → oscillation leaks through → radiation.
Too large → sluggish response → loses track of oscillon motion.

**Optimal τ**: approximately one breathing period T = 2π/ω ≈ 7.2 (at ω=0.87).

**Backprop tuning**: scan τ = {2, 4, 7, 10, 15, 20, 30} and measure:
1. Oscillation amplitude of Φ (should be minimized)
2. Response time to oscillon motion (should be fast)
3. Static Φ(0) (should match Poisson: Φ_Poisson ≈ -0.006 at α=1e-4)
4. Oscillon survival (fc > 0.9)

Select τ* that minimizes oscillation while maintaining quick response.

## Method

### Phase 1: Single oscillon τ scan

1. Equilibrate oscillon t=5000 (no gravity)
2. Turn on □Φ = α⟨ρ⟩_τ with α=-1e-4 (attractive)
3. Evolve t=5000 for each τ in {2, 4, 7, 10, 15, 20, 30}
4. Measure: Φ(0) mean and oscillation amplitude, fc, ω
5. Select τ* = value with smallest Φ oscillation and Φ(0) ≈ Φ_Poisson

### Phase 2: Causality test

6. At τ*: start with equilibrated oscillon + Φ
7. At t=5000: suddenly BOOST the oscillon (give it velocity v=0.01)
8. Track: how fast does Φ respond at distance x=50?
9. For Poisson: Φ(x=50) responds INSTANTLY when the oscillon moves
10. For wave eq: Φ(x=50) responds after delay Δt = 50/c = 50
11. Measure the delay. Is it consistent with c?

### Phase 3: Two-oscillon with causal gravity

12. Two oscillons at D=40, τ=τ*, α=-1e-4
13. Does Φ mediate attraction? At what rate?
14. Compare with Combo 1+2+5+7 (Poisson): is the force the same
    but with a time delay?

## Implementation Notes

The EMA: ⟨ρ⟩_new = (1 - dt/τ)·⟨ρ⟩_old + (dt/τ)·ρ_current

The wave equation for Φ: use Velocity Verlet alongside the scalar fields.
Φ has its own velocity v_Φ and acceleration a_Φ = c²∂²Φ/∂x² + α⟨ρ⟩.
NO mass term (Φ is massless). Absorbing BC for Φ at domain boundaries.

Backreaction on scalar fields: m²_eff = m²(1+2Φ), c²_eff = 1+4Φ
(same as Test F).

## Reference Code

- v24/fundamental/testF_selfref/src/selfref.c (EMA + metric correction)
- v24/fundamental/combo_1257/src/combo1257.c (lattice + gravity)

## Output

- `src/locality_ema.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, α=-1e-4, c_Φ=1.0 (Φ speed = field speed)
τ scan: {2, 4, 7, 10, 15, 20, 30}
Nx=8000, xmax=200, t_equil=5000, t_test=5000

Compile: `gcc -O3 -Wall -o locality_ema src/locality_ema.c -lm`
