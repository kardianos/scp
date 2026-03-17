# Locality Option 3: Telegraph Equation (Damped Wave)

## Thesis

Interpolate between wave (causal) and Poisson (instantaneous) using the
telegraph equation:

    ∂²Φ/∂t² + γ·∂Φ/∂t = c²·∂²Φ/∂x² + α·ρ

γ controls the damping:
- γ=0: pure wave (oscillating Φ, radiation problem)
- γ→∞: overdamped → Poisson (instantaneous)
- γ ~ c/L: causal propagation with rapid relaxation

## Free Parameter: γ (damping rate)

γ must be tuned. The physically motivated value: γ ~ c/L where L is the
oscillon size (~5 code units). So γ ~ 0.2.

**Backprop tuning**: scan γ = {0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 20.0}

At each γ, measure:
1. Φ oscillation amplitude (should decrease with γ)
2. Propagation delay at distance x=50 (should be ~50/c at small γ,
   instantaneous at large γ)
3. Static Φ(0) (should match Poisson at all γ)
4. Oscillon survival (fc > 0.9)

Select γ* that gives CAUSAL propagation (measurable delay) with
SMALL oscillation (Φ_osc/Φ_static < 0.1).

## Method

### Phase 1: γ scan (single oscillon)

1. Equilibrate t=5000
2. Turn on telegraph Φ with α=-1e-4, each γ value
3. Evolve t=5000
4. Measure Φ(0) static value and oscillation amplitude

### Phase 2: Causality test at γ*

5. Boost oscillon at t=5000
6. Measure Φ response delay at x=50
7. Is delay = 50/c (causal) or 0 (instantaneous)?

### Phase 3: Compare with Options 1 and 4

8. At the same α=-1e-4: compare Φ profiles from telegraph, EMA-wave,
   and Poisson. Which is closest to the retarded solution?

## Implementation

The telegraph equation adds a first-derivative damping to Φ's velocity:

    v_Φ_new = v_Φ_old + dt·(c²·∂²Φ/∂x² + α·ρ) - γ·dt·v_Φ_old

This is just Velocity Verlet with velocity damping on the Φ sector.
The scalar fields φ_a are NOT damped (only Φ).

## Reference Code

- v24/fundamental/testF_selfref/src/selfref.c (gravity implementation)

## Output

- `src/locality_telegraph.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, α=-1e-4
γ scan: {0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 20.0}
Nx=8000, xmax=200, t_equil=5000, t_test=5000

Compile: `gcc -O3 -Wall -o locality_telegraph src/locality_telegraph.c -lm`
