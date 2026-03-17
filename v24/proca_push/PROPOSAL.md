# V24-P1: Push λ to Extreme — How Close to m² Can We Get?

## Thesis

V24-ME found oscillons survive at all λ < m² tested (up to 0.995). The
symmetric mode mass m²_S = m² + 2λ INCREASES with λ, stabilizing the
oscillon. The antisymmetric Proca mass m_A = √(m²-λ) → 0 as λ → m².

Test: push λ to {0.999, 0.9995, 0.9999, 0.99999, 0.999999}. If the
oscillon survives, range scales as 1/√(m²-λ) → 32, 45, 100, 316, 1000.

## EOM

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b+φ_c) - μP·dP/dφ_a/(1+κP²)²

Same as V24-ME but with λ much closer to m² = 1.0.

## Method

1. For each λ in {0.999, 0.9995, 0.9999, 0.99995, 0.99999, 0.999999}:
   a. Initialize symmetric oscillon (A=0.8, σ=3.0, all fields equal)
   b. Evolve t=10000
   c. Measure: fc, ω, E, dE/dt, peak amplitude
   d. Check vacuum stability: does any field grow at large |x|?

2. At the most extreme stable λ: run two oscillons at D=50, 80, 120
   to test if the long-range Proca mediates attraction at these distances.

3. Use a LARGE grid: Nx=16000, xmax=500 (range 316 at λ=0.99999 needs
   domain > 600 for proper measurement).

## Key Diagnostics

- Does fc remain > 0.99 at each λ?
- Is there any sign of vacuum instability (field growth at |x| > 50)?
- At what λ does the antisymmetric perturbation propagation distance
  exceed 100 code units?

## Reference Code

- v24/maxwell_e: `/home/d/code/scp/v24/maxwell_e/src/maxwell_e.c`

## Output

- `src/push1d.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
λ scan: {0.999, 0.9995, 0.9999, 0.99995, 0.99999, 0.999999}
Nx=16000, xmax=500, tfinal=10000

Compile: `gcc -O3 -Wall -o push1d src/push1d.c -lm`
