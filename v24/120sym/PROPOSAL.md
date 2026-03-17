# V24-SYM: Symmetric Penalty — Force 120° via (φ₁+φ₂+φ₃)²

## Thesis

Add V_pen = λ_Q(φ₁+φ₂+φ₃)² which ONLY affects the symmetric mode:
m²_S = m² + 6λ_Q. Antisymmetric modes unchanged at m². Combined with
the triple product, this forces the system away from 0° toward 120°.

Unlike the pairwise coupling (V24-PW), this doesn't lower the
antisymmetric mass. The 120° oscillon still has ω close to m.
But the 0° state is strongly disfavored, so 120° is the only option.

## EOM

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - 2λ_Q(φ₁+φ₂+φ₃) - ∂V_triple/∂φ_a

## Method

1. Scan λ_Q ∈ {0.0, 0.1, 0.5, 1.0, 2.0, 5.0} (units of m²)
2. Initialize both 0° and 120° at each λ_Q
3. Evolve t=10000
4. Measure: does the 0° state spontaneously TRANSITION to 120°?
5. Measure: is the 120° state more stable (lower dE/dt)?

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Nx=4000, xmax=100

## Output

- `src/sym120.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -Wall -o sym120 src/sym120.c -lm`
