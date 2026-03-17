# V24-PW: Pairwise Coupling — 120° Phase Binding

## Thesis

Add V_pw = λ(φ₁φ₂ + φ₂φ₃ + φ₃φ₁) to the Lagrangian. This splits the
mass spectrum: m²_anti = m² - λ (lighter, 120°-compatible) and
m²_sym = m² + 2λ (heavier, 0°-penalized). Scan λ from 0 to 0.95m²
to find where the 120° oscillon becomes stable.

## EOM

    ∂²φ_a/∂t² = ∂²φ_a/∂x² - m²φ_a - λ(φ_b + φ_c) - ∂V_triple/∂φ_a

where (b,c) are the other two fields, V_triple = (μ/2)P²/(1+κP²).

## Method

1. Scan λ ∈ {0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.90, 0.95} (units of m²)
2. For each λ: initialize with 120° phases:
   φ₁ = A·g(x), v₁ = 0
   φ₂ = A·g(x)·cos(2π/3), v₂ = -ω·A·g(x)·sin(2π/3)
   φ₃ = A·g(x)·cos(4π/3), v₃ = -ω·A·g(x)·sin(4π/3)
   where ω = √(m² - λ) (the antisymmetric base frequency)
3. Evolve t=10000. Measure: phase survival, ω, dE/dt, fc
4. Also run 0° control at each λ for comparison
5. Find: at what λ does 120° become more stable than 0°?

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Nx=4000, xmax=100

## Output

- `src/pw120.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -Wall -o pw120 src/pw120.c -lm`
