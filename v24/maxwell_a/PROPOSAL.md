# V24-MA: Complex Scalar + U(1) Gauge Field (Minimal EM)

## Background

The V24 180° investigation revealed a Z₂ symmetry (P → -P) in the oscillon.
The simplest way to introduce electromagnetism is to promote two of the three
real fields into a complex field and couple it to a U(1) gauge field (photon).

## Setup

Replace φ₁, φ₂ with one complex field Ψ = φ₁ + iφ₂. Keep φ₃ real.

### Lagrangian

    L = |D_t Ψ|² - |D_x Ψ|² - m²|Ψ|² + ½(∂_t φ₃)² - ½(∂_x φ₃)²
        - ½m²φ₃² - V_coupling - ¼F_{μν}F^{μν}

where D_μ = ∂_μ - ieA_μ is the covariant derivative and F_{μν} = ∂_μA_ν - ∂_νA_μ.

### Coupling Term

The original triple product P = φ₁φ₂φ₃ = Re(Ψ)·Im(Ψ)·φ₃ is NOT gauge
invariant (Ψ → e^{iα}Ψ changes Re and Im separately).

Use a gauge-invariant alternative:

    V_coupling = (μ/2)(|Ψ|²φ₃²)² / (1 + κ(|Ψ|²φ₃²)²)

This replaces P² = (φ₁φ₂φ₃)² with (|Ψ|²φ₃²)² = (φ₁²+φ₂²)²φ₃⁴.
It's gauge-invariant since |Ψ|² is invariant. The saturation is the same form.

Note: this changes the binding from triple-product (P = φ₁φ₂φ₃) to
biquadratic (|Ψ|²φ₃²). The oscillon structure will be different.

### Equations of Motion (1D, temporal gauge A₀=0)

    ∂²Ψ/∂t² = (∂_x - ieA)²Ψ - m²Ψ - ∂V/∂Ψ*
    ∂²φ₃/∂t² = ∂²φ₃/∂x² - m²φ₃ - ∂V/∂φ₃
    ∂²A/∂t² = ∂²A/∂x² - e·Im(Ψ*D_xΨ)    (Maxwell: □A = j)

The current: j = e·Im(Ψ*(∂_x - ieA)Ψ) = e(φ₁∂_xφ₂ - φ₂∂_xφ₁ - eA|Ψ|²)

In component form (Ψ = φ₁ + iφ₂):
    ∂²φ₁/∂t² = ∂²φ₁/∂x² + e(∂_xA)φ₂ + 2eA(∂_xφ₂) + e²A²φ₁ - m²φ₁ - ∂V/∂φ₁
    ∂²φ₂/∂t² = ∂²φ₂/∂x² - e(∂_xA)φ₁ - 2eA(∂_xφ₁) + e²A²φ₂ - m²φ₂ - ∂V/∂φ₂
    ∂²φ₃/∂t² = ∂²φ₃/∂x² - m²φ₃ - ∂V/∂φ₃
    ∂²A/∂t² = ∂²A/∂x² - e(φ₁∂_xφ₂ - φ₂∂_xφ₁) + e²A(φ₁²+φ₂²)

## What to Compute

### Phase 1: Oscillon with gauge coupling

1. Initialize: Ψ = A_init·g(x) (real, so φ₁=A·g, φ₂=0), φ₃ = A·g, A_μ=0
2. Set e=0 first (control, should reproduce standard oscillon with |Ψ|²φ₃² coupling)
3. Then scan e ∈ {0.0, 0.1, 0.3, 0.5, 1.0}
4. Evolve t=10000. Measure: E_total, E_EM (in A field), charge Q, ω

### Phase 2: Two-oscillon EM interaction

5. Place two oscillons at separation D=20
6. Give them OPPOSITE charges (one has Ψ = φ₁+iφ₂, other has Ψ = φ₁-iφ₂)
7. Measure: do they attract (Coulomb)? Force vs D?
8. The EM field A should develop a 1/|x| profile between them (in 1D: linear)

### Phase 3: Check conserved charge

9. Monitor Q(t) = e∫Im(Ψ*∂_tΨ) dx. Should be conserved.
10. Does the oscillon carry a definite charge?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/maxwell_a.c` — solver
- `data/` — TSV output files
- `RESULTS.md` — analysis

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
e scan: {0.0, 0.1, 0.3, 0.5, 1.0}
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o maxwell_a src/maxwell_a.c -lm`
