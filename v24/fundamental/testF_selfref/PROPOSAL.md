# Test F: Self-Consistent Metric from Field Configuration

## Thesis

Compute the metric g_μν implied by the oscillon's stress-energy tensor
via Einstein's equation. Then evolve the fields ON this metric. Iterate
to self-consistency. Does a unique metric emerge?

## Method

### Step 1: Compute metric from T^μν

In 1+1D, Einstein's equation reduces to:
    R = 8πG · T (trace of stress-energy)

For a static, spherically symmetric (1D: reflection-symmetric) source:
    g_00 = -(1 + 2Φ(x)/c²)
    g_11 = (1 - 2Φ(x)/c²) or (1 + 2Ψ(x)/c²)

where Φ is the Newtonian potential:
    ∂²Φ/∂x² = 4πG · ρ(x)

In 1D: Φ(x) = -G ∫|x-x'|·ρ(x') dx' (Green's function is |x|)

### Step 2: Evolve fields on curved metric

The Klein-Gordon equation on curved background:
    (1/√(-g)) ∂_μ(√(-g) g^μν ∂_ν φ) - m²φ = 0

In 1+1D with the metric above:
    -(1+2Φ)⁻¹ ∂²φ/∂t² + (1-2Φ)⁻¹ ∂²φ/∂x² + (corrections from ∂Φ) - m²φ = 0

For weak gravity (Φ << 1):
    ∂²φ/∂t² ≈ (1+4Φ)∂²φ/∂x² + (2∂Φ/∂x)∂φ/∂x - m²(1+2Φ)φ

### Step 3: Iterate

1. Evolve fields on flat metric → compute T^μν → compute Φ
2. Evolve fields on Φ-corrected metric → compute new T^μν → new Φ
3. Repeat until Φ converges

### Practical Implementation

Don't iterate — instead, add Φ as a dynamical field:
    ∂²Φ/∂x² = α · ρ(x)  (Poisson equation, solved at each timestep)

The fields evolve with effective mass m²_eff(x) = m²(1+2Φ) and effective
wave speed c²_eff = 1+4Φ. This is the V22 approach (gravity backreaction)
but with the CORRECT metric coupling.

The coupling constant α = 4πG/c⁴ is a free parameter. Scan to find
self-consistent behavior.

## What to Test

1. Solve Poisson for Φ from the oscillon's ρ(x). Measure Φ(0).
2. Evolve fields with Φ-corrected dynamics. Does the oscillon change?
3. Iterate: does Φ converge to a fixed profile?
4. At self-consistency: what is Φ(0)/c²? This is the gravitational
   redshift at the oscillon center.

## Reference Code

- v21/src/triad1d.c
- v22/src/two_oscillon.c (gravity backreaction)

## Output

- `src/selfref.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
α scan: {0.001, 0.01, 0.1}
Nx=4000, xmax=100, t=5000

Compile: `gcc -O3 -Wall -o selfref src/selfref.c -lm`
