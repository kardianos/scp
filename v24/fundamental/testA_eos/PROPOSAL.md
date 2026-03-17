# Test A: Vacuum Equation of State

## Thesis

Compute the full stress-energy tensor T^μν of the oscillon. Extract the
equation of state P(ρ). Check if P = -ρ holds in any regime (required for
Lorentz-invariant vacuum). Measure pressure anisotropy and shear stress.

This is PURE DIAGNOSTIC — no new dynamics, just measurement of existing
oscillon properties that have never been computed.

## What to Compute

1. Equilibrate the standard oscillon (μ=-20, κ=20, m=1.0) for t=5000
2. At each grid point, compute the stress-energy tensor components:
   - T^00 = energy density ρ = Σ_a [½v_a² + ½(∂_x φ_a)² + ½m²φ_a²] + V(P)
   - T^01 = momentum density = -Σ_a v_a · ∂_x φ_a
   - T^11 = stress (pressure) = Σ_a [½v_a² + ½(∂_x φ_a)² - ½m²φ_a²] - V(P)
   Wait — in 1D the full T^μν is:
   T^00 = Σ_a [½(∂_t φ_a)² + ½(∂_x φ_a)² + ½m²φ_a²] + V
   T^11 = Σ_a [½(∂_t φ_a)² + ½(∂_x φ_a)² - ½m²φ_a²] - V
   T^01 = -Σ_a (∂_t φ_a)(∂_x φ_a)
3. Time-average over several breathing cycles
4. Compute: P(x) = ⟨T^11⟩, ρ(x) = ⟨T^00⟩, momentum flux ⟨T^01⟩
5. Plot P vs ρ at each grid point. Is there a universal curve P(ρ)?
6. Check: in the vacuum (far from oscillon), P/ρ → ?
7. Check: inside the oscillon core, P/ρ → ?
8. Compute the trace T^μ_μ = -ρ + P (in 1+1D). Is it zero (conformal)?

Also do with pairwise coupling λ=0.99.

## Reference Code

- v21/src/triad1d.c

## Output

- `src/eos.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, Nx=4000, xmax=100, t=10000
Also: λ=0.99 variant

Compile: `gcc -O3 -Wall -o eos src/eos.c -lm`
