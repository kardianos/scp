# V26-DynA: Propagating Helical Wave Braid

## Thesis

Replace the STATIC braid initialization with a PROPAGATING helical wave.
The field pattern travels along the braid axis at group velocity v_g = k/ω.
With periodic BC: continuous circulation. The braid is never static.

## Initialization (3D, N=128, L=20)

    φ_a(x,0) = A(r_⊥) · cos(kz + 2πa/3)
    v_a(x,0) = +ω · A(r_⊥) · sin(kz + 2πa/3)

where k = 2π/L, ω = √(k² + m²), A(r_⊥) = A₀·exp(-r_⊥²/(2R²)), r_⊥ = √(x²+y²).

With m=1.0, L=20: k=0.314, ω=1.048, v_g=k/ω=0.30c.

## Constraints

- V > 0: the helical wave propagates at v_g = 0.30c (always moving)
- T > 0: the oscillation at ω=1.048 provides thermal content
- δV > 0: the propagation speed is nonzero and maintained by the initial velocity

## Lagrangian

Same as V26 baseline (simplest, which gave 41.5% l=2):

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - ½m²Σφ_a² - (μ/2)P²/(1+κP²)

NO extra couplings (κ_T=0, κ_S=0, λ_pw=0). The 41.5% l=2 came from this.

## What to Measure

1. Does fc INCREASE compared to static braid (V26 Phase 4: fc=0.21)?
2. Does |P| persist or decay (static: |P| decayed to 0.001)?
3. Does the braid maintain its twist (track the phase structure at t=500)?
4. Is it non-breathing (DFT of ρ_center)?
5. l=2 content on shell at R=8 — still ~41% or better?
6. Propagation: does the pattern move at v_g as expected?

## Parameters

μ=-20, κ=20, m=1.0, A₀=0.8, R_tube=3.0
k=2π/20, ω=√(k²+1)=1.048
N=128, L=20, periodic BC in z, absorbing in x,y, t=500

## Reference

v26/src/v26.c (base 3D braided soliton code)

## Output

`src/dynA.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynA src/dynA.c -lm`
