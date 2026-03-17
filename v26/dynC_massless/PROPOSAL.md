# V26-DynC: Massless Propagating Braid (Emergent Mass from Confinement)

## Thesis

Set m=0 but initialize a PROPAGATING helical wave. With m=0: ω=k, v_g=c
(propagates at light speed). The "mass" of the soliton is the trapped
wave energy confined by the triple product.

V26 mode 2 (m=0, static) survived at fc=0.37 but slowly dispersed.
With propagation: the continuous motion might maintain the structure.

## Initialization (3D, N=128, L=20)

    φ_a(x,0) = A(r_⊥) · cos(kz + 2πa/3)
    v_a(x,0) = +k · A(r_⊥) · sin(kz + 2πa/3)    [ω = k when m=0]

With m=0: ω = k = 2π/20 = 0.314, v_g = 1.0 (light speed).

## Constraints

- V > 0: propagates at c (maximum speed, always moving)
- T > 0: oscillation at ω=0.314
- δV > 0: velocity = c (fixed by masslessness)

## Lagrangian

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - (μ/2)P²/(1+κP²)

NO mass term. The triple product alone provides binding.

## Key Question

Can a PROPAGATING massless helical wave be confined by the triple product?
The V26 Phase 4 showed m=0 + static → collapse (no confinement).
But m=0 + propagating might be different: the wave carries momentum along z,
and the triple product creates a transverse potential well that confines
the wave to the tube. The wave circulates forever (periodic BC), maintaining
the braid structure through continuous motion.

## What to Measure

1. Does the braid survive (fc > 0.3)?
2. Does |P| persist (unlike static m=0 which collapsed)?
3. Is the effective mass M = E/c² well-defined and stable?
4. l=2 content — V26 mode 2 (static, m=0) had 6.76%

## Parameters

μ=-20, κ=20, m=0.0, A₀=0.8, R_tube=3.0
k=2π/20
N=128, L=20, periodic BC in z, t=500

## Output

`src/dynC.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynC src/dynC.c -lm`
