# V26-DynAC: Propagating Massless Braid

## Thesis

Combine axial propagation (A) with masslessness (C). The helical wave
propagates at c (massless) along the braid axis. The "mass" is purely
from the trapped wave energy confined by the triple product.

This is the most radical option: no mass parameter, maximum speed, all
energy from dynamics.

## Initialization (3D, N=128, L=20)

    φ_a(x,0) = A(r_⊥) · cos(kz + 2πa/3)
    v_a(x,0) = k · A(r_⊥) · sin(kz + 2πa/3)    [ω=k, m=0]

## Constraints

- V = c (propagating at light speed)
- T > 0: from ω = k oscillation
- δV > 0: v_g = c > 0

## Parameters

μ=-20, κ=20, m=0.0, k=2π/20
N=128, L=20, periodic BC in z, t=500

## Output

`src/dynAC.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynAC src/dynAC.c -lm`
