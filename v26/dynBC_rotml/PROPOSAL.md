# V26-DynBC: Rotating Massless Braid

## Thesis

Combine rotation (B) with masslessness (C). The braid rotates around
the z-axis with no mass term. Angular momentum from rotation, energy
from dynamics only.

## Initialization (3D, N=128, L=20)

    φ_a(x,0) = A(r_⊥) · cos(θ + kz + 2πa/3)
    v_a(x,0) = Ω · A(r_⊥) · sin(θ + kz + 2πa/3)

## Constraints

- V > 0: angular velocity Ω > 0
- T > 0: from rotation
- δV > 0: angular momentum maintained

## Parameters

μ=-20, κ=20, m=0.0, Ω=0.1, k=2π/20
N=128, L=20, periodic BC in z, t=500

## Output

`src/dynBC.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynBC src/dynBC.c -lm`
