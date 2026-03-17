# V26-DynB: Rotating Braid (Spin from Angular Rotation)

## Thesis

The braid pattern rotates around the z-axis at angular frequency Ω,
creating angular momentum (spin) and gyroscopic stability.

## Initialization (3D, N=128, L=20)

    φ_a(x,0) = A(r_⊥) · cos(θ + kz + 2πa/3)
    v_a(x,0) = +Ω · A(r_⊥) · sin(θ + kz + 2πa/3)

where θ = atan2(y,x) is the azimuthal angle, Ω is the rotation frequency.

The pattern has BOTH axial twist (from kz) AND azimuthal winding (from θ).
The rotation at Ω maintains the azimuthal structure.

## Constraints

- V > 0: the pattern rotates at Ω (angular velocity always nonzero)
- T > 0: the rotation frequency Ω provides thermal content
- δV > 0: angular momentum L = I·Ω > 0 (gyroscopic stability)

## Rotation Frequency

Start with Ω = 0.1 (slow rotation). The angular momentum L ∝ Ω.
Scan Ω = {0.05, 0.1, 0.2, 0.5} to find optimal stability.

## Lagrangian

Same as V26 baseline:
    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)²] - ½m²Σφ_a² - (μ/2)P²/(1+κP²)

## What to Measure

1. Does rotation STABILIZE the braid (fc higher than static)?
2. Is angular momentum conserved (L(t) = const)?
3. Is the braid gyroscopically stable (resists tipping)?
4. l=2 content — does rotation affect the multipole structure?
5. Writhe: does the rotating braid have quantized angular momentum?

## Parameters

μ=-20, κ=20, m=1.0, A₀=0.8, R_tube=3.0
k=2π/20, Ω scan: {0.05, 0.1, 0.2, 0.5}
N=128, L=20, periodic BC in z, t=500

## Output

`src/dynB.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynB src/dynB.c -lm`
