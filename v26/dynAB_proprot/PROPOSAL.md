# V26-DynAB: Propagating + Rotating Braid

## Thesis

Combine axial propagation (A) with azimuthal rotation (B). The braid
both travels along z AND rotates around z. This gives BOTH linear
momentum (from propagation) and angular momentum (from rotation).

## Initialization (3D, N=128, L=20)

    phase_a = θ + kz + 2πa/3
    φ_a(x,0) = A(r_⊥) · cos(phase_a)
    v_a(x,0) = (ω + Ω) · A(r_⊥) · sin(phase_a)

where ω = √(k²+m²) drives axial propagation and Ω drives rotation.
The combined velocity provides both linear and angular momentum.

More precisely: v_a should decompose into axial and rotational parts:
    v_a = ω·A(r)·sin(kz + 2πa/3) + Ω·A(r)·sin(θ + kz + 2πa/3)

But if the phase already includes θ, a single velocity suffices:
    v_a = ω_eff · A(r) · sin(θ + kz + 2πa/3)
    where ω_eff combines both propagation and rotation.

## Constraints

- V > 0: both v_axial > 0 and v_angular > 0
- T > 0: from combined ω + Ω
- δV > 0: nonzero in both directions

## Parameters

μ=-20, κ=20, m=1.0, Ω=0.1, k=2π/20
N=128, L=20, periodic BC in z, t=500

## Output

`src/dynAB.c`, `data/`, `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o dynAB src/dynAB.c -lm`
