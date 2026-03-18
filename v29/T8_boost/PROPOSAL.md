# T8: Lorentz Boost Validation

## Question
Does a boosted bimodal braid Lorentz-contract correctly?

## Motivation
The Lagrangian is Lorentz-invariant. A boosted soliton should contract
by factor γ along the boost direction. This validates the relativistic
structure and confirms the soliton behaves as a particle.

## Method
1. Initialize bimodal braid with Lorentz boost along x:
   φ_a(x,y,z) → φ_a(γ(x-x₀), y, z)
   v_a(x,y,z) → -γv × ∂_x φ_a + γ × v_a_original
2. Test v = 0.1c, 0.3c, 0.5c (γ = 1.005, 1.048, 1.155)
3. N=128, L=30 (extra room for motion), T=200
4. Measure transverse (y) and longitudinal (x) profiles at T=100
5. Compare aspect ratio to expected γ

## IMPORTANT: Dynamics Mass
Use BIMODAL params as-is: mass=1.50 (m²=2.25 in EOM). Pass mass2_override=-1.
The Lagrangian L = ½(∂φ)² - ½m²φ² - V(P) with m=1.50 IS Lorentz-invariant.

## Also check
- Energy: E should equal γ × E_rest (relativistic energy-momentum)
- Momentum: P_x should equal γMv
- Does the braid survive the boost or does it shed/deform?

## Grid & Runtime
- N=128, L=30, T=200: ~5 min per eval
- 3 velocities + 1 control (v=0): 4 runs, ~20 min

## Build
```
cd v29/T8_boost && gcc -O3 -fopenmp -o t8 src/t8.c -lm
./t8
```
