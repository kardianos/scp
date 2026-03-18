# T7: Two-Braid Interaction at the Bimodal Sweet Spot

## Question
Do two bimodal braids attract? What is the force law (1/r² or Yukawa)?
Does the interaction have both "gravitational" and "electromagnetic" components?

## Method
1. Two bimodal braids separated by D=30 along x-axis, centered at (±15,0,0)
2. Three orientation configs:
   a. Same twist (both W=-1)
   b. Opposite twist (W=-1 and W=+1, flip phi₁ → -phi₁ for second braid)
   c. Perpendicular (second braid propagates along y instead of z)
3. N=192, L=40 (dx=0.42, adequate resolution)
4. Run to T=500
5. Track separation D(t) via center-of-mass of |φ|² in each half

## IMPORTANT: Dynamics Mass
Use BIMODAL params as-is: mass=1.50 (m²=2.25 in EOM). Pass mass2_override=-1.

## Key Observable
- ΔD < 0: attraction, ΔD > 0: repulsion
- Same-twist vs opposite: if both attract equally → gravity-like (universal)
  If different → EM-like (charge-dependent)
- Fit D(t): acceleration a = d²D/dt² → force F = M·a
- Force vs distance: F ∝ 1/D² (Coulomb) or F ∝ exp(-mD)/D (Yukawa)?

## Grid & Runtime
- N=192, L=40: 7.1M points, ~3 min per eval with OMP
- 3 orientation configs × 1 run each = ~9 min

## Build
```
cd v29/T7_twobraid && gcc -O3 -fopenmp -o t7 src/t7.c -lm
./t7
```
