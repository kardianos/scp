# T4: Topological Fragility — Critical Perturbation

## Question
How strong a perturbation can the braid survive without losing W=-1?

## Motivation
π₁(R³)=0 means the winding is dynamically conserved, not topologically
protected. Find the critical perturbation amplitude ε_crit. If ε_crit ≫ A₀,
the protection is effectively topological for all practical purposes.

## Method
1. Initialize bimodal braid, evolve for T=200 to settle
2. At T=200, apply localized Gaussian perturbation at the core:
   φ_a(x) → φ_a(x) + ε × G(x) × n_a
   where G(x) = exp(-|x-x_center|²/σ²), σ = R_tube = 3.0
   and n_a = random unit vector (same for all x)
3. Scan ε ∈ {0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0}
4. Continue evolution to T=500
5. Track: winding(t), fc(t), |P|(t) after perturbation

## Also test
- Perturbation aligned with field (n = (1,0,0)) vs random
- Perturbation at core (x_center = origin) vs off-center (x_center = (3,0,0))
- Perturbation to velocity only (ε × G × n added to ∂_t φ, not φ)

## Key Observable
- ε_crit / A₀: ratio of critical perturbation to field amplitude
  If > 5: effectively topological
  If 1-5: moderately robust
  If < 1: fragile

## Grid & Runtime
- N=96, L=20, T=500: ~5 min per eval
- 11 amplitudes × 3 perturbation types = ~33 runs
- But can precompute the T=200 state once and branch from it
- With branching: 1 setup (2 min) + 33 × 3 min = ~100 min
- Or parallelize: 33/16 × 3 min ≈ 7 min (single-threaded per eval)

## Build
```
cd v29/T4_fragility && gcc -O3 -fopenmp -o t4 src/t4.c -lm
./t4
```
