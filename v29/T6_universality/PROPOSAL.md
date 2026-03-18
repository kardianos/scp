# T6: Universality — Bimodal Synergy Across Parameter Space

## Question
Is the bimodal synergy specific to μ=-41.3, κ=50, or does it persist
across parameter space?

## Method
1. Fix the bimodal GEOMETRY from t=0.85 (phases, ellip, amplitudes)
2. Vary the dynamics: scan (μ, κ) pairs
3. For each pair, run Pure A geometry, Pure B geometry, and Bimodal (t=0.85)
4. Check: does the bimodal beat both controls for each (μ, κ)?

## IMPORTANT: Dynamics Mass
All runs use mass=1.50 from BIMODAL params (m²=2.25 in EOM) unless explicitly
varying mass. Use mass2_override=-1 (default from params). When varying (μ,κ),
keep mass=1.50 fixed. The V28 validated result includes this mass.

## Scan Grid
    (μ, κ) ∈ {(-20,20), (-30,30), (-41.3,50), (-60,60), (-80,80), (-100,100)}
    Also: (-41.3,20), (-41.3,100) — vary κ at fixed μ
    And:  (-20,50), (-80,50) — vary μ at fixed κ
    Total: 10 parameter pairs × 3 configs = 30 runs

## Key Observable
- For each (μ,κ): synergy ratio = min(bimodal_trans/A_trans, bimodal_tor/B_tor)
  If > 1.0: synergy exists at this parameter point
- Plot synergy ratio vs μ, vs κ: is there a broad region or a narrow peak?

## Grid & Runtime
- N=80, T=200: ~30s per eval (single-threaded, 16 in parallel)
- 30 runs / 16 cores × 30s ≈ 1 min. Very fast.

## Build
```
cd v29/T6_universality && gcc -O3 -fopenmp -o t6 src/t6.c -lm
./t6
```
