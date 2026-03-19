# V32 High-Resolution Campaign Results

## Overview

Extended the binding-weighted gradient coupling experiments to higher resolution
(N=256, N=512) to validate the attraction signal seen at N=128.

## Key Finding: Model is Unstable at High Resolution

The gradient coupling `w(P) * alpha * (grad_rho/rho) . grad_phi` is
non-conservative and creates runaway energy growth at ALL resolutions.
The instability time depends on resolution:

| N   | L  | dt_factor | Blowup time T_blow |
|-----|----|-----------|--------------------|
| 128 | 30 | 0.12      | > 300 (slow growth) |
| 256 | 30 | 0.12      | ~ 210              |
| 512 | 40 | 0.25      | ~ 80               |
| 512 | 20 | 0.25      | ~ 30 (single braid)|

Higher resolution resolves sharper gradients, which the coupling amplifies.

## E2: Two-Braid Interaction (N=256, L=30, alpha=0.5)

### Attraction Phase (T=0 to 200)

Clear attraction with oscillating D:

    T      D       E/E0
    0    19.99    1.000
   50    19.71    1.098
   80    19.12    1.121   <-- minimum in first oscillation
  100    19.98    1.129   <-- bounces back
  130    19.28    1.173
  160    16.82    1.250   <-- strong attraction
  180    15.85    1.383
  200    13.35    1.787   <-- deepest approach

D decreased from 20.0 to 13.35 (34% reduction) over T=200.
Energy grew 79% during this phase.

### Blowup Phase (T > 200)

    T      D       E/E0
  210    14.59    2.282   <-- begins reversal
  220    25.46    3.207   <-- explosive separation
  230    35.19    4.838
  260    37.71   18.569

## E1: Single Braid (N=256, L=20, alpha=0.5)

Single braid shows same energy instability:

    T     E/E0    fc     max_rho
    0    1.000   0.628   2.69
   30    1.299   0.502   1.33
   60    1.638   0.264   1.48
   80    2.911   0.602   3.04
  100    8.957   0.221   9.35

Braid disperses as energy grows. No steady state achieved.

## L-Dependence of Attraction

| L  | alpha=0.5 | D at T=100 | Verdict  |
|----|-----------|------------|----------|
| 30 | 0.5       | 19.98      | ATTRACTION (long-term) |
| 40 | 0.5       | 23.56      | REPULSION |

The "attraction" at L=30 may involve periodic BC effects:
- At L=30 with D=20, the braids fill 2/3 of the box
- At L=40, braids are more isolated and repel

## E4/E5 Results (from other campaign agent)

E4 sweep at N=128, T~100:
- alpha=0.1-0.2: slight attraction (too short to be conclusive)
- alpha=0.3-0.6: neutral
- alpha=0.7: repulsion

E5 force law at various D:
- D=10 (L=40): delta_D = -3.44 (strong attraction)
- D=15 (L=50): delta_D = -4.80 (strong attraction)
- D=20 (L=60): delta_D = +0.44 (neutral)
- D=25 (L=70): delta_D = +0.01 (neutral)

Attraction only at small D where braids overlap. Short-range, not 1/r.

## Conclusions

1. The binding-weighted gradient coupling creates REAL short-range attraction
   between braids at separations D < 15 (about 5 core radii).

2. The coupling is non-conservative: energy grows exponentially.
   This is NOT a numerical artifact -- it's intrinsic to the model.

3. No steady state exists. All runs eventually blow up.

4. N=512 is too fine for this model -- the resolved gradients
   trigger faster instability. N=256 is the practical maximum.

5. The attraction is L-dependent (boundary effects matter).
   At large L (isolated braids, D=20), the interaction is repulsive.

## Files

- src/v32_run.c: Single-run binary with CLI parameters
- data/E1/: Single braid, N=256, L=20 (unstable past T~80)
- data/E2/: Two braids, N=256, L=30 (attraction T=0-200, blowup T>210)
- data/E2/profile_b*_t*.tsv: Radial profiles at T=0, 100, 200
- data/E4/sweep_summary.tsv: Alpha sweep (from other campaign)
- data/E5/force_law.tsv: Force law scan (from other campaign)
- monitor.sh: Monitoring script for background runs
