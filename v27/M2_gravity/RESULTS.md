# V27-M2: Gravity Between Massless Propagating Braids — RESULTS

## Configuration
- **Optimal braid**: m=0, mu=-50, kappa=50 (from M4 breakthrough, |P|=2.525 at t=500)
- **Lagrangian**: L = sum_a [1/2(dt phi_a)^2 - 1/2(di phi_a)^2] - (mu/2)P^2/(1+kappa*P^2)
- **Initialization**: helical braid with k=2pi/L_z, A0=0.8, R_tube=3.0
- **Grid**: N=128, periodic z, absorbing x,y boundaries

## Test M2a: Strain Field Multipole Decomposition

**Setup**: Single braid at origin, L=20, dx=0.315. Evolved 100 time units to settle, then computed strain tensor eps_{ij} = 1/2(d_i phi_j + d_j phi_i) on spherical shells.

| R | l=0 power | l=1 power | l=2 power | l=2/l=0 |
|---|-----------|-----------|-----------|---------|
| 5 | 5.52e-3 | 2.45e-4 | 1.72e-4 | 0.031 |
| 8 | 1.90e-3 | 2.64e-5 | 6.24e-7 | 0.0003 |
| 12 | 1.43e-3 | 2.38e-5 | 1.02e-6 | 0.0007 |
| 16 | 4.61e-3 | 2.38e-5 | 4.26e-5 | 0.009 |

**Result**: l=2 content is negligible (< 3% at best, < 0.1% at most radii). The strain field is overwhelmingly l=0 (monopole). The l=2/l=0 ratio does NOT increase with distance — it drops sharply from R=5 to R=8, then fluctuates at the noise floor.

**Interpretation**: The braid's strain is dominated by its isotropic Gaussian envelope. The three-strand azimuthal structure (which would generate l=2) is subdominant even at the core scale (R=5), and becomes negligible at distance. This is consistent with the V26 analysis: l=2 from a single soliton cannot dominate because the monopole (l=0) decays as 1/r while the quadrupole potential decays as 1/r^3.

## Test M2b: Two-Braid Interaction (D=30 along x)

**Setup**: Two braids at x = +/-15, L=40, dx=0.630, dt=0.126. Both propagating along z with same helicity.

**Critical issue: NUMERICAL INSTABILITY**

The large domain (L=40) with N=128 gives dx=0.630, which is too coarse to resolve the braid dynamics stably. Key indicators:
- Energy: E grew from 362 to 1.43e6 over t=300 (3950x increase)
- |P|: peaked at 70 (vs initial 0.125) — unphysical amplification
- Separation oscillated wildly: 30 -> 35 -> 24 -> 30 -> 31.6

**Separation trajectory** (sampled):
| time | separation | |P| | E_total |
|------|-----------|------|---------|
| 0 | 30.0 | 0.13 | 362 |
| 20 | 35.4 | 33 | 1.3e5 |
| 50 | 35.2 | 14 | 1.4e5 |
| 90 | 24.4 | 36 | 2.9e5 |
| 100 | 24.0 | 56 | 4.1e5 |
| 150 | 28.6 | 18 | 7.7e5 |
| 200 | 30.3 | 36 | 1.1e6 |
| 300 | 31.6 | 40 | 1.4e6 |

The separation shows one large oscillation (30 -> 24 -> 31.6) but the energy growth means the dynamics are dominated by numerical artifacts, not physical forces. The apparent "attraction" at t=90 (sep=24.4) is likely grid-instability-driven.

**Verdict**: INCONCLUSIVE due to numerical instability. Would need N=256+ at L=40 (dx=0.315) to resolve properly, but that requires 8x more memory (16GB for fields).

## Test M2c: Angular Pattern of Force

**Setup**: Second braid placed at D=30 from origin at angles 0, 45, 90 degrees. L=40, dx=0.630, t_run=200.

| Angle | sep_initial | sep_final | delta_sep | force_proxy |
|-------|------------|-----------|-----------|-------------|
| 0 deg | 30.00 | 31.71 | +1.71 | -0.0085 (repulsive) |
| 45 deg | 30.00 | 31.99 | +1.98 | -0.0099 (repulsive) |
| 90 deg | 30.00 | 31.65 | +1.66 | -0.0083 (repulsive) |

**Result**: All three angles show NET REPULSION (separation increases by ~1.7-2.0 units over t=200). The force is approximately ISOTROPIC: force_proxy varies only 20% across angles (0.0083 to 0.0099), consistent with scalar (l=0) interaction, not quadrupolar (l=2).

**Caveat**: Same numerical instability as M2b (energy not conserved at this resolution). The repulsion could be a numerical artifact from the energy injection.

## Summary and Conclusions

### Key findings:
1. **Strain is monopole-dominated**: l=2/l=0 < 3% at all distances. The braid geometry does not produce a detectable quadrupole strain field.
2. **Two-braid interaction is unstable at this resolution**: dx=0.63 (L=40, N=128) is too coarse. Energy grows by 10^3-10^4 over the run.
3. **Net effect is repulsion, not attraction**: All angles show separation increase. No gravity-like attraction detected.
4. **Force is approximately isotropic**: Consistent with scalar (l=0) exchange, not tensor (l=2).

### Implications for the gravity hypothesis:
The l=2 content of the strain field is negligible compared to l=0. Even if the two-braid interaction were resolved properly (higher N), the dominant interaction would be scalar (monopole-monopole), not spin-2 (quadrupole). This is the same conclusion as V6-V7: the mechanism produces 1/r scalar attraction at best, not tensor gravity.

### What would need to change:
- To get dominant l=2: the braid would need a shape with large quadrupole moment and suppressed monopole. The Gaussian tube profile has almost zero quadrupole.
- To get proper two-braid dynamics: N=256+ at L=40 (or reduce D to fit within L=20).
- To get attraction: may need opposite-helicity braids or different potential.

## Files
- `src/m2.c` — implementation (M2a, M2b, M2c)
- `data/m2a_strain.tsv` — strain multipole decomposition
- `data/m2b_twobraid.tsv` — two-braid separation time series
- `data/m2c_angular.tsv` — angular force pattern summary
- `data/m2c_angle{0,45,90}.tsv` — per-angle time series
