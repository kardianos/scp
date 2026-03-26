# Poynting Flux Multipole Decomposition — Single Braid (Null Model)

**Source**: `/home/d/code/scp/v34/torsion_coupling/data/braid_hires.sfa`
**Grid**: N=80, L=25.0, dx=0.6329, frames=264
**Frame pairs processed**: 263
**Time-averaged centroid**: (0.058, -0.067, 0.082)

## Method

For each consecutive frame pair (n, n+1):
1. E = -dtheta/dt via finite difference
2. B = curl(theta) via centered spatial differences
3. S = E x B (Cosserat Poynting vector)
4. |S| computed at each grid point

Time-averaged <|S|> accumulated over 263 frame pairs.
Decomposed into real spherical harmonics Y_l^m for l=0..6.

Angular grid: 90 x 180 (theta x phi), trapezoidal integration.

## Power Spectrum

C_l = (1/(2l+1)) sum_m |c_lm|^2

| R | l=0 | l=1 | l=2 | l=3 | l=4 | l=5 | l=6 | f(l=0) | f(l=1) | f(l=2) | flux | dominant |
|---|-----|-----|-----|-----|-----|-----|-----|--------|--------|--------|------|----------|
| 3.0 | 2.6589e-04 | 1.4159e-06 | 2.6156e-06 | 2.7682e-07 | 2.6471e-07 | 6.2667e-08 | 1.2214e-07 | 98.2% | 0.5% | 1.0% | 5.7803e-02 | monopole |
| 5.0 | 3.2665e-04 | 7.8747e-07 | 5.5979e-06 | 3.8306e-07 | 2.1389e-07 | 1.0678e-07 | 3.2938e-07 | 97.8% | 0.2% | 1.7% | 6.4069e-02 | monopole |
| 8.0 | 6.0093e-05 | 5.0191e-08 | 4.1811e-06 | 3.7333e-08 | 2.3536e-07 | 1.1791e-08 | 8.2467e-08 | 92.9% | 0.1% | 6.5% | 2.7480e-02 | monopole |
| 10.0 | 3.8477e-05 | 5.3096e-08 | 2.7306e-06 | 2.0116e-08 | 2.5373e-07 | 4.9413e-09 | 4.2954e-08 | 92.5% | 0.1% | 6.6% | 2.1989e-02 | monopole |
| 12.0 | 2.5929e-05 | 1.4484e-08 | 1.9215e-06 | 1.1628e-08 | 2.8522e-07 | 2.0285e-09 | 5.5916e-08 | 91.9% | 0.1% | 6.8% | 1.8051e-02 | monopole |

## Flux vs Radius (1/R^2 Test)

If radiative, total flux F(R) should scale as 1/R^2.
Check: F(R) * R^2 should be constant.

| R | F(R) | F(R)*R^2 |
|---|------|----------|
| 3.0 | 5.7803e-02 | 5.2023e-01 |
| 5.0 | 6.4069e-02 | 1.6017e+00 |
| 8.0 | 2.7480e-02 | 1.7587e+00 |
| 10.0 | 2.1989e-02 | 2.1989e+00 |
| 12.0 | 1.8051e-02 | 2.5993e+00 |

## Anisotropy

| R | max(<|S|>) | min(<|S|>) | max/min |
|---|------------|------------|--------|
| 3.0 | 9.2113e-03 | 2.0833e-03 | 4.42 |
| 5.0 | 1.0891e-02 | 1.7085e-03 | 6.37 |
| 8.0 | 8.6807e-03 | 5.0286e-04 | 17.26 |
| 10.0 | 7.3801e-03 | 3.9603e-04 | 18.64 |
| 12.0 | 7.7661e-03 | 3.9990e-04 | 19.42 |

## Interpretation

This is the **null model**: a single z-aligned braid.

### Key findings:

1. **Monopole fraction is high (92-98%) but this is an artifact of positive-definiteness.**
   |S| >= 0 everywhere, so any non-negative function projected onto Y_00 gets a large
   monopole component. This is the SAME artifact that afflicted the original |theta|
   decomposition (cf. corrected_plans.md). The monopole fraction of a positive-definite
   quantity is NOT a useful discriminant.

2. **The anisotropy ratio is the physically meaningful metric.**
   max/min grows from 4.4 (R=3) to 19.4 (R=12), showing strongly anisotropic radiation.
   A truly isotropic source would have max/min ~ 1. The single braid is FAR from isotropic.

3. **F(R)*R^2 is NOT constant** (0.52 at R=3 to 2.60 at R=12), so the Poynting flux
   does NOT follow 1/R^2 scaling. This indicates we are NOT in the radiation zone.
   The braid oscillations are standing waves in a finite box, not outgoing radiation.
   The increasing F*R^2 may reflect the boundary damping layer or finite-box effects.

4. **The l=2 (quadrupole) fraction grows with R** (1% at R=3 to 6.8% at R=12),
   consistent with a z-aligned oscillating dipole whose radiation pattern goes as
   sin^2(theta) ~ Y_00 + Y_20.

### Implications for UUD comparison:

The useful null model metric is the **anisotropy ratio** (max/min on the shell),
not the monopole fraction f(l=0). For the UUD composite:
- If three orthogonal braids produce isotropic radiation, max/min should approach 1.
- The single braid has max/min ~ 19 at R=12 -- the UUD should be significantly lower.
- A ratio max/min < 2 would indicate effective isotropization.

The monopole fraction will be high (~90%) for ANY positive-definite field regardless
of its angular structure, so it cannot distinguish isotropic from anisotropic sources.
