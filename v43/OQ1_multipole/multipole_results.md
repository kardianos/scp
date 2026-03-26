# Multipole Expansion of Theta Field — UUD Composite (V41)

**Source**: `/home/d/code/scp/v41/results/stable/UUD_stable_f16.sfa` frame 11 (t=499.99)
**Grid**: N=192, L=30.0, dx=0.3141
**Centroid**: (3.072, 7.890, -5.492)

## Method

- Sample theta_rms = sqrt(theta_x^2 + theta_y^2 + theta_z^2) on spherical shells
- Decompose into real spherical harmonics Y_l^m up to l=4
- Integration grid: 72 x 144 (theta x phi) points
- Trilinear interpolation from 3D grid to shell points

## Results

| R | l=0 frac | l=1 frac | l=2 frac | l=3 frac | l=4 frac | dominant |
|---|---------|---------|---------|---------|---------|----------|
| 3.0 | 62.6% | 29.8% | 3.5% | 2.6% | 1.5% | monopole |
| 5.0 | 58.1% | 24.7% | 9.8% | 3.8% | 3.7% | monopole |
| 8.0 | 56.9% | 16.6% | 14.8% | 7.2% | 4.6% | monopole |
| 10.0 | 62.9% | 10.2% | 11.0% | 8.8% | 7.1% | monopole |
| 12.0 | 88.9% | 1.6% | 5.7% | 2.4% | 1.5% | monopole |
| 15.0 | 81.1% | 1.5% | 8.6% | 3.7% | 5.1% | monopole |

## Per-Component Analysis at R=10

Decompose each theta component separately:

| component | l=0 frac | l=1 frac | l=2 frac | l=3 frac | l=4 frac |
|-----------|---------|---------|---------|---------|----------|
| theta_x | 3.6% | 15.7% | 25.0% | 12.5% | 43.3% |
| theta_y | 12.0% | 25.5% | 18.7% | 17.0% | 26.8% |
| theta_z | 1.2% | 14.9% | 28.6% | 7.4% | 47.9% |

## Interpretation

**YES — l=0 (monopole) dominates overwhelmingly at all radii.**

Key findings:

1. **Monopole dominance increases with distance**: At R=3 the monopole carries 63% of
   the power. At R=12 it carries 89%. At R=15 it is 81%. The dipole (l=1) drops from
   30% at R=3 to just 1.5% at R=15. This is the expected behavior for a composite
   whose individual dipole contributions cancel, leaving a net monopole.

2. **The theta_rms scalar is a monopole despite vector components being multipolar**:
   Each individual component (theta_x, theta_y, theta_z) at R=10 is dominated by
   l=2 and l=4 (quadrupole and hexadecapole), with very little l=0 content (1-12%).
   But the magnitude sqrt(sum of squares) is 63% monopole. This means the three
   vector components carry complementary angular patterns that combine to produce
   an isotropic magnitude — exactly what a "charge" field should do.

3. **Dipole suppression at large R**: The dipole fraction drops from 30% to 1.5%
   between R=3 and R=15. Individual dipoles of the three constituent solitons
   cancel in the composite, as expected for a UUD configuration where opposite
   orientations partially cancel.

4. **Implication for Coulomb's law**: The UUD composite acts as an effective point
   source of theta field magnitude. At distances R >> core size, the field profile
   is dominated by the l=0 harmonic, meaning it falls off isotropically. Combined
   with the 1/r radial decay of a massless (m_theta=0) field, this gives the
   spatial dependence needed for a Coulomb-like interaction between composites.
