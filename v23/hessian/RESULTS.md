# V23-A: Hessian Splitting Results

## Parameters

- mu = -20.0, kappa = 20.0, m = 1.000 (v21 3D production values)
- Initial: A = 0.800, sigma = 3.000
- Grid: Nx = 4000, xmax = 80.0, dx = 0.040
- Evolution: tfinal = 5000, dt = 0.020

## Oscillon Profile

- Oscillation frequency: omega = 0.867 (omega/m = 0.867, sub-gap: YES)
- Peak envelope amplitude (t = 4000-5000): f(0) = 0.506
- The 1D oscillon exhibits a breathing modulation; peak amplitude varies over
  long timescales. The envelope represents the maximum |phi| seen in the
  final 1000 time units.

## Hessian Formulas

On the symmetric background phi_1 = phi_2 = phi_3 = f(r):

    lambda_sym(r)  = mu f^4 (5 - 7 kappa f^6) / (1 + kappa f^6)^3
    lambda_anti(r) = -mu f^4 / (1 + kappa f^6)^2

    m^2_sym(r)  = m^2 + lambda_sym(r)
    m^2_anti(r) = m^2 + lambda_anti(r)

With mu = -20 (negative):
- lambda_anti = +|mu| f^4 / (1+kf^6)^2 > 0 always --> m^2_anti INCREASES (hardened)
- lambda_sym = -|mu| f^4 (5-7kf^6) / (1+kf^6)^3 < 0 when 7kf^6 < 5 --> m^2_sym DECREASES (softened)

**Correction to proposal**: Lines 92-97 of PROPOSAL.md state that m^2_anti is
reduced. This is incorrect. The formula lambda_anti = -mu f^4/(1+kf^6)^2 with
mu < 0 gives lambda_anti > 0, meaning m^2_anti = m^2 + lambda_anti > m^2.
The antisymmetric mode is HARDENED, not softened. The symmetric mode is the
one that gets softened.

## Instantaneous Hessian (at peak amplitude)

At the oscillon center (r=0), with f(0) = 0.506, kappa f^6 = 0.336:

| Quantity | Value |
|----------|-------|
| lambda_sym(0) | -1.455 |
| lambda_anti(0) | +0.735 |
| m^2_sym(0) | **-0.455** |
| m^2_anti(0) | +1.735 |
| Delta m^2(0) = m^2_sym - m^2_anti | **-2.190** |
| m^2_sym minimum | **-0.990** at r = 0.96 |

## Time-Averaged Hessian

Numerical average of H(f cos(wt)) over one full oscillation period:

| Quantity | Value |
|----------|-------|
| <lambda_sym>(0) | -1.037 |
| <lambda_anti>(0) | +0.338 |
| <m^2_sym>(0) | **-0.037** |
| <m^2_anti>(0) | +1.338 |
| <Delta m^2>(0) | **-1.375** |
| <m^2_sym> minimum | **-0.047** at r = 0.44 |

## Tachyonic Analysis

**Symmetric mode (compression, eigenvector (1,1,1)/sqrt(3))**:
- Instantaneous: m^2_sym min = -0.990 at r = 0.96 --> **TACHYONIC** at peak amplitude
- Time-averaged: <m^2_sym> min = -0.047 at r = 0.44 --> **MARGINALLY TACHYONIC** even on average
- The symmetric mode develops a tachyonic instability inside the oscillon core

**Antisymmetric mode (shear, 2-fold degenerate)**:
- m^2_anti >= m^2 = 1.0 everywhere --> NOT tachyonic, always hardened
- The antisymmetric mode is stabilized by the triple-product coupling

## Splitting Profile

- Splitting magnitude: |Delta m^2| = 2.19 at center (219% of m^2)
- Splitting half-width: r_1/2 = 1.88
- lambda_sym sign change threshold: f_crit = 0.574 (7 kappa f^6 = 5)
- f(0) = 0.506 < f_crit, so lambda_sym remains negative throughout the core

## Key Findings

1. **Splitting exists and is large**: |Delta m^2| = 2.19 at center, 219% of m^2.
   This is not a perturbative correction -- it is an O(1) effect.

2. **Wrong sector goes tachyonic**: The SYMMETRIC (compression) mode is softened
   and becomes tachyonic (m^2_sym < 0), while the ANTISYMMETRIC (shear) mode is
   hardened (m^2_anti > m^2). This is the OPPOSITE of what the proposal
   anticipated for downstream gravity proposals.

3. **Physical mechanism**: With mu < 0, the potential V(P) = (mu/2)P^2/(1+kP^2)
   is attractive. The Hessian off-diagonal element (a != b) receives an extra
   contribution from d^2P/dphi_a dphi_b = phi != 0 (since P = phi_1 phi_2 phi_3
   has nonzero mixed second derivatives). This off-diagonal piece o > 0 (for
   mu < 0) makes the antisymmetric eigenvalue d - o larger (hardened) and the
   symmetric eigenvalue d + 2o smaller (softened).

4. **Symmetric tachyonic instability**: The compression mode m^2_sym < 0 at peak
   amplitude means the symmetric oscillon sits at a saddle point in field space
   during peak compression. The time-averaged <m^2_sym> is only marginally
   negative (-0.047), so the instability is weak on average. This may contribute
   to the breathing modulation observed in v21 simulations.

5. **Implication for gravity proposals**: The Hessian splitting exists but has the
   wrong sign for scenarios where shear modes need a reduced mass gap inside the
   soliton core. The antisymmetric (shear) sector is STIFFER, not softer.
   Proposals 2+7 that rely on shear softening cannot work with this potential.

6. **Saturation regime**: kappa f^6 = 0.336 at peak, so the oscillon is in the
   mildly saturated regime. Deep saturation (kf^6 >> 1) would flip the sign of
   lambda_sym (via the 5 - 7kf^6 factor), but this does not change the sign of
   lambda_anti which is always positive for mu < 0.

## Data Files

- `data/hessian_instant.tsv`: Instantaneous Hessian eigenvalues vs r (8 columns:
  r, f, P, lambda_sym, lambda_anti, m2_sym, m2_anti, delta_m2)
- `data/hessian_tavg.tsv`: Time-averaged Hessian eigenvalues vs r (7 columns:
  r, f, lambda_sym_avg, lambda_anti_avg, m2_sym_avg, m2_anti_avg, delta_m2_avg)
