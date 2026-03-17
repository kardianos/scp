# V24-S5: What Determines Lambda? -- RESULTS

## Summary

Three mechanisms were tested for determining the pairwise coupling lambda.

**Test A (Gap Margin)**: The gap margin delta = 1 - omega/m_S is NOT constant.
It varies from 0.072 to 0.139 across the lambda range, with no universal value.
The oscillon frequency omega tracks m_S but the ratio is lambda-dependent.

**Test B (Energy Minimum)**: E_osc(lambda) is MONOTONICALLY INCREASING -- no minimum.
E_osc/m_S is also monotonically increasing. There is no energetically preferred lambda.

**Test C (Dynamical Modulus)**: Lambda DOES shift inside the oscillon. At g_Lambda=0.1,
the modulus field settles to Lambda(0) = 0.465 vs vacuum Lambda_0 = 0.500, a shift
of -0.035. The shift is negative (oscillon pushes lambda down) and scales with g_Lambda.

**Conclusion**: Lambda is NOT determined by self-consistency (A) or energy minimization (B).
However, if lambda is a dynamical field (C), the oscillon spontaneously generates a
local environment with a different effective coupling than the vacuum. This is the
dilaton/modulus mechanism: the oscillon "chooses" its own coupling.

## Parameters

mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
Nx=4000, xmax=100, tfinal=10000

## Test A: Gap Margin Universality

| lambda | m_S    | omega  | delta=1-omega/m_S | E_osc  | E/m_S  | peak  | fc     |
|--------|--------|--------|-------------------|--------|--------|-------|--------|
| 0.00   | 1.0000 | 0.8760 | 0.1240            |  1.264 |  1.264 | 0.383 | 0.9996 |
| 0.10   | 1.0954 | 1.0005 | 0.0867            |  1.494 |  1.364 | 0.403 | 0.9992 |
| 0.30   | 1.2649 | 1.0890 | 0.1391            |  5.680 |  4.490 | 0.326 | 0.9996 |
| 0.50   | 1.4142 | 1.2600 | 0.1090            |  6.559 |  4.638 | 0.734 | 0.9999 |
| 0.70   | 1.5492 | 1.4100 | 0.0898            | 10.078 |  6.505 | 0.739 | 0.9999 |
| 0.90   | 1.6733 | 1.5450 | 0.0767            | 11.116 |  6.643 | 0.619 | 0.9997 |
| 0.95   | 1.7029 | 1.5765 | 0.0742            | 11.419 |  6.706 | 0.669 | 0.9997 |
| 0.99   | 1.7263 | 1.6020 | 0.0720            | 11.715 |  6.786 | 0.738 | 0.9995 |

### Analysis

The gap margin delta varies from 0.072 (lambda=0.99) to 0.139 (lambda=0.30).
This is NOT constant -- it changes by nearly a factor of 2 across the range.

The overall trend is decreasing delta with increasing lambda (the oscillon
frequency gets closer to the mass gap at large lambda), but it is non-monotonic:
there is a local maximum at lambda=0.30 (delta=0.139).

The non-monotonicity suggests two competing effects:
- The triple-product coupling (mu, kappa) sets the base frequency
- The pairwise coupling shifts the effective mass and modifies the potential

At lambda=0, omega/m = 0.876 (12.4% below the mass gap), which is typical for
an oscillon with this coupling strength.

At lambda=0.99, omega/m_S = 0.928 (7.2% below gap), closer to threshold.

**Verdict**: The gap margin is NOT universal. It does not fix lambda.

## Test B: Energy Minimum

| lambda | E_osc  | E/m_S  |
|--------|--------|--------|
| 0.00   |  1.264 |  1.264 |
| 0.10   |  1.494 |  1.364 |
| 0.30   |  5.680 |  4.490 |
| 0.50   |  6.559 |  4.638 |
| 0.70   | 10.078 |  6.505 |
| 0.90   | 11.116 |  6.643 |
| 0.95   | 11.419 |  6.706 |
| 0.99   | 11.715 |  6.786 |

Both E_osc and E_osc/m_S are monotonically increasing with lambda.
There is NO minimum. Larger lambda means more pairwise potential energy
stored in the oscillon.

The large jump from lambda=0.1 (E=1.49) to lambda=0.3 (E=5.68) reflects
that the pairwise coupling begins to dominate the energy budget at moderate
lambda, storing significant energy in the phi_a phi_b cross terms.

**Verdict**: No energy minimum. Lambda is NOT determined by minimization.

## Test C: Dynamical Modulus

Parameters: Lambda_0=0.5, m_Lambda=1.0, gamma_L=0.01 (friction)

| g_Lambda | Lambda(center) | Lambda(edge) | Shift   | fc     |
|----------|---------------|--------------|---------|--------|
| 0.01     | 0.4946        | 0.5000       | -0.0054 | 0.9920 |
| 0.05     | 0.4747        | 0.5000       | -0.0253 | 0.9920 |
| 0.10     | 0.4654        | 0.5000       | -0.0346 | 0.9989 |
| 0.20     | 0.4707        | 0.5000       | -0.0293 | 0.9985 |
| 0.50     | (unstable)    | --           | --      | --     |

### Key findings

1. **Lambda shifts downward inside the oscillon.** The modulus field develops a
   "well" centered on the oscillon, with Lambda(0) < Lambda_0. The vacuum value
   Lambda_0 = 0.5 is maintained at the edges.

2. **The shift scales approximately linearly with g_Lambda** at small coupling
   (0.01 to 0.1), then saturates/reverses at g_Lambda=0.2.

3. **Equilibrium shift estimate**: For small g_Lambda, the static equilibrium
   gives delta_Lambda = -g_Lambda * <pw>_avg / m_Lambda^2. With <pw>_avg ~ 0.5
   (time-averaged pairwise sum) and m_Lambda=1.0:
   delta_Lambda ~ -0.5 * g_Lambda, so at g_Lambda=0.1: delta ~ -0.05.
   Measured: -0.035 (same order, slightly less due to oscillation averaging).

4. **Instability at large coupling**: g_Lambda >= 0.5 causes Lambda to go
   sufficiently negative that the effective mass m^2 + 2*Lambda becomes
   tachyonic, destroying the oscillon.

5. **Physical interpretation**: The oscillon creates a local "domain" with a
   different effective coupling constant. Inside the oscillon, the effective
   pairwise coupling is Lambda(0) ~ 0.465 instead of the vacuum value 0.500.
   This is a 7% reduction. The antisymmetric mass changes from
   m_A = sqrt(1-2*0.5) = 0 to m_A = sqrt(1-2*0.465) = 0.265.

### Profile behavior

The Lambda profile is smooth and approximately Gaussian-shaped (matching the
oscillon density profile). It oscillates slightly at the oscillon frequency
(since pw oscillates), but the DC component dominates.

At the boundary, Lambda relaxes cleanly to Lambda_0 = 0.5.

## Conclusions

1. **Lambda is a free parameter** in the static theory. Neither gap-margin
   universality nor energy minimization selects a preferred value.

2. **As a dynamical field, Lambda self-adjusts.** The oscillon generates a
   local environment where the effective lambda differs from the vacuum value.
   This is the modulus mechanism: the coupling constant becomes position-dependent.

3. **Direction of shift**: The oscillon REDUCES the local lambda. This is because
   the source term in the Lambda EOM (proportional to phi_a*phi_b) has a positive
   DC component that pushes Lambda down.

4. **Implications for the Proca mechanism**: If lambda is a dynamical modulus,
   then the antisymmetric mass m_A = sqrt(m^2 - lambda) is LARGER inside the
   oscillon than in vacuum. The Proca range is SHORTER inside the oscillon.
   This is opposite to what one might want for long-range forces.

5. **Self-consistent lambda**: At g_Lambda ~ 0.1, the oscillon settles into a
   state where Lambda(0) ~ 0.465. If the vacuum value Lambda_0 were also 0.465,
   the shift would be smaller, leading to a self-consistent solution. But this
   is just the fixed-point of the iteration, not a new prediction.

## Files

- `src/proca_origin.c` -- source code (Tests A, B, C)
- `data/testA_gap_margin.tsv` -- gap margin scan results
- `data/testA_lam*.tsv` -- time series per lambda value
- `data/testC_modulus_ts.tsv` -- modulus time series
- `data/testC_modulus_profile.tsv` -- spatial profiles at snapshots
