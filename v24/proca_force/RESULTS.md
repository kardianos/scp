# V24-S1: Proca Force Measurement -- RESULTS

## Parameters

- mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
- Equilibration: Nx=4000, xmax=80, t=10000
- Interaction: Nx=8000, xmax=300, t=2000
- Lambda scan: {0.0, 0.5, 0.9, 0.99, 0.999}
- D scan: {15, 20, 30, 40, 60, 80, 100}

## Step 1: Equilibration Summary

All oscillons survive with fc > 0.999 at all lambda values.

| lambda | m_A    | 1/m_A  | omega  | E_total | A_peak |
|--------|--------|--------|--------|---------|--------|
| 0.000  | 1.0000 | 1.00   | 0.8750 | 1.26    | 0.496  |
| 0.500  | 0.7071 | 1.41   | 1.2600 | 6.56    | 0.733  |
| 0.900  | 0.3162 | 3.16   | 1.5450 | 11.12   | 0.753  |
| 0.990  | 0.1000 | 10.00  | 1.6020 | 11.72   | 0.751  |
| 0.999  | 0.0316 | 31.62  | 1.6080 | 11.80   | 0.750  |

Key observations:
- Pairwise coupling raises oscillation frequency: omega goes from 0.875 (lambda=0) to 1.608 (lambda=0.999)
- For lambda >= 0.5, omega > m=1 so the oscillation is ABOVE the single-field mass gap
  but below the symmetric mass m_S = sqrt(m^2 + 2*lambda)
- Energy increases dramatically with lambda (1.26 to 11.80) due to pairwise potential energy
- All lambda values produce stable, long-lived oscillons

## Step 2: Force Measurements F(D)

Raw acceleration values (d^2 sep / dt^2) from quadratic fit over t in [20, 500]:

| lambda | D=15      | D=20      | D=30      | D=40      | D=60      | D=80      | D=100     |
|--------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 0.000  | -3.9e-05  | -1.4e-04  | -6.1e-07  | +9.1e-07  | -8.4e-07  | -6.4e-07  | -1.1e-06  |
| 0.500  | +3.7e-04  | -9.5e-06  | -1.4e-06  | +1.8e-06  | -3.7e-07  | -1.6e-06  | -1.8e-06  |
| 0.900  | +2.4e-04  | -1.3e-05  | -2.7e-06  | +2.0e-06  | +5.7e-06  | -3.0e-06  | -8.6e-06  |
| 0.990  | +2.3e-04  | +2.0e-06  | +8.1e-06  | -1.7e-06  | -1.8e-05  | +5.4e-06  | +1.6e-06  |
| 0.999  | +2.3e-04  | +6.8e-06  | +2.1e-06  | -1.0e-05  | -1.5e-05  | +7.9e-06  | -1.5e-05  |

## Step 3: Yukawa Fits

| lambda | m_A    | 1/m_A (pred) | lam_fit | F0_fit   | ratio |
|--------|--------|--------------|---------|----------|-------|
| 0.000  | 1.0000 | 1.00         | 22.5    | 2.57e-05 | 22.5  |
| 0.500  | 0.7071 | 1.41         | 24.8    | 2.67e-05 | 17.6  |
| 0.900  | 0.3162 | 3.16         | 48.2    | 2.29e-05 | 15.3  |
| 0.990  | 0.1000 | 10.00        | 42.4    | 2.29e-05 | 4.2   |
| 0.999  | 0.0316 | 31.62        | 122.7   | 1.93e-05 | 3.9   |

## Step 4: Critical Comparison

**F at D=30:**
- lambda=0.00: F = -6.1e-07 (range 1/m_A = 1.0)
- lambda=0.99: F = +8.1e-06 (range 1/m_A = 10.0)
- Ratio |F(0.99)|/|F(0)| = 13.4

**F at D=60 for lambda=0.999:**
- F = -1.5e-05 (range 1/m_A = 31.6)
- Detectable: YES (|F| > 1e-10)

## Assessment: INCONCLUSIVE (Noise-Dominated)

The force measurements are **dominated by systematic noise** at the ~1e-5 level.
The evidence for this is:

1. **Sign flips**: F changes sign between adjacent D values (e.g., lambda=0.99:
   D=30 positive, D=40 negative, D=60 negative, D=80 positive). A real Yukawa
   force would decay monotonically.

2. **No exponential decay**: The Yukawa fits return ranges of 22-123, far
   exceeding 1/m_A predictions (1-32). This is because the fit is fitting noise,
   not a real exponential.

3. **Flat noise floor**: For D >= 30, all lambda values show |F| ~ 1e-6 to 1e-5
   with random signs. There is no clear dependence on lambda.

4. **The D=15 signal is real but lambda-independent**: All lambda > 0 give
   F(D=15) ~ +2.3e-4 (identical within 5%). This is direct overlap of oscillon
   tails, not Proca exchange. The lambda=0 case has a different sign (-3.9e-5),
   likely because the lambda=0 oscillon has different shape (different omega,
   different energy).

5. **The D=30 "13x enhancement" is spurious**: lambda=0 at D=30 gives F=-6.1e-07
   which is indistinguishable from the noise floor of ~1e-6. The ratio is
   meaningless.

### Noise source

The centroid oscillation from breathing (~0.5 in amplitude at omega~1) creates
periodic jitter in the energy centroid position. Over the 500-time-unit fit
window, this breathing modulation leaks into the quadratic coefficient at the
~1e-5 level, swamping any real inter-oscillon force.

The cycle-averaging (T_avg=10) mitigates but does not eliminate this because:
- The breathing is not exactly periodic (amplitude modulation)
- Two oscillons breathe at slightly different phases, creating beat patterns
- Energy redistribution between the oscillons modulates centroid positions

### What would be needed

To measure forces below ~1e-5, would need:
- Much longer evolution (t >> 10000) so small acceleration accumulates into detectable displacement
- Or: direct momentum measurement instead of centroid fitting
- Or: use the total energy method (measure interaction energy vs. separation directly)
- Or: adiabatic approximation (compute force from static field overlap, not dynamics)

The Yukawa fit ranges do NOT track 1/m_A. The interaction range cannot be extracted
from this data.
