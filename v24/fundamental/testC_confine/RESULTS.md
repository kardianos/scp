# Test C: Confining Potential — RESULTS

## Answer: NO, the confining potential does NOT stop radiation

The linear confining potential V = -sigma*sqrt(P^2+eps^2) does NOT produce
perfectly stable oscillons. dE/dt -> 0 is never achieved; instead dE/dt ~ t^{-0.65},
a power-law decay that never reaches zero.

## Summary Table

| Mode | Parameters | Stable? | E_loss (t=20000) | dE/dt (t=20000) | fc | omega |
|------|-----------|---------|-------------------|------------------|------|-------|
| 0: Saturating | mu=-10, kappa=0.1 | Runaway | N/A (fields grow) | -4e-5 | 0.000 | N/A |
| 1: Confining | sigma=1 | YES | -14.9% | -7.5e-6 | 0.999 | 0.708 |
| 1: Confining | sigma=2 | BLOWUP (t~3) | - | - | - | - |
| 1: Confining | sigma=3 | BLOWUP (t~2) | - | - | - | - |
| 1: Confining | sigma=5 | BLOWUP (t~1) | - | - | - | - |
| 1: Confining | sigma=10 | BLOWUP (t~1) | - | - | - | - |
| 1: Confining | sigma=20 | BLOWUP (t<1) | - | - | - | - |
| 1: Confining | sigma=50 | BLOWUP (t<1) | - | - | - | - |
| 2: Conf+quartic | sigma=1, kc=0.01 | YES | -10.2% | -8.4e-6 | 1.000 | 0.690 |
| 2: Conf+quartic | sigma=5, kc=0.1 | Runaway | N/A (fields grow) | -1.3e-5 | 0.000 | N/A |
| 2: Conf+quartic | sigma=10, kc=1.0 | Runaway | N/A (fields grow) | -4.1e-5 | 0.000 | N/A |

## Detailed Findings

### 1. Confining potential (mode 1): V = -sigma*sqrt(P^2 + eps^2)

**Only sigma=1 survives.** All sigma >= 2 blow up within t < 3.

The instability mechanism: the confining potential creates an unbounded attractive force.
For the initial Gaussian (A=0.8, sigma_init=3), the product P = phi^3 ~ 0.51 at center.
The confining force scales as sigma * P * dP / sqrt(P^2+eps^2) ~ sigma * dP ~ sigma * phi^2.
At sigma >= 2, this force exceeds the restoring gradient + mass terms, causing
exponential growth: phi -> inf in O(1) time steps.

For sigma=1, the force is weak enough that the mass gap (m=1) and gradient energy
prevent runaway. The oscillon oscillates with omega = 0.708 < m = 1.0 (below mass gap,
confirmed sub-threshold). Core fraction fc = 0.999, indicating excellent spatial
confinement.

**But radiation persists.** The energy loss over t=20000 is 14.9%, with dE/dt following
a power law:

    dE/dt ~ 5.3e-3 * t^{-0.652}

| Time | <dE/dt> |
|------|---------|
| t~150 | -1.6e-4 |
| t~2000 | -5.7e-5 |
| t~5000 | -2.4e-5 |
| t~10000 | -1.3e-5 |
| t~15000 | -8.8e-6 |
| t~20000 | -7.2e-6 |

Extrapolation: at t=10^5, dE/dt ~ 3e-6; at t=10^6, dE/dt ~ 7e-7. Never zero.

### 2. Confining + quartic (mode 2): V = -sigma*sqrt(P^2+eps^2) + (kc/2)*P^4

The quartic term provides a restoring force at large P, preventing blowup.

sigma=1, kc=0.01: Stable oscillon, -10.2% energy loss (slightly better than pure confining).
omega = 0.690, fc = 1.000. Late-time dE/dt ~ -8.4e-6 (comparable to pure confining).

sigma=5, kc=0.1: Does not blow up (quartic prevents that), but the fields grow to large
amplitude (~1.4) and the oscillon delocalizes (fc=0.000). The system runs away to a
spatially extended state.

sigma=10, kc=1.0: Same — no blowup but delocalization. Fields reach ~1.2 everywhere.

### 3. Baseline saturating potential (mode 0): V = (mu/2)*P^2/(1+kappa*P^2)

With mu=-10, kappa=0.1: The strong attractive potential causes the fields to grow from
A=0.8 to ~2.3 and spread across the domain (fc=0). This is NOT a standard oscillon but
rather a spatially-extended condensate. Not a useful comparison.

The saturating potential at these parameters has V_max = |mu|/(2*kappa) = 50, which is
much stronger than the confining potential at sigma=1 (V ~ -sigma*A^3 ~ -0.5 at center).

## Physics Interpretation

The confining potential does NOT solve the radiation problem because:

1. **Confinement != stability.** The linear potential prevents the fields from
   spatially spreading (fc ~ 1), but the oscillon still sheds energy through
   frequency harmonics above the mass gap. The fundamental at omega=0.71 is
   below m=1, but higher harmonics (2*omega=1.42 > m) can propagate.

2. **The radiation mechanism is spectral, not spatial.** The oscillating fields
   generate harmonics nw for n=1,2,3,... Any harmonic with nw > m can radiate
   as a free Klein-Gordon wave. The confining potential doesn't change the
   harmonic structure — it only changes the amplitude.

3. **Narrow stability window.** The confining potential is linearly unstable for
   sigma >= 2: the attractive force overcomes the mass+gradient restoring force.
   The only stable value (sigma=1) gives weak confinement, barely different from
   the quadratic potential.

4. **The quartic stabilizer helps but doesn't eliminate radiation.** Adding
   kc*P^4 prevents blowup but the radiation rate is essentially the same
   (dE/dt ~ 8e-6 vs 7e-6 at t=20000).

## Conclusion

**NEGATIVE RESULT.** Linear confinement does not produce perfectly stable oscillons.
The radiation is fundamentally a spectral phenomenon (harmonics above mass gap),
not a spatial leaking phenomenon. Confining the fields spatially does not prevent
energy loss through propagating wave modes.

The key insight: oscillons radiate because they OSCILLATE, and any periodic motion
with omega < m will have harmonics n*omega > m for sufficiently large n. The only
way to stop radiation completely is to either:
- Remove the mass gap entirely (but then there is no gap at all)
- Make the oscillon non-oscillating (but then it is not an oscillon)
- Find a potential where ALL harmonics are exactly zero (requires exact integrability)

## Parameters Used

- m=1.0, A_init=0.8, sigma_init=3.0 (Gaussian width)
- Nx=4000, xmax=100, dx=0.05, dt=0.025, tfinal=20000
- eps=1e-6 (regularization of |P| cusp)
- Absorbing boundary: outer 25% of domain

## Files

- `src/confine.c` — solver (modes 0,1,2)
- `data/confine_mode{0,1,2}_sc{X}_ts.tsv` — time series
- `data/confine_mode{0,1,2}_sc{X}_spectrum.tsv` — frequency spectra
