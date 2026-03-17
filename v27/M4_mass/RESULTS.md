# V27-M4: Mass from Dynamics — Results

## Summary

The bare mass parameter m is NOT required for braid survival. The braid
survives at m=0 through the triple-product binding alone. However, there
are TWO distinct regimes separated by a gap at m~0.6-0.8 where the braid
unwinds.

## M4a: Mass Reduction Scan

| m    | fc     | |P|    | E_total   | Pz    | Survived |
|------|--------|--------|-----------|-------|----------|
| 1.00 | 0.276  | 0.574  | 281       | 35.8  | YES      |
| 0.80 | 0.157  | 0.003  | 137       | 30.9  | NO       |
| 0.60 | 0.174  | 0.000  | 30.6      | 10.8  | NO       |
| 0.40 | 0.129  | 0.589  | -2871     | 2.8   | YES      |
| 0.20 | 0.108  | 1.105  | -6381     | 36.1  | YES      |
| 0.10 | 0.106  | 1.591  | -7518     | 80.2  | YES      |
| 0.00 | 0.103  | 1.742  | -7915     | 13.3  | YES      |

**Two distinct regimes:**

1. **High-mass regime (m >= 1.0)**: Standard propagating braid. Positive
   energy, moderate |P|~0.57. The mass term provides the confining potential.
   Energy E~281.

2. **Low-mass regime (m <= 0.4)**: Nonlinear collapse regime. Large NEGATIVE
   energy (E ~ -3000 to -8000). The saturated triple-product potential
   V = (mu/2)P^2/(1+kappa*P^2) dominates and drives |P| > 0.5 (above
   initial 0.127). The fields grow until the saturation mechanism (kappa
   term) balances the growth. fc drops to ~0.1-0.13.

3. **Dead zone (0.6 <= m <= 0.8)**: The mass term is too weak to confine,
   but too strong for the nonlinear collapse mechanism to activate. The
   braid disperses. |P| drops below 0.01.

**m_critical**: Formally m_critical = 0.0 (the braid survives at m=0).
But this is a DIFFERENT object — a nonlinear condensate with E < 0, not
a propagating soliton. The minimum mass for the "standard" braid regime
is m_critical ~ 1.0 (m=0.8 already fails).

## M4b: Strong Binding Compensation (m=0)

| mu    | kappa | fc    | |P|   | E_total  | Survived |
|-------|-------|-------|-------|----------|----------|
| -20   | 20    | 0.103 | 1.742 | -7915    | YES      |
| -50   | 50    | 0.095 | 2.525 | -8970    | YES      |
| -100  | 100   | 0.092 | 2.164 | -9690    | YES      |
| -200  | 200   | 0.087 | 1.643 | -10335   | YES      |

All survive at m=0. Stronger binding (larger |mu|) makes the energy more
negative but does NOT improve localization (fc decreases from 0.103 to
0.087). The nonlinear condensate spreads further with stronger binding.

**Key observation**: The saturation V_max = mu^2/(4*kappa) = |mu|/4
(constant when kappa=|mu|) controls the depth. The fc decrease means
the energy is distributed over a larger volume at stronger binding.

## M4c: Pairwise Coupling at Low Mass

| m    | lambda_pw | Result      |
|------|-----------|-------------|
| 0.40 | 0.5       | NaN (blowup)|
| 0.20 | 0.5       | NaN (blowup)|
| 0.10 | 0.5       | NaN (blowup)|
| 0.00 | 0.5       | NaN (blowup)|

**All runs blow up instantly.** The pairwise coupling
lambda_pw*(phi1*phi2 + phi2*phi3 + phi3*phi1) is an unbounded-below
bilinear potential. At low mass, there is insufficient restoring force
to prevent exponential field growth. The fields reach ~10^100 within
t~100 before hitting NaN.

**Conclusion**: Pairwise coupling is INCOMPATIBLE with low mass. It
requires m >= 1.0 for stability (where the mass term provides the
restoring force).

## M4d: Dispersion Relation Inside vs Outside Braid

| k    | omega_in | omega_out | omega_theory |
|------|----------|-----------|--------------|
| 0.20 | 0.980    | 1.080     | 1.020        |
| 0.50 | 0.980    | 1.080     | 1.118        |
| 1.00 | 0.980    | 1.080     | 1.414        |
| 1.50 | 0.980    | 1.080     | 1.803        |
| 2.00 | 0.980    | 1.080     | 2.236        |
| 3.00 | 0.980    | 1.080     | 3.162        |
| 4.00 | 0.980    | 1.080     | 4.123        |
| 5.00 | 0.980    | 1.080     | 5.099        |

**INCONCLUSIVE**: The DFT picks up the braid's own oscillation frequency
(omega ~ 0.98 inside, 1.08 outside) regardless of the perturbation
wavenumber k. The perturbation amplitude (0.01) is too small relative
to the braid background oscillation to be resolved.

The constant frequencies 0.98 and 1.08 are the braid's intrinsic modes:
- omega_in = 0.98 is close to the bare mass m=1.0 (consistent with
  a localized oscillation at roughly the mass frequency)
- omega_out = 1.08 is slightly above m=1.0 (the mass frequency in vacuum)

**To resolve the perturbation dispersion**, one would need either:
(a) much larger perturbation amplitude (risks nonlinear effects), or
(b) subtract the background oscillation before DFT, or
(c) use a static braid (no propagation velocity) as the background.

## Key Findings

1. **Mass is NOT required for braid survival** — the triple-product
   binding alone maintains |P| > 0.01 at m=0. But the m=0 state is
   a nonlinear condensate (E < 0), not a propagating soliton.

2. **The "standard" propagating braid needs m ~ 1.0** — below this,
   the braid either disperses (m=0.6-0.8) or collapses into the
   nonlinear regime (m <= 0.4).

3. **Stronger binding does NOT help localization** — increasing |mu|
   at m=0 makes fc worse (0.103 -> 0.087), spreading the condensate.

4. **Pairwise coupling is unstable at low mass** — the bilinear
   coupling creates unbounded growth without sufficient mass term.

5. **Two-regime structure** is the main discovery: the system has a
   DISCONTINUOUS transition between a standard soliton (m ~ 1) and a
   nonlinear condensate (m << 1), with a dead zone in between.

## Files

- `data/m4a_summary.tsv` — mass scan results
- `data/m4b_summary.tsv` — strong binding results
- `data/m4c_summary.tsv` — pairwise coupling results (all NaN)
- `data/m4_dispersion.tsv` — dispersion relation data
- `data/m4_meff.tsv` — effective mass extraction (inconclusive)
- `data/m4_m4a_m*.tsv` — time series for each mass value
- `data/m4_m4b_mu*.tsv` — time series for each binding strength
- `src/m4.c` — source code
