# V23-E: One Massless Mediator Field -- Results

## Parameters

mu = -20, kappa = 20, m1 = m2 = 1.0, A = 0.8, sigma = 3.0
Grid: Nx = 4000, xmax = 100, dx = 0.050, dt = 0.0255, tfinal = 10000

## Phase 1: Control (m3 = 1.0) vs Massless (m3 = 0.0)

### Control: m3 = 1.0 (all fields massive)

Oscillon SURVIVES to t = 10000 with:
- omega = 0.876 (below mass gap m = 1.0) -- confirmed sub-threshold oscillon
- Peak amplitudes: A1 = A2 = A3 = 0.434 (symmetric, fields locked together)
- E_total = 1.264 (slowly decreasing from initial 3.365 as transients radiate)
- f_core = 1.000 (all energy localized in core)
- Very stable: energy loss < 0.1%/1000 time units after t > 2000

### Massless: m3 = 0.0

Oscillon DOES NOT SURVIVE in the intended sense. Instead:
- Fields 1,2 retain moderate oscillation (peak ~ 0.41) at the center
- Field 3 undergoes TACHYONIC RUNAWAY: amplitude grows to ~32 at center
- phi3 develops a broad dome-shaped profile spanning the full domain
- Energy goes deeply negative: E = -46.8 (from initial +1.66)
- f_core = 0.000 (field 3 energy spread across entire domain)

**Mechanism**: With mu = -20, the coupling dV/dphi3 = mu*phi1*phi2*P/(1+kP^2)^2
gives phi3 an effective tachyonic mass m^2_eff = mu*<phi1^2*phi2^2>/... < 0 inside
the oscillon core. Since phi3 is massless, nothing prevents it from growing without
bound. The field condenses to a large static value phi3 ~ 30, on top of which
fields 1,2 continue to oscillate.

**Field 3 spatial profile at t = 10000**: Approximately linear decay from center,
with slope increasing from ~0.29 at |x|=20-30 to ~0.58 at |x|=60-70. This is
NOT an equilibrated profile -- the massless field continues to grow with time
(1D retarded Green's function: phi3 ~ Q*t + linear-in-x). The absorbing
boundaries partially damp the growth at large |x|.

## Phase 2: Mass Scan

```
m3      survived  omega   A1_peak  A3_peak  E_total     fc     Notes
------  --------  ------  -------  -------  ----------  -----  -------------------------
1.0000  YES       0.8760  0.434    0.434    +1.264      1.000  Stable oscillon, symmetric
0.8000  YES(weak) 1.0020  0.011    0.058    +0.281      0.212  Mostly dispersed, marginal
0.6000  NO        1.0020  0.040    0.051    +0.561      0.090  Dispersed
0.4000  NO*       1.0020  0.027    0.089    +0.149      0.126  Dispersed (false positive)
0.2000  NO        1.9620  0.552    2.143    -14.12      0.000  Tachyonic condensation
0.1000  NO        ---     0.438    3.761    -30.34      0.000  Tachyonic condensation
0.0500  NO        ---     0.427    8.067    -40.76      0.000  Tachyonic condensation
0.0000  NO        ---     0.406    31.952   -46.80      0.000  Tachyonic condensation
```

*m3=0.4 marked survived only because A1 barely exceeds 0.01 threshold; fc=0.126
indicates essentially dispersed.

### Two Distinct Regimes

**Regime I (m3 >= 0.8): Oscillon regime**
- The mass gap of field 3 is large enough to prevent tachyonic condensation
- At m3 = 1.0: robust oscillon with all three fields locked
- At m3 = 0.8: greatly weakened but still localized (fc = 0.21)
- The S3 symmetry breaking disrupts the resonance; fields 1,2 lose their
  coupling partner and disperse faster

**Regime II (m3 <= 0.2): Tachyonic condensation**
- Field 3 undergoes tachyonic instability: m^2_eff,3 = m3^2 + mu*<phi1*phi2> < 0
- phi3 condenses to a large static (or slowly oscillating) value
- The condensate is NOT localized -- it spreads across the domain
- A3_peak grows as m3 decreases: 2.14 (m3=0.2), 3.76 (0.1), 8.07 (0.05), 32.0 (0.0)
- Energy goes deeply negative due to the potential V = (mu/2)P^2/(1+kP^2)
- Fields 1,2 continue oscillating inside the condensate (pk ~ 0.4-0.55)

### Critical Mass

m3_c is between 0.4 and 0.8. The transition is not sharp:
- m3 = 0.8: weakly surviving (fc = 0.21, very small amplitudes)
- m3 = 0.6: dispersed (fc = 0.09)
- m3 = 0.4: dispersed / tachyonic onset
- m3 = 0.2: full tachyonic condensation

The oscillon requires m3 close to m1 = 1.0 for stability. Even m3 = 0.8
severely weakens the oscillon because the S3 symmetry breaking disrupts
the triple-product resonance that stabilizes the bound state.

## Key Physics Conclusions

1. **The mediator path FAILS for this model**: A massless field 3 does not
   produce a stable oscillon with a 1/r tail. Instead, the tachyonic
   instability drives phi3 to condensation, breaking the oscillon.

2. **Root cause**: mu < 0 (attractive coupling) creates a tachyonic effective
   mass for phi3 inside the oscillon. A massive field can resist this because
   m3^2 > |mu*<phi1*phi2>|, but a massless field cannot.

3. **The condensate is NOT a localized charge**: phi3 ~ 30 everywhere, not
   just inside the core. There is no sharp boundary or Yukawa falloff.
   The 1D massless Green's function (linear growth) means phi3 extends to
   infinity.

4. **Interestingly, fields 1,2 survive**: Despite the condensation, the massive
   fields maintain oscillation at pk ~ 0.4. They oscillate in the background
   of a large phi3 condensate. But this is NOT an oscillon -- it's two massive
   fields in a modified vacuum (phi3 != 0).

5. **Could mu > 0 (repulsive) work?** With positive mu, the effective mass
   m^2_eff,3 = m3^2 + mu*<phi1*phi2> > 0 even at m3 = 0. But then the
   triple-product coupling is repulsive, and the oscillon itself may not
   form (the v21 results specifically require mu < 0 for oscillon binding).

## Files

- `src/mediator1d.c` -- solver with per-field masses
- `data/mediator_m3_X.XXXX_ts.tsv` -- time series for each m3
- `data/mediator_m3_X.XXXX_profile.tsv` -- spatial profiles at t=10000
- `data/mediator_m3_X.XXXX_spectrum.tsv` -- DFT power spectra
- `data/scan_summary.tsv` -- summary table
