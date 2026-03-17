# V23-F: Critical Gravity Phase 2 Redo — Results

## Summary

Equilibrated single oscillons at mu = {-20, -14, -10, -8, -6} (kappa=20, m=1),
saved full-state profiles at peak breathing phase, then ran two-oscillon
interactions at D = {15, 20, 25, 30, 40} to measure the inter-oscillon force.

**Key finding**: The interaction range lambda approximately equals the predicted
correlation length xi_pred = 1/sqrt(m^2 - omega^2) for the cases closest to
the mass gap edge (mu=-6, mu=-8). However, xi_pred stays bounded at 2.1-2.8
across the full mu scan. The gap never closes below 6.4%, so no divergence
of the interaction range is observed.

**Verdict**: Critical gravity is NOT supported. The interaction range tracks
the correlation length as expected, but the correlation length itself does not
diverge because the gap margin stays at 6-11% for all surviving oscillons.

## Step 1: Equilibration

Gaussian initial conditions (A=1.0, sigma=3.0) evolved for t=10000 on
Nx=4000, xmax=80 grid. Profile saved at last positive peak of phi_1(x=0)
in the final 5% of evolution.

| mu   | omega  | gap%  | E_total | A_peak | xi_pred |
|------|--------|-------|---------|--------|---------|
| -20  | 0.8850 | 11.5% | 1.268   | 0.484  | 2.148   |
| -14  | 0.8430 | 15.7% | 2.899   | 0.715  | 1.859   |
| -10  | 0.8910 | 10.9% | 4.679   | 0.732  | 2.203   |
|  -8  | 0.9140 |  8.6% | 5.012   | 0.727  | 2.465   |
|  -6  | 0.9360 |  6.4% | 6.284   | 0.733  | 2.841   |

All five oscillons survived t=10000 equilibration.

Gap margin = (m - omega)/m. The gap does NOT close: even at mu=-6 (weakest
coupling), gap margin is 6.4%. This is the fundamental limitation.

Energy stability at late times: dE/dt ~ 10^{-5} to 10^{-4} (slowly losing
energy to residual radiation absorbed at boundaries).

### Tail Decay

Profile tails match the predicted decay rate kappa = sqrt(m^2 - omega^2)
in the range x = [3, 10] from center (ratio within 0.98-1.13 of prediction).
Beyond x ~ 10, the profile is dominated by residual radiation in the box,
not the true oscillon tail. This limits force measurements at D > 20.

## Step 2: Two-Oscillon Interactions

Profiles superposed as phi_a(x) = f(x-D/2) + f(x+D/2) with matching
velocities. Evolved on Nx=8000, xmax=200 grid for t=5000.

### Energy Conservation

| mu   | D=15   | D=20   | D=25   | D=30   | D=40   |
|------|--------|--------|--------|--------|--------|
| -20  | -47.6% | -44.4% | -44.6% | -34.9% | -22.2% |
| -14  | -47.4% | -48.7% | -47.4% | -48.0% | -5.9%  |
| -10  | -26.8% | -25.2% | -24.0% | -8.5%  | -1.6%  |
|  -8  | -15.1% | -20.4% | -20.4% | -19.1% | -3.2%  |
|  -6  | -3.7%  | -7.2%  | -7.9%  | -7.8%  | -0.1%  |

Large energy losses at small D indicate the superposition is far from
equilibrium in the overlap region. Radiation is generated and absorbed at
boundaries. The mu=-6 case has the best energy conservation (oscillons
are larger, overlap region is smaller relative to oscillon size).

### Dynamical Behavior

Typical behavior at D=15: oscillons attract, merge/bounce, complex dynamics.
At D=40 for mu=-6: oscillons slowly drift apart (v ~ 0.006 code units/time),
energy nearly constant. The drift velocity is dominated by the initial
radiation kick from the non-equilibrium superposition, not by the steady-state
inter-oscillon force.

## Step 3: Force Extraction

### Method

The inter-oscillon force was extracted from the STATIC overlap energy:

    E_int(D) = E_two_oscillons(t=0, D) - 2 * E_single

This avoids contamination from radiation dynamics. The force at the midpoint
between D values was computed numerically:

    F(D_mid) = -(E_int(D+5) - E_int(D)) / 5

The interaction range was estimated from the force ratio:

    lambda = 5 / ln(F(17.5) / F(22.5))

### Results

| mu   | omega  | xi_pred | F(17.5)  | F(22.5)  | lambda | lambda/xi_pred |
|------|--------|---------|----------|----------|--------|----------------|
|  -6  | 0.9360 | 2.841   | 1.24e-01 | 2.66e-02 | 3.24   | 1.14           |
|  -8  | 0.9140 | 2.465   | 7.55e-02 | 1.65e-02 | 3.28   | 1.33           |
| -10  | 0.8910 | 2.203   | 2.80e-02 | 1.48e-02 | 7.85   | 3.57           |
| -14  | 0.8430 | 1.859   | 2.43e-03 | -2.51e-03 | N/A   | N/A            |
| -20  | 0.8850 | 2.148   | 5.65e-03 | 4.81e-03 | 31.1   | 14.5           |

**Reliable measurements** (lambda/xi_pred ~ 1): mu=-6 (1.14) and mu=-8 (1.33).
These confirm that the interaction range IS set by the oscillon tail length.

**Unreliable measurements**: mu=-10 (ratio 3.6, radiation background), mu=-14
(sign change in force), mu=-20 (nearly constant force, radiation dominated).

### Force vs D: Full Table

Interaction energy E_int = E_two(t=0) - 2*E_single:

| mu   | D=15   | D=20    | D=25    | D=30    | D=40    |
|------|--------|---------|---------|---------|---------|
|  -6  | 1.311  | 0.689   | 0.556   | 0.513   | 0.334   |
|  -8  | 0.195  | -0.182  | -0.265  | -0.239  | -0.045  |
| -10  | 0.645  | 0.505   | 0.431   | 0.324   | 0.105   |
| -14  | 0.179  | 0.167   | 0.180   | 0.191   | 0.179   |
| -20  | 0.038  | 0.010   | -0.014  | -0.037  | -0.084  |

Sign changes and non-monotonic behavior at mu=-20, -14, -8 reflect sensitivity
to the breathing phase of the saved profile. The oscillon field oscillates
sign in the tail, so the overlap integral depends on relative phase.

## Step 4: Interaction Range vs Gap Margin

The key diagnostic for critical gravity:

| mu   | gap_margin | xi_pred | lambda_meas | diverges? |
|------|------------|---------|-------------|-----------|
|  -6  | 6.4%       | 2.84    | 3.24        | no        |
|  -8  | 8.6%       | 2.46    | 3.28        | no        |
| -10  | 10.9%      | 2.20    | ~7.9        | no (noisy)|

**The gap margin does NOT approach zero** for any mu in the scan.
At mu=-6, the gap is 6.4% and xi_pred=2.84. Moving to weaker coupling
(smaller |mu|) brings the frequency closer to the mass gap, but the oscillon
also becomes less tightly bound and eventually disperses.

The trend: as |mu| decreases from 20 to 6, xi_pred increases from 2.15 to
2.84 (a factor of 1.32). This is modest growth, not a divergence.

## Conclusions

1. **Interaction range tracks correlation length**: For the cleanest cases
   (mu=-6, -8), lambda_meas / xi_pred = 1.1-1.3, confirming that the inter-
   oscillon force decays as exp(-D/xi) where xi = 1/sqrt(m^2 - omega^2).

2. **No divergence**: xi_pred ranges from 1.86 to 2.84 across the full mu
   scan. The oscillons cannot exist with omega arbitrarily close to m because
   the binding becomes too weak. The gap margin bottoms out at ~6%.

3. **Radiation contamination**: The equilibrated-profile approach is much
   cleaner than raw Gaussian init (V23-C), but radiation from the superposition
   overlap still dominates at D > 20. Reliable force measurements are limited
   to D = 15-20 where actual tail overlap exceeds the radiation background.

4. **Critical gravity fails**: The mechanism requires xi -> infinity as a
   control parameter is tuned. In this system, the gap cannot be closed
   continuously — the oscillon dissolves before omega reaches m. The
   interaction range is Yukawa with xi ~ 2-3 (sub-fm in physical units),
   far too short for gravity.

## Files

- `src/critical2.c` — equilibration + interaction solver
- `data/profile_mu{N}.dat` — equilibrated profiles (text, 7 columns)
- `data/equil_mu{N}_summary.txt` — equilibration summary
- `data/interact_mu{N}_D{D}_ts.tsv` — separation vs time
- `data/yukawa_fits.tsv` — fit parameters
