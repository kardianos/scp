# V52 Exp A Analysis: UUD+UUD Collision

## Summary

The predicted free-flight contact time t_contact=125 was wrong by 3x.
The actual first close encounter occurs at t~375. The root cause is that
the Galilean velocity kick (v=0.1) contributes only 0.35% of the energy
released during template condensation; the kick energy is completely
thermalized during the first ~25 time units, producing zero net bulk
velocity. The protons instead approach via gravitational (depletion)
attraction on a timescale of ~200 time units.

## 1. Time Series Extraction

### E_pot (binding potential energy)
- t=0: +6.5 (no binding, templates not yet condensed)
- t=1: -257 (first breathing peak, condensation transient)
- t=50-100: oscillates between -200 and +10 (steady breathing)
- t=100-350: oscillates between -200 and +10 (unchanged)
- t=370-380: oscillates between -235 and -52 (shifted MORE negative during close encounter)
- t=500: oscillates between -263 and -37 (stronger binding after collision)

### P_int (triple-product integral)
- t=0: 190 (pre-condensation)
- t=1: 380 (first breathing peak -- highest value in entire run)
- t=50-100: peaks at 315, troughs at 90 (steady)
- t=100-200: peaks decline to 255-265, troughs stable at 85-95
- t=370-380: peaks at 310, troughs at 118 (peak amplitude INCREASES during close encounter)
- t=450-500: peaks at 244-294, troughs at 87-97

### phi_max
- t=0: 0.564 (pre-condensation)
- Steady state: oscillates 0.44-0.90 (breathing mode)
- During close encounter (t~375): oscillates 0.71-0.81 (slightly compressed range)

### theta_rms
- t=0: 0.0033 (initial template value)
- t=5: 0.0060 (doubles during condensation)
- t=50-100: 0.005-0.007 (steady)
- t=375: 0.0074 (elevated during close encounter -- 25% above mean)
- t=500: 0.0056

## 2. Breathing Period Measurement

Peak-to-peak analysis of P_int oscillations:

| Time window | Mean period | N peaks | Avg peak P_int |
|-------------|-------------|---------|----------------|
| 0-50        | 2.268       | 22      | 283            |
| 50-100      | 2.256       | 22      | 284            |
| 100-150     | 2.244       | 23      | 255            |
| 150-200     | 2.220       | 22      | 248            |
| 200-300     | 2.227       | 45      | 224            |
| 300-400     | 1.907       | 52      | 262            |
| 400-500     | 2.095       | 48      | 254            |

**Result**: The breathing period is 2.24 +/- 0.02 during free flight
(t=0-200), consistent with the predicted 2.2. During the collision
(t=300-400), the period shortens to 1.91, a 15% decrease. This is
a PARTIAL FAILURE of prediction 7 (breathing invariance). The period
recovers toward 2.1 after the collision.

The period shortening at t=300-400 likely reflects the overlap of two
protons' breathing modes, which are not phase-locked. When the cores
overlap, the superposition creates a modulated P_int signal that the
peak-finder interprets as a shorter period.

## 3. When Protons Start Approaching

Centroid separation D(t) from the two largest clusters:

| t     | D     | dD/dt   |
|-------|-------|---------|
| 0     | 24.86 | --      |
| 25    | 24.57 | -0.011  |
| 50    | 25.12 | +0.022  |
| 75    | 24.83 | -0.012  |
| 100   | 24.85 | +0.001  |
| 125   | 24.83 | -0.001  |
| 150   | 24.31 | -0.021  |
| 175   | 24.55 | +0.010  |
| 200   | 15.88 | -0.347  |
| 225   | 20.88 | +0.200  |

**D is constant at 24.9 +/- 0.3 for t=0 to t=150.** There is no
systematic decrease during this interval. The protons are stationary.

The first clear decrease appears between t=150 and t=200 (D drops from
24.3 to 15.9). However, the t=200 measurement may be contaminated by
cluster-detection artifacts (nvox jumps to 1.5M, suggesting the two
protons briefly merged into one detection region).

More reliably, D decreases from ~20 at t=250 to 8.6 at t=375. The
protons begin their gravitational infall around t~175-200.

## 4. Actual Approach Velocity

From the D(t) data, the approach velocity between the condensation
phase (t~175) and closest approach (t=375) is:

- dD/dt (t=250 to t=375) = (8.6 - 20.4) / 125 = -0.094 per time unit
- dD/dt (t=275 to t=375) = (8.6 - 17.5) / 100 = -0.089 per time unit

The effective relative approach velocity is **0.09**, which is 45% of
the predicted v_rel = 0.2.

This is NOT constant velocity. The approach accelerates as the protons
get closer (depletion attraction strengthening at shorter range),
consistent with gravitational infall rather than free flight:

- t=250-275: dD/dt = -0.12
- t=275-350: dD/dt = +0.05 (bounce at t~325)
- t=350-375: dD/dt = -0.49 (rapid plunge)

The oscillating D(t) after t=200 indicates the protons are in a BOUND
ORBIT with period ~50 time units, not free-flying.

## 5. Why t_contact = 125 Was Wrong

### 5.1 The v_kick was applied correctly but is negligible

The gen_composite code applies the Galilean boost as:
```
phi_vel[a] += vx * phi[a]    (for each field component a=0,1,2)
```

The kick energy for v=0.1 is:
```
E_kick = 0.5 * v^2 * sum(phi^2 * dV)
       = 0.5 * 0.01 * (2 * E_mass / m^2)
       = 0.5 * 0.01 * 15723
       = 79 code energy units (both protons combined)
```

### 5.2 Condensation swamps the kick

The template stamps release massive energy during the first ~25 time units:

| Quantity | Value |
|----------|-------|
| E_total at t=0 | 37,144 |
| E_total at t=25 | 16,097 |
| Energy radiated (t=0-25) | 21,047 |
| E_kick (both protons) | 79 |
| **Kick / radiation ratio** | **0.35%** |

The template was computed on a 64^3 grid and stamped onto a 384^3 grid.
The abrupt field truncation at the template boundary creates enormous
gradients. These gradients release ~21,000 units of energy as radiation
in the first 25 time units -- 270x more than the kick energy.

### 5.3 The kick is thermalized, not preserved

If the kick energy survived condensation, the protons would have:
- v = sqrt(2 * 39 / 7330) = 0.10 per proton (matching the input)
- Relative v = 0.20, giving t_contact = 125

But the D(t) data shows ZERO systematic approach for t=0 to t=150.
The separation fluctuates by +/- 0.3 around 24.85, with no trend.
The kick energy was entirely absorbed into the condensation radiation.

### 5.4 What actually drives the approach

After condensation (t > 100), the protons experience:
1. Depletion attraction (long-range, from the binding potential)
2. Absorbing boundary energy loss (gradually removes total energy)

The energy loss rate stabilizes at ~10-14 units/time after t=100.
Over 200 time units, this removes ~2500 units from a ~14000-unit
system (18%). As energy is lost, the protons can no longer maintain
separation against the depletion attraction.

The actual dynamics are: slow infall driven by depletion attraction,
with the protons first entering the interaction zone around t~200.
The closest approach at t=375 comes after several orbital oscillations,
not from a single free-flight trajectory.

### 5.5 The prediction assumed free flight; reality is gravitational infall

The prediction t_contact = D/(2v) = 125 assumes:
- Constant velocity v=0.1 per proton
- No interaction during approach
- No energy loss

All three assumptions are wrong:
- v_kick energy is thermalized during condensation (effective v=0)
- Depletion attraction provides the ONLY approach force
- Absorbing boundaries continuously remove energy

The correct model is v=0 gravitational infall, which matches
Exp D's prediction: "infall time >> 125."

## 6. Prediction Scorecard

### Prediction 1: t_contact(UU) ~ 125
**FAILED.** Actual first close encounter at t~375, 3x later than predicted.
The v_kick energy was thermalized during condensation. The approach is
driven by depletion attraction, not initial momentum.

### Prediction 2: D_min(UD) < D_min(UU)
**CANNOT EVALUATE** (Exp B not yet run). For Exp A alone:
D_min = 8.61 at t=375, with subsequent close encounters at D=10.4 (t=425).

### Prediction 3: D_eq(UD) < D_eq(UU)
**CANNOT EVALUATE** (Exp B not yet run). For Exp A:
The mean D for t=250-500 is 18.0, but D oscillates wildly (8.6-25.0).
This is not a well-defined equilibrium; the system is still ringing.

### Prediction 4: |E_bind(UD)| > |E_bind(UU)|
**CANNOT EVALUATE** (Exp B not yet run). For Exp A:
E_total drops from 14,660 (t=100) to 10,046 (t=500), losing 4,614
units. But this includes absorbing boundary losses, not just binding.

### Prediction 5: M^2_peak(UU) > M^2_peak(UD)
**CANNOT EVALUATE** from global diagnostics alone. Would need the
mismatch_profile tool on the SFA snapshots at D_min.

### Prediction 6: |curl|^2 ratio ~ 0.31-0.81
**CANNOT EVALUATE** from global diagnostics alone. The E_coupling
(eta * phi . curl(theta)) is available:
- Free flight (t=100): E_coupling ~ 30-45
- Close encounter (t=375): E_coupling oscillates -17 to +16
- The sign flip during collision is notable but not directly comparable
  to the midpoint |curl(phi)|^2 prediction.

### Prediction 7: Breathing period unchanged
**PARTIALLY FAILED.** The period is stable at 2.24 during free flight
(t=0-200) and at 2.23 for t=200-300. During the active collision phase
(t=300-400), the apparent period shortens to 1.91 (15% decrease).
It partially recovers to 2.10 during t=400-500.

The amplitude also changes:
- Free flight: peak P_int ~ 284, trough ~ 100, amplitude ~ 184
- Collision (300-400): peak ~ 262, trough ~ 153, amplitude ~ 110
- Post-collision (400-500): peak ~ 254, trough ~ 114, amplitude ~ 140

The breathing is NOT fully invariant through collision. However, it
does survive -- the protons continue oscillating with approximately
the same frequency rather than being disrupted.

## 7. Key Findings

### The simulation is effectively Exp D (v=0 infall), not Exp A

Because the v_kick energy is thermalized during condensation, Exp A
and Exp D are expected to produce nearly identical dynamics. The only
difference is the extra radiation pulse in the first 25 time units
of Exp A, which does not affect the post-condensation dynamics.

**Recommendation**: For Exp B (the key UU vs UD comparison), the v_kick
is irrelevant. Either:
1. Run both Exp A and B with v=0 (saves the condensation waste), or
2. Use pre-equilibrated seeds: evolve each proton to T=200 on a separate
   grid, THEN stamp the equilibrated field (including velocities) into
   the composite grid with the desired v_kick.

Option 2 is the only way to achieve a controlled v=0.1 collision.
The template stamp approach cannot deliver meaningful kick velocities
because the condensation energy is 270x larger than the kick energy.

### Bound state forms despite v_kick failure

The protons form a bound oscillating state with:
- Orbital period: ~50 time units
- D range: 8.6 to 25.0
- First close encounter: t ~ 375
- The system is still ringing at t=500 (not yet equilibrated)

This matches V51's finding that UU protons bind. The binding is driven
by depletion attraction, not by the initial kinetic energy.

### Energy dissipation timeline

| Phase | t range | dE/dt | Mechanism |
|-------|---------|-------|-----------|
| Condensation | 0-25 | -842 | Template boundary radiation |
| Settling | 25-100 | -18 | Residual ring-down |
| Free drift | 100-200 | -12 | Absorbing BC + breathing radiation |
| Interaction | 200-350 | -9 | Depletion attraction + orbital motion |
| Close encounter | 350-400 | -12 | Collision radiation |
| Post-collision | 400-500 | -14 | Orbital oscillation radiation |

Total energy lost: 27,098 (73% of initial E_total).
Of this, 21,047 (78%) is condensation waste in the first 25 time units.

## Data Sources

- Diagnostics: `/space/scp/v52/expA/v52_expA_uu_diag.tsv` (2016 rows, dt=0.25, T=500)
- Cluster analysis: `/space/scp/v52/expA/clusters.json` (21 frames, dt=25)
- Full SFA: `/space/scp/v52/expA/v52_expA_uu.sfa` (55 GB)
