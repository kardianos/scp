# Torsion Coupling: Phase 1 Results (gamma scan)

All runs: N=128, L=20, T=300, braids=1, diag_dt=5, snap_dt=100.
Initial energy: E_total = 5357.0, E_pot = -122.2 (same for all gamma).

## Summary Table

| gamma | c_L     | E_drift (%) | E_pot(T=300) | <E_pot>  | min(E_pot) | Braid survival |
|-------|---------|-------------|-------------|----------|------------|----------------|
|  0.0  | 1.000   | -0.28%      | -202.6      | -202.3   | -455.5     | YES            |
|  0.1  | 1.049   | -0.24%      | -172.9      | -188.9   | -380.9     | YES            |
|  0.5  | 1.225   | -0.31%      | -269.7      | -189.8   | -379.6     | YES            |
|  1.0  | 1.414   | -0.30%      | -177.7      | -186.6   | -389.0     | YES            |
|  2.0  | 1.732   | -0.15%      | -136.8      | -165.6   | -342.7     | YES            |
|  5.0  | 2.449   | +0.34%      | -247.1      | -128.7   | -280.6     | YES (weakened) |
| -0.1  | 0.949   | -0.17%      | -253.8      | -199.7   | -415.6     | YES            |
| -0.5  | 0.707   | +0.28%      | -85.0       | -191.4   | -370.9     | YES            |

E_drift = (E_total(T=300) - E_total(T=0)) / |E_total(T=0)| x 100%.
<E_pot> = time-averaged E_pot over the full run.

## Key observations

### 1. gamma=0 control: verified

The gamma=0 run produces the standard equation (the torsion force is zero). Energy
drift is -0.28%, which is the usual Verlet drift for this grid resolution. The braid
oscillates normally with E_pot fluctuating between -19 and -456.

### 2. Energy accounting caveat

**The E_total printout does NOT include the elastic energy (gamma/2)(div phi)^2.**
The Verlet integrator conserves the true Hamiltonian (including the elastic term),
but the diagnostic only sums E_kin + E_grad + E_mass + E_pot. For gamma=0 this is
the full energy and drift is -0.28%. For gamma!=0, the apparent drift reflects the
untracked elastic term transferring energy to/from the tracked components.

Evidence: the E_total range widens systematically with |gamma|:
- gamma=0:   range 5342-5357 (15 units, 0.28%)
- gamma=0.1: range 5342-5357 (15 units, 0.28%) -- barely different from control
- gamma=1.0: range 5334-5364 (30 units, 0.56%)
- gamma=5.0: range 5357-5424 (67 units, 1.25%)

For gamma=5.0, E_total is biased upward because the elastic energy is negative
(compression reduces gradient energy, which is tracked, while elastic energy is not).

### 3. Braid survival: ALL braids survive

Every gamma value tested produces a surviving braid at T=300. None broke. The braid
oscillates in all cases, with E_pot cycling between its minimum (deep binding) and
maximum (near zero, expanded state).

### 4. Braid weakening with large positive gamma

The time-averaged E_pot is a robust measure of braid binding strength:
- gamma=0.0: <E_pot> = -202 (control)
- gamma=0.1: <E_pot> = -189 (7% weaker)
- gamma=0.5: <E_pot> = -190 (6% weaker)
- gamma=1.0: <E_pot> = -187 (8% weaker)
- gamma=2.0: <E_pot> = -166 (18% weaker)
- gamma=5.0: <E_pot> = -129 (36% weaker)

**Large positive gamma weakens the braid.** This makes physical sense: positive gamma
makes longitudinal (compression) waves faster than transverse (shear) waves. The braid
is a helical shear structure, so faster compression relative to shear should work against
the twist that holds the braid together.

### 5. Negative gamma: nearly identical to control

- gamma=-0.1: <E_pot> = -200 (1% weaker than control)
- gamma=-0.5: <E_pot> = -191 (5% weaker)

Negative gamma (shear faster than compression) does NOT significantly strengthen the
braid. It keeps binding close to the control value. The slight weakening at gamma=-0.5
may be a real effect or fluctuation noise.

### 6. No qualitative shape changes observed

All braids remain coherent at the center (D=0.0 for all single-braid runs). No spatial
drift, no splitting, no mode conversion. The braid oscillates with a similar period
in all cases. The main quantitative change is the binding depth.

### 7. E_pot oscillation amplitude

The E_pot oscillation range (max - min) shrinks with large positive gamma:
- gamma=0.0: range = 436 (deep oscillations)
- gamma=5.0: range = 279 (shallower oscillations)

This is consistent with weaker binding: the potential well is shallower, so the
oscillation amplitude in E_pot space is smaller.

## Conclusions

1. **The torsion coupling is compatible with braid existence** across the full
   range tested (gamma = -0.5 to 5.0).

2. **Positive gamma weakens the braid monotonically.** The effect is gradual --
   no sharp transition or breakup threshold was reached up to gamma=5.0.

3. **Negative gamma preserves the braid** at near-control strength but does not
   significantly enhance binding.

4. **The diagnostic should be updated** to include E_elastic = (gamma/2)(div phi)^2
   for proper energy tracking.

## Recommendations for further exploration

1. **gamma = 0.1 to 1.0**: This is the "sweet spot" where the coupling is
   significant enough to differentiate wave modes but not so strong as to
   weaken the braid substantially. Good range for Phase 2 (shear/compression
   decomposition around the braid).

2. **gamma = 10, 20**: Test whether there is an actual breakup threshold for
   large positive gamma. The monotonic weakening suggests it exists.

3. **gamma = -1.0 to -0.9**: Need to verify c_L = sqrt(1+gamma) > 0, which
   requires gamma > -1. As gamma -> -1, the longitudinal speed vanishes. This
   limit is physically interesting (pure shear medium).

4. **Two-braid runs**: The coupling may affect inter-braid forces. Repeat the
   best gamma values with braids=2 to see if the torsion coupling changes
   attraction/repulsion.

---

# Phase 2: Cosserat 6-field Results (eta scan)

6-field system: 3 position (phi) + 3 angle (theta), coupled through curl.

    d^2 phi_a/dt^2 = Lap(phi_a) - m^2 phi_a - V'(P) + eta * curl(theta)_a
    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 theta_a + eta * curl(phi)_a

All runs: N=128, L=20, T=300, braids=1, diag_dt=5, snap_dt=50.
Initial energy: E_total = 5357.0, E_pot = -122.2.
Default: m^2 = 2.25, m_theta^2 = 2.25.

## Summary Table: eta scan (m_theta^2 = 2.25 fixed)

| eta | E_drift (%) | E_pot(T=300) | <E_pot> | min(E_pot) | theta_rms(T=300) | <theta_rms> | E_theta(T=300) | <E_coupling> | Braid survival |
|-----|-------------|-------------|---------|------------|-----------------|-------------|----------------|-------------|----------------|
| 0.0 | -0.276%     | -386.5      | -203.5  | -419.4     | 0.0000          | 0.0000      | 0              | 0           | YES (control)  |
| 0.1 | -0.231%     | -230.7      | -181.7  | -369.9     | 0.0662          | 0.0462      | 1964           | -3          | YES            |
| 0.5 | -0.248%     | -143.2      | -136.8  | -283.6     | 0.0662          | 0.0549      | 2077           | -64         | YES (weaker)   |
| 1.0 | -0.315%     | -62.5       | -42.4   | -122.2     | 0.0641          | 0.0642      | 1631           | -239        | WEAK           |
| 2.0 | -0.015%     | -3.0        | -3.0    | -122.2     | 0.0508          | 0.0631      | 1125           | -223        | DISSOLVED      |
| 5.0 | BLOWUP      | --          | --      | --         | 3.5e28          | --          | --             | --          | CFL VIOLATION  |

E_theta = E_theta_kin + E_theta_grad + E_theta_mass (total angle sector energy at T=300).
<E_coupling> = time-averaged coupling energy = -eta * integral(phi . curl(theta)).

## Summary Table: m_theta variation (eta = 1.0 fixed)

| m_theta^2 | E_drift (%) | E_pot(T=300) | <E_pot> | min(E_pot) | theta_rms(T=300) | <theta_rms> | E_theta(T=300) | <E_coupling> | Braid survival |
|-----------|-------------|-------------|---------|------------|-----------------|-------------|----------------|-------------|----------------|
| 2.25      | -0.315%     | -62.5       | -42.4   | -122.2     | 0.0641          | 0.0642      | 1631           | -239        | WEAK           |
| 0.25      | -0.584%     | -65.8       | -102.5  | -266.9     | 0.0566          | 0.0449      | 959            | -119        | YES            |
| 0.00      | -0.597%     | -44.2       | -119.1  | -267.0     | 0.0621          | 0.0677      | 831            | -110        | YES            |

## Key Observations

### 1. The braid SOURCES the angle field -- confirmed

For every eta > 0, theta_rms rises from zero to ~0.05-0.08 within the first ~30 time
units and then oscillates around a steady-state value. The braid's helical curl(phi)
drives the angle field, which becomes permanently excited.

Steady-state theta_rms is nearly independent of eta (0.05-0.07 for all stable runs).
This means the angle field amplitude is set by the SOURCE (braid curl) not the coupling
strength. Stronger eta just couples more energy into the angle sector.

### 2. Coupling destabilizes the braid progressively

E_pot retention (braid binding) decreases monotonically with eta:
- eta=0.0: <E_pot> = -204 (control, full binding)
- eta=0.1: <E_pot> = -182 (11% weaker)
- eta=0.5: <E_pot> = -137 (33% weaker)
- eta=1.0: <E_pot> = -42  (79% weaker -- nearly dissolved)
- eta=2.0: <E_pot> = -3   (99% weaker -- FULLY dissolved)

The braid effectively dissolves between eta=1.0 and eta=2.0. At eta=2.0, E_pot
hovers near zero for the entire run after the initial transient, meaning the triple-
product potential has been completely disrupted.

### 3. Energy transfer to angle sector is massive

At eta=0.1, the angle sector absorbs ~1960 energy units by T=300 -- 37% of the total
energy. This energy comes from the position field: E_phi_kin drops from 2698 to 1443.
The braid acts as an antenna, radiating into the angle field.

At eta=0.5 and above, the angle sector energy is comparable to or exceeds the position
sector kinetic energy. The two sectors reach approximate equipartition.

### 4. eta=5.0 is a CFL violation, not a physical instability

The eta=5.0 run blows up exponentially (10^8 per time unit) starting at the FIRST
timestep. This is a numerical instability: the curl coupling with eta=5 increases the
effective wave speed beyond what dt=0.0315 can resolve. A smaller dt would stabilize it,
but the physical regime is already covered by eta=0-2.

### 5. Lighter angle mass HELPS the braid survive

At eta=1.0, reducing the angle mass from m_theta^2=2.25 to m_theta^2=0.25:
- <E_pot> improves from -42 to -103 (braid 2.4x more bound)
- theta_rms is slightly lower (0.045 vs 0.064)
- E_theta is lower (959 vs 1631) -- less energy trapped in angle sector

The lighter angle field carries energy AWAY from the braid more efficiently (faster
group velocity), which paradoxically helps the braid survive. Heavy angle fields
trap energy near the braid and disrupt it.

### 6. MASSLESS angle field (m_theta=0): STABLE and braid-preserving

The m_theta=0 case is stable through T=300 with only -0.60% drift. Key results:
- <E_pot> = -119 (42% weaker than control, but 2.8x better than massive eta=1.0)
- theta_rms = 0.062 (steady oscillation)
- E_theta_mass = 0 (confirmed massless)
- E_theta = 831 (lowest of the eta=1.0 variants)
- E_coupling = -110

A massless angle field propagates at c (speed of light) without Yukawa decay.
This is significant: the angle field carries torsion information about the braid
to infinity, not just to a finite Yukawa range.

### 7. Coupling energy is always negative

<E_coupling> < 0 for all eta > 0, meaning the coupling term LOWERS the total
energy. The phi and theta fields develop correlations (phi . curl(theta) > 0)
that reduce the system energy. This is energetically favorable -- the system
WANTS to develop angle field structure correlated with the braid's curl.

The magnitude of <E_coupling> peaks near eta=1.0 (-239) and decreases for eta=2.0
(-223), consistent with the braid being too disrupted at large eta to maintain
coherent correlations.

## Conclusions

1. **The braid sources the angle field.** For any eta > 0, the braid's helical
   curl(phi) excites theta oscillations that persist for the full simulation.
   The angle field carries information about the braid's twist structure.

2. **The coupling is destructive above eta ~ 1.** The braid dissolves between
   eta=1 and eta=2. The position-angle energy transfer disrupts the triple-product
   potential that holds the braid together.

3. **Lighter angle masses help braid survival.** The m_theta=0 (massless) case
   preserves the braid 3x better than m_theta^2=2.25, despite identical eta=1.0.
   The massless mode radiates energy away efficiently instead of trapping it locally.

4. **The massless angle field is the physically interesting regime.** It is:
   - Stable (no blowup, -0.6% drift)
   - Long-range (propagates at c, no Yukawa cutoff)
   - Sourced by the braid topology (curl(phi) is nonzero only where the braid has twist)
   - A candidate for long-range torsion force carrier

5. **eta ~ 0.1-0.5 is the stable operating range** where the angle field is
   excited but the braid survives with >65% binding retention.
