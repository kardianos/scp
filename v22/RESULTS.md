# V22 Results: Two-Oscillon Gravitational Interaction

## Summary

Two v21 oscillons placed at separation D=20 on the z-axis. Control (no gravity)
shows strong direct field-tail attraction leading to merger at t~409. Gravity run
(alpha=0.1, beta=0.1) shows monotonic separation increase — the gravity-enhanced
mass gap compactifies the oscillons, suppressing the tail overlap that drives
the direct attraction. The scalar gravity mediator's 1/r^2 force is negligible
compared to the exponential tail-coupling force at D=20.

## Run Parameters

| Parameter | Gravity Run | Control Run |
|-----------|-------------|-------------|
| N | 400 | 350 |
| L | 60.0 | 60.0 |
| dx | 0.301 | 0.344 |
| dt | 0.075 | 0.086 |
| pts/sigma | 10.0 | 8.7 |
| Memory | 7.2 GB | 4.8 GB |
| alpha, beta | 0.1, 0.1 | 0.0, 0.0 |
| D (initial sep) | 20.0 | 20.0 |
| tfinal | 500 | 500 |
| Wall time | 371 min | 234 min |
| CFL | 0.25 | 0.25 |

Common: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0.

## Separation Trajectories

```
t       Control(sep)  Gravity(sep)  Delta
  0       20.00         19.99       -0.01
 25       20.20         20.38       +0.18
 50       21.39         21.60       +0.21   (shedding peak)
 75       21.09         21.74       +0.65
100       20.64         22.02       +1.38
125       20.18         22.40       +2.22   (post-shedding)
150       19.86         22.76       +2.90
175       19.60         23.47       +3.87
200       19.29         24.30       +5.01
225       18.85         —           —
250       18.32         25.10       +6.78
275       17.58         25.94       +8.36
300       16.72         26.94      +10.22
325       15.61         27.98      +12.37
350       13.89         29.03      +15.14
375       11.22         30.22      +19.00
400        8.41         31.51      +23.10
410        4.96         —           —       (control MERGER)
415        3.81         —           —       (control minimum sep)
425        4.59         —           —       (control BOUNCE)
450       10.84         33.73      +22.89
475       14.03         34.99      +20.96
500       14.16         36.21      +22.05   (FINAL)
```

## Control Run (No Gravity): Detailed Analysis

### Phases

1. **Shedding (t=0-50)**: Both oscillons shed radiation. Separation increases
   from 20.0 to 21.4 (radiation pressure pushes centroids outward).

2. **Direct attraction (t=50-400)**: Separation decreases monotonically from
   21.4 to 8.4. Rate accelerates from -0.03/t.u. to -0.5/t.u. as tail overlap
   grows exponentially with decreasing distance.

3. **Merger (t=400-420)**: Rapid infall. Sep reaches minimum 3.81 at t~415.
   Cores deeply overlapping. fc drops to 0.43 (heavily disrupted).

4. **Bounce (t=420-440)**: Oscillons re-separate to sep~10. The merger is
   partially elastic — the combined configuration is not a stable single oscillon.

5. **Post-bounce (t=440-500)**: Two disrupted blobs at sep~14, fc~0.33.
   Energy dropping rapidly (176 -> 145, losing 18% to radiation).
   Not clear if they will re-merge or stabilize.

### Final state
- sep = 14.16, fc = 0.334 (severely disrupted)
- E = 145.3 (46% of initial 268.6 radiated)
- omega = 0.948 (still below mass gap)

### Interpretation: Direct Field-Tail Coupling

The triple-product interaction V = (mu/2)P^2/(1+kappa*P^2) creates an
attractive force between overlapping exponential tails. Each oscillon's
field falls as phi ~ exp(-sqrt(m^2 - omega^2) * r) at large r. At D=20
with m=1.0, omega~0.95, the decay length is 1/sqrt(1-0.9025) = 3.2.
The tail overlap at the midpoint is ~ exp(-10/3.2) = 0.044 — small but
nonzero. The coupling V(P) amplifies this through the triple product,
creating a force that grows exponentially as separation decreases.

## Gravity Run (alpha=0.1): Detailed Analysis

### Phases

1. **Shedding (t=0-55)**: Similar to control. Sep peaks at 21.6.

2. **Steady expansion (t=55-500)**: Sep increases monotonically from 21.6
   to 36.2. Rate is ~0.04/t.u., slowly increasing.

3. **No merger**: Oscillons remain well-separated throughout.

### Final state
- sep = 36.21 (16.2 above initial D=20)
- fc = 0.590 (declining — oscillons approaching absorbing boundary)
- E = 180.1 per pair (90.0 each)
- omega = 0.960 (below mass gap)
- Phi_mid = +0.06 (oscillating)

### Why Gravity Causes Expansion (Not Attraction)

The expected gravitational attraction was F ~ alpha*E/(4*pi*r^2) ~ 0.007,
giving acceleration a ~ 8e-5 and displacement dz ~ 5 over 350 t.u.

Instead, the dominant effect is **tail compactification**:
- m_eff^2 = m0^2 - beta*Phi, with Phi < 0 near each oscillon
- m_eff > m0 -> deeper mass gap -> steeper exponential tail decay
- Shorter tails -> less overlap -> weaker direct attraction
- The direct attraction that drives merger in the control is suppressed

The gravitational 1/r^2 attraction (F~0.007) is negligible compared to the
exponential tail force that it suppresses. The net effect is repulsive.

Additionally, the expansion rate (~0.04/t.u.) corresponds to a constant
velocity v ~ 0.04, not acceleration. This suggests an initial "kick" from
asymmetric shedding (gravity modifies the shedding radiation pattern) rather
than a sustained force. The oscillons may be coasting apart on their
shedding momentum.

### Late-time fc decline (t>400)

fc drops from 0.96 to 0.59 at t=500. With sep=36.2, each oscillon's
centroid is at z=+/-18.1. The absorbing boundary starts at 70% of L/2 = 21.
The oscillons are approaching the damping zone (only 3 units away), which
explains the declining fc. The box is becoming too small for this separation.

## Key Findings

### 1. Direct field attraction dominates at D=20

The triple-product tail coupling is the dominant inter-oscillon force at
this separation. It is attractive and grows exponentially as separation
decreases. Without gravity, this leads to merger in ~400 t.u.

### 2. Scalar gravity SUPPRESSES inter-oscillon attraction

At alpha=0.1, the gravity-enhanced mass gap makes oscillons more compact.
The reduced tail overlap more than compensates for the gravitational
attraction. The net effect is separation, not attraction.

### 3. Merger dynamics are partially elastic

The control merger produces a bounce at t~415 (minimum sep=3.8) followed
by re-separation. The post-merger objects are severely disrupted (fc=0.33)
and losing energy rapidly. This is qualitatively similar to soliton-
antisoliton scattering: inelastic pass-through with significant radiation.

### 4. The gravitational signal is unobservable at alpha=0.1

The predicted gravitational displacement (~5 units) is swamped by the
direct field coupling (which moves oscillons by ~7+ units in the control).
To isolate gravity, need either:
- Much larger initial separation (D >> 20, so tails don't overlap)
- Much smaller coupling (mu, kappa -> 0 while keeping oscillon alive)
- Heavier oscillons with shorter tails (larger m-omega gap)

## Comparison with PLAN Predictions

| Quantity | Predicted | Observed |
|----------|-----------|----------|
| Grav force | 0.0068 | Unobservable (swamped) |
| Grav displacement (350 t.u.) | ~5 | N/A (competing effects) |
| Shedding duration | ~150 t.u. | ~50-100 t.u. |
| Post-shedding E/oscillon | ~85 | ~90 (gravity), ~89 (control) |
| Direct attraction | Not predicted | DOMINANT effect |
| Merger | Not predicted | Occurs at t~409 (control) |

## Data Files

- `data/production_grav.txt` — gravity run log (N=400, alpha=0.1)
- `data/control_nograv_n350.txt` — control run log (N=350, alpha=0)
- `data/control_nograv.txt` — killed N=200 run (INVALID, insufficient resolution)
- `data/gravity_ts.tsv` — gravity per-step time series
- `data/control_ts.tsv` — CORRUPTED (mixed gravity+control data from shared file)
- `data/two_osc_profile_t*.tsv` — z-axis profiles (from last run to write each t)

Note: Both runs wrote to the same `two_oscillon_ts.tsv` file simultaneously,
corrupting it. The log files are the authoritative data source.

## Phase Offset Test (N=200, L=40, D=20, t=200, no gravity)

Two independently formed oscillons generically have different breathing phases.
The sign and magnitude of the direct tail-coupling force depends on relative phase.

### Setup

Upper oscillon at phase 0: phi = A*f(r), v = 0
Lower oscillon at phase delta: phi = A*cos(delta)*f(r), v = -A*omega*sin(delta)*f(r)

### Results at t=200

| Phase offset | Final sep | Delta_sep | Direction | Rate |
|---|---|---|---|---|
| delta=0 (in-phase) | 19.40 | -0.60 | Attracting | -0.004/t.u. |
| delta=pi (anti-phase) | 21.06 | +1.06 | **Repelling** | +0.006/t.u. |
| delta=pi/2 (quadrature) | 20.59 | +0.59 | **Repelling** | +0.004/t.u. |

### Interpretation

- **In-phase**: Tail overlap is coherent (P contributions add). Time-averaged
  V(P) creates attractive well. This was the (unrealistic) assumption of the
  original production runs.

- **Anti-phase**: At any given instant, when one oscillon's fields are at
  maximum, the other's are at minimum or negative. The triple product P in
  the overlap region oscillates in sign. The time-averaged force is repulsive
  (net outward radiation pressure from the oscillating overlap).

- **Quadrature**: Intermediate but closer to repulsive. The overlap region
  has P oscillating with half the period, creating a weaker repulsion.

- **Generic case**: Two independently formed oscillons have random relative
  phase. For most phase offsets (|delta| > ~pi/4), the interaction is
  REPULSIVE. Only nearly in-phase oscillons attract. This means:
  1. The merger observed in the control run is an artifact of identical-phase initialization
  2. Generic oscillon pairs will repel at D=20
  3. Gravitational attraction must overcome this phase-dependent repulsion

### Asymmetry in quadrature case

The delta=pi/2 case shows asymmetric centroid motion: the upper (full amplitude)
barely moves (z_upper=+9.99), while the lower (zero initial amplitude, given
velocity kick) moves outward (z_lower=-10.60). The lower oscillon was "kicked"
by the phase mismatch during initial transient.

## D=40 Wide-Box Test (Anti-Phase, N=400/350, L=70, t=500)

### Motivation

At D=40, tail overlap ~ exp(-kappa*20) ~ 2e-3 is negligible. The direct
coupling that dominated at D=20 is absent. This isolates the gravitational
signal.

### Setup

Both runs use anti-phase (delta=pi) initialization.
- Gravity: N=400, dx=0.351, alpha=0.1, beta=0.1
- Control: N=350, dx=0.401, alpha=0, beta=0

### Separation Trajectories

```
t       Control(sep)  Gravity(sep)  Delta
  0       40.00         40.00       +0.00
 50       39.84         39.96       +0.12
100       39.87         40.26       +0.39
150       39.94         40.77       +0.83
200       39.89         41.22       +1.33
250       39.81         41.95       +2.14
300       39.68         42.64       +2.96
350       39.55         43.54       +3.99
400       39.42         44.47       +5.05
450       39.27         45.58       +6.31
500       39.13         46.38       +7.25
```

### Analysis

Control shows slow contraction: 40.00 -> 39.13 = -0.87 over 500 t.u.
(residual tail coupling at D=40, rate ~ -0.002/t.u.)

Gravity shows monotonic expansion: 40.00 -> 46.38 = +6.38 over 500 t.u.
The expansion is ACCELERATING (d2sep/dt2 ~ +5e-5), not decelerating.

### Why gravity doesn't attract (even at D=40)

The expected inter-oscillon gravitational deceleration is:
  a_grav = alpha*E/(4*pi*(D/2)^2 * E) = alpha/(4*pi*(D/2)^2) ~ 2e-5

But the observed acceleration is +5e-5 (outward). The self-gravitational
modification of each oscillon dominates:

1. m_eff = sqrt(m0^2 - beta*Phi) > m0 near each oscillon's own core
2. Higher m_eff -> higher breathing frequency (omega=0.960 vs 0.948)
3. Modified shedding radiation pattern -> persistent outward velocity
4. The "rocket effect" from asymmetric shedding overwhelms 1/r^2 attraction

The self-modification is a LOCAL effect (each oscillon modifies itself),
while the gravitational attraction is an INTER-oscillon effect scaling as
1/D^2. The local effect does not decrease with D, so increasing separation
does not help.

### Resolution: need alpha << 0.1

To see gravitational attraction, need alpha small enough that:
- Self-modification (delta_m_eff ~ beta*Phi_self ~ alpha*E*beta/(4*pi*R))
  is negligible compared to the gravitational force signal
- This requires alpha << (R/D)^2 ~ (3/20)^2 ~ 0.02
- Try alpha = 0.001-0.01: self-effect reduced 10-100x, gravity reduced
  same factor but accumulates over longer runs (t=2000+)

Alternatively: compare TWO gravity-on runs at different D. The self-modification
is a common mode that cancels in the difference, leaving only the 1/r^2 signal.

## Next Steps

1. **Weak-gravity test (alpha=0.001, D=40, t=2000)**: Reduce alpha by 100x
   to make self-modification negligible. Run longer to accumulate signal.
   Predicted displacement: alpha*t^2/(8*pi*D^2) ~ 0.001*4e6/(8*pi*1600) ~ 0.1
   Marginal but detectable against ~0 control drift.

2. **Differential D test (D=30 vs D=50, alpha=0.1)**: DONE — NEGATIVE.
   D=30 expands +8.91, D=50 expands +5.37. D=30 expands MORE (opposite of
   gravity prediction). Self-modification is NOT a common mode: at smaller D,
   the Phi wells overlap more (Phi(0)=-0.141 vs -0.085), creating stronger
   m_eff enhancement and faster shedding-driven expansion. Ratio 1.66 vs
   predicted 2.78 for 1/D^2. Approach fails.

3. **Phase-averaged interaction**: Scan delta from 0 to pi in steps of pi/6.
   Map out full phase-dependent potential at D=20 (where tail coupling is strong).

4. **Oscillon merger physics**: Re-run in-phase merger at higher resolution
   to study whether remnant re-forms a stable oscillon.
