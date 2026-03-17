# V26-DynABC: Full Dynamic Braid — RESULTS

## Summary

**DEFINITIVE NEGATIVE.** The full dynamic braid (propagating + rotating + massless)
does NOT survive. Both the dynamic and static configurations disperse completely
by t=500. The massless braided helical wave is not self-sustaining.

## Parameters

| Parameter | Value |
|-----------|-------|
| mu        | -20   |
| kappa     | 20    |
| m         | 0.0 (massless) |
| A0        | 0.8   |
| R_tube    | 3.0   |
| Omega     | 0.1   |
| k_axial   | 0.1571 (= 2*pi/40) |
| v_init    | 0.2571 (= k + Omega) |
| BC        | Periodic z, absorbing x,y |

Runs at N=48 (L=20, dx=0.85, t=500) and N=128 (L=20, dx=0.315, t=5 for initial).

## 1. SURVIVAL: NO

Both configurations fail the survival criteria (fc > 0.3, |P| > 0.01):

### Dynamic braid (propagating + rotating):
| Time | E_total | fc    | peak_phi | |P|_max  |
|------|---------|-------|----------|---------:|
| 0    | 250.7   | 0.393 | 0.784    | 0.120    |
| 50   | 38.0    | 0.016 | 0.115    | 0.000396 |
| 100  | 9.1     | 0.221 | 0.056    | 0.000041 |
| 200  | 2.7     | 0.057 | 0.031    | 0.000007 |
| 300  | 1.0     | 0.101 | 0.019    | 0.000002 |
| 500  | 0.47    | 0.085 | 0.009    | 0.000000 |

Energy drops from 251 to 0.47 (99.8% radiated). Peak |P| drops 6 orders of magnitude.

### Static braid (same phase, vel=0):
| Time | E_total | fc    | peak_phi | |P|_max  |
|------|---------|-------|----------|---------:|
| 0    | 217.1   | 0.393 | 0.784    | 0.120    |
| 50   | 30.8    | 0.014 | 0.098    | 0.000255 |
| 100  | 6.8     | 0.235 | 0.048    | 0.000029 |
| 200  | 2.0     | 0.056 | 0.028    | 0.000005 |
| 300  | 0.66    | 0.108 | 0.013    | 0.000000 |
| 500  | 0.44    | 0.084 | 0.004    | 0.000000 |

Energy drops from 217 to 0.44 (99.8% radiated). Static disperses at similar rate.

### Comparison: Dynamic vs Static
Both reach essentially identical final states (E ~ 0.4-0.5, fc ~ 0.08). The
initial velocity provides no stabilization advantage. The dispersal is dominated
by the massless wave equation — without mass (m=0), there is no restoring force
to confine the waves.

## 2. NON-BREATHING: Irrelevant (dispersed)

Dynamic (second half of t=500 run):
- rho(center) mean = 5.75e-5, relative variance = 0.66
- Peak omega(rho) = 0.02 (DC drift, not oscillation)
- Peak omega(phi_0) = 2.37

The configuration is "non-breathing" in the trivial sense: there is essentially
nothing left to breathe. The signal is noise on a near-zero background.

## 3. QUADRUPOLE (l=2 fraction)

### At t=5 (N=128, still intact):
| l | Dynamic | Static |
|---|---------|--------|
| 0 | 0.690   | 0.675  |
| 1 | 0.016   | 0.018  |
| 2 | **0.294** | **0.307** |

At early times, the braid geometry does produce significant l=2 content (~30%).
This is from the helical structure imprinting aspherical strain. However, this
is a transient feature of the initial condition, not a stable equilibrium.

### At t=500 (N=48, dispersed):
| l | Dynamic | Static |
|---|---------|--------|
| 0 | 0.860   | 0.962  |
| 1 | 0.032   | 0.008  |
| 2 | 0.108   | 0.030  |

The l=2 fraction decays as the configuration disperses. The residual l=2 in the
dynamic case (10.8%) is slightly higher than static (3.0%), suggesting the
angular momentum preserves some asphericity in the remnant radiation, but the
signal is 10,000x weaker than initial.

## 4. PROPAGATION SPEED

The phase tracking shows v_prop ~ 3-4 (not c=1 as expected). This is an artifact:
the initial braid pattern disperses rapidly, so the "maximum" of phi_0 along z
jumps between different wave packets rather than tracking a single propagating mode.
The R^2 of the linear fit is 0.98, but this reflects overall drift of the radiation
pattern, not coherent soliton propagation.

At early times (N=128, t=0..5): v_prop ~ 4.1, R^2 = 0.46 (poor fit, transient).
The configuration never establishes coherent propagation.

## 5. ANGULAR MOMENTUM CONSERVATION

| Quantity | t=0     | t=500   | Ratio  |
|----------|---------|---------|--------|
| L_z      | 246.68  | 0.41    | 0.0017 |
| P_z      | 40.92   | 0.072   | 0.0018 |
| E_total  | 250.70  | 0.47    | 0.0019 |

All three conserved quantities (energy, linear momentum, angular momentum) decay
by the same factor (~0.2%). This is consistent: the absorbing boundary conditions
in x,y remove energy, momentum, and angular momentum from outgoing radiation.
The ratio L_z/P_z/E remains roughly constant, confirming that the loss mechanism
is radiation through the absorbing boundaries, not numerical dissipation.

At early times (N=128, t=0..5): L_z drops only 0.4%, confirming conservation
in the interior. The loss is entirely from radiation reaching the boundary.

## 6. DYNAMIC vs STATIC: No Advantage

The key finding: the dynamic braid (with initial velocities for propagation +
rotation) disperses at essentially the same rate as the static braid (same
spatial pattern, zero velocity). By t=500, both have radiated 99.8% of their
energy.

This means the initial velocity (k+Omega = 0.257) does not provide any
stabilization. The configuration is fundamentally unstable because:

1. **No mass gap**: m=0 means all excitations propagate at c with no cost.
   There is no potential well to trap the wave.

2. **Triple-product potential too weak**: V(P) = mu*P^2/(1+kappa*P^2) only
   acts where all three fields overlap. As the fields spread, the overlap
   decreases as |P| ~ A^3, making the nonlinear force vanish cubically.

3. **No Derrick-type balance**: For massless fields in 3D, Derrick's theorem
   requires higher-derivative terms (Skyrme term) for static soliton stability.
   The triple-product potential is a standard 0-derivative potential and cannot
   stabilize against dispersal.

## Physics Interpretation

The "particle as process" hypothesis — that a dynamical, massless, rotating,
propagating braided configuration maintains itself through continuous motion —
is **falsified** by this simulation. The braid disperses because:

- Massless waves in 3+1D always spread (1/r decay)
- The triple-product coupling V(P) cannot overcome this because it requires
  3-field overlap which shrinks faster than the linear spreading
- Neither propagation nor rotation provide confinement
- The initial angular momentum and linear momentum are simply carried away
  by the outgoing radiation

**For stable solitons, you need either:**
1. A mass gap (m > 0) providing a potential well
2. Topological protection (Skyrmion, hedgehog map)
3. Higher-derivative terms (Skyrme L_4 term) for Derrick balance
4. Or some combination of the above

The braided geometry alone, however interesting its l=2 content at early times,
is not sufficient for self-sustaining configurations.

## Data Files

- `data/dynABC_dynamic_evolution.tsv` — full time series for dynamic run
- `data/dynABC_static_evolution.tsv` — full time series for static run
- `data/dynABC_dynamic_dft.tsv` — DFT of center density (dynamic)
- `data/dynABC_static_dft.tsv` — DFT of center density (static)
- `data/dynABC_dynamic_rho_history.tsv` — rho(center,t) time series
- `data/dynABC_dynamic_strain.tsv` — strain on R=8 shell (dynamic)
- `data/dynABC_dynamic_multipoles.tsv` — multipole decomposition (dynamic)
- `data/dynABC_dynamic_phase_track.tsv` — z-position of phi_0 max vs time
- `data/dynABC_static_strain.tsv` — strain on R=8 shell (static)
- `data/dynABC_static_multipoles.tsv` — multipole decomposition (static)

## Compile and Run

```
gcc -O3 -fopenmp -Wall -o dynABC src/dynABC.c -lm
./dynABC -N 128 -tfinal 500 -o data     # Full spec (requires ~50 min per run)
./dynABC -N 48 -tfinal 500 -o data      # Fast version (~7 min total)
```
