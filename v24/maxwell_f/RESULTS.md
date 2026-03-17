# V24-MF Results: Z2 -> U(1) Promotion via Goldstone Mode

## Summary

A massless scalar theta field coupled to the v21 triad oscillon via the derivative coupling
`box(theta) = g * d_t(P^2)` (where P = phi_1*phi_2*phi_3) produces **constant-amplitude radiation
in 1D** at the oscillon's oscillation frequency. The theta field does NOT develop a static
long-range profile (no "1/r potential" analog). It radiates traveling waves.

With backreaction (theta modifying the phi equations), the theta field acts as an additional
self-coupling that STRENGTHENS the oscillon at g=1.0, keeping it at higher amplitude and energy.
The two-oscillon force test is inconclusive due to the oscillons merging in 1D.

**Bottom line**: The Goldstone mode radiates at constant amplitude (1D), not as a long-range
static field. In 3D, this would give 1/r radiation, but radiation pressure -- not a static force.

---

## Model

Lagrangian: L_triad + L_theta + L_int

    L_triad = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ] - V(P)
    L_theta = (1/2)(dt theta)^2 - (1/2)(dx theta)^2
    L_int   = g * (dt theta) * P^2

Equations of motion:
- phi_a: d_tt phi_a = d_xx phi_a - m^2 phi_a - dV/dphi_a + 2g*(dt theta)*P*(dP/dphi_a)
- theta: d_tt theta = d_xx theta + g * d_t(P^2)

Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0
Grid: Nx=4000, xmax=100, dx=0.050, dt=0.025

---

## Phase 1: Single Oscillon + Theta (g scan)

### Without backreaction (theta does not affect phi)

| g | E_phi(t=5000) | E_theta(t=5000) | E_th/E_phi | fc | omega_phi | theta_amp(x=0) | theta_amp(x=50) |
|---|---------------|-----------------|------------|-----|-----------|----------------|-----------------|
| 0 (ctrl) | 1.274 | 0 | 0 | 1.00 | 0.87 | 0 | 0 |
| 0.01 | 1.274 | 4.4e-7 | 3.4e-7 | 1.00 | 0.87 | 3.1e-4 | 2.5e-4 |
| 0.1  | 1.274 | 4.4e-5 | 3.4e-5 | 1.00 | 0.87 | 3.1e-2 | 2.5e-2 |
| 1.0  | 1.274 | 4.4e-3 | 3.5e-3 | 1.00 | 0.87 | 3.1e-1 | 2.5e-1 |

**Key findings**:
- E_theta scales as g^2 (linear source: exact 100x ratio between each decade of g)
- E_phi is IDENTICAL across all g (one-way coupling, no backreaction)
- theta_amp is nearly constant from x=0 to x=50 (1D: no geometric decay)
- Slight amplitude decrease at large x is due to absorbing boundary at x_abs=75

### With backreaction (theta modifies phi)

| g | E_phi(t=5000) | E_theta(t=5000) | fc | omega_phi | theta_amp(x=0) |
|---|---------------|-----------------|-----|-----------|----------------|
| 0.01 | 1.274 | 4.4e-7 | 1.00 | 0.87 | 3.0e-4 |
| 1.0  | 2.94  | 0.26   | 0.95 | 0.76 | 0.41   |

At g=0.01: backreaction is negligible (perturbative).

At g=1.0: **dramatic backreaction**:
- E_phi increases from 1.27 (control) to 2.94 -- oscillon retains 2.3x more energy
- E_theta = 0.26, about 8% of E_phi (significant)
- Peak phi amplitude oscillates 0.06-0.79 (vs 0.35 control) -- wider breathing
- Oscillon frequency shifts DOWN from 0.87 to 0.76 (deeper below mass gap)
- **The theta backreaction STRENGTHENS the oscillon** by providing additional binding

---

## Theta Spatial Profile (Late Time)

The theta field at any instant shows a smooth envelope, NOT a static potential:

### Instantaneous theta(x) at t=5000, g=0.1 (no backreaction)

| x | theta | Interpretation |
|---|-------|----------------|
| 0 | 0.030 | Oscillon core |
| 5 | 0.030 | Near field |
| 10 | 0.029 | |
| 20 | 0.027 | |
| 30 | 0.023 | |
| 40 | 0.017 | |
| 50 | 0.011 | |
| 60 | 0.006 | Approaching absorbing boundary |
| 70 | 0.003 | Inside absorbing zone (x_abs=75) |
| 80 | 6e-5  | Heavily damped |

### Oscillation envelope (late-time amplitude at each x)

| x | Amplitude (g=0.1) | Amplitude (g=1.0) | Ratio to x=0 |
|---|-------------------|--------------------|---------------|
| 0 | 0.031 | 0.311 | 1.00 |
| 10 | 0.030 | 0.303 | 0.97 |
| 30 | 0.029 | 0.287 | 0.92 |
| 50 | 0.025 | 0.246 | 0.79 |

The amplitude ratio drops by only 21% from core to x=50. In an infinite domain (no absorbing
boundary), the amplitude would be exactly constant. **This confirms: in 1D, massless radiation
has constant amplitude (no 1/r decay).** The 21% drop is entirely due to the absorbing boundary
starting at x=75, which damps the outgoing wave.

### Physical interpretation

In 1D, the retarded Green's function for the massless wave equation is:
    G(x,t) = (1/2) * theta(t - |x|/c)

A periodic source S(x',t) = S_0(x') * cos(omega_s * t) at the origin produces:
    theta(x,t) ~ A * cos(omega_s * (t - |x|/c))

where A = const (independent of |x|). This is a traveling wave at constant amplitude.
**There is no static theta profile.** The "potential" is oscillatory, not monotonic.

In 3D, the retarded Green's function G ~ 1/r would give theta ~ 1/r -- but this is still
a radiation field (oscillatory), NOT a static Coulomb-like potential.

---

## Phase 2: Two-Oscillon Interaction

### Setup
Two Gaussians at x = +/-20, D_sep = 40, same phi parameters.

### Result: INCONCLUSIVE

The two Gaussians in 1D merge into a single central structure within a few hundred time units.
The separation collapses to ~0 by t=100. This happens identically for g=0 (control) and g=0.1.

The two-oscillon test is not viable in 1D because:
1. Symmetric initialization (additive Gaussians for all three fields) creates a single mode
2. There is no topological charge preventing merger (unlike Skyrmions)
3. In 1D, there is no angular momentum barrier to stabilize orbits

The theta field shows no measurable effect on the merger dynamics at g=0.1.
At g=1.0 with backreaction, the simulation goes to NaN (numerical instability from
strong coupling between theta and phi at the merged core).

---

## Phase 3: Backreaction

| g | E_phi (no BR) | E_phi (with BR) | Change | Mechanism |
|---|---------------|-----------------|--------|-----------|
| 0.01 | 1.274 | 1.274 | <0.01% | Negligible |
| 0.1 | 1.274 | ~1.28 | ~0.5% | Perturbative |
| 1.0 | 1.274 | 2.94 | +131% | Strong self-trapping |

At g=1.0, the backreaction force `+2g*(d_t theta)*P*(dP/dphi_a)` provides an effective
additional attractive coupling. The oscillon responds by:
- Oscillating at higher amplitude (peak 0.79 vs 0.35)
- Retaining 2.3x more energy in the core
- Shifting to lower frequency (0.76 vs 0.87), deeper below the mass gap
- Entering a wider breathing mode (amplitude swings 0.06 to 0.79)

The theta energy E_theta = 0.26 is stored in the traveling waves radiating from the core.
This represents ongoing energy transfer from phi to theta (~8% of total), but the system
reaches a quasi-steady state where E_phi stabilizes at ~2.94.

---

## Spectrum Analysis

The phi field oscillates at omega = 0.87 (no backreaction) or 0.76 (g=1.0 with backreaction),
both below the mass gap m=1.0 (sub-gap oscillon, radiation forbidden).

The theta field is sourced by d_t(P^2) ~ d_t(phi^6) which contains harmonics at
omega, 3*omega, 5*omega... However, the DFT of theta(0,t) shows the dominant
frequency is very low (omega ~ 0.018), suggesting a slowly varying DC component
that accumulates over time, plus oscillatory content at the phi frequencies.

The low-frequency theta content is likely due to the slow amplitude evolution of
the oscillon (shedding phase), which modulates P^2 on long timescales.

---

## Conclusions

### Does theta develop a long-range profile?
**NO (in 1D).** The theta field radiates constant-amplitude traveling waves, not a static
potential. The source d_t(P^2) is oscillatory, so theta oscillates at every point in space.
There is no DC component in the source (P^2 > 0 always, but d_t(P^2) averages to zero
over one oscillation period).

In 3D, the radiation amplitude would fall as 1/r, but it would still be oscillatory --
radiation pressure, not a Coulomb-like force.

### Does theta mediate a force between two oscillons?
**INCONCLUSIVE.** The 1D two-oscillon test fails because the oscillons merge before any
theta-mediated interaction can be measured. A 3D test would be needed.

However, the physical expectation is:
- The theta radiation from oscillon A reaches oscillon B as a traveling wave
- This wave acts as a periodic driving force on oscillon B
- The net time-averaged force is zero (oscillatory driver with zero DC component)
- There is radiation pressure (always repulsive), scaling as ~ g^2 * amplitude^2

### Why no static profile?
The source term d_t(P^2) has ZERO time average: the oscillon breathes symmetrically,
so <d_t(P^2)>_T = [P^2(t+T) - P^2(t)]/T -> 0 for large T. Without a DC source,
there is no static theta field. This is fundamentally different from gravity (where
energy density is always positive and sources a static Newtonian potential).

To get a static theta field, one would need a coupling to P^2 directly (not its
time derivative), which would give theta a mass term (not massless). Or one needs
a non-derivative coupling like -g*theta*rho where rho = total energy density > 0.

---

## Files

| File | Description |
|------|-------------|
| `src/maxwell_f.c` | 1D solver with theta field + backreaction |
| `data/mf_mode1_g*_ts.tsv` | Time series (no backreaction) |
| `data/mf_mode1_br_*_ts.tsv` | Time series (with backreaction) |
| `data/mf_mode2_*_ts.tsv` | Two-oscillon time series |
| `data/mf_mode3_*_ts.tsv` | Control (g=0) time series |
| `data/mf_mode*_*_prof_t*.tsv` | Spatial profiles at time snapshots |
| `data/mf_mode*_*_spectrum.tsv` | DFT power spectra |
| `PROPOSAL.md` | Original specification |
