# Combo 2+5: Inertia + Condensate — Lattice Deformation — RESULTS

## Setup

8-oscillon periodic chain with pairwise coupling lambda=0.5 (from Test E).
Parameters: mu=-20, kappa=20, mass=1.0, N_osc=8, d=16, Nx=2560.

Protocol:
1. Equilibrate single oscillon (t=5000) -> M_osc=9.45, omega_breath=1.26
2. Build 8-oscillon chain, pre-equilibrate with damping (t=500)
3. Conservative settle (t=5000) — chain fully relaxed, A12=0 everywhere
4. Apply localized antisymmetric perturbation at oscillon #4:
   phi_1 += eps*Gaussian, phi_2 -= eps*Gaussian (eps=0.01)
5. Evolve conservatively for t=10000
6. Decompose into symmetric (compression) and antisymmetric (shear) modes
7. Measure dispersion, propagation speed, amplitude decay

## Key Results

### 1. The antisymmetric perturbation PROPAGATES

The perturbation launched at oscillon #4 reaches all 8 oscillons.
Peak |A12| at distant oscillons is comparable to (or exceeds) the source,
indicating the perturbation is NOT damped — it freely propagates.

| Oscillon | Distance | Peak |A12| | Ratio to source | Arrival t |
|----------|----------|-------------|-----------------|-----------|
| 4 (src)  | 0        | 9.57e-2     | 1.00            | 0         |
| 3        | 1        | 1.30e-1     | 1.36            | 34.9      |
| 5        | 1        | 1.14e-1     | 1.19            | 24.9      |
| 2        | 2        | 1.10e-1     | 1.15            | 59.9      |
| 6        | 2        | 1.03e-1     | 1.07            | 79.8      |
| 1        | 3        | 1.33e-1     | 1.39            | 109.8     |
| 7        | 3        | 1.10e-1     | 1.15            | 104.8     |
| 0        | 4        | 1.39e-1     | 1.45            | 124.7     |

The peak amplitudes at distant sites EXCEED the source — the antisymmetric
mode grows as it propagates (likely through parametric coupling with the
large-amplitude symmetric breathing oscillation).

### 2. Dispersion: OPTICAL mode, not acoustic

The antisymmetric mode has a FLAT dispersion relation:

| q | k       | omega_asym | c_s_asym | omega_sym | c_s_sym |
|---|---------|------------|----------|-----------|---------|
| 0 | 0.000   | 0.522      | -        | 0.213     | -       |
| 1 | 0.049   | 0.515      | 10.48    | 0.047     | 0.947   |
| 2 | 0.098   | 0.516      | 5.25     | 0.093     | 0.947   |
| 3 | 0.147   | 0.508      | 3.45     | 0.140     | 0.947   |
| 4 | 0.196   | 0.501      | 2.55     | 0.190     | 0.967   |

**Symmetric mode**: Linear dispersion omega ~ 0.95*k. This is an ACOUSTIC
phonon with sound speed c_s = 0.947 (matching Test E result c_s ~ 0.9).

**Antisymmetric mode**: Nearly FLAT at omega ~ 0.51 for all k. This is an
OPTICAL mode — a gapped excitation with omega_gap ~ 0.5. The apparent
c_s values (10.5, 5.3, ...) are artifacts of dividing a nearly-constant
omega by varying k. The true group velocity d(omega)/dk is near zero.

### 3. Physical interpretation

The antisymmetric mode is an **internal (optical) oscillation** of the
oscillon's three-component structure, not a propagating lattice phonon.

- **Symmetric phonon** (omega ~ c_s * k): oscillons move in space as rigid
  objects. This is a true acoustic mode — the lattice compresses and rarefies.

- **Antisymmetric "phonon"** (omega ~ 0.51 = const): the INTERNAL structure
  of each oscillon oscillates (phi_1 vs phi_2 detuning). This is an optical
  mode analogous to the optical branch in a diatomic lattice, where the
  gap comes from the internal restoring force of the pairwise coupling.

The gap frequency omega_gap ~ 0.51 is set by the pairwise coupling lambda=0.5
through the antisymmetric mass: m_A = sqrt(m^2 - lambda) = sqrt(0.5) = 0.707.
The optical mode frequency is below the antisymmetric mass gap, consistent
with a bound internal oscillation.

### 4. Propagation speed

From arrival times, the antisymmetric signal travels at roughly
c_prop ~ 0.4-0.6, consistent with the GROUP velocity of the optical band
(which has small but nonzero dispersion: omega drops from 0.52 at q=0 to
0.50 at q=4, giving v_g ~ -0.02/0.20 ~ -0.1, suggesting weak backward
propagation mixed with the overall signal front).

The actual propagation mechanism is likely through the weak COUPLING between
the optical and acoustic branches at the oscillating cores, not through the
optical branch dispersion itself.

### 5. Energy conservation

E_initial = 96.616, E_final = 96.114 (0.52% loss over t=10000).
The perturbation energy Delta_E = 3.1e-4 is negligible compared to the
total chain energy — this is truly a small perturbation.

### 6. Amplification effect

The antisymmetric amplitude at distant oscillons EXCEEDS the source.
Initial perturbation: A12 = 0.015 at osc #4.
Late-time A12: ~0.03-0.14 at all oscillons (2x-10x initial).

This means the symmetric breathing mode (which has amplitude ~0.5) acts as
an **energy reservoir** that feeds the antisymmetric mode via parametric
resonance. The antisymmetric perturbation is AMPLIFIED, not damped.

## Conclusions

1. **YES: The antisymmetric perturbation propagates coherently.**
   All 8 oscillons develop antisymmetric content. The mode is long-lived
   and even amplified by coupling to the breathing mode.

2. **The antisymmetric mode is an OPTICAL phonon (gapped), not acoustic.**
   omega ~ 0.51 independent of k. Gap set by pairwise coupling lambda.
   This is an internal oscillation mode, not a propagating lattice wave.

3. **As a spin-2 candidate: PARTIALLY positive.**
   - Good: a new propagating degree of freedom exists (optical branch).
   - Bad: it is GAPPED (massive), not gapless (massless).
   - Bad: in 1D, the distinction between "shear" and "compression" is limited.
     True spin-2 character requires transverse directions (2D or 3D lattice).
   - The gap means this mode is more like a MASSIVE spin-2 (Proca-like) than
     a massless graviton analog.

4. **Sound speeds:**
   - Symmetric (compression): c_s = 0.947
   - Antisymmetric (optical): omega_0 = 0.51 (gapped, no well-defined c_s)
   - Ratio: not applicable (different mode types)

5. **The amplification effect is notable.** It suggests the antisymmetric
   mode is parametrically unstable — the breathing oscillation acts as a
   pump that amplifies any symmetry-breaking perturbation. This could have
   implications for the stability of the symmetric oscillon configuration.

## Files

| File | Description |
|------|-------------|
| `src/combo25.c` | Simulation code |
| `data/combo25_asym_ts.tsv` | A12, A13, S vs time at each oscillon |
| `data/combo25_energy_ts.tsv` | Total energy vs time |
| `data/combo25_positions.tsv` | Oscillon positions vs time |
| `data/combo25_modes.tsv` | Dispersion relation (sym vs asym) |
| `data/combo25_asym_spectrum_2d.tsv` | Full (k, omega) spectrum for A12 |
| `data/combo25_sym_spectrum_2d.tsv` | Full (k, omega) spectrum for u_n |
| `data/combo25_propagation.tsv` | Propagation: arrival times, peak amplitudes |
