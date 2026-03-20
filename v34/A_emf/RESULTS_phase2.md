# Phase 2 Results: Torsion Waves as EM Carrier

## Summary

Torsion waves are NOT a viable electromagnetic force carrier in the SCP theory.
All four experiments produced negative results for the torsion-as-EM hypothesis:

1. Torsion (det(J)) and depletion (delta_phi2) have different radial profiles,
   confirming they are distinct quantities — but both are short-range.
2. Torsion perturbations do NOT propagate as coherent wave packets. They
   disperse within ~5 time units and convert to amplitude perturbations.
3. Braid response to torsion vs amplitude pulses is qualitatively the same
   (both give small +x kicks), differing only by a factor ~3x in magnitude.
4. Opposite-winding braids (W=+1 vs W=-1) do NOT respond in opposite
   directions to the same torsion pulse. No charge-dependent force observed.

---

## Experiment 2a: Torsion Profile of Settled Braid

**Code**: `src/torsion_profile.c`
**Data**: `data/torsion_profile.tsv`
**Input**: Settled braid snapshot from phonon test (N=256, L=60, t=200)

### Method
Read the settled braid snapshot and computed at each grid point:
- Jacobian J_ai = dphi_a/dx_i (3x3 matrix, central differences)
- det(J) = epsilon_ijk (dphi_0/dx_i)(dphi_1/dx_j)(dphi_2/dx_k)
- Sigma_phi2 = phi_0^2 + phi_1^2 + phi_2^2
- |grad phi|^2 = sum of all gradient components squared

Averaged in spherical shells (dr=0.5) centered on the braid.

### Results

Braid center: (0.000, 0.000, 6.487)
Background Sigma_phi2 = 1.581e-2

| Quantity | Core (r<3) | Shell (r=3-8) | Far field (r>13) |
|----------|-----------|---------------|-------------------|
| <det(J)> | ~10^-20 | ~10^-21 | ~10^-22 |
| <det(J)^2> | ~10^-13 | ~10^-12 | ~10^-13 |
| delta_phi2 | 0.28 | 0.03-0.09 | <0.001 |
| <grad2> | 0.07 | 0.005-0.05 | 0.0004 |

**Key findings**:

1. **<det(J)> is consistent with zero at all radii.** The mean Jacobian
   determinant averages to ~10^-20, which is numerical noise. This is expected:
   the braid is approximately axially symmetric, so det(J) has equal positive
   and negative contributions in each shell.

2. **<det(J)^2> peaks at r ≈ 4-5, not at the core.** The torsion fluctuation
   variance peaks in the braid's outer shell, while the depletion (delta_phi2)
   peaks at the core (r < 1). These are clearly different radial profiles.

3. **Different decay rates**: delta_phi2 drops from 0.28 at core to ~0.001
   by r=13 (factor 280). det(J)^2 drops from peak 8.5e-12 at r=4.75 to
   ~1e-13 at r=13 (factor 85). The torsion signal is more extended than
   the depletion but both are localized (no long-range tails).

4. **No long-range torsion halo**: det(J)^2 does not show power-law decay.
   It drops to background levels by r ~ 15, consistent with the braid's
   Yukawa envelope.

**Interpretation**: Torsion and amplitude depletion probe different aspects
of the braid structure, but both are confined to the braid's core region.
Neither produces a long-range field that could mediate a 1/r force.

---

## Experiment 2b: Torsion Wave Packet Propagation

**Code**: `src/torsion_wave.c`
**Data**: `data/torsion_wave/` (torsion), `data/amplitude_wave/` (control)

### Method
Initialized a traveling wave packet in the uniform background (no braid):
- **Torsion pulse**: Rotate (phi_0, phi_1) by delta_theta(x) = 0.05 * exp(-(x+15)^2/8) * cos(pi*x/2)
  This preserves Sigma_phi2 exactly but changes P = phi_0*phi_1*phi_2.
- **Amplitude pulse (control)**: phi_a *= (1 + 0.05 * env * cos(k*x))
  This changes both Sigma_phi2 and P.

Parameters: N=128, L=30, T=60. Both initialized with v_group = 1.0.

### Results

**Energy conservation**: Both runs conserve energy to <0.001% over T=60.
- Torsion: E_total = 7500.10 +/- 0.02 (0.0003% variation)
- Amplitude: E_total = 7501.13 +/- 0.02 (0.0003% variation)

**Torsion pulse does NOT propagate coherently**:

| Time | P centroid x | P_rms | phi2_rms | phi2_max |
|------|-------------|-------|----------|----------|
| 0 | -15.0 | 3.0e-8 | 0.0 | 0.0 |
| 5 | -0.0 | 3.6e-6 | 1.4e-4 | 6.5e-4 |
| 10 | -0.0 | 1.1e-6 | 1.2e-4 | 4.9e-4 |
| 20 | -0.0 | 3.2e-6 | 1.6e-4 | 7.4e-4 |
| 30 | -0.0 | 3.9e-6 | 1.4e-4 | 6.1e-4 |

By t=5 the P perturbation centroid has jumped from x=-15 to x~0, which
is NOT coherent propagation at v=1 (would give x=-10). Instead, the
perturbation has dispersed across the entire box.

**Torsion-to-amplitude conversion**: At t=0, phi2_rms = 0 exactly (the
rotation preserves Sigma_phi2). By t=5, phi2_rms = 1.4e-4, with
phi2_max = 6.5e-4. The torsion perturbation converts ~30% of its
energy into amplitude perturbations within one crossing time.

**Amplitude pulse (control)**: Shows similar dispersive behavior.
Both modes oscillate at the background frequency omega ~ 1.5 and
bounce around the periodic box. No coherent propagation for either mode.

**Interpretation**: In this field theory, all three linearized modes are
massive (m_eff ~ 1.5) in the uniform background. A localized perturbation
of any type disperses as massive wave packets with group velocity
v_g = k/sqrt(k^2 + m^2) < 1. The perturbation wavelength lambda=4
gives k = pi/2 ~ 1.57, so v_g ~ 0.72, and the packet would need
lambda/v_g ~ 5.5 time units to cross one wavelength — consistent with
the rapid dispersal we observe.

There is no massless mode in the far field, so there is no candidate
for a long-range force carrier (whether EM or otherwise).

---

## Experiment 2c: Torsion vs Amplitude Wave Hitting a Braid

**Code**: `src/braid_response.c`
**Data**: `data/braid_torsion/` (torsion), `data/braid_amplitude/` (amplitude)

### Method
1. Initialize W=+1 braid at origin, settle for T=50
2. Apply either torsion or amplitude pulse at x=-15 (same eps=0.05)
3. Track braid x-position for T=200

### Results

**Both pulses give tiny +x kicks (radiation pressure)**:

| Mode | Net x-drift | Drift velocity | Late <x> | SNR |
|------|------------|----------------|----------|-----|
| Torsion | 5.0e-4 | 3.2e-6/t | 5.4e-4 | 3.5 |
| Amplitude | 1.8e-3 | 1.0e-5/t | 1.8e-3 | 5.6 |
| Ratio amp/tor | 3.4x | 3.1x | 3.4x | — |

The amplitude pulse kicks the braid ~3.4x harder than the torsion pulse.
Both kicks are in the +x direction (away from the pulse source at x=-15).

**Energy comparison**: The amplitude pulse deposits 1.48 more energy units
into the system (12146.6 vs 12145.1), which is only 0.012% more. Yet it
produces 3.4x more displacement, suggesting the amplitude perturbation
couples more efficiently to the braid's translational mode.

**Peak amplitude**: Both braids show similar peak_phi2 oscillations over
time (~1.0-1.8), indicating internal breathing modes are not significantly
different between the two pulse types.

**Interpretation**: Both pulse types produce the same qualitative effect:
a tiny radiation-pressure kick in the direction of wave propagation.
The amplitude pulse is ~3x more effective because it directly modulates
Sigma_phi2, which couples to the braid's density structure, while the
torsion pulse must first convert to amplitude perturbations before
interacting. There is no qualitative difference — both behave as
radiation pressure from massive scalar waves.

---

## Experiment 2d: Opposite Winding Braids

**Code**: `src/braid_response.c`
**Data**: `data/braid_torsion/` (W=+1), `data/braid_torsion_neg/` (W=-1)

### Method
Apply identical torsion pulse to braids with opposite winding:
- W=+1: delta = (0, +3.0005, +4.4325)
- W=-1: delta = (0, -3.0005, -4.4325)

### Results

**CAVEAT**: The W=-1 initialization produces a braid with 8.5% less energy
(11117 vs 12145). This is NOT simply a mirror of W=+1 — reversing the
phase offsets changes the relative alignment with the background, producing
a structurally different braid. The energy ratio E(-1)/E(+1) = 0.9153
remains constant throughout, confirming both are stable but different.

**No opposite response observed**:

| Quantity | W=+1 torsion | W=-1 torsion |
|----------|-------------|-------------|
| Late <x> | 5.4e-4 | 4.7e-4 |
| Net drift | 5.0e-4 | 4.8e-4 |
| Drift v_x | 3.2e-6 | 6.4e-6 |
| SNR | 3.5 | 1.3 |

Both braids drift in the SAME direction (+x, away from the pulse).
The magnitudes are comparable. There is no sign of charge-dependent
interaction: W=+1 and W=-1 do not experience opposite forces.

The W=-1 braid's lower SNR (1.3) means its drift is not even
statistically significant above noise.

**Interpretation**: The torsion pulse acts as radiation pressure on both
braids regardless of winding number. The winding number does not serve
as "electric charge" under torsion wave interactions.

---

## Overall Conclusions

### What we tested
Whether torsion (antisymmetric field gradients / Jacobian determinant)
could serve as an electromagnetic force carrier, distinct from the
amplitude depletion channel that mediates gravity.

### What we found

1. **Torsion and depletion have different profiles** (positive):
   det(J)^2 peaks at the braid shell (r~4-5), while delta_phi2 peaks
   at the core (r<1). They probe different structural features.

2. **Torsion waves do not propagate** (fatal):
   All modes are massive (m_eff ~ 1.5) in the far field. Torsion
   perturbations disperse within ~5 time units and convert to amplitude
   perturbations. There is no coherent torsion wave to serve as a photon.

3. **No differential scattering** (fatal):
   Torsion and amplitude pulses produce the same qualitative effect on
   braids: radiation pressure in the propagation direction. The amplitude
   pulse is 3.4x more effective, but this is quantitative, not qualitative.

4. **No charge-dependent force** (fatal):
   Opposite-winding braids respond the same way (both kicked in +x) to
   the same torsion pulse. Winding number does not function as charge.

### Why torsion fails as EM

The fundamental obstacle is that ALL linearized modes around the uniform
background are massive with m_eff ~ 1.5. This includes both amplitude
and torsion perturbations. The mass gap prevents long-range propagation
and makes all perturbations behave as Yukawa (exponentially decaying)
rather than Coulomb (power-law).

For EM to emerge, we would need a MASSLESS mode that:
- Preserves Sigma_phi2 (no gravity coupling)
- Couples to winding number (charge-dependent)
- Propagates at speed c

None of these conditions are met by torsion perturbations in this theory.

### Implications for SCP

The SCP theory has one long-range interaction: gravity, mediated by
the nonlinear amplitude depletion mechanism (confirmed in earlier work
as 1/r^1.2 power-law). It does NOT have a second long-range force.

To obtain EM, the theory would need one of:
1. **A Goldstone mode**: Requires spontaneous symmetry breaking with a
   continuous symmetry. The current V(P) doesn't provide this.
2. **A gauge field**: Requires promoting a global symmetry to local,
   adding connection/gauge degrees of freedom.
3. **A different background**: If a background exists where one mode
   becomes massless, that could provide the photon.
4. **Complex fields / U(1) extension**: The V29 approach with complex
   fields and explicit U(1) gauge coupling.

## Files

- `src/torsion_profile.c` — Experiment 2a: torsion profile analysis
- `src/torsion_wave.c` — Experiment 2b: wave propagation (torsion + amplitude)
- `src/braid_response.c` — Experiments 2c+2d: braid response + winding test
- `data/torsion_profile.tsv` — Radial profile data
- `data/torsion_wave/` — Torsion wave propagation data
- `data/amplitude_wave/` — Amplitude wave propagation data
- `data/braid_torsion/` — W=+1 braid + torsion pulse tracking
- `data/braid_amplitude/` — W=+1 braid + amplitude pulse tracking
- `data/braid_torsion_neg/` — W=-1 braid + torsion pulse tracking
