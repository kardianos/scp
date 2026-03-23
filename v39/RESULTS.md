# V39 Results

## Original Braid Reconstruction

The original V34 plain braid was successfully reconstructed and confirmed operational
with the full 6-field Cosserat equation.

**SFA**: `data/original_braid.sfa` (4.6 GB, 61 frames, N=128 L=15 T=300)
**View**: `/home/d/code/scp/sfa/viewer/volview /home/d/code/scp/v39/data/original_braid.sfa`

### Parameters (CORRECT — matching V34)
- Equation: 6-field Cosserat (3 phi + 3 theta)
- m²=2.25, m_θ²=0 (massless theta), η=0.5
- μ=-41.345, κ=50, A_bg=0.1
- δ = {0, 3.0005, 4.4325}
- R_tube=3.0, ellip=0.3325
- **NO z-truncation** — braid fills entire z-axis, periodic BC
- Single elliptical tube, traveling wave

### Results
- E_pot oscillates -65 to -157 (strong breathing, binding maintained)
- theta_rms grows from 0 to 0.070 (active theta field, EMF-like pattern visible)
- E_total drift: 0.37% over T=300 (excellent conservation)
- Survival: YES

---

## Density-Dependent κ Retest (6-field Cosserat) — DEFINITIVE NEGATIVE

### Background

The 3-field-only `collapse_3d.c` showed a "stable collapsed core" at γ=10 (E_pot = -218,000).
This was flagged as suspect because the simulator was missing the theta/curl coupling.

`collapse_3d.c` was **rewritten** with full 6-field Cosserat physics:
- 18 arrays (6 fields × {val, vel, acc})
- curl(θ) coupling in phi forces, curl(φ) coupling in theta forces
- Correct Lagrangian derivative for mode 3 (Term I + Term II density feedback)
- 6-column SFA output (phi_x/y/z + theta_x/y/z, compressed f32, BSS+zstd)
- theta_rms diagnostic for verification

### Oscillon Init (Gaussian blob, no helical structure)

All runs: N=128, L=10, spherical absorbing BC, dt_factor=0.025, η=0.5, m_θ²=0.

| γ | T | E_pot(final) | phi_max(final) | theta_rms(final) | Outcome |
|---|---|---|---|---|---|
| 0.6 | 50 | -1.9 | 0.48 | 0.012 | Breathes then disperses |
| 2.0 | 50 | ~0 | 0.17 | 0.002 | Disperses completely by t=40 |
| 5.0 | ~35 | NaN | NaN | NaN | **BLOWUP** at t~25 |
| 10 | ~28 | NaN | NaN | NaN | **BLOWUP** at t~4, NaN by t=28 |

The stability boundary for oscillons is between γ=2 and γ=5. Neither stable run
shows collapse — they just disperse through the absorbing BC.

### Braid Init (helical traveling wave, standard V34 parameters)

All runs: N=128, L=15, spherical absorbing BC, dt_factor=0.025, η=0.5, m_θ²=0,
A_bg=0.1, braid R_tube=3.0, ellip=0.3325.

| Run | γ | E_pot(0) | E_pot(100) | E_pot(200) | E_tot(200) | phi_max(200) | θ_rms(200) |
|-----|---|---------|-----------|-----------|-----------|-------------|-----------|
| **Control** | 0 | -92 | -85 | **-67** | **554** | **0.75** | **0.011** |
| Marginal | 0.6 | -111 | -3 | -13 | 305 | 0.60 | 0.006 |
| Moderate | 2.0 | -136 | -0.5 | **~0** | **79** | **0.19** | **0.001** |

**Key finding: higher γ causes FASTER dispersal, not collapse.**

- **Control (γ=0)**: Braid maintains strong binding throughout. E_pot stays -67 to -85.
  phi_max holds at 0.75. theta_rms stable at 0.011. **Braid survives T=200.**

- **γ=0.6**: Initial binding is deeper (-111 vs -92, matching theory prediction).
  But binding weakens much faster — E_pot drops from -111 to -3 by t=100.
  Total energy 305 at T=200 (vs 554 for control). Braid weakened but alive.

- **γ=2.0**: Deepest initial binding (-136), but braid **disperses completely** by t=100.
  E_pot = -0.001 at T=200, phi_max = 0.19 (background noise). Structure is GONE.
  The deeper well causes more violent oscillations that radiate away through theta.

### Mechanism of Failure

The density-dependent kappa is **counterproductive** with 6-field physics:

1. Higher γ → deeper V(P) well at core → more violent oscillations
2. More violent oscillations → stronger curl(φ) → stronger θ field sourcing
3. θ waves radiate energy outward through the absorbing boundary
4. Net effect: the braid burns through its binding energy faster
5. **The θ curl coupling acts as an efficient energy extraction channel**

The 3-field "stable collapsed core" at γ=10 was an artifact of:
- Missing theta coupling (no radiation channel for collapse energy)
- κ saturation providing a hard floor that trapped energy in the core
- Without theta, the only energy loss was through slow phi radiation

With 6-field physics, theta curl coupling extracts energy efficiently and
prevents any self-reinforcing collapse. **Variant F is ruled out.**

### SFA Files (6-field, compressed f32)

| File | Description | Size | Status |
|------|-------------|------|--------|
| `data/braid_control.sfa` | Control braid, mode 0, T=200 | 1.5 GB | GOOD |
| `data/braid_g0.6.sfa` | Braid + κ-dep γ=0.6, T=200 | 1.6 GB | GOOD |
| `data/braid_g2.0.sfa` | Braid + κ-dep γ=2.0, T=200 | 1.5 GB | GOOD |
| `data/kappa6f_g0.6.sfa` | Oscillon + κ-dep γ=0.6, T=50 | 222 MB | GOOD |
| `data/kappa6f.sfa` | Oscillon + κ-dep γ=10 (blowup) | 184 MB | BLOWUP |

View: `/home/d/code/scp/sfa/viewer/volview <path>`

All SFA files contain 6 columns: phi_x, phi_y, phi_z, theta_x, theta_y, theta_z.

---

## BLV Effective Metric — Semi-Analytical Estimate

### Motivation

With the density-dependent κ collapse ruled out, we asked: can we estimate the
internal metric effects of a concentrated-ρ braid analytically, without a full
simulation? The BLV (Boillat-Lanczos-Vernieri) effective metric formalism gives
the propagation speed of perturbation waves through a background field
configuration. Regions where waves slow down act as attractive gravitational
potential wells.

### Method

The solver (`src/blv_cosserat.c`) takes the measured radial energy density
profile from the V34 phonon test (N=256, L=60, T=200) and computes:

1. **Background reconstruction**: Convert ρ(r) to field amplitude A(r) via
   scaling from the known braid core amplitude (A_core=0.8, ρ_core=0.086).

2. **Time-averaged V(P) coupling**: At each radius, numerically average the
   second derivative of V(P) over one oscillation period of the braid's helical
   phase structure (δ = {0, 3.0005, 4.4325}). This gives the effective mass
   matrix W_ab(r) that perturbation waves see.

3. **Level 1 — Scalar perturbation** (phi only, no theta): Effective mass
   m_eff²(r) = m² + W_diag(r). Group velocity v_g(r) at reference wavenumber
   k = m (Compton wavelength). Potential Φ(r)/c² = -Δm_eff²/(2k²).

4. **Level 2 — Coupled phi-theta**: For z-propagating transverse modes, the
   curl coupling gives the dispersion relation:
   ```
   (α - m² - W)(α - m_θ²) = η²k²
   ```
   where α = ω² - k². With m_θ²=0 (massless theta), the phi-dominant branch
   has higher effective mass than Level 1 (theta adds a stiffening channel).

5. **Deflection angle**: Born approximation integral θ(b) = -2∫ (dn/dr) dt
   using cosh substitution r = b·cosh(t).

### Results

**Level 1 (scalar, phi-only):**

| Quantity | Value |
|----------|-------|
| Φ_min/c² | **-0.0215** (attractive) |
| Location | r ≈ 2.75 (braid core) |
| Potential depth | **20.2 MeV** |
| Core speed reduction | 1.5% slower than vacuum |
| Half-width | ~5 code = 2.8 fm |

**Level 2 (coupled phi-theta, η=0.5):**

| Quantity | Value |
|----------|-------|
| Φ_min/c² | **-0.0198** (attractive) |
| Location | r ≈ 2.75 |
| Potential depth | **18.6 MeV** |
| Core speed reduction | 1.3% slower than vacuum |
| Theta correction | **8.1% shallower** than Level 1 |

The theta coupling weakens the potential because massless theta provides a fast
propagation channel (speed c) that partially short-circuits the slowing effect.
Perturbation energy can convert from phi (slow in the core) to theta (always
fast) via the curl coupling, raising the effective group velocity.

**Potential well structure:**

```
   Phi/c²
   +0.08 |         * * * * * * * * *  (repulsive barrier: depleted region)
         |       *                   *
   +0.04 |     *                       *
         |   *                           *   *   *   *
    0.00 |--*-----------------------------*-----------*-------> r
         | *                                           (approaches 0)
   -0.02 |*  (attractive core)
         0  2  4  6  8  10  15  20  25  30  35  40
```

The well is attractive only inside r ≈ 5 (braid core). Beyond the core, Φ flips
**positive** — the depleted region acts as a repulsive barrier. This creates a
confining potential: waves are attracted INTO the core but repelled from the
depletion zone. The sign flip occurs because V''(P) changes sign as P crosses
through the κ-saturation regime.

### Physical Comparison

| System | Φ/c² at surface |
|--------|----------------|
| GR proton (Schwarzschild) | ~10⁻³⁹ |
| BLV Cosserat braid (Level 2) | ~10⁻² |
| **Ratio** | **~10³⁷** |

The braid's internal metric effect is 10³⁷ times stronger than Newtonian gravity
at the proton scale. This is consistent with:
- Prior Skyrme model results (v3 null-rotor: ratio ~10³⁰–10³⁵)
- The L₆ sextic term analysis (v2: ratio ~10³⁷ at λ₆=10)
- Known physics: nuclear binding (~20 MeV) is ~10³⁸ times stronger than
  gravitational binding at nuclear scales

The result is **nuclear physics, not gravity**. The V(P) coupling creates a
~20 MeV potential well at ~3 fm scale — comparable to the nuclear potential
well depth (~50 MeV) and range (~1.5 fm). This is the braid's self-binding
energy, not a gravitational effect.

### Theta Coupling Effect

The 8.1% reduction from theta is modest because η²/(m² + W) ~ 0.25/2.25 ~ 0.11.
The correction scales as:

```
Phi_L2/Phi_L1 ≈ 1 - η²/(m² + W + η²) ≈ 0.92
```

This matches the 0.919 ratio from the numerical computation. The theta channel
acts as a parallel propagation path that dilutes the metric perturbation. For
larger η (stronger coupling), the correction would be larger — at η=2.0 (where
simulations show braid dissolution), the theta correction would be ~50%.

### Deflection Angles

| Impact parameter b | θ_L1 (rad) | θ_L2 (rad) |
|--------------------|-----------|-----------|
| 1.0 | -0.016 | -0.015 |
| 3.0 | +0.030 | +0.026 |
| 5.0 | +0.038 | +0.034 |
| 10.0 | -0.004 | -0.003 |
| 20.0 | -0.006 | -0.005 |

The sign alternation reflects the potential structure: negative θ (inward
deflection) at b < 2 and b > 8 (attractive core and far-field depletion),
positive θ (outward deflection) at b ≈ 3–6 (repulsive barrier at the braid
surface). Peak deflection is ~0.03 rad at b ≈ 5 — detectable in a simulation
but too small for dramatic lensing.

### Interpretation

These results are **not conclusive** — they rely on a spherically averaged,
time-averaged background reconstruction from the measured density profile. The
actual braid is anisotropic (elliptical, z-aligned) and time-dependent
(oscillating). A full 3D Level 3 computation (reading a single SFA frame and
computing the BLV metric on the grid) would capture these effects.

However, the estimates are consistent with known physics at every level:

1. **The potential depth (~20 MeV) is nuclear-scale**, matching the braid's
   binding energy from the V(P) coupling. This is self-consistent — the metric
   perturbation IS the binding.

2. **The ratio to GR (~10³⁷) matches** the known ratio of nuclear to
   gravitational forces at the proton scale.

3. **The theta correction (8%) is modest**, consistent with the V39 simulation
   result that theta weakens but does not destroy braid stability at η=0.5.

4. **No small-scale black hole is possible.** The potential Φ/c² ~ 0.02 is
   far from the trapping threshold (Φ/c² → -0.5 for horizon formation). Even
   concentrating N braids, the superposition gives Φ_N ~ N × 0.02/r^1.2, which
   reaches -0.5 only at N ~ 25 at the core — but at that density, the braids
   overlap and the superposition approximation breaks down. The gradient energy
   (which grows as R⁻⁴ under compression, per THEORY_density_kappa.md §3.3)
   always prevents collapse.

5. **Gravity requires collective, macroscopic accumulation** — not single-particle
   self-trapping. This is exactly how real gravity works: ~10¹⁹ protons make a
   Planck mass, ~10⁵⁷ make a solar mass. The BLV metric of a single braid is
   nuclear binding, not gravity.

### Files

| File | Description |
|------|-------------|
| `src/blv_cosserat.c` | Level 1+2 BLV metric solver |
| `../v34/phonon_test/data/depletion_t0100.tsv` | Input: measured ρ(r) profile |

---

## Issues Found During V39

### RESOLVED: collapse_3d.c was 3-field only
The collapse_3d.c simulator was **rewritten** with full 6-field Cosserat physics.
The rewrite includes:
- 18 arrays (was 9)
- curl coupling in forces (was absent)
- theta in Verlet step, boundary damping, energy diagnostics
- Correct mode 3 force (Term I + Term II density feedback)
- 6-column SFA output with theta fields

Verification: all runs show theta_rms growing from zero (0.001-0.011 at T=200),
confirming the theta coupling is active.

### CRITICAL: "truncated" is NOT the original braid
The V37 "truncated" geometry wraps the braid in a Gaussian z-envelope (σ_z=3.0),
making it a compact ~6 code unit blob. The ORIGINAL V34 braid has NO z-truncation
and extends the full box length via periodic BC.

| | Original V34 | V37 "truncated" | V37 "braid3" |
|---|---|---|---|
| Z-envelope | NONE | σ_z=3.0 | σ_z=3.0 |
| Strands | 1 (single tube) | 1 (single tube) | 3 (helical) |
| R_tube | 3.0 | 3.0 | 2.0 |
| R_helix | N/A | N/A | 1.0 |
| Ellipticity | 0.3325 | 0.3325 | none |
| P_int(0) | ~160 | ~23 | ~142 |
| E_pot(0) | -92 | -17 | -74 |

---

## Files

| File | Description | Status |
|------|-------------|--------|
| `data/braid_control.sfa` | Control braid, 6-field, mode 0, T=200 | **GOOD** |
| `data/braid_g0.6.sfa` | Braid + κ-dep γ=0.6, 6-field, T=200 | **GOOD** |
| `data/braid_g2.0.sfa` | Braid + κ-dep γ=2.0, 6-field, T=200 | **GOOD** |
| `data/kappa6f_g0.6.sfa` | Oscillon + κ-dep γ=0.6, 6-field, T=50 | **GOOD** |
| `data/kappa6f.sfa` | Oscillon + κ-dep γ=10, 6-field (blowup) | BLOWUP |
| `data/original_braid.sfa` | V34 plain braid, 6-field, confirmed | GOOD |
| `data/kappa_3d.sfa` | Old 3-field κ-dep collapse | SUPERSEDED |
| `data/kappa_3d_periodic.sfa` | Old 3-field κ-dep periodic | SUPERSEDED |
| `data/kappa_3d_sphere.sfa` | Old 3-field κ-dep sphere | SUPERSEDED |
| `src/collapse_3d.c` | 6-field Cosserat simulator (rewritten) | **GOOD** |
| `src/oscillon_1d.c` | 1D simulator | 3-field only |
| `src/critical_density_1d.c` | 1D phase scanner | 3-field only |
| `src/blv_cosserat.c` | BLV effective metric solver (Level 1+2) | **GOOD** |
| `THEORY_density_kappa.md` | Theoretical analysis | Valid (math) |
