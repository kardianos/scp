# Breathing Oscillator Characterization

## The Superluminal Problem

**BOTH the UUD proton and deuterium have superluminal field velocities.**

| Structure | Global v_max | Where | Fraction > 0.5c |
|-----------|-------------|-------|-----------------|
| UUD proton (T=500) | **2.22c** | r ≈ 15-20 (outer shell) | 0.4-1.6% |
| Deuterium (T=500) | **2.02c** | r ≈ 3-8 (inner core) | 0.05% of volume |

These are NOT numerical artifacts — they persist across frames and are
localized at specific radii. The field velocities (∂φ/∂t) exceed c.

**However**: this may not violate causality. The Cosserat equation IS
Lorentz-invariant (∂²φ/∂t² = ∇²φ - ...), so the GROUP velocity of wave
packets is always ≤ c. The PHASE velocity of the carrier wave can exceed c
(this is normal for dispersive media — the phase velocity of de Broglie
waves in QM is c²/v > c for massive particles).

The superluminal |v| values are the FIELD velocity (time derivative of φ),
not the velocity of any information or energy packet. The braid as a whole
moves at V < c. The internal field oscillation can exceed c in the same way
that the intersection point of scissors blades can exceed c.

**What this means physically**: the breathing oscillator has internal field
oscillations faster than c. This is consistent with the field being a
standing wave (two counter-propagating waves each at c produce a standing
wave with |∂φ/∂t| that can exceed c at the antinodes).

**What we need to verify**: that the ENERGY transport velocity and the
GROUP velocity of perturbations remain ≤ c. The BLV metric analysis (V39)
already showed this — the effective group velocity is always subluminal
(v_g = k/ω < c for ω² = k² + m²).

## UUD Proton Breathing Character

### Oscillation Structure

| Metric | Value |
|--------|-------|
| Breathing period (E_pot peaks) | **~150 t** |
| E_kin → E_pot correlation | **-0.15 (uncorrelated)** |
| Global v_max | 2.22c (outer shell, r≈15) |
| v_mean (T=500) | 0.064c |
| T_kin (velocity dispersion) | 0.007 |
| Fraction v > 0.5c | 0.4% |
| R_rms range | 10.2 → 20.2 (oscillating) |

**Key finding**: E_kin and E_pot are UNCORRELATED (r = -0.15). This means
the breathing is NOT a simple virial oscillation (where E_kin and E_pot would
be strongly anti-correlated). Instead, the kinetic and potential energy evolve
semi-independently — the structure has multiple internal degrees of freedom
(three braids, each breathing on its own timescale, plus the composite mode).

### Temperature Profile (T=500)

```
Region          T_kin    |v|_mean  v_radial    Interpretation
Core (r<5)      0.001    0.05     -0.010      COLD, slowly contracting
Mid (r=5-10)    0.003    0.10     -0.011      WARM, contracting
Shell (r=10-20) 0.015    0.11     -0.015      HOT, contracting (fastest)
Halo (r>20)     0.003    0.04     -0.005      COLD, slowly contracting
```

**The entire structure is contracting** (v_radial < 0 at all radii). The
hottest region is the outer shell (r=15-20) where T_kin peaks at 0.015.
The core is cold (T_kin = 0.001).

This inverted temperature profile (hot shell, cold core) is characteristic
of a confined, self-gravitating system — the outer material falls inward,
heats up through compression, and the core remains cool because binding
energy is radiated via theta.

### The T > 0 Constraint

The kinetic temperature is everywhere > 0:
- Minimum T_kin = 0.0001 at r = 0.4 (core center)
- The particle is ALWAYS in thermal motion
- Even the coldest region has |v| = 0.041c — never at rest

This satisfies the T > 0 requirement. The breathing oscillator is a
genuinely thermal dynamical system, not a frozen static configuration.

## Deuterium Breathing Character

### Oscillation Structure

| Metric | Value |
|--------|-------|
| Breathing period (E_pot peaks) | **~300 t** |
| E_kin → E_pot correlation | **+0.12 (multi-mode, CORRECTED)** |
| Global v_max | 2.02c (inner core, r≈3) |
| v_mean (T=500) | 0.12c |
| T_kin at core | **0.16** (VERY HOT) |
| Fraction v > 0.5c | 0.05% |
| R_rms range | 56 → 74 |

**CORRECTION (from review-02)**: The E_kin/E_pot correlation of +0.61 was a
subsampling artifact from the 7-frame SFA. On the full 101-point diag series
(dt=5), the correlations are:
- Raw (signed): **+0.12** (nearly zero)
- Detrended: **+0.30** (weakly positive, below "in-phase" threshold)
- Detrended |E_pot| vs E_kin: **-0.30** (weakly anti-correlated)

**The deuterium breathing is MULTI-MODE (uncorrelated)**, more similar to
the UUD proton's free mode (r=-0.15) than to a driven coherent oscillation.
The inter-baryon attraction does NOT synchronize the breathing into a single
coherent mode. Instead, the two baryons breathe semi-independently within
the bound state, each with its own multi-braid internal modes.

The breathing period is 300t (2× the UUD's 150t) because the deuterium is
a larger, more massive system with more internal degrees of freedom.

### Temperature Profile (T=500)

```
Region          T_kin    |v|_mean  v_radial    Interpretation
Core (r<5)      0.16     0.53     -0.08       EXTREMELY HOT, contracting
Inner (r=5-15)  0.01     0.30     -0.03       Warm, contracting
Mid (r=15-40)   0.004    0.20     -0.02       Cool, contracting
Outer (r=40-60) 0.001    0.29     -0.04       Cold, contracting
Halo (r=60-80)  0.0004   0.18     -0.001      Very cold, ~static
Far (r>80)      0.0005   0.05     +0.005      Cold, slightly expanding (BC)

```

**The deuterium core is EXTREMELY hot**: T_kin = 0.16 at r < 5, with
42% of core voxels having |v| > 0.5c. This is the nuclear interaction
zone where the two baryons' fields overlap. The high temperature reflects
the intense field dynamics at the nuclear bond.

**Contracting everywhere except the very edge**: v_radial < 0 out to r≈70.
Only the absorbing boundary zone (r > 80) shows slight expansion. This
confirms the 43% R_rms contraction is real dynamics, not just boundary drain.

### Comparison: UUD vs Deuterium

| Property | UUD Proton | Deuterium |
|----------|-----------|-----------|
| Breathing period | 150 t | 300 t |
| E_kin/E_pot correlation | -0.15 (free) | +0.12 (multi-mode, corrected) |
| Core T_kin | 0.001 | **0.16** (160× hotter) |
| Core |v| | 0.05c | **0.56c** (11× faster) |
| v_radial (core) | -0.01 | **-0.08** (8× faster contraction) |
| Global v_max | 2.22c | 2.02c |
| Fraction > 0.5c | 0.4% | 0.05% |

**The deuterium core is 160× hotter than the single proton.** The nuclear
interaction creates an intense thermal region at the bond zone. The single
proton's core is cool and slow; the deuterium's core is hot and fast.

This temperature difference is the nuclear binding energy being expressed as
kinetic energy at the interaction zone. The binding energy (E_pot = -54 to -277)
is comparable to the core kinetic energy (Ekin_core ~ several hundred).

## Implications for Particle Behavior

### Mass = Total oscillation energy

The particle's inertial mass is its total energy E_total = E_kin + E_grad +
E_mass + E_pot. The breathing oscillation IS the mass — a faster-breathing
particle has more kinetic energy and is therefore heavier.

### Inertia = Resistance to acceleration

When the particle is pushed (external gradient), the breathing oscillation
must rearrange. A more massive (faster-breathing) particle has more internal
KE to rearrange → more inertia. This naturally gives m_inertial = m_gravitational.

### Charge response = θ coupling modulation

When an external θ field impinges on the particle, it modifies the breathing
via the η×curl(θ) term. The response depends on the breathing phase — the
particle is not a static charge but an oscillating one. The time-averaged
response gives the effective charge.

### Temperature determines interaction rate

A hotter particle (higher T_kin) has faster internal dynamics → couples more
strongly to external fields → interacts faster. The deuterium core at T_kin=0.16
would respond to perturbations much faster than the single proton at T_kin=0.001.

This has implications for nuclear reaction rates: hot nuclei interact faster
(consistent with thermonuclear physics).

## Files

- `sfa/analysis/breathing.c` — the breathing characterization tool
- `v41/results/stable/breathing.json` — UUD proton data
- `v42/results/analysis/breathing.json` — deuterium data
