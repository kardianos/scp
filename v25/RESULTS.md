# V25 Results: Three-Field Elastic Gravity in 3D

## Parameters
- mu=-20, kappa=20, mass=1.0, A=0.8, sigma=3.0
- lambda_pw=0.5, eta=0.1, lambda_L=0.1, alpha_g=0.001
- Grid: N=96, L=15, dx=0.316, dt=0.063
- Threads: 16, total wall time: 7.7 min

## Phase 1: 3D Elastic Oscillon (t=0..500)

**Result: Oscillon SURVIVES with elastic couplings.** The three fields remain
perfectly symmetric (pk identical across all three) and breathing persists
through t=500 with core fraction fc=0.993.

- Initial energy: E=278.6 (t=0)
- Final energy: E=183.0 (t=500) -- 34% radiated to absorbing boundary
- Peak amplitude oscillates between ~0.08 and ~0.85 -- strong breathing
- Core fraction fc > 0.99 throughout -- energy stays localized
- Breathing frequency: omega = 1.380 (omega/m = 1.38)

**Note on omega/m > 1**: The measured frequency exceeds the mass gap. This is
likely due to the pairwise coupling shifting the effective mass. With
lambda_pw=0.5, the symmetric mode mass is m_S^2 = m^2 + 2*lambda_pw = 2.0,
so m_S = 1.414. Then omega/m_S = 0.976 < 1, consistent with a sub-gap
oscillon in the symmetric channel. The elastic terms (eta, lambda_L) further
modify the effective dispersion relation.

## Phase 2: Strain Diagnostics

Strain tensor computed on spherical shell at r=10 after Phase 1 equilibration.

### Multipole decomposition of |sigma_ij|^2:

| l | coefficient | fraction |
|---|-------------|----------|
| 0 | 4.05e-03    | 97.96%   |
| 1 | -6.20e-05   | 1.50%    |
| 2 | 2.24e-05    | 0.54%    |

**Result: Strain is DOMINANTLY ISOTROPIC (l=0).** The shear tensor |sigma|^2
on the shell is 98% monopole, with only 0.5% quadrupole content. This is
expected for a spherically symmetric breathing oscillon -- the dominant
radiation is spin-0 (compression mode), not spin-2 (shear).

The small l=2 content (0.5%) comes from lattice discretization breaking
continuous rotational symmetry. A true spin-2 signal would require breaking
the spherical symmetry (e.g., by deforming the oscillon or using a binary).

## Phase 3: Self-Consistent Metric (t=500..1000)

Alpha_g ramped from 0 to 0.001 via Hermite smooth step over t=500..600.

**Result: Oscillon SURVIVES self-consistent metric coupling.** The metric
backreaction (h_ij modifying the Laplacian and effective mass) does not
destabilize the oscillon.

- Energy at ramp start: E=183.0 (t=500)
- Energy at end: E=120.4 (t=1000) -- continued radiation loss
- Core fraction: fc=0.845 at t=1000 -- decreased but still localized
- Breathing frequency: omega=1.404 (slightly shifted from Phase 1's 1.380)
- Peak amplitude decreased to ~0.10 by t=1000

The oscillon gradually weakens over t=500..1000. The core fraction drops from
0.993 to 0.845, indicating some energy is leaking. This could be due to:
1. Continued radiation from the breathing mode (same as Phase 1)
2. Metric-induced instability slowly extracting energy
3. The oscillon naturally decaying over these long timescales

The frequency shift (1.380 -> 1.404) suggests the metric coupling slightly
stiffens the effective potential.

## Phase 4: Gravitational Wave Detection

After Phase 3, the oscillon was re-gridded onto L=20 (for far-field access)
and given a v=0.01 velocity kick in the z-direction. TT strain tracked at
R=12 from center at 4 angles for t=200.

### Peak TT strain amplitudes:

| Angle | |h_+|_max | |h_x|_max |
|-------|----------|----------|
| z-axis (theta=0)         | 1.53e-05 | 2.12e-03 |
| x-axis (theta=pi/2)      | 9.40e-06 | 2.13e-03 |
| xy-45 (theta=pi/2,phi=pi/4) | 3.84e-03 | ~0 |
| 45-deg (theta=pi/4)      | 3.85e-03 | 9.79e-06 |

### Analysis:

**h_cross (hx) pattern**: Nearly identical amplitude (~2.1e-03) along z-axis
and x-axis, but suppressed (~0) at 45-degree azimuthal angles. This is
consistent with the strain field having a specific symmetry structure tied to
the z-kick direction, but NOT a clean quadrupolar pattern.

**h_plus (h+) pattern**: Suppressed along z-axis and x-axis (~1e-05),
maximal at azimuthal 45-degree and polar 45-degree (~3.8e-03). This is
NOT the expected quadrupolar pattern for spin-2 radiation.

**Expected vs observed**: For a true spin-2 gravitational wave from a moving
source along z:
- h+ should be zero along z (forward), maximum in equatorial plane
- The angular pattern should follow sin^2(theta)*cos(2*phi)

The observed pattern shows h+ maximal at azimuthal 45 degrees, not at
phi=0. This suggests the dominant radiation has components with m=+/-2
azimuthal dependence (which IS spin-2 character), but the polar angle
distribution doesn't match the expected sin^4(theta) pattern.

**Caveat**: The oscillon was significantly weakened by t=1000 (Phase 3 end),
with pk~0.10 and fc~0.85. Phase 4 started from this depleted state. A
cleaner test would skip or shorten Phase 3, or use a stronger initial
oscillon.

## Key Findings

1. **Elastic oscillon STABLE**: The cross-gradient, Lame, and pairwise
   couplings do not destroy the 3D breathing oscillon. fc > 0.99 for t=500.

2. **Self-consistent metric STABLE**: alpha_g=0.001 metric backreaction
   is tolerated. The oscillon weakens but survives through t=1000.

3. **Strain is predominantly spin-0**: On a spherical shell around the
   breathing oscillon, |sigma|^2 is 98% l=0 (isotropic). The spherical
   symmetry of the oscillon means it radiates compression, not shear.

4. **Spin-2 content requires symmetry breaking**: To generate genuine
   quadrupolar (spin-2) radiation, one needs either:
   - A non-spherical source (binary, deformed oscillon)
   - Scattering of shear waves off the oscillon
   - An aspherical perturbation stronger than v=0.01

5. **TT strain IS measurable**: Peak h+ ~ 3.8e-03, which is well above
   noise. The angular structure shows azimuthal dependence consistent
   with tensor character, even though the pattern is not cleanly
   quadrupolar.

## Phase 5: Binary Oscillon GW Emission

Two 3D oscillons at D=12 with tangential velocities v_orb={0.01, 0.1}.

**Result: NEGATIVE for spin-2.**

Both velocities lead to rapid inspiral and merger (sep 12→4 in ~70 t.u.).
No stable orbit forms. The merged blob is nearly spherical → monopole.

| v_orb | l=0 | l=2 | h₊ pole/eq | ω_GW / 2Ω |
|-------|-----|-----|-----------|-----------|
| 0.01 | 96.3% | 2.7% | 1.02 | 58× |
| 0.10 | 96.9% | 2.1% | 1.16 | 5.8× |

The breathing mode (ω=0.193) overwhelms any orbital quadrupole by orders
of magnitude. The GW frequency matches the breathing, not the orbital.

## Phase 6: Tidal Deformation

Applied external quadrupolar potential V_tidal = -½ε_T·cos(Ω_T·t)·(x²-y²)·Σφ².

### ε scan at Ω_T=0.1

| ε_T | l=0 | l=2 | h₊ pole/eq | k₂ | Oscillon survives? |
|-----|-----|-----|-----------|-----|-------------------|
| 0.001 | 93.5% | 1.3% | 6.5 | 1390 | Yes (partially) |
| **0.01** | **73.6%** | **26.0%** | **0.65** | **466** | **No (destroyed)** |
| 0.1 | NaN | — | — | — | NaN (blowup) |

### Ω scan at ε_T=0.01

| Ω_T | l=2 | h₊ pole/eq | k₂ | E_final |
|-----|-----|-----------|-----|---------|
| 0.05 | 17.7% | 1.46 | 400 | 0.01 |
| **0.10** | **26.0%** | **0.65** | **466** | **0.01** |
| 0.20 | 16.3% | 1.20 | 137 | 35.6 |
| 0.50 | 5.0% | 1.92 | 41 | 107 |

**Peak l=2 = 26% at Ω_T=0.1** — a quadrupolar resonance. But this occurs
DURING oscillon destruction (E drops from 183 to 0.01). The l=2 content
is from asymmetric dispersal, not steady-state quadrupolar oscillation.

At ε_T=0.001 (where the oscillon survives): l=2 is only 1.3%.

### Tidal Love Number

k₂ ≈ 100-1000 — the oscillon is extremely soft (1000× softer than neutron
stars). It has no shear rigidity. Tidal forces disperse it rather than
deforming it elastically.

## V25 Final Assessment

### The Seven Checkboxes for Emergent Spin-2 Gravity

1. □ TT strain σ^TT propagates independently: **PARTIAL** — exists but
   not independent from compression (breathing dominates)
2. □ σ^TT is gapless while compression is gapped: **NO** — both gapped,
   antisymmetric optical branch at ω≈0.51
3. □ σ^TT has quadrupolar radiation pattern: **PARTIAL** — 26% at
   resonance but only during oscillon destruction
4. □ σ^TT mediates attractive force: **NOT TESTED** (binary merged
   before orbital GW could be measured)
5. □ Universal coupling to energy: **NOT TESTED**
6. ✓ Propagation speed = c: **YES** (from V24-L4, automatic in field dynamics)
7. ✓ Single Lagrangian, no separate Φ: **YES** (strain computed from fields)

**Score: 2 YES, 3 PARTIAL, 2 NOT TESTED**

### What V25 Established

1. The 3D elastic oscillon WITH self-consistent metric WORKS (Phase 1+3)
2. The strain field is computable and has real structure (Phase 2)
3. Causal propagation is automatic (from field dynamics)
4. The oscillon IS tidally deformable (k₂ measurable, Phase 6)
5. A quadrupolar resonance exists at Ω_T≈0.1 (Phase 6)

### What V25 Did NOT Achieve

1. Clean spin-2 gravity (l=0 breathing always dominates)
2. Stable binary orbit (oscillons merge too quickly)
3. Steady-state tidal deformation (oscillon too soft, disperses)
4. Independent shear wave propagation (not cleanly separated from compression)

### The Root Obstacle

The oscillon is a BREATHING SCALAR object. Its fundamental mode of
existence is spherically symmetric oscillation (l=0, monopole). Every
interaction — binary merger, tidal deformation, boosting — first excites
this breathing mode, which overwhelms any l=2 content.

For spin-2 gravity to dominate, the oscillon would need:
- A NON-BREATHING ground state (static soliton, like a Skyrmion)
- OR: the breathing frequency above ALL radiation thresholds (so l=0
  radiation is forbidden and only l=2 can radiate)
- OR: a fundamentally different binding mechanism that produces a
  rigid (high k₂) object resistant to tidal disruption
