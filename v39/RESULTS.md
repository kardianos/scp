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

## Issues Found During V39

### CRITICAL: collapse_3d.c was 3-field only
The collapse_3d.c simulator built by the agent used ONLY 3 phi fields — no theta,
no curl coupling, no eta parameter. This means:
- **ALL V39 collapse results are suspect** (inverse coupling, density-dependent κ)
- The "stable collapsed core" at γ=10 may not exist with 6-field physics
- The plain braid run also showed no theta waves (because theta wasn't simulated)
- The torus + κ-dependent run (γ=7) is also suspect

This was caught when the user noticed the braid SFA showed no theta field and
the structure looked wrong compared to previous versions.

**Root cause**: The agent prompt for collapse_3d.c focused on "configurable mass
coupling" and didn't explicitly require the full 6-field equation. The code was
based on a simplified template, not on the reference v33_cosserat.c.

**Fix applied**: Added to CLAUDE.md and memory: ALL simulations MUST use the full
6-field Cosserat equation. Verification checklist: 18 arrays, curl terms, ETA
parameter, theta_rms growing from zero.

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

### Density-Dependent κ Collapse (NEEDS RE-VERIFICATION)
The 1D results showed:
- γ=10: 220,000× E_pot growth (massive collapse)
- γ=5: dispersal
- Threshold between 5 and 10

The 3D results showed:
- γ=10: stable collapsed core (E_pot = -218,000, nearly static)
- γ=7 + torus: collapsed to a featureless blob (aspect → 1.0)
- γ=5: dispersed

**These 3D results are ALL from the 3-field simulator and must be re-done with
the full 6-field equation before they can be trusted.** The theta coupling may:
(a) prevent the collapse entirely
(b) allow collapse but add internal structure
(c) not affect the collapse (unlikely given η=0.5 is a significant coupling)

## Files

| File | Description | Status |
|------|-------------|--------|
| `data/original_braid.sfa` | V34 plain braid, 6-field, confirmed | GOOD |
| `data/kappa_3d.sfa` | κ-dependent collapse, spherical BC | 3-FIELD ONLY — SUSPECT |
| `data/kappa_3d_periodic.sfa` | κ-dependent collapse, periodic BC | 3-FIELD ONLY — SUSPECT |
| `data/kappa_3d_sphere.sfa` | κ-dependent collapse, spherical abs BC | 3-FIELD ONLY — SUSPECT |
| `data/plain_braid.sfa` | Plain braid, 3-field only | WRONG (no theta) |
| `data/braid_6field.sfa` | Truncated braid, 6-field | WRONG GEOMETRY (truncated ≠ original) |
| `data/torus_kappa7.sfa` | Torus + κ-dep, 3-field | 3-FIELD ONLY — SUSPECT |
| `src/collapse_3d.c` | 3D simulator | MISSING THETA — needs full rewrite |
| `src/oscillon_1d.c` | 1D simulator | 3-field only (acceptable for 1D exploration) |
| `src/critical_density_1d.c` | 1D phase scanner | 3-field only |
| `THEORY_density_kappa.md` | Theoretical analysis | Valid (math, not simulation) |

## Next Steps (for future versions)

1. **Rewrite collapse_3d.c** with full 6-field Cosserat (use v33_cosserat.c as base)
2. **Re-run κ-dependent collapse at γ=10** with 6-field equation
3. **Re-run torus + κ at γ=7** with 6-field equation
4. All new simulators MUST be verified against the original braid (theta_rms > 0)
