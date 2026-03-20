# Phase 1a Results: Analytical Mode Decomposition

## Summary

The linearized equation around the uniform background has THREE modes,
all massive (m_eff ≈ 1.5). There is NO massless Goldstone mode at the
standard background amplitude A_bg = 0.1. The mode splitting from V(P)
is only Δm² = 0.004, negligible compared to m² = 2.25.

However, the splitting grows DRAMATICALLY with local field amplitude.
Near the braid surface (A ≈ 0.5), Mode 1 becomes nearly massless
(m_eff² = 0.55). At the braid core (A ≈ 1.0), Mode 1 goes tachyonic.
The three modes experience fundamentally different physics near braids.

## Mode Identification

At any point in the background, the coupling matrix M_ab = ∂²V/∂φ_a∂φ_b
has three eigenmodes:

**Mode 1 (Amplitude)**: Perturbation aligned with the local field vector.
Changes Σφ² (field density) and P (triple product). This mode DEPLETES.
- Far field (A=0.1): m_eff² = 2.247 (massive, same as bare m²)
- Braid surface (A=0.5): m_eff² = 0.55 (nearly massless!)
- Braid core (A=1.0): m_eff² = -0.59 (tachyonic → drives braid binding)

**Mode 3 (Phase)**: Perturbation perpendicular to the local field vector.
Does NOT change Σφ² or P. This mode is "neutral" — no depletion.
- Far field: m_eff² = 2.251 (massive, essentially = bare m²)
- Braid surface: m_eff² = 2.59 (HEAVIER than bare m²)
- Braid core: m_eff² = 5.05 (much heavier)

**Mode 2 (Mixed)**: Intermediate character. Slightly above bare m².

## Implications

### For the phonon / long-range gravity:
The power-law depletion (δρ ∝ 1/r^1.2) measured in the phonon test is
NOT a linearized Goldstone mode — all linear modes are massive in the
far field. The depletion is a NONLINEAR collective effect originating
from the braid's surface region where Mode 1 becomes nearly massless.

The braid's surface (A ≈ 0.3-0.5) acts as an "antenna" that radiates
Mode 1 perturbations with very low effective mass. These perturbations
propagate outward with a Yukawa range of 1/√0.55 ≈ 1.35 (not 0.67).
The cumulative effect of this radiation from the extended braid surface
creates the observed power-law halo.

### For EMF:
The phase mode (Mode 3) is distinct from the amplitude mode (Mode 1)
but has the SAME mass in the far field. Near the braid, Mode 3 is
HEAVIER than Mode 1 — the phase perturbation is more confined, not less.

This means: if EM is mediated by the phase mode, it would be a
SHORT-RANGE force near braids (stronger confinement than gravity),
not a long-range force. This is more like the strong nuclear force
than electromagnetism.

For long-range EM, we need either:
1. A mode that becomes LIGHTER (not heavier) near braids — none found
2. A nonlinear mechanism analogous to how gravity works
3. Complex fields / gauge coupling (the original U(1) path from V29-T10G)

### The amplitude-dependent mode splitting:

| A_local | m²(Mode 1) | m²(Mode 3) | Splitting | Physics |
|---------|-----------|-----------|-----------|---------|
| 0.01 | 2.250 | 2.250 | 0.000 | All modes identical |
| 0.10 | 2.247 | 2.251 | 0.004 | Negligible splitting |
| 0.20 | 2.201 | 2.259 | 0.058 | Small splitting |
| 0.50 | 0.550 | 2.587 | 2.037 | Mode 1 nearly massless |
| 1.00 | -0.587 | 5.053 | 5.640 | Mode 1 tachyonic |

The braid creates a "lens" where Mode 1 propagates easily (low mass)
but Mode 3 is confined (high mass). This is gravitational lensing
from the mode perspective — the braid bends Mode 1 perturbations
toward itself while deflecting Mode 3 perturbations away.

## Files

- `mode_analysis.py` — analytical computation of M_ab eigenvalues
- This file — results and interpretation
