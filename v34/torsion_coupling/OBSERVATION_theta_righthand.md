# Observation: θ Field Follows Right-Hand Rule Around Braid

**Date**: 2026-03-20
**Simulation**: cosserat_sfa, N=80, L=25, T=50, η=0.5, m_θ²=0

## What Was Observed (Volumetric Viewer)

1. The braid travels along z (helical axis), wrapping in periodic BC
2. The θ (angle) field forms circular wave patterns PERPENDICULAR to
   the braid's travel direction
3. The pattern follows the RIGHT-HAND RULE: curl fingers in the
   direction of θ circulation, thumb points in the braid's travel direction
4. Discernible waves propagate outward from the braid or are maintained
   as a standing pattern

## Why This Matters

This is the spatial pattern of a MAGNETIC FIELD around a current-carrying
wire. The braid = current (moving topological charge). The θ field = B-field
(circular, perpendicular, right-hand rule).

This pattern was NOT imposed. It EMERGES from:
- The curl coupling: ∂²θ/∂t² = ∇²θ + η×curl(φ)
- The braid's helical structure (nonzero curl(φ) along the helix axis)
- The massless propagation of θ (m_θ=0 → no Yukawa suppression)

The curl operator inherently produces rotational patterns perpendicular
to the source — this is the mathematical origin of the right-hand rule.
Adding it as a physical coupling between position and angle fields makes
the right-hand rule a DYNAMICAL consequence, not an axiom.

## Caveats

- The braid is initialized with a z-preferred helical structure
- Periodic BC in z means the braid wraps → the "current" is continuous
- Not yet quantified: does the θ field decay as 1/r from the braid axis?
- Not yet tested: does the pattern reverse for opposite winding (W=-1)?
- The observation is visual from the volume viewer, not yet numerically
  confirmed with angular decomposition

## Next Steps

1. Quantitative measurement: decompose θ into cylindrical components
   around the braid axis, verify circular symmetry
2. Radial decay: θ_circular(r) should fall as 1/r (Biot-Savart analog)
3. Winding test: W=-1 braid should produce opposite circulation
4. Two braids: parallel braids (same winding) should attract via θ overlap
   (like parallel currents), antiparallel should repel

## Data

- SFA archive: data/sfa_hires.sfa (264 frames, N=80, L=25)
- Viewer: ./volview_sfa data/sfa_hires.sfa
- Best viewed: press O twice for translucency, toggle 3 (blue only)
  to isolate θ, rotate to view perpendicular to z-axis
