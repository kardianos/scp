# V37: Compact Braid Geometries

## Problem

The current braid is a helical TUBE — infinite along z, wrapping through
periodic boundaries. It has length >> width. This is unphysical: a proton
is roughly spherical (aspect ratio ~1).

The tube geometry also affects the physics:
- **Gravity**: the depletion profile is cylindrical, not spherical.
  The measured 1/r^1.2 exponent may be an artifact of the tube geometry
  rather than a fundamental prediction. Spherical depletion should give
  closer to 1/r or 1/r² depending on dimension.
- **EM**: the θ field follows the right-hand rule around the z-axis
  because the helix IS along z. A compact braid would have isotropic
  θ radiation.
- **Mass**: the braid's energy depends on the box length L (longer tube
  = more energy). A compact braid has definite mass independent of L.
- **Scattering**: two tube-braids interact over their full length.
  Compact braids would interact as point-like objects at large separation.

## Goal

Find a COMPACT (spherical or near-spherical) braid topology that retains
all discovered properties:
1. Stable, self-sustaining (survives T=500+)
2. Triple product |P| > 0 (nonlinear binding)
3. Topological winding (conserved charge)
4. Depletes the surrounding field (gravity mechanism)
5. Sources the θ field through curl coupling (EM mechanism)
6. Distinct winding handedness W=±1 (charge conjugation)

## Candidate Geometries

### 1. Borromean Rings (V26, most promising)

Three interlocking rings, one per field, in orthogonal planes.
V26 results: fc=0.92, |P| growing 70×, survived T=500.

**Advantages**: Already tested and stable. Naturally compact (no preferred
axis). Three-fold symmetry matches the three fields.

**Issues**: V26 used a different equation (no mass term, no θ fields).
Need to verify survival in the full Cosserat equation (Eq. 10) with
m²=2.25 and θ coupling.

**Test**: Initialize Borromean rings in the Cosserat code. Does the
compact topology survive? Is the depletion spherical?

### 2. Torus Knot

A helix that closes on itself by wrapping around a torus.
The (p,q) torus knot wraps p times through the hole and q times
around the tube. The (2,3) trefoil knot is the simplest.

**Advantages**: Compact, finite energy, definite winding number.
The knot topology provides topological protection.

**Issues**: Complex initialization. The three fields need to be
arranged along the knot with proper phase offsets.

### 3. Hopfion

A field configuration where the preimage of every point on the target
sphere is a circle (Hopf fibration). Used in condensed matter, BECs.

**Advantages**: Compact, topologically protected (Hopf charge ≠ 0),
inherently 3D (no preferred axis).

**Issues**: The triple-product potential V(P) may not support Hopfions.
V2/hopfion_composition explored this but with different physics.

### 4. Skyrmion (Hedgehog)

The original B=1 Skyrmion from the CHPT work (V2-V12). A radial
hedgehog configuration where each field component maps to a spatial
direction.

**Advantages**: Well-studied, compact, stable (with mass term).
We have extensive Skyrmion solver code from V2.

**Issues**: The Skyrmion uses a different Lagrangian (Skyrme model
with 4-derivative term). Our V(P) triple-product potential is
different. The hedgehog may not be a solution of Eq. (1).

### 5. Truncated Helix

Keep the helical structure but TRUNCATE it with a smooth envelope:
instead of extending to ±L, confine the helix to a compact region
of length ~R_tube (comparable to the tube radius).

    φ_a = A × env(z) × env_xy(x,y) × cos(kz + δ_a)
    env(z) = exp(-z²/(2σ_z²))  with σ_z ≈ R_tube

This creates a "braid sausage" — finite, compact, with definite ends.

**Advantages**: Simplest modification of current init. Smooth, no
topological complications. Easy to test.

**Issues**: No topological protection at the ends (the helix just
fades out). May unwind from the ends inward. The winding number
may not be conserved without periodic BC in z.

### 6. Double Helix (Closed Loop)

Two half-helices joined at both ends, forming a closed loop.
Like a DNA double helix bent into a circle.

**Advantages**: Closed, compact, definite winding. Combines
helical binding with closed topology.

**Issues**: Requires a large enough radius to fit the helix
wavelength around the loop. Complex initialization.

## Evaluation Criteria

For each candidate:

| Test | What to measure | Target |
|------|----------------|--------|
| Survival | E_pot retention at T=500 | > 50% |
| Compactness | Aspect ratio (max/min extent) | < 2 |
| Depletion | δρ(r) profile symmetry | Isotropic |
| Winding | Phase winding conservation | W = ±1 stable |
| θ coupling | θ_rms in Cosserat equation | > 0.02 |
| Mass | Total energy (L-independent) | Definite |

## Priority Order

1. **Truncated Helix** (simplest, tests if compactness is feasible)
2. **Borromean Rings** (already proven stable, needs Cosserat test)
3. **Torus Knot** (topologically protected compact braid)
4. **Hopfion** (if others fail, most exotic)
5. **Skyrmion** (different physics, may not apply)

## Expected Impact on Physics

### Gravity exponent

A spherical braid produces spherical depletion → δρ(r) ∝ 1/r^n
where n depends on 3D geometry:
- Monopole source in 3D: n=1 (Coulomb/Newton)
- The current n=1.2 may be from the cylindrical tube geometry
- Spherical braid might give n closer to 1.0 or 2.0

This would directly address F1 (force law exponent) and potentially
bring gravity into exact Newtonian agreement.

### EM pattern

A compact braid has no preferred axis → the θ field would be
radiated isotropically (or with multipole structure determined by
the braid's internal geometry, not the box orientation).

### Mass spectrum

Compact braids of different topology (trefoil, figure-eight, etc.)
would have different masses → the beginning of a particle spectrum.

## Data Format

All V37 simulation output MUST use the SFA format.
