# V27 Milestone 3: Electromagnetic Structure from Torsion — RESULTS

## Configuration
- m=0, mu=-50, kappa=50 (V27-M4 optimal massless propagating braid)
- A0=0.8, R_tube=3.0, N=128, L=20 (single braid), L=40 (two-braid)
- Periodic z, absorbing x,y boundaries
- Settling time: t=200 before measurements (M3a, M3b)

## M3a: Torsion Field Mapping

### Key Numbers
- |Omega|_max on z=0 plane: 0.818
- |Omega_z|_max on z=0 plane: 0.661
- |Omega|_max on x=0 plane: 1.108

### Spatial Structure: NOT LOCALIZED

The vorticity |Omega| shows a **flat radial profile** from r=0 to r~10,
with magnitudes ~0.25-0.30 throughout. It then RISES to a peak of ~0.63
at r~15 (near the absorbing boundary inner edge at r=14) before decaying
exponentially beyond r~16.

| r range | |Omega|_avg |
|---------|------------|
| 0-2     | 0.28-0.30  |
| 2-6     | 0.16-0.25  |
| 6-10    | 0.12-0.15  |
| 10-14   | 0.12-0.30  |
| 14-16   | 0.30-0.63  |
| 16-20   | <0.05      |

**Localization fraction**: Only 7.7% of |Omega|^2 lies within r < 2*R_tube = 6.
The vorticity is **delocalized** — it extends throughout the entire domain,
including radiation that has propagated to the absorbing boundary region.

### Omega_z is ZERO

The azimuthally-averaged Omega_z component is identically zero to machine
precision (~1e-14) at ALL radii. This means the torsion field has **no net
z-component** when averaged over angle. The individual Omega_x, Omega_y, Omega_z
components have similar magnitudes locally (~0.3), but Omega_z averages to zero.

**Physical interpretation**: The helical braid has a vorticity field that rotates
in the xy-plane as a function of z (following the braid twist), but has no
net axial vorticity. This is consistent with the braid's azimuthal symmetry —
the three strands produce local vorticity that cancels when integrated over
a full cross-section.

## M3b: Torsion Flux Through Cross-Sections

### Key Result: FLUX IS ZERO

| Slice  | z      | Phi_T (total)  | Phi_T (core)   |
|--------|--------|----------------|----------------|
| z=-L   | -20.00 | -1.66e-13      | 4.47e-14       |
| z=-L/2 | -9.92  | -5.20e-14      | 3.39e-14       |
| z=0    | 0.16   | -4.01e-14      | 9.91e-14       |
| z=L/2  | 10.24  | -1.53e-14      | -2.84e-14      |

- Mean Phi_T = -1.74e-14 (machine zero)
- Std Phi_T = 1.12e-13 (pure noise)
- **Phi_T / (2*pi) = 0.000000**

The torsion flux is identically zero at ALL z-slices. There is no quantized
flux, no flux conservation to test — there is simply no flux.

**Why?** The torsion omega_{ij} = (1/2)(d_i phi_j - d_j phi_i) requires the
SPATIAL indices i,j to map to DIFFERENT field components a,b. For the helical
ansatz phi_a = A(r)*cos(kz + 2*pi*a/3):

    d_x phi_1 - d_y phi_0

Both terms involve transverse derivatives of the Gaussian envelope times
oscillating phases. When integrated over the full xy-plane, the cos/sin phases
from different strands (offset by 2*pi/3) cancel EXACTLY by symmetry.

This is a **fundamental symmetry obstruction**: a cylindrically symmetric braid
produces zero net torsion flux through any cross-section perpendicular to its axis.

## M3c: Two Braids — Torsion-Mediated Force

### Setup
- Two braids at D=30 separation along x-axis
- L=40, N=128, dx=0.630 (coarser grid for larger domain)
- Config A: same twist (both +k_z)
- Config B: opposite twist (+k_z and -k_z)
- Evolution: t=200

### Results

| Config         | sep(0) | sep(200) | delta_sep | Interpretation |
|----------------|--------|----------|-----------|----------------|
| Same twist     | 30.000 | 22.831   | -7.169    | ATTRACTION     |
| Opposite twist | 30.000 | 24.973   | -5.027    | ATTRACTION     |

**Both configurations attract.** Same-twist attracts MORE strongly (-7.2 vs -5.0).

### Caveats

1. **Energy instability**: Total energy grows dramatically during evolution
   (E: 265 -> 5.0e5 for same-twist, 265 -> 1.1e5 for opposite). The coarse
   resolution dx=0.63 (only ~10 points across R_tube=3) is likely causing
   numerical instability. Results should be interpreted cautiously.

2. **Not torsion-mediated**: Since torsion flux is zero (M3b), the observed
   attraction cannot be torsion-mediated. It is more likely gradient energy
   overlap — two nearby field configurations sharing the same potential well
   attracts via the triple-product potential V(P). This is analogous to
   Skyrmion-Skyrmion binding, not EM.

3. **Same > Opposite**: Same-twist braids attract more because their field
   profiles are coherent (constructive interference in the overlap region),
   while opposite-twist braids have partial destructive interference. This
   is the opposite of EM behavior (like charges repel, unlike attract).

## Summary

| Test | Target | Result | Status |
|------|--------|--------|--------|
| M3a: Vorticity localized? | Yes | No — delocalized, 7.7% in core | FAIL |
| M3b: Quantized flux? | Phi_T = n * 2*pi | Phi_T = 0 (machine zero) | FAIL |
| M3b: Flux conserved along z? | Phi_T(z) = const | Phi_T = 0 everywhere | N/A |
| M3c: Opposite twist attracts? | Yes (unlike charges) | Both attract, same > opposite | FAIL |

## Conclusion

**The torsion field of a helical braid does NOT produce electromagnetic structure.**

The fundamental problem is that the cylindrical symmetry of the braid ansatz
(three strands at 2*pi/3 phase offset, Gaussian radial envelope) guarantees
ZERO net torsion flux by angular cancellation. This is not a numerical issue —
it is an exact symmetry property.

For EM to emerge from torsion, the braid would need:
1. **Broken azimuthal symmetry** — e.g., an asymmetric radial profile, or
   strands that are not uniformly spaced at 2*pi/3.
2. **Non-abelian structure** — the simple triple-product coupling treats all
   three fields democratically. A coupling that distinguishes between strands
   could break the cancellation.
3. **A different torsion definition** — perhaps using the triple product P
   itself as the "charge density" rather than the antisymmetric gradient tensor.

The inter-braid force observed in M3c is real but is NOT electromagnetic in
character. It is binding via gradient/potential overlap, analogous to nuclear
binding, and has the wrong sign pattern (same attracts more than opposite).

## Parameters
- mu=-50, kappa=50, m=0, A0=0.8, R_tube=3.0
- N=128, L=20 (M3a,b), L=40 (M3c), cfl=0.20
- BC: periodic z, absorbing x,y
