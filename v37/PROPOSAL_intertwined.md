# Intertwined Compact Braids: First Principles

## The Co-Location Principle

V(P) needs φ₀, φ₁, φ₂ simultaneously nonzero at the same points.
The helix achieves this: three phases of ONE wave on ONE tube.
Rings fail: three fields on THREE separate tubes → P ≈ 0.

**Solution**: ALL three fields must travel along the SAME closed path,
intertwined as phase offsets — not as separate geometric objects.

## Candidate 1: Toroidal Braid

A single torus tube carrying all three fields, wound along the torus
circumference with phase offsets:

    φ_a(x) = A × env(d_tube) × cos(k × θ_tor + δ_a) + A_bg × cos(...)

where θ_tor = atan2(y,x) is the toroidal angle, d_tube is distance
from the torus center-circle, and k = n/R_major for n full oscillations.

Key difference from the failed "loop" geometry: the failed version used
k_twist=3 (three twists of the CROSS-SECTION). This version uses k=1
(one oscillation ALONG the torus), identical to the helix but bent
into a ring. The physics is the same as the helix — three phase-shifted
cos() on one tube — just curved into a circle.

For closure: k × 2πR = 2πn → k = n/R. With R=4, n=2: k = 0.5.
The three fields complete 2 full oscillations around the ring.

## Candidate 2: Triple-Strand Helix (Truncated)

Three helical sub-tubes wound around a common axis, each carrying
primarily one field but overlapping with the others:

    Tube a center: (R_h cos(kz + 2πa/3), R_h sin(kz + 2πa/3), z)

    φ_a(x) = A × Σ_b env(d(x, tube_b)) × weight(a,b) × cos(kz + δ_a)

With R_h small (< tube radius), the three tubes overlap significantly.
The triple product P is large in the overlapping core.
Truncated with z-envelope for compactness.

This is a physical braid — three strands twisted together. The
co-location comes from the overlap of the strands.

## Candidate 3: Spherical Harmonic Modes

Three l=1 spherical harmonics as the three fields:

    φ₀ = A × f(r) × cos(θ)              (Y₁₀ — dipole along z)
    φ₁ = A × f(r) × sin(θ)cos(φ)        (Y₁₁ — dipole along x)
    φ₂ = A × f(r) × sin(θ)sin(φ)        (Y₁₋₁ — dipole along y)

    f(r) = exp(-r²/2σ²) × cos(kr)       (radial envelope × oscillation)

Properties:
- Naturally spherical and compact
- All three fields nonzero almost everywhere (except on symmetry planes)
- P = φ₀φ₁φ₂ ∝ cos(θ)sin(θ)cos(φ)sin(θ)sin(φ) = sin³θ cosθ sinφ cosφ
  This is nonzero in most of the sphere — co-location achieved!
- No topological winding (unlike the helix) — may lack protection

The radial oscillation cos(kr) replaces the axial cos(kz). Combined
with the angular structure, this creates a 3D pattern with no
preferred axis.

## Candidate 4: Toroidal Braid (Three Strands on One Torus)

Combine candidates 1 and 2: a torus carrying three intertwined strands.
Each strand winds helically around the torus tube cross-section AS the
tube goes around the torus. Like a braided rope bent into a circle.

    strand_a center (in tube cross-section):
        (r_s cos(m×θ_tor + 2πa/3), r_s sin(m×θ_tor + 2πa/3))

    where r_s is the strand offset, m is the number of braiding twists.

    φ_a strong near strand_a, but nonzero near all strands due to overlap.

This has THREE kinds of structure:
1. Toroidal closure (compact, no z-drift)
2. Braiding of strands (topological protection)
3. Phase oscillation along the torus (V(P) binding)

## Candidate 5: Single Ring, Single Oscillation

The SIMPLEST possible closed braid: one torus tube, one wavelength,
three phase-shifted fields. Like the helix but exactly one period
bent into a circle.

    R_major = λ/(2π) where λ = 2π/k_helix
    φ_a = A × env(d_tube) × cos(θ_tor + δ_a) + background

With k_helix = π/L ≈ 0.21 (current value): R_major ≈ 4.8.
With tube_radius = 2, this fits in L=12 easily.

The failed "loop" used k_twist (twisting the cross-section) instead
of oscillating ALONG the torus. This version oscillates along, like
the helix — just closed.

## Priority Order

1. **Candidate 5** (single ring, single osc) — simplest, most direct
   translation of the helix to closed form
2. **Candidate 1** (toroidal braid) — same idea, two oscillations
3. **Candidate 3** (spherical harmonics) — radically different, naturally
   spherical, but may lack topological protection
4. **Candidate 2** (triple-strand truncated) — physical braid, but
   still has z-axis like the truncated helix
5. **Candidate 4** (toroidal braid, three strands) — most complex
