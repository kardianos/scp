# V50/C4: V44 + Cosserat Strain + Curl-Squared Hardening

## Base: C3 (V44 + Cosserat strain α|curl(φ)/2 - θ|²)

C3 showed that the Cosserat strain constraint keeps theta geometrically
aligned with the braid's twist. Energy conservation verified at α=0.1.
Theta develops spatial structure that follows curl(φ) instead of filling
the box uniformly.

## Problem: No shell hardening

The Cosserat constraint prevents theta saturation but doesn't create a
rigid boundary that divides interior from exterior. Two approaching
braids can still merge because there's no energetic barrier at the
surface where their fields overlap.

## Addition: Curl-squared hardening

    L_harden = -(β/2) |θ|² |∇×φ|²

Physical interpretation: where the fabric is TWISTED (|∇×φ|² large)
AND there's ROTATION present (|θ|² large), the medium becomes rigid.
Neither twist alone nor rotation alone creates stiffness — both must
be present together. This naturally forms a shell at the braid surface.

### Why |∇×φ|² and not P² or |φ|²

| Order parameter | In vacuum | At core | At surface | Strobes? |
|----------------|-----------|---------|------------|----------|
| P² | ~3×10⁻⁸ | large | moderate | Yes (6ω) |
| \|φ\|² | 0.015 (39% ripple) | 0.1-0.3 | moderate | Yes (2ω) |
| **\|∇×φ\|²** | **0 exactly** | moderate | **MAXIMUM** | **No** |

|∇×φ|² is zero for the 1D carrier wave (∂_x φ = ∂_y φ = 0). It's
maximum at the braid SURFACE where the helical twist has the steepest
spatial gradients. This creates hardening exactly at the interface —
not at the core (P²) or everywhere (|φ|²).

### What it does to each field

Theta force: -β|∇×φ|²θ_a
  → Theta becomes heavy where twist is strong
  → Theta confined to a shell at the braid surface
  → Massless in vacuum (|∇×φ|=0) → photon propagates at c

Phi force: from EL variation of |∇×φ|², involves ∇×(β|θ|²∇×φ)
  → Twist stiffened where theta is present
  → Creates the delta-rho displacement at the shell
  → Two overlapping shells repel (hardened interface resists merger)

### The shell structure

    Interior (r < R_core):
      P large, |∇×φ| moderate, θ ≈ curl(φ)/2 (from Cosserat)
      → Strong binding, theta follows geometry

    Surface (r ≈ R_core):
      P transitioning, |∇×φ| MAXIMUM, θ peaked
      → Maximum hardening from β|θ|²|∇×φ|²
      → This IS the shell — rigid, divides in/out
      → Creates delta-rho from the phi stiffening

    Exterior (r > R_core):
      P ≈ 0, |∇×φ| ≈ 0, θ → 0 (Yukawa from Cosserat mass)
      → Normal vacuum, massless photon, 1/r EM

### Parameters

    α:  Cosserat strain (from C3, 0.1)
    β:  Curl-squared hardening (NEW, start 0.1-1.0)
    κ_h: Chiral helicity (from C2, set to 0 initially — test one term at a time)
    η, m², μ, κ: V44 standard

## EL Derivation

The theta force is algebraic (simple):
    F_θ_a = -β|∇×φ|²θ_a

The phi force requires variation of |∇×φ|² w.r.t. ∂_j φ_a.
This has the same structure as the Cosserat strain but with |θ|²
as a coefficient. Specifically:

    ∂(|∇×φ|²)/∂(∂_j φ_a) = 2 × curl(φ)_b × ∂(curl_b)/∂(∂_j φ_a)

The result is: F_φ_a = -∇×(β|θ|² ∇×φ)_a

This is the NEGATIVE curl of (β|θ|² times the curl of phi).
On the grid, this requires computing β|θ|²curl(φ) at each voxel,
then taking its curl using the standard stencil — a two-pass
computation similar to C3's Cosserat mismatch.

TODO: Verify sign with Maxima before running.
