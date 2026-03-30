# V50/C3: V44 + Chiral Helicity + Cosserat Strain Energy

## Base: V44 equations (clean, geometric theta from curl source)

```
∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η curl(θ)_a
∂²θ_a/∂t² = ∇²θ_a                     + η curl(φ)_a
```

## Problem: theta accumulates without limit

The curl coupling is asymmetric in practice:
- phi → theta: braid core has LARGE, concentrated curl → strong source
- theta → phi: box-filling theta has SMALL, diffuse curl → weak return

Energy flows out of the braid efficiently but returns inefficiently.
Theta fills the box and overwhelms the simulation.

Root cause: V44's curl coupling provides a SOURCE for theta but no
GEOMETRIC CONSTRAINT. Theta isn't required to match the actual twist
of the phi field. Once sourced, it propagates freely.

## Fix: Cosserat strain energy

In a real Cosserat elastic solid, the rotation field θ is constrained
to match the displacement gradient. The strain energy penalizes the
mismatch between the actual geometric twist (curl(φ)/2) and the
rotation field (θ):

    E_strain = α |curl(φ)/2 - θ|²
             = α Σ_a (curl(φ)_a/2 - θ_a)²

This single term does three things:

1. AT THE BRAID (curl(φ) large):
   θ is driven toward curl(φ)/2 → theta follows the geometry.
   The "correct" amount of theta is determined by the twist.

2. IN THE FAR FIELD (curl(φ) ≈ 0):
   E_strain = α|θ|² → theta is PENALIZED for existing without
   corresponding twist. Effectively a position-dependent mass that's
   large where there's no twist and zero where theta matches the twist.

3. RETURN PATH:
   If θ deviates from curl(φ)/2, the energy cost drives it back.
   Energy that escaped into "free theta" is recaptured by the
   geometric constraint. This is the missing dissipation mechanism.

## Forces from the Cosserat strain term

    E_strain = α Σ_a (curl(φ)_a/2 - θ_a)²

Define the mismatch vector:
    M_a = curl(φ)_a/2 - θ_a

### Force on theta (simple, algebraic):

    F_θ_a = -∂E_strain/∂θ_a = +2α M_a = 2α(curl(φ)_a/2 - θ_a)
           = α curl(φ)_a - 2α θ_a

This is a spring pulling theta toward curl(φ)/2:
- If θ < curl(φ)/2: F > 0, drives theta up (toward the geometric value)
- If θ > curl(φ)/2: F < 0, drives theta down (back to geometric value)
- At equilibrium (θ = curl(φ)/2): F = 0 (no force, geometrically aligned)

The term -2α θ_a acts as an effective mass m_θ² = 2α. But unlike a bare
mass, it has a COMPENSATING source term (α curl(φ)_a) that cancels at
the geometrically correct value.

### Force on phi (involves curl of the mismatch):

E_strain depends on curl(φ), which involves derivatives of φ. The
Euler-Lagrange variation gives:

    F_φ_a = +∂_j[∂E_strain/∂(∂_j φ_a)]

Since E_strain = α Σ_b (curl(φ)_b/2 - θ_b)² and curl(φ) depends on
∂_j φ_a, the variation produces:

    F_φ_a = α × curl(M)_a = α × curl(curl(φ)/2 - θ)_a

This can be computed numerically as: take the mismatch M_a at each
voxel, then compute its curl using the standard stencil.

Expanding: F_φ_a = (α/2) curl(curl(φ))_a - α curl(θ)_a

The first term (α/2 curl(curl(φ))) is the twist stiffness — it resists
changes to the rotational content of phi. This is the "displacement
stiffening" that prevents blobs.

The second term (-α curl(θ)) looks like the existing η curl(θ) coupling
but with coefficient -α instead of +η. At equilibrium (θ = curl(φ)/2),
these combine: the net curl(θ) force on phi is (η - α) curl(θ). If
α > η, the sign flips and theta actually OPPOSES phi twist changes
(stability). If α < η, the existing coupling dominates (radiation).

## Combined equations (V44 + chiral + Cosserat strain)

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a
                 + η curl(θ)_a                   [V44 curl coupling]
                 + κ_h chiral_force(P,φ)_a       [C2 chiral helicity]
                 + α curl(M)_a                   [Cosserat strain, M = curl(φ)/2 - θ]

    ∂²θ_a/∂t² = ∇²θ_a
                 + η curl(φ)_a                   [V44 curl coupling]
                 + 2α M_a                        [Cosserat strain, pulls toward curl(φ)/2]

Where M_a = curl(φ)_a/2 - θ_a (the geometric mismatch).

## Implementation

The mismatch M_a = curl(φ)_a/2 - θ_a is computed at each voxel from
existing quantities (curl already computed for the eta term).

The theta force is trivial: +2α × (curl(φ)/2 - θ) = α curl(φ) - 2α θ.

The phi force requires curl(M), which means computing curl of the
mismatch field. This needs M at the 6 neighbors. Two approaches:

Option A (two-pass): First pass computes M at all voxels. Second pass
computes curl(M) and applies the force.

Option B (inline): Compute curl(M) directly from neighbor values of
curl(φ) and θ. Since curl(M) = curl(curl(φ)/2 - θ) = curl(curl(φ))/2
- curl(θ), and curl(θ) is already computed (for eta coupling), we just
need curl(curl(φ)) which is the "vector Laplacian correction."

curl(curl(φ))_a = ∂_a(∇·φ) - ∇²φ_a

But ∇·φ requires the divergence of phi, which is another stencil pass.
Actually simpler: just compute curl of the mismatch directly from
neighbor values.

Best approach: Option A. Store M_a in a temporary array (3 × N³ doubles),
compute curl(M) from the stencil. Cost: one extra array + one stencil
pass. Modest overhead.

## Parameters

    α: Cosserat strain coupling strength (NEW, start 0.1-1.0)
    κ_h: chiral helicity coupling (from C2, start 0.5-1.5)
    η: curl coupling (V44, 0.5)
    m², μ, κ: standard V44 params

Total: 6 physics params (V44's 4 + κ_h + α). Clean.

## What this should fix

1. Theta accumulation: theta that doesn't match curl(φ)/2 is penalized
   by 2α|θ|². Energy returns to phi instead of accumulating.

2. Theta follows geometry: at equilibrium, θ = curl(φ)/2 everywhere.
   The theta field IS the twist, not an independent wave.

3. Anti-blob (from C2 chiral term): braids resist merger because the
   helicity topology is stiffened.

4. EM propagation: the coupled theta-phi wave still propagates, but
   theta is now geometrically constrained. "Light" is a coupled
   displacement-rotation wave where theta tracks curl(phi)/2. This IS
   how EM works in a Cosserat medium.

5. Massless photon: in the far field, curl(φ) ≈ 0 and θ ≈ 0 (at
   equilibrium). The effective theta mass from the strain term is 2α,
   so theta IS massive in the background. BUT: this mass is exactly
   cancelled by the source term (α curl(φ)) at the geometric value.
   A propagating curl(φ) wave automatically generates the matching θ
   wave at no energy cost — the photon propagates as a coupled mode,
   not as a free theta wave.

   Note: the "photon" is no longer a pure theta wave. It's a coupled
   (curl(φ), θ) wave where θ = curl(φ)/2. The speed might differ from
   c — needs to be measured from the dispersion relation.
