# Spherical Seed Postmortem

## What Happened

The spherical phase soliton seed (gen_spherical_seed.c) initially showed
promising signs — P_int oscillated 6-34 with periodic BC, phi_max grew
to 1.69 — but closer analysis reveals the seed is DISPERSING, not
condensing.

### Evidence of Dispersal

**P_peak collapses**:
- t=0: P_peak = 0.004 (seed value)
- t=20: P_peak = 0.33 (initial resonance — the "amazing" frame)
- t=40: P_peak = 0.022 (15× drop from peak)
- t=60: P_peak = 0.013
- t=80: P_peak = 0.011

**Cluster structure dilutes**:
- t=0: 10 shells centered at origin (seed structure)
- t=20: 18 clusters, dominant at center (P_peak=0.037)
- t=40: 45 clusters scattered across grid (fragmented)
- t=60: 17 clusters, one fills 727K/2M voxels (homogenized)
- t=80: 5 clusters, one fills 772K voxels (nearly uniform)

**E_pot stays positive and GROWS**:
- t=0: +1095 (no binding well)
- t=4: +541 (temporary dip during resonance)
- t=80: +467 (settles to a positive background value)

E_pot should be NEGATIVE for a bound structure. The V(P) potential with
μ=-41.345 gives E_pot < 0 when P is large enough. Our seed never
reaches sufficient P concentration.

**E_total increases** (periodic BC):
- t=0: 32,679
- t=14: 35,222 (peak, +7.8%)
- t=80: 34,702 (+6.2% above initial)

Energy is being CREATED, not conserved. This 7.8% increase in the first
14 time units indicates the seed is far from equilibrium — the initial
θ configuration is inconsistent with the dynamics and releases energy
as it adjusts.

### The phi_max growth was misleading

phi_max going from 0.19 to 1.69 looked like condensation. But phi_max
measures the MAXIMUM field value at any point, not the concentrated
binding. What happened:
- The initial phase gradient creates curl → θ force → φ acceleration
- The φ field overshoots (phi_max=1.69 at t=24) due to the initial θ impulse
- The overshoot propagates outward as radiation
- The radiation bounces off periodic boundaries and fills the box
- phi_max stays elevated because the box is full of radiation, not
  because a particle condensed

### P_int oscillation was also misleading

P_int oscillated 6-34, which looked like breathing. But for a real proton,
P_int peaks at 200-300 (V43 template after condensation). Our P_int=34
peak is 10× too low — it's not a proton, it's background fluctuation
amplified by the initial transient.

## Why the Ansatz Failed

### Problem 1: No confinement mechanism

The spherical ansatz φ_a = A(r) cos(k·r + δ_a + χ_a·ψ(r)) places all
three field components at the same radial position. At any point r, all
three φ_a are nonzero simultaneously. This means P = φ₀φ₁φ₂ is nonzero
EVERYWHERE within the envelope — there's no spatial structure that
creates P=0 boundaries to confine the binding.

In contrast, the 3-braid structure creates P=0 along each axis (where
only one braid contributes) and P≠0 only in the central overlap region.
The carrier phase cancellation Σcos(δ+Δ_b) = 0 creates a topological
boundary where P vanishes — this is the confinement.

The spherical ansatz has no such boundary. P is simply proportional to
A(r)³ × C_P everywhere within the Gaussian envelope. The "particle" has
no defined edge.

### Problem 2: Wrong topology

A real proton has ANGULAR structure — the three braids are along x,y,z
axes, each contributing a different carrier phase. The carrier phase
cancellation depends on the angular direction (looking along x sees one
phase, along y sees another). This angular structure IS the topology.

The spherical ansatz has φ_a(r) depending only on radius r. There is no
angular dependence. The chirality enters only through the phase argument
χ_a × ψ(r), which modulates the radial oscillation but doesn't create
angular structure. The result is a spherically symmetric blob with no
internal topology.

### Problem 3: Self-consistent θ is necessary but not sufficient

The θ initialization (θ = G × curl(φ)) correctly accounts for the
Cosserat equilibrium. This is why the seed has 260× more θ than the old
braid seeds and shows immediate activity. But self-consistent θ prevents
disruption — it doesn't create confinement. The seed needs BOTH correct
θ AND correct spatial topology.

## What a Correct Seed Needs

1. **Carrier phase cancellation in angular directions**: The three field
   components must have DIFFERENT angular dependence, not just different
   radial phase shifts. At any point on a sphere of radius r, the three
   phases must partially cancel. Full cancellation (Σcos = 0) at the
   boundary defines the confinement surface.

2. **Angular-dependent curl structure**: The curl of φ must have the
   right multipole structure. For UUD: two "Up" braids contribute
   constructive curl_x, one "Down" braid contributes opposite. This
   creates the charge asymmetry. A spherically symmetric ansatz has
   no such multipole structure.

3. **Self-consistent θ** (already implemented): θ = G × curl(φ) with
   G = 0.6/(0.2 + 0.5S). This part works.

4. **Finite breathing amplitude**: The proton is a limit cycle, not a
   fixed point. The seed should have nonzero velocities (phi_vel, theta_vel)
   that match the breathing mode. Zero velocities (our current approach)
   mean the seed starts at a turning point of the breathing cycle — the
   first half-cycle releases energy that radiates away before the structure
   can stabilize.

## Proposed Next Ansatz: Angular Braid Encoding

Instead of φ_a(r) = A(r) cos(k·r + δ_a + χ_a·ψ(r)), use:

  φ_a(r,θ,φ) = A(r) × cos(k·r + δ_a + Δ_braid(θ,φ))

where Δ_braid(θ,φ) encodes the three carrier phases as angular sectors:

  Δ_braid(θ,φ) = Σ_b Δ_b × w_b(θ,φ)

with w_b being smooth angular weight functions peaked along each axis:
  w_0(θ,φ) peaked along x-axis (braid 0)
  w_1(θ,φ) peaked along y-axis (braid 1)
  w_2(θ,φ) peaked along z-axis (braid 2)

Normalized: Σ_b w_b = 1 everywhere on the sphere.

This creates angular-dependent carrier phases: looking from the x-direction,
you predominantly see Δ₀; from y, you see Δ₁; from z, you see Δ₂. In
between, you see a weighted average.

The chirality would enter through the SIGN of the carrier wave in each
angular sector: w_b carries χ_b as a phase modification.

This is the 3-braid structure ENCODED as a sphere — merged from the
start, with the angular topology built in. The key question is whether
the smooth angular weighting produces carrier phase cancellation at the
boundaries.

### Alternative: Hedgehog-like encoding

Use the Skyrme hedgehog structure adapted to three real scalars:

  φ_a(r,θ,φ) = A(r) × cos(k × n̂_a · r⃗ + δ_a)

where n̂_a are the three axis directions: n̂₀ = x̂, n̂₁ = ŷ, n̂₂ = ẑ.

Each field component has its carrier wave propagating along its own axis.
At the center (r=0), all three overlap with the right phase structure.
Along any axis, only one component oscillates — the other two see
n̂_a · r⃗ = 0 and contribute the background.

This IS the 3-braid structure, but written as a smooth field rather than
three separate tubes. The "tube" shape emerges from the cos(k × n̂_a · r⃗)
factor: field component a has maximum amplitude along axis a and decays
as cos(k × r × cos_angle) = cos(0) = 1 along axis a, oscillating and
averaging to zero perpendicular to it.

The difference from gen_phase_confined: envelope. The old generator used
Gaussian tubes; this uses a common spherical envelope A(r) for all three
components, with the angular structure entirely from the cosine argument.

Chirality: cos(k × n̂_a · r⃗) for Up, cos(-k × n̂_a · r⃗) for Down.

## Priority

Test the hedgehog-like encoding first — it's closest to what the
working V43 template actually looks like (three braids merged into a
sphere), but as a single coherent formula rather than three separate
tubes.
