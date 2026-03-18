# V31: Characterize the M7+c(ρ_B) Gravity Model

## Starting Point

V30 established: M7 two-component + c_eff(ρ_B) is stable at all tested
α_c (0-0.5). The B-field shows simultaneous accretion (core c>1) and
depletion (shell c<1). The braid survives with fc=0.84-0.98.

The model:
    S: ∂²S_a/∂t² = c_eff²(ρ_B)∇²S_a - m²S_a - ∂V/∂S_a - g(ΣB²)S_a
    B: ∂²B_a/∂t² = ∇²B_a - m²B_a - g(ΣS²)B_a
    c_eff² = 1 - α_c(1 - ρ_B/ρ₀)

## V31 Goals (three focused tests)

### T1: Radial c_eff Profile — Is the Depletion 1/r?

Measure c_eff(r) and ρ_B(r) as functions of distance from braid center.
Use the differential method (braid+bg vs bg-only control) to isolate
the braid's effect.

- N=128, L=30, T=500, α_c=0.2
- Measure ρ_B(r) and c_eff²(r) at 100 radial bins, time-averaged T>300
- Fit: δc_eff²(r) = c_eff²(r) - 1 ~ A/r^α
- Is α closer to 1.0 than T12's 1.2? (The c coupling may improve it)
- Also measure the ACCRETION profile at the core (r<5)

### T2: Two Braids — Does c-Depletion Create Attraction?

Two bimodal braids at separation D=20, both in M7+c(ρ_B) model.
Each braid creates a c-depletion well. The OTHER braid propagates
slower in the depleted zone → geodesic curves toward the first braid.

- N=192, L=40, T=500, α_c=0.2
- Two braids at (±10, 0, 0)
- Track separation D(t) — if D decreases: ATTRACTION from c-depletion
- Compare to M7 WITHOUT c coupling (α_c=0, control)
- Also test: same-twist vs opposite-twist (universal or charge-dependent?)

This is THE gravity test. If two braids attract through c-depletion,
the model produces gravity.

### T3: Pinned Boundaries — Steady-State Accretion/Depletion

Pin the domain edges at ρ_B = ρ₀ (fixed background density). The edges
act as an infinite reservoir. The braid at the center accretes B from
the reservoir, creating a steady-state depletion profile.

- N=128, L=40, T=1000 (long run for steady state)
- Pinned BC: B_phi = A_bg at edges (Dirichlet for B, absorbing for S)
- α_c=0.2
- Measure ρ_B(r) at intervals — does it stabilize?
- At steady state: measure the accretion RATE (energy flowing inward)
  This rate IS the gravitational mass.

## NOT in V31 (future work)

The following are tantalizing but beyond scope:

- **Gravitational lensing**: Send test waves through the c-depletion well.
  Measure bending angle vs impact parameter. Compare to GR prediction.
  (Needs T1 results first to know the c_eff profile.)

- **Dark matter halos**: The accretion+depletion profile is a halo analog.
  Does it match NFW or Burkert profiles? Does it explain flat rotation curves?
  (Needs T3 steady-state profile first.)

- **Gravitational waves**: Accelerate a braid, measure c_eff ripples at
  far field. Is the radiation quadrupolar (spin-2)?
  (Needs T2 dynamics first.)

- **Multi-braid cosmology**: Many braids in a large domain with c(ρ_B).
  Does the system form clusters, filaments, voids? Large-scale structure?
  (Needs all three T1-T3 first.)

These are all natural extensions but each requires confirmed T1-T3 results.

## Parameters

Base model (from V30):
    m = 1.50, μ = -41.3, κ = 50
    g = 0.01 (S-B coupling)
    α_c = 0.2 (moderate c variation, 10% speed range)
    A_bg = 0.1 (background amplitude)

## Implementation

All three tests use the M7+c(ρ_B) solver from V30 (v30_m7cfield.c).
Extend it with:
- Radial profile measurement (T1)
- Two-braid initialization (T2)
- Pinned BC for B fields (T3)

Shared core: copy v30_m7cfield.c as the base, add per-test modifications.

## Success Criteria

- **T1**: If δc_eff² ~ 1/r^α with α < 1.5 → improving toward gravity
- **T2**: If D(t) decreases → FIRST GRAVITATIONAL ATTRACTION from field depletion
- **T3**: If ρ_B(r) stabilizes → steady-state mass = accretion rate

If T2 shows attraction, this is a major milestone: gravity from field
dynamics, not from an imposed metric or auxiliary graviton field.
