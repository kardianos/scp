# V32 Extended Experiments — Overnight Campaign

## The Working Model

    acc_a = ∇²φ_a + w(P)×α×(∇ρ/ρ)·∇φ_a - m²φ_a - V'(φ_a)
    w(P) = 1/(1 + |P|/P_thresh)
    α = 0.5, m²=2.25, μ=-41.3, κ=50, A_bg=0.1
    Single field, periodic BC, N=1024 where possible

## Experiment List (priority order)

### E1: Single Braid Steady State (N=1024, L=20, T=2000)
Goal: Does the braid reach a true steady state with the gradient coupling?
Measure: E(t), fc(t), D_profile(r,t), ρ(r,t), w(r,t), flux(r,t)
Look for: energy stabilization, radial profile convergence, the dark-matter
bump (accretion at core, depletion at intermediate r)
Pinned border at edges (B fields held at A_bg)

### E2: Two Braids Long Run (N=1024, L=40, T=2000)
Goal: Full attraction dynamics with high resolution
Measure: D(t) trajectory, force = d²D/dt², radial profiles around each braid
Look for: force law (1/r²?), attraction rate, whether braids merge or orbit
Same-twist and opposite-twist variants

### E3: Five Braids (N=1024, L=60, T=2000)
Goal: Multi-body dynamics, clustering, binding
Init: 5 braids at random positions, D>15 between each pair
Measure: all pairwise separations D_ij(t), total energy, cluster formation
Look for: do they form bound states? pairs? a cluster? hierarchy?

### E4: Parameter Sweep (N=128, L=30, T=300 each)
Goal: Map the α sweet spot and P_thresh sensitivity
Sweep: α ∈ {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}
Sweep: P_thresh ∈ {5%, 10%, 20%, 50%} of peak |P|
Quick runs to find optimal parameters before long runs

### E5: Separation Scan (N=128, L=varies, T=300)
Goal: Force law F(D)
Two braids at D ∈ {10, 15, 20, 25, 30, 40}
Measure ΔD(T=300) for each → F ∝ ΔD/T²
Fit: F ~ 1/D^n → what is n?

### E6: Composite Braids (N=256, L=30, T=500)
Goal: Can two braids MERGE into a bound state?
Init: two braids very close (D=8, inside each other's interaction surface)
Also: two braids with different phase (180° offset = "neutron" analog)
Look for: stable bound state, binding energy, new topology

### E7: Mass Spectrum (N=128, L=20, T=500)
Goal: Do different braid topologies have different masses?
Test: winding W=1, W=2, W=3 braids (vary k_factor)
Also: elliptical vs circular, different phase offsets
Measure: gravitational mass = ∫w(P)×|∇φ|×|∇ρ| at the surface

## Characterization Targets

After experiments, analyze for:

1. **Dark matter bump**: Is there ρ_B > ρ_0 at intermediate r (accretion shell)?
   How does it compare to NFW/Burkert profiles?

2. **Attraction modes**: Is the force universal (all braids attract equally)?
   Or charge-dependent (same-twist vs opposite-twist)?

3. **Force law**: F ∝ 1/D^n — what is n? Does it approach 2 (Newton)?

4. **Energy conservation**: How much does E drift? Can we reduce it
   by tuning α or using a symplectic correction?

5. **Proton/neutron analogs**: Can two braids with different internal
   phase form a bound state (deuteron analog)?

6. **Binding energy**: If braids bind, E_bound < E_1 + E_2? By how much?

7. **Steady-state profiles**: What does the equilibrium ρ(r), w(r), flux(r)
   look like? How do they compare to the initial anatomy analysis?

8. **Correlated effects**: Do attraction strength, dark matter bump size,
   and braid mass correlate? Is there a single parameter that controls all?

## Runtime Estimates

| Exp | N | L | T | Est. time | Priority |
|-----|---|---|---|-----------|----------|
| E4 | 128 | 30 | 300 | 7×4 = 28 × 5 min = ~2.5 hr | 1st |
| E5 | 128 | var | 300 | 6 × 5 min = 30 min | 2nd |
| E1 | 1024 | 20 | 2000 | ~8 hr (N³ scaling) | 3rd (overnight) |
| E2 | 1024 | 40 | 2000 | ~24 hr (large grid) | 4th (overnight) |
| E6 | 256 | 30 | 500 | ~1 hr | 5th |
| E7 | 128 | 20 | 500 | 3 × 8 min = 24 min | 6th |
| E3 | 1024 | 60 | 2000 | ~48 hr (huge) | 7th (if time) |

Total for overnight: E4 + E5 + E1 + E6 + E7 ≈ 12 hr
E2 and E3 may need to run longer or at lower resolution.
