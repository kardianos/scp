# Locality Option 5: Field-Derived Metric (No Separate Φ)

## Thesis

Instead of adding a separate gravitational field Φ and solving a PDE for it,
DEFINE the effective metric directly from the existing scalar fields:

    g_eff(x,t) = 1 + β · ⟨ρ(x,t)⟩ / ρ_0

where ρ is the total energy density (computed from φ_a at each point),
⟨⟩ is an EMA with time constant τ, and ρ_0 is a reference density.

The metric IS the fields. Changes in ρ propagate at the field's speed
(≤ c) → automatic causality. No separate equation for Φ.

## Free Parameters: β (coupling) and τ (EMA timescale)

β determines the strength of the metric correction.
τ determines the filtering of oscillations.

**Backprop tuning**: scan β and τ jointly.

Fix τ = 7 (one breathing period, from Option 1 analysis).
Scan β = {-0.001, -0.003, -0.01, -0.03, -0.1} (negative for attractive).

The effective dynamics:
    m²_eff(x) = m² · g_eff(x) = m²(1 + β⟨ρ⟩/ρ_0)
    c²_eff(x) = g_eff(x)² (or 1 + 4β⟨ρ⟩/ρ_0 in weak-field limit)

## Why This Differs from Test F

Test F used Poisson to compute Φ from ρ — the Φ-ρ relationship was
NON-LOCAL (Poisson kernel ∝ |x|). Here, g_eff is computed LOCALLY from
ρ at the SAME point — no spatial integration.

The gravitational effect is LOCAL: only the energy density at x affects
the metric at x. There's no 1/r potential, no Poisson kernel.

BUT: the fields propagate, carrying energy density information from one
point to another. So the INDIRECT effect IS nonlocal — energy from
oscillon A propagates via field dynamics to oscillon B, where it modifies
g_eff. The nonlocality comes from the field propagation, not from Φ.

The "gravitational interaction range" = the range over which ρ is nonzero
= the oscillon's tail length. With the pairwise coupling (Proca mediator),
this range is 1/m_A — the SAME as the Proca range.

## Method

### Phase 1: Single oscillon, β scan

1. Equilibrate t=5000 (no metric correction)
2. Turn on g_eff with τ=7, each β value
3. Evolve t=5000
4. Measure: ω shift, E change, fc
5. Compare with Test F (Poisson, α=1e-4):
   - Does the local metric produce SIMILAR effects?
   - Is the oscillon modified in the same way?

### Phase 2: Lattice with local metric

6. Build 8-oscillon chain with λ=0.5 and g_eff (best β from Phase 1)
7. Evolve t=10000
8. Measure phonon spectrum
9. Compare with Combo 1+2+5+7 (Poisson gravity):
   - Same phonon speed shift?
   - Better stability?

### Phase 3: Causality test

10. Boost one oscillon in the lattice
11. Track g_eff at distant oscillon positions
12. The g_eff change should propagate at the FIELD speed (not instantaneous)
13. Measure the delay — is it consistent with c?

### Phase 4: Self-consistency check

14. At equilibrium: compute ρ_0 = average energy density in the lattice
15. Compute g_eff everywhere
16. Are the field dynamics on g_eff consistent with the ρ that produces g_eff?
17. This is automatic (by construction), but verify numerically

## Choosing ρ_0

ρ_0 is the reference density that normalizes the metric correction.
Physical choice: ρ_0 = average energy density of the lattice = E_total / L.
With E_total ≈ 9.5 (8 oscillons at λ=0.5) and L = 128: ρ_0 ≈ 0.074.

Then β·⟨ρ⟩/ρ_0 at the oscillon center (⟨ρ⟩ ≈ 0.5):
    g_eff - 1 = β · 0.5/0.074 ≈ 6.8β

For |g_eff - 1| ≈ 0.01 (1% correction, similar to Test F):
    β ≈ -0.0015

Start scan at β = -0.001.

## Reference Code

- v24/fundamental/testF_selfref/src/selfref.c (metric correction)
- v24/fundamental/testE_lattice/src/lattice.c (chain code)

## Output

- `src/locality_fieldmetric.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ=0.5
τ=7.0, ρ_0 computed from lattice average
β scan: {-0.001, -0.003, -0.01, -0.03, -0.1}
Single: Nx=4000, xmax=100
Lattice: Nx=2560, N_osc=8, d=16
t_equil=5000, t_test=5000

Compile: `gcc -O3 -Wall -o locality_fieldmetric src/locality_fieldmetric.c -lm`
