# Derived Parameters for SCP Field Theory — Revision 02

## Context

V45 tested deuterium (UUD+UDD) at 8 separations (D=0,5,10,15,20,40,60,80)
with identical seeds and grid. Result: energy monotonically decreases toward
D=80. No binding well at any separation. Phase confinement repulsion works
(+1253 at D=0), but no attraction emerges at nuclear distances.

This document explores whether the current parameters (m²=2.25, μ=-41.345,
κ=50, η=0.5, A_bg=0.1) can be derived from first principles, and whether
their relationships explain the absence of binding.

## First Principles (carried forward)

The SCP medium is a single 3-component vector field φ with dynamics governed
by a relativistic wave equation plus a nonlinear potential V(P) where
P = φ₀φ₁φ₂. All structures are attractors of this equation. Free parameters
break predictive power. A mature theory should have parameters fixed by
internal consistency.

## Parameter Analysis

### 1. Mass term m² — KEEP current derivation approach

m² = 2.25 (m = 1.5) was optimized for braid survival. V34 Track G showed
m² ≥ 1.25 required for braid binding, m² ≥ 0.25 for vacuum stability.
The current value sits in the sweet spot.

**Revised derivation**: m² is constrained by a window, not a single value:
- Lower bound: braid stability (m² ≥ 1.25)
- Upper bound: gravitational range (larger m² → shorter Yukawa tail → weaker
  gravity). V34 phonon test showed depletion decays as r^{-1.2}, not Yukawa,
  so this constraint is weaker than expected.
- The V28 CMA-ES optimization selected m=1.5 from a continuous search — it's
  the stability optimum, not arbitrary.

**Not valid (from -01)**: Relating m² to k_bg² = (π/L)² conflates a physical
mass with a simulation box artifact. The background dispersion ω² = k² + m²
is always satisfied; it doesn't constrain m².

### 2. Potential parameters μ and κ — REVISE approach

**From -01**: The proposal μ = -(m²P_core)/(dV/dP) is circular — P_core
depends on μ and κ. Can't derive a parameter from the equilibrium it creates.

**Revised**: Use the virial theorem as the constraint. For the Cosserat
equation, the equilibrium condition for a localized structure is:

    E₂ - E₄ + 3E_V = 0  (virial)
    Mc² = 2E₄ - 2E_V    (mass-energy)

where E₂ is kinetic/gradient energy, E₄ is quartic coupling, and E_V is the
potential energy. This was derived and verified in the v2 Skyrme model work
(finite-lambda solver).

The virial relation connects μ, κ, and the braid profile. Given m² and the
measured braid profile f(r), the virial equation determines the μ/κ ratio:

    μ/κ = function of (∫f²r²dr, ∫f⁴r²dr, ∫(f')²r²dr)

This eliminates one free parameter. The remaining one (say κ alone) sets
the overall energy scale. This is a genuine first-principles derivation
because the virial theorem follows from the Lagrangian, not from fitting.

**Key question**: Does the virial-derived μ/κ ratio match the empirical
μ=-41.345, κ=50 (ratio = -0.827)? If yes, the current values are already
self-consistent. If no, the mismatch might explain the missing binding.

### 3. Coupling η — REVISE substantially

**From -01**: Set η so θ/φ = α ≈ 1/137. This is conceptually appealing but
premature. The V34 measurement (θ/φ = 28% at η=0.5) is at the single-braid
level, not the composite baryon level. The V44 OQ1 result shows the proton's
effective charge emerges from the composite radiation pattern (3 orthogonal
dipoles → monopole). The mapping η → α requires:
1. Compute the composite radiation pattern (done — OQ1)
2. Compute the time-averaged radiation pressure between two composites
3. Extract the effective coupling constant from the radiation pressure

This hasn't been done. Setting η from α without this calculation is
hand-waving.

**Revised proposal**: η should be constrained by the V42/V45 force balance.
V42 showed strong:EM → 1:1 at equilibrium (D≈40). V45 showed no binding
at any D. The question is: at what η would binding appear?

Experiment: Run the V45 separation sweep at η = {0.1, 0.3, 0.7, 1.0, 2.0}.
If binding appears at some η_crit, that constrains η from the requirement
"nuclear binding must exist." This is a physical derivation — the parameter
is fixed by the demand that the theory reproduce observed physics.

**Not valid (from -01)**: Making η field-dependent (η = f(|curl φ|)) adds
complexity without clear physical motivation. V39 tested density-dependent κ
and found faster dispersal, not binding. Field-dependent couplings also
break the Lagrangian structure unless derived from a gauge principle.

### 4. Background amplitude A_bg and phases δ — KEEP with caveats

A_bg = 0.1 sets the background energy density ρ_bg ≈ 0.03. This is the
"fabric density" that determines gravitational coupling. Deriving it from
a cosmological initial condition (total energy → expansion → equilibrium
density) is physically motivated but beyond current simulation capability.

**Phases δ = {0, 3.0005, 4.4325}**: These are NOT arbitrary — V28 CMA-ES
optimization found them by maximizing braid stability over a continuous
parameter space. 3.0005 and 4.4325 are close to but not exactly 2π/3 and
4π/3. The deviation from exact thirds may be physically meaningful — it
breaks a degeneracy that would otherwise allow the three field components
to rotate into each other.

**Not valid (from -01)**: Replacing empirical δ with exact 2π/3 multiples
needs validation. The CMA-ES found the empirical values for a reason.

## What V45 Actually Tells Us

The null result has two possible interpretations:

### Interpretation A: Missing physics
The current equation lacks the mechanism for nuclear binding. Real nuclear
binding comes from meson exchange — a massive mediator field that creates
a Yukawa attractive potential at ~1-2 fm. The θ field is massless (m_θ=0),
so it creates 1/r radiation pressure, not a Yukawa well. To get binding:
- Need a massive θ mode (m_θ > 0 at some scale), OR
- Need a new coupling term not present in the current equation, OR
- Need η large enough that the radiation pressure creates a dynamic binding

### Interpretation B: Seed quality
The gen_deuterium seeds produce fragmented baryons (6-8 cores visible in
volview instead of 2 compact objects). The V45 result may be measuring the
interaction between fragmented debris, not between two equilibrated baryons.
The V42 deuterium used the same seeds but showed binding-like behavior at
T=500. This inconsistency needs resolution.

**Test**: Run V45 D=40 with pre-converged template seeds (extract from V43
proton formation output). If properly equilibrated baryons show binding
but gen_deuterium seeds don't, the null result is a seed artifact.

### Interpretation C: The binding IS the confinement
Perhaps nuclear binding in this theory is not a potential well at finite D,
but the phase confinement itself. Three braids form a baryon because of
phase confinement (P→0 at overlap). Two baryons don't "bind" in the
potential energy sense — they are CONFINED together by the same mechanism,
just at a larger scale. The V42 observation of stable deuterium may reflect
confinement at the composite level, requiring longer timescales to manifest
than T=500 in the V45 format.

## Proposed Experiments (v46)

1. **η sweep**: Run D=40 deuterium at η = {0.1, 0.3, 0.7, 1.0, 2.0}.
   If binding appears at some η, that's the constraint.

2. **Template-seeded D=40**: Use pre-converged baryon templates instead of
   gen_deuterium. Compare with V45 D=40 to isolate seed effects.

3. **Virial check**: Compute the virial ratio E₂/E₄ from existing V44/V45
   data. If it deviates from the theoretical value, the parameters aren't
   self-consistent.

4. **Massive θ test**: Run deuterium with m_θ = 0.1, 0.5, 1.0. A small
   θ mass creates a Yukawa attractive potential that could produce binding.
   This tests Interpretation A directly.

5. **Long-run D=40**: Run gen_deuterium D=40 to T=5000 to test whether
   binding emerges on longer timescales (Interpretation C).
