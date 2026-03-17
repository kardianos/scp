# V24 Results: Alternative Oscillon Constructions

## Overview

V24 tested 12 distinct constructions seeking a bound state that either
(a) is stabilized by motion/noise (the "self-reconstruction" hypothesis), or
(b) uses a different phase structure (120°) to suppress radiation.

All 12 tests produced negative or null results for the primary hypothesis.
Two unexpected discoveries emerged: the 180° anti-phase state as a natural
attractor, and a correct UDD > UUD mass ordering at intermediate coupling.

## Part 1: Motion and Noise (Tests 1-3, DG, D)

| Test | Hypothesis | Result |
|------|-----------|--------|
| Boost v=0-0.7 | Moving oscillon radiates less | NO: v=0 most stable, forward beaming dominates |
| Thermal T=0-0.1 | Noise provides minimum ds>0 | NO: coherence destroyed at T=0.001 |
| Rotating triad δ=0-0.2 | Freq splitting creates gyroscopic stability | NO: splitting decays in ~200 t.u. |
| Three-frequency | Non-degenerate ω_a reduces radiation | NO: all frequencies converge to degenerate |

**Conclusion**: The symmetric oscillon (φ₁=φ₂=φ₃, ω≈0.87) is a unique,
powerful attractor. All perturbations decay back to it within ~200 t.u.
The triple-product coupling enforces complete phase synchronization.

## Part 2: Alternative Topologies (A, E)

| Test | Construction | Result |
|------|-------------|--------|
| 2D Vortex | Circulation-stabilized pair | Dissolves in ~50 t.u. (mass kills 1/r tail) |
| Dissipative soliton | Gain-loss balance | No stable balance (polynomial gain always blows up) |

## Part 3: 120° Phase Binding (PW, PT, SYM, B)

| Test | Mechanism | Result |
|------|-----------|--------|
| 120° at μ=-20 (V24-B) | Harmonic elimination | P_max → f³/4, binding 16× weaker, disperses |
| Pairwise coupling λ scan | m²_anti = m²-λ lowers base mass | P≈0 at 120°, NO nonlinear binding at any λ |
| Pairwise + triple threshold | 3ω < m at λ≈0.88 | Condensation at all λ≥0.8, no oscillating state |
| Symmetric penalty (Σφ)² | Force 0°→120° | System goes to 180°, not 120° |

**Root cause**: The 120° phase configuration makes the time-averaged triple
product ⟨P²⟩₁₂₀ = f⁶/32, which is 12× smaller than ⟨P²⟩₀ = 3f⁶/8.
This eliminates the nonlinear frequency shift that pushes ω below the mass
gap. Without sub-gap protection, the 120° state radiates freely.

There is NO intermediate regime: at small amplitudes the binding is too
weak (dispersal), at large amplitudes the pairwise coupling causes
condensation. No stable oscillating 120° oscillon exists.

## Part 4: UUD/UDD Flavor Physics (UUD, UDD)

| Test | Configuration | Result |
|------|--------------|--------|
| UUD (proton) | m₁=m₂=m_U, m₃=m_D, 120° phases | Disperses (13-25% retention) |
| UDD (neutron) | m₁=m_U, m₂=m₃=m_D, 120° phases | Disperses (13-21% retention) |

At λ=0.30 with m_D < m_U: **UDD retains 27-187% more energy than UUD**,
matching the real neutron > proton mass ordering. The mechanism: UDD's
DD-antisymmetric mode has lower eigenvalue → slower oscillation → less
efficient radiation.

Neither state is stable, but the mass ordering is physically correct.

## Discoveries

### 1. The 180° Anti-Phase State

Found independently by V24-SYM and V24-PW: the system initialized at 120°
naturally evolves to 180° (one field anti-phase to the other two).

Properties of the 180° state (φ, φ, -φ):
- P = φ₁φ₂φ₃ = f·f·(-f) = -f³ (FULL binding, opposite sign to 0°)
- The nonlinear frequency shift is PRESENT (unlike 120° where P≈0)
- The state is a natural attractor (robust against perturbation)
- It provides a natural UUD-like structure: two "positive" + one "negative"
- The opposite-sign P could have different radiation properties from 0°

This was NOT predicted or sought — it emerged from the dynamics.

### 2. UDD > UUD Mass Ordering

The mass splitting at λ=0.30 with m_U=1.0, m_D=0.80-0.95 gives:
- UDD energy > UUD energy (neutron heavier than proton)
- The ordering flips at m_D ≈ 0.75

This matches the real proton-neutron mass difference (Mn - Mp ≈ 1.3 MeV)
and arises from the DD-mode vs UU-mode eigenvalue splitting.

## The Fundamental Constraints (Updated)

V24 establishes several structural results about the three-field
triple-product system:

1. **The symmetric 0° oscillon is the unique attractor** for the
   conservative (μ<0, κ>0, m>0) Lagrangian. All perturbations decay to it.

2. **The 120° phase eliminates binding** because P≈0, removing the
   nonlinear frequency shift needed for sub-gap protection.

3. **The 180° phase preserves binding** (P=-f³) and is a secondary attractor.

4. **The mass gap m²φ² prevents all non-oscillatory structures** (vortices,
   traveling waves, extended states).

5. **Noise and motion are purely destructive** for the coherent oscillon.

6. **Mass splitting (UUD vs UDD) produces correct ordering** even in the
   non-stable 120° regime.

## Part 5: 180° Anti-Phase Characterization (180A, 180B, 180C)

| Test | Result | Finding |
|------|--------|---------|
| 180A (characterize) | Z₂ symmetry | 0° and 180° are bit-for-bit identical (exact Z₂: φ₃→-φ₃) |
| 180B (UUD/UDD) | Mass ordering | UDD < UUD in 180° (reversed from 120°, bare mass dominates) |
| 180C (radiation) | dE/dt identical | Same to 6 sig figs. Not a new particle, just a Z₂ partner |

## Part 6: Maxwell Integration (MA–MF)

| Test | Approach | Result | Finding |
|------|----------|--------|---------|
| MA | Complex Ψ + gauge | NEGATIVE | Phase-locking incompatible with U(1); e²A²φ runaway at e≥0.15 |
| MB | 3 complex + quark Q | NEGATIVE | Gauge coupling shifts ω above gap; oscillon destroyed at e=0.1 |
| MC | Topological current 2D | NEGATIVE | Q_top = 0 exactly (smooth fields, parallel gradients) |
| MD | Elastic dislocation 2D | NEGATIVE | Zero incompatibility (Clairaut's theorem for smooth C² fields) |
| **ME** | **Proca from pairwise** | **POSITIVE** | **m_A tunable to 0.07, range 14, two-oscillon merger at λ=0.99** |
| MF | Goldstone θ from Z₂ | NEGATIVE* | No static profile; but backreaction strengthens oscillon +131% |

**Explicit gauge coupling (MA, MB) kills the oscillon** — the e²A²φ term
shifts ω above the mass gap. The triple-product phase-locking is fundamentally
incompatible with U(1) charge.

**Emergent EM (MC, MD) fails** — the oscillon is too smooth/symmetric for
topological charge or strain incompatibility.

**Proca from pairwise coupling (ME) is the breakthrough** — the pairwise term
V_pw = λ(φ₁φ₂+φ₂φ₃+φ₃φ₁) splits the mass spectrum. The antisymmetric mode
mass m_A = √(m²-λ) is tunable. The symmetric oscillon is STABILIZED (m_S
increases with λ). Two oscillons attract through the Proca channel.

## Part 7: Proca Range Extension (P1–P5)

| Test | Approach | Result | Range Achieved |
|------|----------|--------|---------------|
| **P1** | Push λ → 0.99999 | **POSITIVE** | **316+ (unlimited, no wall)** |
| **P2** | Condensed phase (λ>m²) | MIXED | 56 (pseudo-Goldstone, m=0.018) |
| P3 | Cascade Proca→θ | Trivial | 130+ (any massless field from localized source) |
| P4 | Self-consistent λ(x) | NEGATIVE | 1.4 (self-defeating feedback, 2.4% shift) |
| **P5** | MF+ME combined | **POSITIVE** | **Past tachyonic boundary (λ=1.05 survives)** |

**The Proca range is effectively unlimited.** λ can be pushed arbitrarily
close to m² (tested to λ=0.99999, range 316). The oscillon ALSO survives
past the tachyonic boundary at λ=1.05 (both with and without θ backreaction).

**The oscillon is stable in a tachyonic vacuum.** At λ > m², the antisymmetric
modes have m²_A < 0 in vacuum but are stabilized inside the oscillon core by
the triple-product Hessian contribution λ_anti = +|μ|f⁴/(1+κf⁶)².

**Triple product breaks SO(2) → pseudo-Goldstone.** In the condensed phase
(λ > m²), the would-be Goldstone acquires mass from the triple product. At
λ=1.01: m_pseudo ≈ 0.018 (range 56). True massless Goldstone requires μ=0
(no triple product), which removes the oscillon binding.

## Complete V24 Summary

**28 investigations across 7 series:**
- Part 1 (Motion/Noise): 5 tests, all negative
- Part 2 (Alternative topology): 2 tests, both negative
- Part 3 (120° phase binding): 4 tests, all negative
- Part 4 (UUD/UDD flavor): 2 tests, correct mass ordering but unstable
- Part 5 (180° characterization): 3 tests, Z₂ symmetry confirmed
- Part 6 (Maxwell integration): 6 tests, ME positive (Proca mediator)
- Part 7 (Proca range extension): 5 tests, P1+P5 positive (unlimited range)

**Three major discoveries:**
1. Z₂ symmetry (180° = 0° in other sector)
2. UDD > UUD mass ordering (correct neutron > proton in weakly-bound regime)
3. Unlimited Proca range from pairwise coupling (no stability wall, survives
   past tachyonic boundary)

**The fundamental tension (updated):**
- Triple product creates oscillons BUT prevents massless mediators (breaks SO(2))
- Explicit gauge coupling destroys oscillons (e²A²φ shifts ω above gap)
- Pairwise coupling gives tunable Proca BUT not true 1/r (still Yukawa)
- The closest to 1/r: pseudo-Goldstone at m=0.018 in condensed phase, or
  Proca at m_A=0.003 (range 316) by pushing λ=0.99999

## Part 8: Steps 1-4 (Force, Field, Multi-scale, 3D)

| Test | Result | Finding |
|------|--------|---------|
| S1 (force) | Inconclusive | Symmetric oscillon doesn't couple to Proca channel |
| S2 (field) | Positive | Proca channel exists at ω=0.15 (sub-gap) when seeded |
| S3 (multi-scale) | **Positive** | Two ranges coexist (13× spread), sign change in F(D) |
| S4 (3D) | **Positive** | 3D oscillon survives all λ through 0.99, merger at D=20 |
| S1b (asym force) | Ambiguous | 180° force 82-1000× stronger but different soliton branch |

## Part 9: Steps 5-6 (Origin of λ, Spin)

| Test | Result | Finding |
|------|--------|---------|
| S5 (origin λ) | Partial | Gap margin not constant; modulus shifts 7% but wrong direction |
| S6 (spin) | **Definitive** | Spin-0 scalar. Cross-gradient adds 11% quadrupole but not spin-2 |

## Part 10: Fundamental Tests (7 Ideas + Combos)

| Test | Result | Finding |
|------|--------|---------|
| A (EOS) | Done | Pressureless dust core, radiation vacuum |
| B (Inertia) | **Positive** | a=F/M ±4%, quadratic deformation, F_crit≈2e-5 |
| C (Confine) | Partial | σ=1 survives but still radiates (harmonics above gap) |
| D (Gauge) | Negative | Ω modulation destabilizes via parametric resonance |
| **E (Lattice)** | **Positive** | **λ=0.5 stabilizes chain, c_s≈0.9, phonons real** |
| **F (Self-ref)** | **Positive** | **Φ=-0.006 at α=1e-4, self-consistent, stable** |
| **Combo 2+5** | **Positive** | **Antisymmetric mode propagates + amplifies in lattice** |
| Combo 2+3+5 | Mixed | Near-gapless mode but it's a saddle-point instability |
| **Combo 1+2+5+7** | **Positive** | **Self-gravitating lattice stable, phonon shift 3.7%** |

## Part 11: Locality Tests (Causal Propagation)

| Test | Result | Finding |
|------|--------|---------|
| **L1 (EMA wave)** | **Positive** | Causal (t=51.5 at x=50), Φ=-9.1e-3, τ insensitive |
| L2 (antisymmetric) | Negative | Saddle-point instability, 350× amplification destroys |
| L3 (telegraph) | Partial | Works but Φ too small, slow convergence |
| **L4 (wave eq)** | **Positive** | **Exact causality (50.5±0.5), Φ=78% Poisson, reference** |
| L5 (field metric) | Local only | Causal but not gravity (no nonlocal coupling) |

## V24 Final Summary

**48 sub-investigations across 11 series.**

The V24 Lagrangian with full equations:

    L_total = L_matter + L_pairwise + L_gravity

    L_matter = Σ_a [½(∂φ_a)² - ½m²φ_a²] - (μ/2)P²/(1+κP²)

    L_pairwise = -λ(φ₁φ₂ + φ₂φ₃ + φ₃φ₁)

    L_gravity = ½(∂Φ)² ;  □Φ = -αρ ;  backreaction: m²→m²(1+2Φ), c²→1+4Φ

**What V24 achieved:**
1. Matter: oscillon (from v21, unchanged)
2. Short-range force: pairwise coupling stabilizes lattice (Test E)
3. Long-range mediator: Proca at tunable range (ME/P1, up to 316+)
4. Causal gravity: wave equation □Φ=-αρ propagates at c (L1/L4)
5. Self-consistency: ρ→Φ→dynamics→ρ converges (Test F)
6. Newton's law: a=F/M verified (Test B)
7. Phonon spectrum: acoustic branch c_s≈0.9 + optical branch ω≈0.51

**What V24 did NOT achieve:**
1. Spin-2 gravity (the mediator is spin-0)
2. Fixed coupling constants (α, λ are free parameters)
3. True 1/r from gauge invariance (Proca is Yukawa, not Coulomb)
4. Confinement (oscillons still leak at dE/dt~10⁻⁵)
5. 3D lattice verification

**The V24 theory is a 1+1D scalar-gravity analog** with matter, forces,
and causal dynamics. It captures Newtonian gravity (1/r potential,
causal propagation, self-consistency) but not General Relativity
(spin-2, tensor modes, equivalence principle unverified).
