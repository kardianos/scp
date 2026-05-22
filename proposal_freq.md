# Proposal: Emergent Frequency Confinement

**Status**: Proposal
**Date**: April 2026
**Motivation**: Oscillons disperse because they never reach amplitude levels
where the nonlinear potential creates self-consistent frequency confinement.

---

## 1. The Ball-and-Wire Metaphor

Picture a track that oscillates back and forth at a fixed frequency. On this
track sit three pendulums — rigid wires of different lengths, each with a ball
at the end. The wires have different natural frequencies because they have
different lengths.

When the track oscillates, only the pendulum whose natural frequency matches
the track frequency will swing with large amplitude. The other two balls barely
move. Resonance selects the active channel.

Now extend this: each pendulum's "length" isn't fixed. It depends on how much
energy has concentrated at that ball. A heavy, concentrated ball pulls its wire
taut — changing its effective length and therefore its natural frequency. The
system is self-referential: the amplitude determines the frequency, and the
frequency determines which channel absorbs energy.

In this picture:

- **Stability** comes from resonance lock — once a ball is swinging at its
  natural frequency, it efficiently absorbs energy from the track, maintaining
  its motion.
- **Separation** comes from frequency mismatch — balls at different frequencies
  don't exchange energy. They're dynamically isolated.
- **Cohesion** comes from matched frequencies — when two systems share the same
  frequency, they couple constructively and bind.

---

## 2. Displacement-Based Frequency

The key insight: the frequency is not a parameter. It is an emergent property
of the local field state.

In the current Cosserat equations, the effective spring constant for each field
component is:

    k_eff = m² + nonlinear correction from V(P)

where `P = φ₀ φ₁ φ₂` and `V(P) = (μ/2) P² / (1 + κP²)`.

The potential is saturating and μ is negative, so inside a concentrated region:

- **High P** (particle interior): V'(P) is large and negative → k_eff drops
  below m² → effective frequency **decreases** → the "pendulum" is longer.
- **Low P** (background): V'(P) ≈ 0 → k_eff ≈ m² → frequency stays at ω₀.

The particle interior oscillates at a **different frequency** than the
surrounding medium. Waves generated inside the particle encounter a frequency
mismatch at the boundary. They reflect rather than propagating freely into the
background.

This is self-consistent confinement: the particle digs its own frequency well
by depleting the field around it. The depletion itself creates the barrier that
prevents further radiation. The mass of the particle IS the field depletion —
the two are the same thing.

The mechanism provides:

1. **Resistance to mixing**: Particles at different concentration levels have
   different internal frequencies. They can't merge because their oscillation
   patterns are incommensurate. Only structures at the same concentration level
   (same effective frequency) can coherently overlap and bind.

2. **Particle identity**: The frequency fingerprint of a structure encodes its
   mass and internal state. Two protons have the same frequency profile and can
   interact. A proton and a random field fluctuation have mismatched frequencies
   and pass through each other.

3. **Natural quantization**: Only certain concentration levels produce
   self-consistent frequency wells deep enough for confinement. This could
   create discrete mass levels — a spectrum of allowed particle states.

---

## 3. Implementation Analysis

### 3.1 The Mechanism Already Exists in the Equations

The saturating potential `V(P) = (μ/2)P²/(1+κP²)` already produces
amplitude-dependent effective frequency. No equation changes are needed.

The second derivative of V with respect to φ_a gives the nonlinear correction
to the effective spring constant. At a local amplitude A, the effective
frequency squared is approximately:

    ω_eff² ≈ m² + μ · g(A, κ)

where g(A, κ) encodes how the saturation changes with local triple-product
amplitude. Since μ < 0, the correction is negative — higher amplitude means
lower effective frequency.

### 3.2 But It's Too Weak at Current Amplitudes

Estimated frequency shift for current parameters (m²=2.25, μ=-41.345, κ=50):

| Local P | κP²  | |V'(P)·P⁻¹| | Frequency shift Δω/ω₀ |
|---------|------|--------------|------------------------|
| 0.05    | 0.13 | 0.40         | ~9%                    |
| 0.1     | 0.50 | 1.55         | ~17%                   |
| 0.3     | 4.50 | 8.74         | ~46%                   |
| 0.5     | 12.5 | 12.1         | ~66%                   |
| 1.0     | 50.0 | 0.80         | ~19% (saturated)       |

The frequency shift peaks at intermediate P (around 0.3–0.5) where the
potential is steepest, then decreases at very high P where the saturation
flattens V.

**Current oscillons operate at P ~ 0.05–0.1**, giving only a 9–17% frequency
shift. This is not enough to create a strong confinement barrier. The oscillon
breathes, disperses to low amplitude where the barrier vanishes, and radiation
escapes.

**For genuine confinement, the structure needs to sustain P ~ 0.3–0.5**, where
the frequency shift reaches 40–66%. At that level, the interior and exterior
frequencies differ by nearly a factor of 2 — waves at one frequency are almost
totally reflected at the boundary.

### 3.3 Two Paths to Stronger Self-Confinement

**Path A: Parameter tuning** — Adjust μ, κ, m² to push the confinement
threshold down into the regime the oscillons naturally reach.

The frequency shift scales roughly as |μ|/(m² · κ). To make the barrier
stronger at low P:

- **Increase |μ|**: Stronger nonlinearity. Currently |μ|/m² ≈ 18.4, well
  above the |μ|/m² < 7 dispersal limit from V53. But that limit was derived
  for the *absence* of frequency confinement — if confinement works, higher |μ|
  might be tolerable because the particle protects itself.
- **Decrease κ**: Less saturation → the potential stays steep to higher P.
  Currently κ=50. Values of κ=10–20 would keep the potential steep through
  P ~ 0.3–0.5.
- **Decrease m²**: Lower bare mass → the nonlinear correction is a larger
  fraction of the total. m² = 1.0–1.5 would make the same |μ| produce a
  proportionally larger frequency shift.

Candidate parameter sets to test:

| Label | m²    | μ       | κ   | Rationale                              |
|-------|-------|---------|-----|----------------------------------------|
| P0    | 2.25  | -41.345 | 50  | Baseline (current)                     |
| P1    | 2.25  | -41.345 | 15  | Lower saturation, steeper potential     |
| P2    | 1.50  | -30.0   | 15  | Lower mass + lower saturation           |
| P3    | 1.50  | -20.0   | 10  | Aggressive: maximize |μ|/(m²κ) ratio   |
| P4    | 2.25  | -60.0   | 50  | Higher |μ| at current κ (stress test)   |

**Path B: Seed design** — Initialize structures that start in the
high-concentration regime (P ~ 0.3–0.5) and test whether the emergent
frequency barrier sustains them.

The analytical proton seeds (gen_proton_analytical) produce P_int ~ 1270 per
baryon, which is far above the saturation regime. But the pre-converged
templates (proton_template.sfa at 64³) operate in a more moderate range. A
seed generator that targets P ~ 0.3–0.5 at initialization would test whether
the self-confinement mechanism can hold from the start.

### 3.4 Diagnostics for Frequency Confinement

Standard diagnostics (mass, breathing depth, fragmentation) will show whether
a structure survives. To specifically test the frequency confinement hypothesis,
add:

1. **Local frequency measurement**: At each diagnostic snapshot, compute the
   instantaneous oscillation frequency ω_local(r) from the time derivative
   ∂φ/∂t and the field φ. At a local extremum of φ, ω² ≈ |∂²φ/∂t²|/|φ|.
   Map this as a radial profile from the particle center.

2. **Frequency barrier ratio**: ω_interior / ω_background. If this ratio
   exceeds ~1.5, the barrier should be strong enough for significant wave
   reflection.

3. **Radiation leakage rate**: Measure energy flux through a shell at radius
   R from the particle center. Compare to the total particle energy. A strong
   frequency barrier should make this rate exponentially small rather than
   algebraic.

4. **P profile**: Track the radial profile of P(r) over time. If
   self-confinement works, P at the center should stabilize rather than
   oscillating through zero.

---

## 4. Experimental Plan

### Phase 1: Parameter Scan (GPU, ~2 hours total)

Run the standard oscillon seed at N=128 under absorbing BC for T=500 with
each parameter set P0–P4. Measure:

- Oscillon lifetime (time until mass drops below 50% of initial)
- Breathing depth (max/min mass ratio)
- Peak P value reached during evolution
- ω_interior/ω_background ratio at peak concentration

**Success criterion**: Any parameter set that produces lifetime > 2× baseline
(P0) with breathing depth < 5:1 (vs current ~20:1).

### Phase 2: Frequency Profile Measurement (GPU, ~1 hour)

For the best parameter set from Phase 1, run at N=192, T=1000 with high
diagnostic cadence (diag_dt=0.5). Extract:

- Radial ω(r) profile at multiple times
- P(r) profile evolution
- Energy flux through shells at r = 2, 5, 10, 20 code units

**Success criterion**: ω_interior/ω_background > 1.3 sustained for T > 200.

### Phase 3: Composite Test (GPU, ~2 hours)

If Phase 1-2 show self-confinement, test composite structures (UUD, UDD) with
the best parameters. Does frequency confinement help the composite hold
together? Do two composites at the same frequency bind more strongly than
two at different frequencies?

**Success criterion**: UUD composite survives T=1000 under absorbing BC with
<50% mass loss.

---

## 5. Risks and Fallbacks

**Risk 1**: Lower κ or higher |μ| destabilizes the background oscillation,
preventing any structure from forming.

*Mitigation*: Start from P1 (only κ change) which is the most conservative.
If background destabilizes, increase A_bg to provide stronger anchoring.

**Risk 2**: The frequency barrier exists but oscillon breathing still
disperses the structure below the confinement threshold before the barrier
can form.

*Mitigation*: Use Path B seeds that start above the threshold. If the barrier
holds an initialized high-P structure, the mechanism works even if
self-assembly from low-P seeds fails.

**Risk 3**: Higher |μ|/m² ratio violates the V53 dispersal condition and
structures radiate faster, not slower.

*Mitigation*: The V53 condition assumed no frequency confinement. If
confinement is real, it overrides the dispersal condition. But we should
verify empirically — if P3 or P4 disperse faster than P0, the mechanism is
too weak at those parameters.

**Fallback**: If parameter tuning alone doesn't produce self-confinement,
consider an explicit amplitude-dependent mass term: m_eff² = m² + γ|φ|².
This is a standard nonlinear Klein-Gordon extension that directly implements
"displacement-based frequency." It requires a kernel change but is a
well-studied soliton mechanism.
