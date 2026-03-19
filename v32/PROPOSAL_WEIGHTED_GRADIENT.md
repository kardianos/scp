# Binding-Weighted Gradient Coupling: The Braid's Interaction Surface

## The Core Insight

The braid is NOT a uniform blob. It has internal geometry:

- **Tightly bound core** (r < 3): All three fields aligned, |P| = |φ₀φ₁φ₂|
  large. The triple-product coupling holds these regions together. They are
  self-contained — they don't interact with the surrounding fabric.

- **Weakly bound surface** (r ≈ 4-6): Fields become misaligned, |P| → 0.
  The binding loosens. This is WHERE the braid exchanges energy with the
  fabric. The intake/outtake happens here.

- **Unbound fabric** (r > 6): |P| ≈ 0. Pure ambient field. No self-structure.

The critical realization: **the braid interacts with the fabric ONLY at its
weakly-bound surface.** The tightly-bound interior is invisible to the
ambient field. Any coupling term (gradient, c(ρ), etc.) that operates
uniformly across the whole field is WRONG — it's dominated by the core's
own massive gradients, causing self-interaction blowup.

## The Mechanical Analogy

Think of the braid like a structural object under stress-strain analysis:

- The **bonds** between the three fields (triple product P) are like
  molecular bonds in a crystal
- The **core** is the bulk crystal — strong bonds, rigid structure
- The **surface** is the crystal face — broken/dangling bonds, reactive
- The crystal interacts with its environment ONLY at the surface
- Internal stress doesn't couple to external forces — only surface stress does

In our field theory:
- The triple product P measures bond strength (field alignment)
- Where P is large: fields are locked together (bonded)
- Where P is small: fields are free to interact with the environment
- The gradient coupling should operate proportional to the WEAKNESS of
  the binding: w(P) = 1/(1 + |P|/P_threshold)

## The Braid's Anatomy (from V32 analysis)

Time-averaged radial profile of the equilibrated bimodal braid:

| r | |P| | ρ | w(P) | flux_r | Role |
|---|-----|---|------|--------|------|
| 0-3 | 0.11-0.15 | 1.2-1.6 | 0.22-0.28 | INTAKE (-0.01) | Tightly bound core |
| 3-4 | 0.06-0.11 | 0.8-1.2 | 0.28-0.43 | Mixed | Transition zone |
| **4-6** | **0.01-0.06** | **0.3-0.8** | **0.43-0.90** | **OUTTAKE (+0.01)** | **Interaction surface** |
| 6-12 | <0.005 | 0.03-0.17 | 0.96-1.0 | Weak outtake | Nearly pure fabric |
| 12+ | ~0 | 0.025 | ~1.0 | Weak intake | Background fabric |

### Key findings from the anatomy:

1. **The interaction surface is at r ≈ 4-6** where the binding weakness
   w(P) crosses 0.5.

2. **The energy flux REVERSES at the surface**: the core INTAKES (pulls
   field inward), the surface OUTTAKES (releases field outward). The braid
   is a PUMP — it processes field from the fabric through its helical
   structure and releases it at the surface.

3. **The density gradient |∇ρ| peaks at the surface** (r ≈ 4-5), not at
   the core center. This is exactly where the gradient coupling matters
   most — and where the binding weight allows it to operate.

4. **Far-field return flow** (r > 12): weak intake. The fabric slowly
   flows back toward the braid. This return flow IS the gravitational
   attraction in the depletion picture.

## The Equation

Standard wave equation (uniform medium, what we've been using):

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a

Correct wave equation on non-uniform medium (divergence form):

    ∂²φ_a/∂t² = (1/ρ)∇·(ρ∇φ_a) - m²φ_a - ∂V/∂φ_a
               = ∇²φ_a + (∇ρ/ρ)·∇φ_a - m²φ_a - ∂V/∂φ_a

The gradient coupling term (∇ρ/ρ)·∇φ_a makes the field drift along the
density gradient. But applied uniformly, it's dominated by the braid's own
gradient (self-interaction → blowup, as tested in V32 gradient test).

**The binding-weighted version:**

    ∂²φ_a/∂t² = ∇²φ_a + w(P)×α×(∇ρ/ρ)·∇φ_a - m²φ_a - ∂V/∂φ_a

    w(P) = 1 / (1 + |P| / P_threshold)

    P = φ₀φ₁φ₂ (triple product = binding strength)
    P_threshold = 10% of initial peak |P|

This is a single equation on a single field. No split, no smoothing, no c(ρ).
The binding weight w(P) naturally separates:

- **Core** (|P| >> P_thresh): w ≈ 0 → gradient coupling OFF → no self-interaction
- **Surface** (|P| ~ P_thresh): w ≈ 0.5 → partial coupling → interaction zone
- **Fabric** (|P| ≈ 0): w ≈ 1 → full coupling → ambient dynamics

## Why This Should Produce Gravity

1. **Braid 1** creates a density profile ρ(r) — high at core, low in depletion
   zone, background at far field.

2. **∇ρ points toward Braid 1** in the surrounding fabric.

3. **Braid 2** propagates through this fabric. At Braid 2's weakly-bound
   surface (r ≈ 4-6 from Braid 2's center), the binding weight w ≈ 0.5.

4. The gradient coupling w(P)×(∇ρ/ρ)·∇φ at Braid 2's surface is:
   - ∇ρ from Braid 1's density gradient (pointing toward Braid 1)
   - ∇φ from Braid 2's own field (pointing outward from Braid 2's surface)
   - The dot product biases Braid 2's dynamics toward Braid 1

5. Braid 2 drifts toward Braid 1 = **gravitational attraction**.

6. The strength is proportional to:
   - Braid 1's density gradient at Braid 2's location (∝ M₁/r²)
   - Braid 2's surface field strength (∝ M₂)
   - → Force ∝ M₁M₂/r² = Newtonian gravity

## What's Different from Previous Approaches

| Approach | Problem | This approach |
|----------|---------|---------------|
| M7 split (V29-V31) | Artificial S/B separation | Single field, no split |
| c(ρ) on grid (V30-V31) | Core freezing or blowup | No c(ρ) modification |
| Uniform gradient (V32) | Self-interaction blowup | Weighted by binding → core suppressed |
| SPH (V32) | Signal too weak (0.9%) | Explicit gradient term at surface |

The binding weight is the key innovation. It's not arbitrary — it comes
directly from the braid's own structure (the triple product P that defines
the braid is also the mask that controls where it interacts with the fabric).

## Implementation

Source: `v32/src/v32_weighted_grad.c`

Per timestep:
1. Compute ρ(x) = energy density (kinetic + gradient + mass + potential)
2. Compute |P(x)| = |φ₀φ₁φ₂| (binding strength)
3. Compute w(x) = 1/(1 + |P|/P_thresh) (binding weakness)
4. Compute ∇ρ(x) via central differences (periodic)
5. Compute ∇φ_a(x) via central differences
6. Force: acc_a = ∇²φ_a + w×α×(∇ρ/ρ)·∇φ_a - m²φ_a - V'(φ_a)
7. Velocity Verlet (symplectic, energy-conserving)

Cost: ~1.5× the standard equation (extra ∇ρ and ∇φ passes).
Grid: N=128, L=30, periodic BC, dt=0.12×dx.

## Test Configuration

Two braids at D=20, scanning α = {0, 0.5, 1.0, 2.0}:
- α=0: control (standard equation, no gradient coupling)
- α=0.5: half-strength
- α=1.0: physical (correct divergence form)
- α=2.0: double-strength (stress test)

Success criteria:
- Braids survive (no blowup — the binding weight should prevent it)
- D(t) decreases for α>0 but not α=0 (gradient coupling creates attraction)
- Energy approximately conserved (symplectic Verlet)

## Connection to Mass and Gravity

If the gradient coupling produces attraction:

- **Mass** = the braid's surface integral of w(P)×|∇φ| (how much the
  braid's surface interacts with the fabric gradient)
- **Gravitational force** ∝ (surface area × interaction strength) × (∇ρ from
  other braid) ∝ M₁M₂/r²
- **Mass is emergent**: it's a property of the braid's GEOMETRY (how much
  weakly-bound surface it has), not a Lagrangian parameter

Different braid topologies → different surface areas → different masses.
This is the mass spectrum from topology that we've been looking for.
