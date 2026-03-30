# V48: Emergent θ Mass from Local Field Structure

## The Requirement

m_θ cannot be a constant. A constant m_θ gives ALL θ modes a mass
everywhere in space — making the "photon" massive in free space, which
contradicts the measured v=0.906c propagation (V43 OQ5) and the
massless 1/r Coulomb behavior at long range.

m_θ must be a function of the local field state — specifically, of
quantities accessible within the local causal neighborhood (the
light cone, or practically, within dx of each grid point).

## What m_θ(local) Should Do

- **In empty space** (far from braids): m_θ → 0. θ waves propagate
  at c. Electromagnetism works normally. Coulomb force is 1/r².

- **Near/between braids** (nuclear distances): m_θ > 0. θ propagation
  becomes Yukawa e^(-m_θr)/r. Creates an attractive potential well
  at distance ~1/m_θ. Nuclear binding.

- **Inside a braid core**: m_θ could be larger still — the θ field
  is "heavy" inside matter. This naturally confines θ energy to the
  braid (consistent with V40 stability signature: θ confinement
  outer/inner < 0.7).

## What Local Quantities Are Available

At each grid point, without looking at neighbors (within ds=0):
- φ_a (field values) → amplitude |φ|², triple product P
- θ_a (angle values) → amplitude |θ|²
- φ̇_a, θ̇_a (velocities) → kinetic energy density

Within one dx (using derivatives already computed for the Laplacian):
- ∇φ, ∇θ (gradients) → gradient energy density
- ∇×φ, ∇×θ (curls) → current density, magnetic field
- |P| = |φ₀φ₁φ₂| → binding density (topological content detector)

## The Natural Choice: m_θ² ∝ |φ|²

The simplest Lorentz-invariant, Lagrangian-compatible, local choice:

    m_θ²(x) = λ_θ × |φ(x)|²  =  λ_θ × (φ₀² + φ₁² + φ₂²)

This comes from a Lagrangian term:

    L_mass_θ = -½ λ_θ |φ|² |θ|²

This is a standard **quartic interaction** between φ and θ — the same
structure as the Higgs mechanism. It's:

- Lorentz invariant ✓ (scalar contraction of field amplitudes)
- Local ✓ (depends only on field values at one point)
- Lagrangian ✓ (derives from -½λ_θ φ²θ², variation gives mass terms)
- No derivatives ✓ (doesn't change canonical momenta)
- No history ✓ (instantaneous, no memory)
- Physically motivated ✓ (see below)

### Physical interpretation

Where |φ|² is large (inside braids, A_core ≈ 0.2-0.3):
    m_θ² = λ_θ × 0.04-0.09 → m_θ ≈ √(λ_θ) × 0.2-0.3

Where |φ|² is small (far from braids, |φ| ≈ A_bg = 0.1):
    m_θ² = λ_θ × 0.01 → m_θ ≈ √(λ_θ) × 0.1

Wait — the background |φ|² ≈ 3 × A_bg² × ½ = 0.015 (three components,
time-averaged cos²). This is nonzero everywhere because the background
oscillates. So m_θ is never exactly zero — even in free space.

This is actually CORRECT physically. The vacuum IS the background field.
The θ field propagates through the background medium, which gives it
an effective mass. The "photon" has a tiny mass from the background,
but it's much smaller than the mass near braids.

### The mass ratio

Background: |φ|²_bg ≈ 0.015 → m_θ_bg² = 0.015 λ_θ
Braid core: |φ|²_core ≈ 0.1-0.3 → m_θ_core² = (0.1-0.3) λ_θ

Ratio: m_θ_core / m_θ_bg ≈ √(0.2/0.015) ≈ 3.6

The θ mass at the braid core is ~4× the background mass. This means:
- Long-range θ (far from braids): nearly massless, slight Yukawa
  suppression at very long range
- Short-range θ (between braids): 4× heavier, Yukawa range 4× shorter
- This creates a NATURAL two-scale force: weak long-range EM and
  stronger short-range nuclear

### Setting λ_θ

For the nuclear range to be ~5-10 code units (where V45 showed the
transition region D≈10-15):

    1/m_θ_core ≈ 5-10 → m_θ_core ≈ 0.1-0.2
    m_θ_core² = λ_θ × |φ|²_core ≈ λ_θ × 0.15
    λ_θ ≈ m_θ_core² / 0.15 ≈ 0.01-0.04 / 0.15 ≈ 0.07-0.27

So λ_θ ≈ 0.1-0.3 gives the right nuclear range.

The background photon mass: m_θ_bg = √(0.015 × 0.15) ≈ 0.047.
Yukawa range in free space: 1/0.047 ≈ 21 code units.
This means EM is effectively 1/r for D < 20 and Yukawa-suppressed
beyond D ≈ 20. At simulation scales (L=100), this is fine — the
boundary effects dominate before the Yukawa cutoff matters.

## Alternative: m_θ² ∝ P²

Instead of |φ|², use P² = (φ₀φ₁φ₂)²:

    m_θ²(x) = λ_θ × P²

Lagrangian: L = -½ λ_θ P² |θ|² = -½ λ_θ (φ₀φ₁φ₂)² (θ₀² + θ₁² + θ₂²)

This activates ONLY where there's binding structure (P ≠ 0), not
where the field merely has amplitude. In the background, P ≈ 0
(the three phase-offset carrier waves cancel in the triple product).
So m_θ_bg ≈ 0 exactly — the photon IS truly massless in free space.

At the braid core (P ≈ 0.08): m_θ² = λ_θ × 0.0064
For m_θ_core ≈ 0.15: λ_θ ≈ 0.0225/0.0064 ≈ 3.5

This is cleaner:
- **Free space**: P = 0 → m_θ = 0 → exactly massless photon ✓
- **Braid core**: P ≈ 0.08 → m_θ ≈ 0.15 → Yukawa range ≈ 7 ✓
- **Between two braids**: P > 0 where tails overlap → m_θ > 0 → binding ✓
- **The "pion" IS the massive θ mode inside the P ≠ 0 region**

The Lagrangian term -½λ_θ P² |θ|² is:
- Sixth order in φ (P² = φ₀²φ₁²φ₂²), second order in θ
- Lorentz invariant ✓ (no derivatives)
- Local ✓
- Energy-conserving ✓ (Lagrangian → Noether)

## The Equation of Motion

With L_mass_θ = -½ λ_θ P² |θ|²:

### θ equation:
∂²θ_a/∂t² = ∇²θ_a - λ_θ P² θ_a + η(∇×φ)_a

The mass term λ_θ P² θ_a acts as a POSITION-DEPENDENT mass. Where P
is large (braid cores), θ oscillations are heavy and decay quickly.
Where P is zero (free space), θ is massless and propagates at c.

### φ equation:
∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P)∂P/∂φ_a
             + η(∇×θ)_a
             - λ_θ P (∂P/∂φ_a) |θ|²

The last term is new: it couples the φ field to the θ energy density
|θ|². Where θ is large AND P is nonzero, this creates an additional
force on φ. The sign depends on P and ∂P/∂φ — it can be attractive
or repulsive depending on the field configuration.

## Why This Might Create Binding

Consider two baryons at separation D ≈ 10:

1. Each baryon has a braid core where P ≈ 0.08 → m_θ ≈ 0.15
2. The θ radiation from each baryon has Yukawa range ~7 code units
3. At D=10, the baryons' θ Yukawa tails OVERLAP
4. The overlap creates a SHARED region of nonzero ⟨θ²⟩ between them
5. The shared θ energy couples to the φ field through the new term
6. If the coupling is attractive, it creates a binding well at D ≈ 10

The binding range is set by 2/m_θ_core ≈ 14 code units — matching the
V45 transition zone at D ≈ 10-15.

At D > 2/m_θ: the Yukawa tails don't overlap → no binding → pure EM (1/r)
At D < 2: phase confinement (P → 0 at overlap) → m_θ → 0 → no Yukawa
          → repulsion from gradient penalty

This gives: REPULSIVE at D < 2, ATTRACTIVE at D ≈ 5-14, COULOMB at D > 14.
The nuclear force curve — from a single Lagrangian term.

## The Derrick Scaling Check

The new term -½λ_θ P² |θ|² has no derivatives. Under x → λx:
- P² → P²(x/λ) (unchanged pointwise)
- |θ|² → |θ|²(x/λ) (unchanged pointwise)
- d³x → λ³

Total scaling: λ³ (same as mass term). This is a λ³ term.

For Derrick's virial: it adds to the E₃ sector (same as m²φ² and V(P)).
It doesn't change the virial structure — it's "just another potential."

BUT: unlike the pure φ potential V(P), this term couples φ AND θ.
The binding mechanism isn't through Derrick scaling — it's through
the YUKAWA RANGE of the θ field. The θ mass creates an exponentially
decaying attractive potential that doesn't require field deformation
at the interaction point. This EVADES Derrick's theorem because the
attraction is mediated by a SEPARATE field (θ), not by deformation
of the same field (φ).

This is exactly how nuclear binding works in real physics: the pion
(a separate field from the nucleon) mediates attraction without
requiring the nucleons to deform.

## Implementation

### Kernel changes

In the force computation, replace:

    // Old: constant m_theta
    force_theta[a] = laplacian_theta - m_theta2 * theta[a] + eta * curl_phi[a]

With:

    // New: field-dependent m_theta
    double P = phi[0] * phi[1] * phi[2];
    double m_theta2_eff = m_theta2 + lambda_theta * P * P;
    force_theta[a] = laplacian_theta - m_theta2_eff * theta[a] + eta * curl_phi[a]

And add the back-reaction on φ:

    // New: theta energy coupling to phi
    double theta2 = theta[0]*theta[0] + theta[1]*theta[1] + theta[2]*theta[2];
    double dPdphi_a = P / phi[a];  // = phi_b * phi_c
    force_phi[a] += -lambda_theta * P * dPdphi_a * theta2;

### Config

New parameter: `lambda_theta` (default 0 for backward compatibility)
Suggested starting value: λ_θ = 3.5 (gives m_θ ≈ 0.15 at braid core)

### Diagnostics

Report m_θ_eff_max = max over grid of √(λ_θ P²) to monitor the
effective mass. Should be ~0.1-0.2 at braid cores, ~0 in background.

## Summary

| Property | m_θ = const | m_θ² = λ_θ|φ|² | m_θ² = λ_θP² |
|----------|------------|----------------|--------------|
| Free-space photon | Massive ✗ | Slightly massive | **Exactly massless ✓** |
| Braid core θ | Massive | Massive | **Massive ✓** |
| Lorentz invariant | ✓ | ✓ | **✓** |
| Lagrangian | ✓ | ✓ | **✓** |
| Energy conserved | ✓ | ✓ | **✓** |
| No new fields | ✓ | ✓ | **✓** |
| Nuclear range | Fixed | Varies | **Varies with P ✓** |
| Free-space EM | Yukawa | Slight Yukawa | **Exact 1/r ✓** |
| Pion analog | Global | Local | **Topology-gated ✓** |

**m_θ² = λ_θP² is the recommended form.** It gives exactly massless
photons in free space (P=0) and massive mediators only where there's
topological content (P≠0). One new parameter λ_θ, one new Lagrangian
term, full Lorentz/energy/Lagrangian consistency.

The "pion" is not a separate field — it's the θ field itself, made
heavy by the local binding density P². The pion exists only where
baryons exist. Between two baryons, the overlapping P tails create
a shared massive-θ channel that mediates nuclear attraction.
