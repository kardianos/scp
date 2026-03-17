# V23-E: One Massless Mediator Field (Path A)

## Thesis

Break the S₃ permutation symmetry of the three fields: give mass m to
fields 1 and 2, but leave field 3 MASSLESS. The triple product P = φ₁φ₂φ₃
couples the massive sector (oscillon) to the massless mediator (field 3).

The time-averaged source ⟨φ₁φ₂⟩ ≠ 0 in the oscillon drives field 3 via:

    □φ₃ = -∂V/∂φ₃ = -μφ₁φ₂·P/(1+κP²)²

Since φ₃ is massless, the static solution is φ₃ ~ Q/(4πr) at long range,
giving a 1/r potential between oscillons.

This is the simplest possible path to 1/r gravity from the existing model.
The trade-off: it's scalar (spin-0, not spin-2) and breaks the S₃ symmetry.

## Mathematical Setup

### Modified Lagrangian

    L = Σ_{a=1,2} [½(∂φ_a)² - ½m²φ_a²] + ½(∂φ₃)² - V(φ₁φ₂φ₃)

Note: field 3 has NO mass term. Fields 1,2 have mass m=1.0.

### Equation of Motion

    ∂²φ₁/∂t² = ∇²φ₁ - m²φ₁ - μφ₂φ₃·P/(1+κP²)²       (massive)
    ∂²φ₂/∂t² = ∇²φ₂ - m²φ₂ - μφ₁φ₃·P/(1+κP²)²       (massive)
    ∂²φ₃/∂t² = ∇²φ₃         - μφ₁φ₂·P/(1+κP²)²       (MASSLESS)

### Expected Behavior

In the oscillon: φ₁ ≈ φ₂ ≈ A·f(r)·cos(ωt) (massive, breathing).
Field 3 is driven by the source S₃ = -μφ₁φ₂·P/(1+κP²)².

With P = φ₁φ₂φ₃ = A²f²cos²(ωt)·φ₃:

For small φ₃: S₃ ≈ -μA²f²cos²(ωt)·φ₃ (linear in φ₃)

This gives φ₃ an EFFECTIVE MASS inside the oscillon:

    m²_eff,3 = μA²f²⟨cos²⟩ = μA²f²/2

With μ=-20, A≈0.5: m²_eff,3 = -20·0.25/2 = -2.5 → TACHYONIC inside core.

This means φ₃ will be DRIVEN to a nonzero value inside the oscillon by the
tachyonic instability, creating a localized "charge" for the massless field.
Outside the oscillon, φ₃ satisfies □φ₃ = 0, so φ₃ ~ Q/(4πr) → 1/r tail.

### The Key Question

Does the oscillon survive when one of the three fields is massless?

Concern: the mass gap protects the oscillon from radiating. Field 3 has no
gap — it can radiate at any frequency. The triple product coupling could
cause the oscillon to drain its energy through field 3 radiation.

Counter-argument: the coupling goes through P = φ₁φ₂φ₃, which requires ALL
THREE fields to be nonzero. If φ₃ is small far from the core, the coupling
is weak. The oscillon might still be long-lived even if field 3 slowly leaks.

## What to Compute

### Phase 1: Single Oscillon with Massless Field 3 (1D)

1. Modify the v21 1D triad code: set m₃ = 0 while keeping m₁ = m₂ = 1.0.

2. Initialize all three fields with the same Gaussian (A=0.8, σ=3.0).

3. Evolve for t=10000. Monitor:
   a. Does the oscillon survive? (energy, f_core, peak amplitude)
   b. What is the breathing frequency ω?
   c. What happens to field 3? Does it grow, oscillate, or disperse?
   d. What is the spatial profile of field 3 at late times?
   e. Does field 3 develop a 1/r tail?

4. Compare with the standard case (m₃ = 1.0) as control.

### Phase 2: Scan m₃ from 1.0 to 0

5. Scan m₃ ∈ {1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.0}.

6. For each m₃: measure oscillon lifetime, ω, field 3 profile, radiation rate.

7. Find the critical m₃_c below which the oscillon dies.

8. If m₃ = 0 oscillons survive: measure the field 3 tail profile.
   Fit to φ₃(r) ~ Q/r at large r. Extract the "charge" Q.

### Phase 3: Two-Oscillon Interaction via Massless Mediator

9. If Phase 2 is positive: place two oscillons at separation D with m₃=0.

10. The equilibrated field 3 profiles overlap → interaction energy.

11. Measure the force vs D. Check if it scales as 1/D² (expected for 1/r
    potential in 1D... actually in 1D: φ₃ ~ Q·|x| gives F = const, not
    1/D². In 3D: φ₃ ~ Q/r gives F ~ 1/r². So the 1D test gives a constant
    force, which is still a distinctive signal.)

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (base code)
- Minimal change: make mass per-field (m[3] = {m, m, m3} instead of scalar m)

## Output Structure

- `src/mediator1d.c` — 1D three-field solver with per-field masses
- `data/` — output data (time series, profiles per m₃)
- `RESULTS.md` — results and analysis

## Parameters

μ=-20, κ=20, m₁=m₂=1.0, m₃ scanned from 1.0 to 0.0
Grid: Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o mediator1d src/mediator1d.c -lm`
