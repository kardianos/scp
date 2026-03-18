# V30: Rotating Expansion — Unified Attack on Pins 2,3,4,6,7

## The Idea

Start with a DENSE, ROTATING three-field configuration. Expand it.
As the rotating dense field expands, angular momentum conservation
forces the formation of vortex lines. With three coupled fields and
triple-product binding, these vortices naturally BRAID.

This is how quantized vortices form in spinning superfluids — but here
the three-field coupling creates braided defects instead of single vortices.

## Why This Addresses Five Pins at Once

| Pin | Question | How expansion+rotation answers it |
|-----|----------|-----------------------------------|
| 7 | Do braids form spontaneously? | Rotation creates vortex lines → triple coupling braids them |
| 2 | Where does mass come from? | Mass = field energy trapped in each braid during expansion |
| 3 | Where does charge come from? | Charge = angular momentum carried by each braid |
| 4 | Do braids interact? | Multiple braids form simultaneously → observe interaction |
| 6 | What is the substrate? | The expanding field IS the substrate; its structure is observable |

## Physical Picture

```
t=0: Dense rotating blob        t=T₁: Expansion + vortex formation
┌─────────┐                     ┌─────────────────────┐
│ ░░░░░░░ │ ← high density      │  ·    ⊗    ·       │
│ ░░⟳░░░ │ ← rotating          │     ·    ⊙   ·     │ ← vortex lines
│ ░░░░░░░ │                     │  ⊗    ·    ·   ⊗   │    (braids?)
└─────────┘                     └─────────────────────┘

t=T₂: Braids formed             t=T₃: Braids interact
┌─────────────────────┐         ┌─────────────────────┐
│         ·           │         │                     │
│    ⊗←——→⊙          │         │    ⊗——→⊙            │ ← attract/repel?
│         ·           │         │                     │
└─────────────────────┘         └─────────────────────┘
```

## Variants (run as a sequence)

### V30-A: Dense Rotating Blob → Expand → Count Defects (Pin 7)

Initialize: all three fields at high amplitude A₀ = 2.0 (dense) in a
sphere of radius R₀ = 5, with rigid-body rotation Ω around z-axis.

    φ_a(x,0) = A₀ × exp(-r²/(2R₀²)) × cos(k·z + 2πa/3 + Ω×θ_perp)
    v_a(x,0) = matching velocity from rotation

where θ_perp = atan2(y,x) is the azimuthal angle in the xy-plane.
The Ω×θ_perp term gives the azimuthal phase twist (rotation).

Let it evolve freely. The dense blob expands under pressure (gradient
energy), and the rotation creates vortex lines where the phase wraps.

Measure at intervals: number of localized energy peaks (braids),
their positions, winding numbers, and the background field between them.

Scan: Ω ∈ {0.0, 0.1, 0.5, 1.0, 2.0} — how much rotation is needed
for braid formation? Ω=0 is the control (T10D showed no braids without
rotation).

### V30-B: Complex Fields + Rotation → Charge Spectrum (Pin 3)

Same as V30-A but with complex fields (6 real DOF → 12).
The U(1) symmetry gives conserved charge Q = ∫ Im(ψ* ∂_t ψ) d³x.

After expansion, measure the charge of each formed braid.
Is it quantized? Does it correlate with winding number?
Do opposite-charge braids form in pairs (charge conservation)?

### V30-C: Post-Expansion Interaction (Pin 4)

After V30-A produces braids, continue the simulation and track:
- Separation between braids D(t) — attract or repel?
- Is the force charge-dependent (same vs opposite rotation)?
- Does the depleted background between braids create a gradient?

This might be the FIRST test where braids attract via depletion,
because they formed IN the background (not inserted by hand), so
the background has the right structure.

### V30-D: Expansion Rate → Mass Spectrum (Pin 2)

Vary the initial density A₀ ∈ {1.0, 1.5, 2.0, 3.0, 5.0} with fixed Ω.
Higher density → more energy to partition → heavier braids?

Or vary R₀ (initial radius) with fixed A₀:
R₀ ∈ {3, 5, 8, 12} — larger initial region → more braids?

The mass of each braid = its energy after separation from the others.
Does mass depend on expansion dynamics (initial conditions) or only
on the final braid topology (universal)?

### V30-E: Background Characterization (Pin 6)

After braids form, characterize the FIELD between them:
- Energy density profile ρ(r) between braids
- Is there a condensate (uniform |φ|)? Plane waves? Thermal noise?
- What is the "vacuum state" that the braids sit in?
- Does the background structure match the M7 "B field" that we
  imposed by hand in T12?

This tells us what the substrate NATURALLY looks like, rather than
what we assumed (plane waves with A_bg=0.1).

## Implementation

### Shared Core
Use braid_core.h for the PDE solver. The initialization is NEW
(dense rotating blob, not helical braid). The force computation
is standard (triple-product potential, same as all V28/V29 tests).

For V30-B: extend to complex fields (double the arrays, complex
triple product V(|ψ₀ψ₁ψ₂|²)).

### Grid
- V30-A: N=128, L=40, T=500 (enough room for expansion + defect formation)
- V30-B: N=96, L=40, T=500 (complex fields = 2× memory, reduce N)
- V30-C: Continue V30-A to T=1000
- V30-D: N=128, L=40, T=500 × 5 density values
- V30-E: Analysis of V30-A/C output, minimal new simulation

### Diagnostics
New: defect detection (find localized energy peaks, measure their
winding, position, and charge). Plus all standard metrics from
braid_core.h (fc, l2, torsion, winding).

## M7 Two-Component Option

Alternative: use M7 (S+B separation) from the start.
Initialize S=0 (no braids) and B=dense rotating blob.
If braids form, they form in the S fields (through B→S conversion).
The B fields remain as the background.

This would directly connect to the T12 depletion results and test
whether the M7 framework supports spontaneous braid formation.

## Parameters

Lagrangian: same as V28/V29 bimodal sweet spot
    μ = -41.3, κ = 50, m = 1.50 (mass provides confinement)

Initial condition (V30-A):
    A₀ = 2.0 (dense, about 2.5× the equilibrium braid amplitude)
    R₀ = 5.0 (initial blob radius)
    Ω = 1.0 (rotation rate, scan 0-2)
    k = π/L (axial wavenumber, same as bimodal)

Boundary: absorbing at r > 0.85L (expansion products absorbed at edges)

## Expected Outcomes

**Best case**: Dense rotating blob expands and fragments into multiple
braids with quantized winding and charge, sitting in a naturally
depleted background. Braids interact through the depletion gradient.
Mass proportional to trapped angular momentum. This would confirm
the entire framework in one experiment.

**Likely case**: Some vortex-like structures form but are unstable
or don't clearly braid. The three-field coupling creates correlated
defects but not clean braids. Background structure is observable but
messy. Still useful — constrains what initialization is needed.

**Worst case**: The dense blob just disperses uniformly (like T10D
without cooling). Rotation doesn't help. No defects form. This would
mean braids require MORE specific conditions than just density +
rotation, pushing toward the "braids are metastable artifacts" conclusion.

## Priority

V30-A first (does rotation create braids?). This is the gate.
If YES: proceed with B, C, D, E.
If NO: try M7 variant, different Ω, different A₀, before giving up.
