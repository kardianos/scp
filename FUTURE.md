# SCP Future Directions and Open Questions

**Purpose**: This document tracks future research directions, unresolved
questions, and exploratory paths that have been identified through
development (V21-V33), internal analysis, and external review. Items
here are BEYOND the current confirmed results but represent viable
directions for investigation. Maintain this document as directions are
explored, resolved, or abandoned. Keep current with CONCEPT.md.

---

## Critical Priority (blocks fundamental claims)

### F1: Mass Parameter and Long-Range Gravity (Largely Resolved)

**The original problem**: m²=2.25 creates Yukawa decay (range ~0.67).
Individual field excitations die as e^{-mr}/r. How can gravity be long-range?

**The resolution (V34 phonon test)**: The depletion IS NOT Yukawa. Direct
measurement (N=256, L=60, T=200) shows δρ ∝ 1/r^1.2 (power law, R²=0.98).
Yukawa m=1.5 is excluded at 500,000× discrepancy (R²=-0.44). When the
Yukawa mass is fit freely, it converges to m≈0.02 — effectively zero.

The braid's collective perturbation of the background propagates as a
PHONON (massless Goldstone mode of the background's broken translational
symmetry). The m² parameter confines the braid but does NOT limit the
gravitational range. Mass and range are decoupled.

**What V34 Track G showed**: m²≥1.25 needed for braid binding. m²≥0.25
for vacuum stability. Braid survival is the constraint, not range.
Field-dependent mass (Track GB) does not help — constant m² is optimal.

**Remaining questions**:
- Why is the depletion exponent n≈1.2 (not 2.0)? Need longer runs (T=1000+)
  with equilibrated braids and isotropic backgrounds
- Analytically identify the massless phonon mode in the linearized equation
- Does the force exponent approach 2.0 with better equilibration?
- Is the n=1.2 depletion exponent a feature of the z-oriented background
  (anisotropy artifact) or fundamental?

**Source**: V34 phonon_test, V34 G_metastability, V34 GB_field_mass.

### F2: Gravity Mechanism — Resolved (Dynamic Footprint, NOT Energy Minimization)

**The original tension**: Is the force from DYNAMIC processing (intake
asymmetry → momentum) or from ENERGY MINIMIZATION (overlapping
depletions → lower total energy)?

**ANSWER**: Neither simple picture is correct. The force is a DYNAMIC
response to the local density gradient, linear in ∇ρ. It is NOT
derivable from the static energy landscape (which is purely repulsive).

**Complete experimental findings** (V33, March 2026):

1. **Footprint asymmetry CONFIRMED**: In a ρ gradient, the braid's
   perturbation profile (δρ above background) is asymmetric. The
   half-width on the low-ρ side is 1.09–1.57× larger than on the
   high-ρ side (5/6 time snapshots, N=128 gradient test).

   Mechanism: the effective mass m_eff² = m² + V''(P_bg) is lower
   where ρ is lower → longer Yukawa range → perturbation extends
   further into the depleted side.

2. **Asymmetric drag REFUTED**: A kicked braid in high-ρ (A_bg=0.15)
   retains/amplifies momentum (1.40× over T=125), while the same kick
   in low-ρ (A_bg=0.05) decays to 0.30×. Coupling is STRONGER in
   high ρ, not weaker. The braid has more "grip" in dense field.
   Code: v33/src/v33_drag_test.c. Data: v33/data/drag_{high,mid,low}/.

3. **F ∝ ∇ρ confirmed (R² = 0.9998)**: Gradient sweep at 4 strengths
   (ρ ratios 1.22, 1.49, 3.45, 9.00) shows drift perfectly proportional
   to ∇ρ. The ratio drift/∇ρ = -186 ± 3 (constant across all gradients).
   This is exactly F = -C × ∇ρ, the form of Newtonian gravity.
   Code: v33/run_f2_tests.sh. Data: v33/data/f2_gradient/.

4. **F ≠ -dE/dD (energy minimization REFUTED)**: Total energy E(D) for
   two braids at separation D=8–50 is monotonically DECREASING with D.
   Interaction energy E_int(D) = E_pair(D) - E_pair(D→∞) is:
     D=8: +2939, D=12: +808, D=15: +236, D=20: +16, D=30: +0.5
   ALL POSITIVE (repulsive). No attractive well at any D. But C1
   measured ATTRACTION at D=12–80. Therefore the attractive force
   CANNOT come from static energy minimization.
   Code: v33/src/v33_energy_vs_D.c. Data: v33/data/f2_energy/.

5. **Footprint shift scales linearly with ∇ρ**: At t=0, the field-weighted
   center shift is proportional to gradient strength (ratios 1:2:6
   matching gradient ratios 1:2:6 exactly). The geometric footprint
   asymmetry provides the mechanism; the gradient sets the magnitude.

**Summary of mechanism**:
- The braid's oscillation cycle creates an asymmetric spatial perturbation
  in a ρ gradient (footprint extends further into depleted side)
- This asymmetry produces a net force proportional to ∇ρ
- The force is DYNAMIC (from the braid's ongoing interaction with the
  gradient field), not from static energy minimization
- The coupling strength is proportional to ρ (stronger in dense field)
- The force direction is toward depletion = toward other braids = GRAVITY

**Still open**:
- Precise analytical derivation of the proportionality constant C ≈ 186
- Connection between C and braid parameters (m, μ, κ, braid size)
- N=512 gradient test in progress — higher resolution confirmation

**Source**: Gemini-02 §1, V33 tests (footprint, drag, energy_vs_D, gradient sweep).

---

## High Priority (needed for physical viability)

### F3: Lorentz Contraction Verification

**Test**: Boost a braid at v=0.1c, 0.3c, 0.5c. Does the shape contract
by γ? Does the oscillation frequency slow (time dilation)?

**Method**: Initialize φ_a(x) → φ_a(γ(x-vt)) with corresponding velocity.
Check aspect ratio and internal frequency at T=100.

**Expected**: YES — the equation is Lorentz-invariant, so boosted solutions
must transform correctly. This would confirm emergent Special Relativity.

**Challenges**: Grid resolution during contraction (at v=0.5c, γ=1.15,
the braid is 15% shorter in one direction — needs adequate dx).

**Source**: Gemini-01 §2, Gemini-02 §3, Grok-02 §2.

### F4: Isotropic Background

**Test**: Initialize the background with random phases across all directions
(not just z-aligned). Do braids survive? Is the force isotropic?

**Why it matters**: The current z-preferred background is a simulation
convenience. If the model requires it, there's a hidden preferred-frame
problem. If braids work in an isotropic background, the model is
genuinely Lorentz-invariant in practice, not just in principle.

**Source**: Gemini-02 §3, Grok-02 §4.

### F5: Equivalence Principle — Do All Energy Forms Gravitate?

**Question**: In GR, gravity couples to ALL energy (including radiation,
EM fields, kinetic energy). In our model, gravity comes from depletion
of the φ field. If EM is added (complex fields + gauge), do EM waves
also deplete the background?

**Test**: Once EM is integrated (see F8), send a localized EM wave packet
through the field. Does it create a depletion zone? Does a braid respond
to it gravitationally?

**If NO**: the model violates the equivalence principle. Only "matter"
(braids) gravitates, not radiation. This would be a fundamental failure.

**Source**: Gemini-02 §4.

---

## Medium Priority (extensions and refinements)

### F6: Electromagnetism — Cosserat Angle Fields (CONFIRMED CHARGE-DEPENDENT FORCE)

**The path**: 3 angle fields θ_a coupled to 3 position fields φ_a via curl.
Massless θ (m_θ=0) sourced by the braid's helical twist.

**Confirmed results** (V34 θ characterization):

1. **Charge-dependent force**: Same-winding braids attract 27% more than
   the 3-field (gravity-only) baseline. Opposite-winding braids attract
   57% less. The θ field mediates a force whose sign depends on winding.

2. **Winding = charge**: Time-averaged θ_φ reverses sign when winding
   reverses (far-field ratio ≈ -1.0). The DC component carries charge
   information.

3. **Right-hand rule**: Volumetric visualization shows circular θ patterns
   perpendicular to the braid axis. Confirmed dt-converged (2% variation
   across 4× dt range).

4. **Wave-mediated, not static**: The θ field oscillates (period ~4t) with
   99.8% wave and 0.2% DC bias. The force is mediated by wave exchange
   (QFT-like), not by a static 1/r field (Biot-Savart).

**Remaining questions**:
- Is there a static (Coulomb) regime for non-moving braids?
- What determines the θ/φ force ratio (analog of fine structure constant)?
- Does θ support bound orbital modes (electron analog)?
- Can the wave-exchange picture be connected to virtual photon exchange?

**Source**: V34 θ characterization, v34/torsion_coupling/theta_characterize/RESULTS.md.

### F7: Multi-Braid Formation (Nucleosynthesis)

**Concept**: Heavier "atoms" are braids sharing helical structure, formed
by condensation from hot dense field (stellar nucleosynthesis analog),
NOT by collision of existing braids (which scatter, V33-C3).

**Test**: High-density field initialization with controlled cooling.
V30's expansion tests failed (no spontaneous braids from FRW expansion).
Need higher density, slower cooling, or different phase structure.

**Questions**:
- What initial conditions produce multi-braid bound states?
- Is there a critical temperature/density for braid condensation?
- Do different cooling rates produce different "elements"?

**Source**: User concept (session discussion), V30 results.

### F8: Binding Energy Accounting

**Question**: When a braid forms, does the total field energy DECREASE?
(True binding energy, like nuclear binding.) Or is it just rearranged?

**Test**: Compare E_total of (braid + depleted background) vs
(uniform background at same total field content). If E_braid < E_uniform,
there is genuine binding energy and the braid is energetically favored.

**Source**: Grok-02 §3.

### F9: Analytical Effective Potential

**Goal**: Derive the effective braid-braid interaction potential V_eff(D)
analytically from linearized perturbation theory around the braid.

**Approach**: Treat the second braid as a test perturbation in the first
braid's depletion field. Solve the linearized wave equation for the
perturbation. The overlap integral gives V_eff(D).

This would predict the force law WITHOUT simulation — a theoretical
confirmation of the numerical results.

**Source**: Grok-01 suggestions.

---

## Lower Priority (long-term / speculative)

### F10: Spin and Helical Handedness

**Question**: The braid has a helical twist (winding number W = ±1).
Does this correspond to spin? Left-handed vs right-handed braids would
be spin-up vs spin-down analogs.

**Test**: Do left-handed and right-handed braids interact differently?
Is there a spin-statistics connection?

### F11: Gravitational Waves

**Question**: Does an accelerated braid produce radiation with spin-2
(quadrupolar) tensor structure? V29-T5 showed 63% quadrupolar strain
at R=10, but this was static geometry, not radiation.

**Test**: Accelerate a braid (give it a kick), measure the tensor
structure of the energy radiation at far field.

### F12: Dark Matter Profiles

**Question**: Does the depletion profile around a braid match observed
dark matter halo profiles (NFW, Burkert)? V30 showed accretion +
depletion structure. Quantitative comparison needed.

### F13: Quantization

**Question**: The theory is purely classical. What happens when quantized?
Does the triple-product potential have a well-defined quantum field theory?
Are the braids stable quantum mechanically?

### F14: Cosmological Constant

**Question**: The background field has nonzero energy density ρ_bg ≈ 0.03.
This acts like a cosmological constant. Can its value be related to
the observed dark energy density?

### F15: Background Origin

**Question**: What sets m², μ, κ, and A_bg? Are these fundamental
constants, or do they emerge from a deeper theory (T9 substrate)?

### F16: Gradient Test Long-Term Conservation

**Question**: In v33_gradient_test.c, the pinned x-boundaries impose a
fixed ρ gradient. Over long runs (T >> 100), does total system energy
stay flat, or is there slow injection/drift from the artificial BC?

**Source**: Grok-02 §7.

---

## Resolved / Abandoned Directions

### X1: S/B Two-Component Split (V29-V31)
Abandoned. The split is artificial. The standard equation on a single
field produces the same physics without the split.

### X2: c(ρ) Speed-of-Light Modification (V30-V31)
Abandoned. c faster in dense → blowup. c slower in dense → freezing.
M7+inverted c works but is artificial. Not needed — standard equation
produces gravity without c modification.

### X3: Binding-Weighted Gradient Coupling (V32)
Abandoned. The gradient coupling is net repulsive when properly
controlled. The intrinsic attraction from the standard equation is
the real mechanism. The binding anatomy (core/surface/fabric) is real
and useful for understanding, but shouldn't modify the equation.

### X4: SPH Particle-Based Simulation (V32)
Shelved. Viable (braid survives, particles cluster) but gravity signal
too weak (0.9%) and energy conservation poor (+78% drift). May revisit
for visualization or adaptive resolution.

### X5: Rotating FRW Expansion (V30)
Failed. Boundary artifacts, shell fragmentation, no bulk braid formation.
Braids require specific initialization, not cosmological expansion alone.

### X6: Asymmetric Drag Hypothesis (V33 Drag Test)
Refuted. The hypothesis was: motion cost is proportional to ρ (not length),
so braids decelerate slower in low-ρ → preferential drift toward depletion.
Drag test (v_kick=0.1 in A_bg=0.05/0.10/0.15) showed the OPPOSITE:
momentum retention is 0.30× in low-ρ vs 1.40× in high-ρ. Coupling
(both force and drag) is stronger in high-ρ, not weaker. The gravity
mechanism is geometric (asymmetric footprint), not dynamic (friction).
Code: v33/src/v33_drag_test.c. Data: v33/data/drag_{high,mid,low}/.
