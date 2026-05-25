# Shared Vocabulary: Constraints, Troughs, Wells, and Stability Bounds

**Purpose**: When we say a "free variable" or algebraic choice is not free, we mean it is fixed by a physical bound or extremum in the medium. This short note collects the recurring physical concepts we will use when making hypotheses quantitative. Every future calculation or Lean statement should be explicit about which of these (or new ones) is doing the forcing.

---

## Core Concepts

### 1. Density Well (ρ_M well)
- Local maximum (or sufficiently deep local minimum of an effective potential) in the medium's density deviation `ρ_M` (or equivalent measure: norm, connectivity, number of closed relations).
- A configuration is in a density well if small perturbations do not cause it to radiate or disperse on timescales short compared to the age of the system.
- Depth of the well = how much higher the local density is than the background vacuum.
- Width of the well = how much "room" the configuration has before it becomes unstable.

The observed particles/generations/sectors are the known stable occupants of such wells.

### 2. Protection Budget / Protection Cost
- The fraction of available degrees of freedom (or algebraic components, or internal relations) that must be "spent" to keep the high-density state from leaking.
- Examples in the current work:
  - The S³ radius `|ξ| = 1/√2` (fraction of quaternionic freedom committed to protection).
  - Restricted bivector support in the v58 protected-chirality mechanism.
  - The choice of which graded pieces (L vs. F) to activate.
- Every protection technology has a cost. Stacking L + F for u-quarks is only viable because the algebra permits additive combination without double-paying the cost in the long-range channels.

### 3. Relationship Trough (or Causal Link Maximization)
- In a pre-geometric or relational substrate: the configuration that maximizes the number (or usefulness) of closed relations, causal links, or algebraic compositions per emergent macroscopic volume, without over-constraining the light-cone for small fluctuations.
- A "trough" because moving away from the optimum in either direction (too few relations = low density; too many = rigidity or dispersion) raises the effective energy or reduces stability.
- The 2/9 phase and triality may be the geometric expression of three equivalent ways to maximize closed relations inside the same bounded algebraic cell.

### 4. Geometric / Causal Stability Margin
- The amount by which the emergent causal structure (light-cone, propagation speed, absence of dispersion) remains robust when the local medium state (density, ξ, protection level) is varied within the range occupied by real particles and gradients.
- Non-dispersive universality (GW170817-level) is not an input; it is a *stability requirement* on any viable protection + density mechanism. Configurations that would introduce frequency-dependent effects or species-dependent refraction are disallowed.

### 5. Force-Channel Separation Budget
- The "room" left in the algebraic structure, after protection costs are paid, for the Newtonian (grade-1, density-sourced) and Maxwell (grade-2, bivector) responses to remain linearly independent and to couple to density gradients in a controlled way (`f(ρ)` modulation without destroying the other channel).
- The v58 living candidate lives inside a safe band where this separation holds. Any algebraic choice that pushes the system outside that band (too much cross-talk) is physically disallowed for long-lived objects.

### 6. Energy / Action Trough (Effective Potential Minimum)
- When an effective description exists (or can be derived), the stable configurations sit at local minima of some functional that balances the drive for higher density against the costs of protection, gradients, and radiation.
- The multi-well structure for ξ (different minima for different sectors) is an example.

---

## How to Use This Vocabulary in Future Work

When proposing or refining a hypothesis, answer explicitly:

- Which of the above (or a new named bound) is doing the *forcing*?
- What would happen (quantitatively or structurally) if the bound were violated (e.g., if the protection budget were 10% lower, or the relationship trough were shallower)?
- What measurable signature (in the multivector code, in a small Lean model, or phenomenologically) would reveal the operation of that particular bound?

Example (for Hypothesis 2 above):

"The radius 1/√2 is forced by the requirement that the protection budget (S³ radius) simultaneously satisfy the density-well depth needed for observed lepton masses *and* the force-channel separation budget needed for the v58 living candidate to remain inside its safe band. If the radius were a free parameter, we would expect a continuum of lepton spectra; the physics selects only the value that keeps both the density well and the separation budget non-empty."

---

This document is deliberately short and will grow only by addition of new named concepts when they prove necessary. The goal is shared, precise language so that quantitative work and relationship-forcing analysis stay coupled.