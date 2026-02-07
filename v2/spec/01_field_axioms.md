# 01 — Field Axioms

This chapter establishes the foundational postulates of Chiral Harmonic Particle Theory (CHPT). Everything else in the theory must be derivable from these axioms or explicitly added as additional postulates. Axioms are numbered for reference throughout the spec.

---

## Axiom 1 — Existence of the Field

There exists a single, continuous field defined over three spatial dimensions and one temporal dimension (3+1D). This field is the sole fundamental entity. All matter, energy, forces, and observable phenomena are configurations or dynamics of this field.

There is nothing else. No separate particles, no independent force carriers, no additional fields. One field, one substance.

### Field Type (Resolved)

The field $\Psi$ is valued in the **even subalgebra** $Cl^+(3,0,1)$ of Projective Geometric Algebra (see [math/01_algebra.md](math/01_algebra.md)). This is isomorphic to the Dual Quaternions and decomposes as:

$$ \Psi = \rho + \vec{J} e_0 + \mathbf{F} + I \tau $$

with scalar (mass density), degenerate bivector (flux), spatial bivector (field strength), and pseudoscalar (helicity) components.

### Remaining Unknowns

- **Topology of the field domain**: Is the spatial domain R^3 (infinite flat space), a compact manifold (e.g., 3-torus), or something else? This matters for boundary conditions and total energy accounting.
- **Whether spacetime is fundamental or emergent**: The axiom as stated assumes a pre-existing 3+1D spacetime backdrop. An alternative (more ambitious) formulation would have spacetime itself emerge from the field. The current proposal assumes the former, which means Lorentz structure must either be imposed or derived.
- **Time evolution**: The degenerate basis $e_0$ ($e_0^2 = 0$) is projective, not timelike. Time must enter either via a separate parameter $\partial_t$ or by embedding in a larger algebra. See [15_open_problems.md](15_open_problems.md), A7.

---

## Axiom 2 — Field Density

The field possesses a well-defined, non-negative, real-valued local density at every point in space and time. This density is the single measurable quantity of the field.

All physical quantities — energy, mass, momentum, charge — ultimately reduce to descriptions of this density distribution and its dynamics.

### Clarification

"Density" here is not mass-density or energy-density in the conventional sense. It is the primitive quantity of the field itself. It has no prior physical definition; physical meaning is assigned to it through the consequences of the axioms.

### Unknowns

- **Units and scale**: What sets the absolute scale of density? Is there a natural unit, or is only relative density physically meaningful?
- **Minimum and maximum**: Is density bounded? Can it reach zero (true vacuum) or infinity (singularity)? The knot mechanism suggests local density can be much higher than background, but whether there is a hard upper bound is unspecified.
- **Relationship to field type (Clarified)**: With the field valued in $Cl^+(3,0,1)$, "density" $\rho$ is the **scalar (grade-0) component** of $\Psi$. The full field norm $|\Psi|^2 = \Psi\widetilde{\Psi}$ includes contributions from all components (scalar, bivector, pseudoscalar). The Axiom 2 "density" maps to $\rho$; the potential $V = (\lambda/4)(|\Psi|^2 - \rho_0^2)^2$ constrains the full norm. See [math/01_algebra.md](math/01_algebra.md).

---

## Axiom 3 — Conservation of Total Density

The total integrated density over all of space is constant in time. Density is neither created nor destroyed; it is only redistributed.

Formally: if rho(x,t) is the density at point x and time t, then the integral of rho over all space is a constant independent of t.

### Consequences

- This is the foundational conservation law from which energy conservation will be derived.
- Any process that increases local density (knot formation) must decrease density elsewhere (field depletion).
- There is no mechanism for energy input or output to/from the field. The field is a closed system.

### Tension with Multivector Field

With $\Psi \in Cl^+(3,0,1)$, the scalar component $\rho$ is only one part of the field. Two distinct conservation laws are in play:

1. **Density conservation** (this axiom): $\int \rho \, d^3x = \text{const}$. This conserves the scalar "stuff."
2. **Energy conservation** (from the Lagrangian): $\int T_{00} \, d^3x = \text{const}$, where $T_{00} = \frac{1}{2}|\nabla\Psi|^2 + V(\Psi)$ includes contributions from all field components.

These are **not the same quantity**. A process that converts scalar density $\rho$ into bivector field strength $\mathbf{F}$ (e.g., a knot emitting radiation) conserves energy but changes the distribution of $\rho$. Whether total $\int \rho \, d^3x$ is separately conserved depends on whether the field equation admits a corresponding Noether current. This must be checked against the proposed Lagrangian in [math/03_dynamics.md](math/03_dynamics.md).

### Unknowns

- **What enforces conservation**: The proposed Lagrangian ([math/03_dynamics.md](math/03_dynamics.md)) provides energy-momentum conservation via Noether's theorem. Whether $\int \rho \, d^3x$ is independently conserved (as this axiom claims) is an open question — it requires an internal symmetry of the Lagrangian that rotates the scalar into other components.

---

## Axiom 4 — No Intrinsic Dissipation

The field dynamics are fundamentally non-dissipative. There is no friction, viscosity, or drag intrinsic to the field. A disturbance, once created, propagates indefinitely without loss of amplitude from dissipation.

### Clarification

This does not mean energy cannot be redistributed, radiated, or transformed. It means there is no mechanism that converts field energy into some other form (heat, entropy sink, etc.) because there IS no other form. All energy remains in the field.

Apparent dissipation — such as a knot losing energy by emitting null-rotor radiation — is redistribution of field density, not loss. The emitted radiation continues to exist as field configuration.

### Tension with Observations

- **Thermodynamic irreversibility**: Macroscopic irreversibility (entropy increase) must emerge from the microscopic dynamics without being built in. This is plausible (statistical mechanics achieves this for classical mechanics) but non-trivial to demonstrate.
- **Decoherence**: Quantum decoherence looks like dissipation. CHPT must account for this as information dispersal into field degrees of freedom, not true loss.

---

## Axiom 5 — Finite Maximum Process Speed (c)

There exists a finite, invariant maximum speed c at which the *local state* of the field can update or reconfigure. No information or energy transfer (signalling) can occur faster than c.

### Clarification on Non-Locality
This axiom enforces **Signal Locality** (no superluminal communication). It does **not** forbid **Topological Non-Locality**.
*   A knot is a topological structure designated by the field configuration.
*   If a single topological structure extends across space (e.g., a flux tube or entangled pair), changes to the *global topology* can be instantaneous constraints, provided they do not allow superluminal signaling.
*   This distinction allows CHPT to be consistent with quantum non-locality (Bell's theorem) without violating the field's differential structure.

### Critical Unknowns
- **What determines c**: Is c set by a parameter in the field equation, or does it emerge from deeper structure? In standard physics, c arises from the permittivity and permeability of the vacuum (c = 1/sqrt(epsilon_0 * mu_0)). CHPT needs an analogous derivation or must take c as a free parameter.
- **Invariance**: Stating c is "invariant" implies Lorentz invariance. The field equation must be constructed such that the "process speed" is frame-independent.

---

## Axiom 6 — The Field Supports Stable Process Structures (Knots)

The field dynamics admit solutions that are spatially localized, self-sustaining, and persistent **processes**. These are called **knots**.

### Clarification: Process Ontology
A knot is not a static object moving *through* the field (like a marble in water). A knot is a dynamic **event**—a standing wave of "knotting and unknotting" encoded in the field.
*   **Motion**: When a knot moves, the field does not flow. Rather, the *pattern of knotting* propagates. The field reconfigures itself in the direction of motion, while relaxing in the wake.
*   **Mass**: Mass is not "stuff" carried along; it is the inertia of this reconfiguring process. It takes energy to drive the "unknotting/knotting" cycle forward, which we perceive as kinetic energy.

### Stability Mechanisms
The stability of these processes is primarily **Topological**. The field configuration belongs to a non-trivial homotopy class (due to the multivector structure of the field algebra). A knot cannot simply "disperse" because there is no continuous deformation that unwinds it to the vacuum state. This connects CHPT to Modern Knot Theory (Vassiliev invariants, Khovanov homology) as the basis for the particle spectrum.

---

## Summary of Axioms

| # | Axiom | Status |
|---|-------|--------|
| 1 | Single continuous 3+1D field | Postulated; field type adopted: $Cl^+(3,0,1)$ (see [math/01_algebra.md](math/01_algebra.md)) |
| 2 | Well-defined local density | $\rho$ = scalar component of $\Psi$; units/bounds TBD |
| 3 | Total density conserved | Tension: $\int\rho$ vs $\int T_{00}$ — separate conservation TBD |
| 4 | No intrinsic dissipation | Postulated; emergent irreversibility TBD |
| 5 | Finite max propagation speed c | Postulated; invariance group TBD |
| 6 | Stable localized configurations (knots) | Postulated; existence proof TBD |

### Optional Methodology: Minimal vs. Extended Axiom Sets

Two paths forward:

1. **Minimal**: Accept only Axioms 1-5 and attempt to DERIVE Axiom 6 (knot existence) from the field equation. This is more powerful but requires finding a specific equation that naturally produces soliton-like solutions. Precedent exists: sine-Gordon equation, Skyrme model, etc.

2. **Extended**: Accept all six axioms and focus on deriving consequences. Axiom 6 becomes a constraint on admissible field equations rather than a prediction.

Both are valid approaches. The minimal path is scientifically stronger but mathematically harder. The extended path lets conceptual development proceed while the existence question is resolved in parallel.
