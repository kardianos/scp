# 01 — Field Axioms

This chapter establishes the foundational postulates of Chiral Harmonic Particle Theory (CHPT). Everything else in the theory must be derivable from these axioms or explicitly added as additional postulates. Axioms are numbered for reference throughout the spec.

---

## Axiom 1 — Existence of the Field

There exists a single, continuous field defined over three spatial dimensions and one temporal dimension (3+1D). This field is the sole fundamental entity. All matter, energy, forces, and observable phenomena are configurations or dynamics of this field.

There is nothing else. No separate particles, no independent force carriers, no additional fields. One field, one substance.

### Unknowns

- **Mathematical type of the field**: The proposal does not yet specify whether the field is scalar-valued (a single number at each point), vector-valued, tensor-valued, or spinor-valued. A scalar field is the simplest candidate, but may be insufficient to support chirality and the full rotational structure needed for particle physics. A field valued in a geometric algebra (specifically Cl(3,0,1) for 3D PGA) is the most natural candidate given the null-rotor machinery invoked later (see [06_null_rotors.md](06_null_rotors.md)), but this is not yet established.
- **Topology of the field domain**: Is the spatial domain R^3 (infinite flat space), a compact manifold (e.g., 3-torus), or something else? This matters for boundary conditions and total energy accounting.
- **Whether spacetime is fundamental or emergent**: The axiom as stated assumes a pre-existing 3+1D spacetime backdrop. An alternative (more ambitious) formulation would have spacetime itself emerge from the field. The current proposal assumes the former, which means Lorentz structure must either be imposed or derived.

---

## Axiom 2 — Field Density

The field possesses a well-defined, non-negative, real-valued local density at every point in space and time. This density is the single measurable quantity of the field.

All physical quantities — energy, mass, momentum, charge — ultimately reduce to descriptions of this density distribution and its dynamics.

### Clarification

"Density" here is not mass-density or energy-density in the conventional sense. It is the primitive quantity of the field itself. It has no prior physical definition; physical meaning is assigned to it through the consequences of the axioms.

### Unknowns

- **Units and scale**: What sets the absolute scale of density? Is there a natural unit, or is only relative density physically meaningful?
- **Minimum and maximum**: Is density bounded? Can it reach zero (true vacuum) or infinity (singularity)? The knot mechanism suggests local density can be much higher than background, but whether there is a hard upper bound is unspecified.
- **Relationship to field type**: If the field is not purely scalar, "density" may refer to the magnitude or norm of a multi-component field value. This needs to be made precise once the field type (Axiom 1) is resolved.

---

## Axiom 3 — Conservation of Total Density

The total integrated density over all of space is constant in time. Density is neither created nor destroyed; it is only redistributed.

Formally: if rho(x,t) is the density at point x and time t, then the integral of rho over all space is a constant independent of t.

### Consequences

- This is the foundational conservation law from which energy conservation will be derived.
- Any process that increases local density (knot formation) must decrease density elsewhere (field depletion).
- There is no mechanism for energy input or output to/from the field. The field is a closed system.

### Unknowns

- **What enforces conservation**: Is conservation a consequence of the field's equation of motion (like a continuity equation), or is it a separate constraint imposed on top? A continuity equation (drho/dt + div(J) = 0 for some current J) is the standard mathematical mechanism and would be the natural choice, but the form of J depends on the field dynamics, which are not yet specified.

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
| 1 | Single continuous 3+1D field | Postulated; field type TBD |
| 2 | Well-defined local density | Postulated; units/bounds TBD |
| 3 | Total density conserved | Postulated; enforcement mechanism TBD |
| 4 | No intrinsic dissipation | Postulated; emergent irreversibility TBD |
| 5 | Finite max propagation speed c | Postulated; invariance group TBD |
| 6 | Stable localized configurations (knots) | Postulated; existence proof TBD |

### Optional Methodology: Minimal vs. Extended Axiom Sets

Two paths forward:

1. **Minimal**: Accept only Axioms 1-5 and attempt to DERIVE Axiom 6 (knot existence) from the field equation. This is more powerful but requires finding a specific equation that naturally produces soliton-like solutions. Precedent exists: sine-Gordon equation, Skyrme model, etc.

2. **Extended**: Accept all six axioms and focus on deriving consequences. Axiom 6 becomes a constraint on admissible field equations rather than a prediction.

Both are valid approaches. The minimal path is scientifically stronger but mathematically harder. The extended path lets conceptual development proceed while the existence question is resolved in parallel.
