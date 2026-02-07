# 03 — Propagation

This chapter describes how disturbances move through the field, establishes the speed limit, and explains why CHPT is not an aether theory.

Depends on: [01_field_axioms.md](01_field_axioms.md), [02_energy_and_density.md](02_energy_and_density.md)

---

## Disturbances in the Field

Any localized change in density will propagate outward. This is a direct consequence of the field being continuous and having dynamics — perturb it somewhere, and the perturbation spreads. The details of how it spreads depend on the (unspecified) field equation, but the general behavior is:

1. **Linear waves**: Small perturbations around the vacuum $\rho_0$ propagate as waves with speed c. The perturbation $\psi = S + \mathbf{F} + IP$ has scalar ($S$), bivector ($\mathbf{F}$), and pseudoscalar ($P$) components (see [math/04_electromagnetism.md](math/04_electromagnetism.md)). The bivector component $\mathbf{F}$ is identified with electromagnetic radiation (photons / null-rotors, see [06_null_rotors.md](06_null_rotors.md)). The scalar and pseudoscalar modes are additional propagating degrees of freedom whose observational status is an open problem (see [15_open_problems.md](15_open_problems.md), B5).

2. **Nonlinear structures**: Large perturbations (knots) are localized and do not simply disperse. They propagate as coherent structures — solitons or topological defects — that maintain their shape. Their propagation speed can be anywhere from 0 to c (exclusive).

3. **Radiation**: Accelerating knots and transitioning knot states emit propagating density disturbances. This radiation always moves at c, regardless of the emitter's motion.

---

## The Speed Limit c

Axiom 5 states that c is the maximum propagation speed. What does this mean physically?

### Interpretation

The field has a characteristic "stiffness" — how density changes at one point influence neighboring points. This stiffness, combined with the field's effective inertia (resistance to density change), determines the propagation speed, just as tension and mass density determine the speed of a wave on a string:

    v_wave = sqrt(restoring_force / inertia)

For the CHPT field, c = sqrt(stiffness / inertia) where both "stiffness" and "inertia" are properties of the field equation. c is the speed of infinitesimal perturbations around the uniform background.

### Why c is a Maximum

Nonlinear structures (knots) always move slower than c. Intuition: a knot is a large-amplitude density feature that interacts with its own field. This self-interaction provides effective mass/inertia that increases with speed, requiring ever more energy to accelerate further. In the limit v -> c, the required energy diverges. This is exactly the relativistic mass-energy relationship, but here it emerges from field dynamics rather than spacetime geometry.

### Unknowns

- **Dispersion**: Do waves of different frequencies travel at exactly the same speed c? In most media, wave speed depends on frequency (dispersion). If the CHPT field is dispersive, then "the speed of light" would depend on photon energy — which is very tightly constrained by observations (gamma-ray bursts show no energy-dependent speed to ~10^-20 relative precision). The field equation must either be non-dispersive or have dispersion below this threshold.
- **Does c vary with background density?**: If c depends on rho_0, then regions with different background densities (e.g., near massive objects) would have different light speeds. This could mimic gravitational lensing and time dilation, which is interesting. But it could also produce unwanted effects. See [11_relativity.md](11_relativity.md).

---

## Why CHPT Is Not an Aether Theory

CHPT posits a field that fills all of space. This superficially resembles the "luminiferous aether," but the resemblance is false. CHPT is a **Process Theory**, not a material medium theory.

### The Aether Fallacy: Object vs. Medium
Classical aether theory assumed two things:
1.  **The Medium** (Aether): A background substance.
2.  **The Object** (Matter): A separate thing moving *through* the medium, creating drag.

### The CHPT Reality: Process Ontology
In CHPT, there is no "object" distinct from the "medium."
*   **Knotting and Unknotting**: An electron is not a ball swimming through fluid. It is a standing wave of **knotting** (field organization).
*   **Motion**: Motion is the propagation of this organizational state. As a knot moves, the space ahead "knots" and the space behind "unknots."
*   **No Drag**: There is no friction because there is no "thing" rubbing against the field. The field IS the thing. The energy of motion is simply the energy required to drive this reconfiguring process.

This is why the Michelson-Morley null result is automatic. You cannot measure your speed relative to the medium when you *are* the medium's activity.

### Mathematical Consequence: Lorentz Covariance
This physical insight requires that the field equation be **Lorentz Covariant**. The laws of "knotting/unknotting" must look the same in any inertial frame, ensuring that no experiment performed by a knot (an observer) can detect the state of the background field.

---

## Three Modes of Propagation

It is crucial to distinguish three distinct "speeds" in CHPT:

1.  **Wave Propagation (Speed = c)**: Small perturbations (null-rotors, radiation) propagate at the characteristic speed of the field's differential structure. This is the "process speed limit" for causal signaling.
2.  **Knot Propagation (Speed < c)**: Stable process structures (particles) propagate at $v < c$. They possess inertia because the "knotting/unknotting" cycle requires energy to sustain and transport. As $v \to c$, the energy required to update the configuration diverges.
3.  **Topological Updates (Instantaneous)**: Global topological constraints (e.g., the connectedness of a flux tube) apply to the entire field simultaneously. If a knot is topologically connected to another distant knot, this connection is a non-local state. Changes to the topology (like a string breaking) are global events that do not "propagate" at $c$ but simply *are*. This allows for quantum entanglement without violating signal locality.

---

## Propagation of Knots vs. Propagation of Waves

A crucial distinction:

- **Waves** (small perturbations, radiation, null-rotors) propagate AT c. They are massless. They do not maintain a fixed spatial extent; they spread as they propagate (unless confined by boundary conditions or nonlinear focusing).

- **Knots** (particles) propagate at v < c. They are massive. They maintain a fixed spatial structure (modulo Lorentz contraction). They carry localized density excess.

The distinction between wave-like and knot-like propagation is not sharp at all scales — a knot in motion continuously emits and absorbs radiation (see [06_null_rotors.md](06_null_rotors.md)), creating a dynamic equilibrium between the localized and extended parts of the field configuration. This is the CHPT version of the "dressed particle" concept in quantum field theory.

---

## Group Velocity, Phase Velocity, and Information

For waves in the field:

- **Phase velocity**: The speed at which wave crests move. Can in principle exceed c without violating causality.
- **Group velocity**: The speed at which a wave packet's envelope moves. This is bounded by c and represents the speed of energy/information transfer.

The distinction matters if the field is dispersive (different frequencies, different speeds). For a non-dispersive field, phase and group velocity are equal to c.

### Unknown

- **Nonlinear dispersion**: Even if the linear wave equation is non-dispersive, nonlinear effects (which must exist to support knots) could modify the effective dispersion for large-amplitude waves. Whether this affects observable propagation speeds is TBD.

---

## Summary

| Concept | CHPT Status | Critical Issue |
|---------|------------|----------------|
| Wave propagation at c | Required by axiom | Dispersion must be negligible |
| Knot propagation at v < c | Required for massive particles | Soliton existence TBD |
| Lorentz invariance | Required by the physics (no object-medium duality) | Must be verified in specific equation |
| Aether problem | Does not arise (knots ARE the field) | Not an aether theory |
| Preferred frame | None (no object-medium distinction) | Field equation must formalize this |
