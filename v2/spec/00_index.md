# Chiral Harmonic Particle Theory (CHPT) — Specification

## Purpose

This specification documents the conceptual foundations of CHPT: a proposal to derive all of physics from a single continuous density field in 3+1 dimensional spacetime. Particles are stable field configurations ("knots"), forces are density gradients, and radiation consists of oscillating field patterns ("null-rotors").

The goal of this document set is to be as detailed, honest, and precise as possible about what the theory proposes, what it can explain, what it cannot, and what mathematical work is required to turn the proposal into a testable physical theory.

## One-Sentence Summary

All physics is configurations of a single frictionless density field: particles are knots, forces are gradients, light is oscillation.

---

## Document Map

### Foundations (Axioms and Primitives)

| File | Title | Summary |
|------|-------|---------|
| [01_field_axioms.md](01_field_axioms.md) | Field Axioms | The six foundational postulates: a continuous 3+1D field with conserved density, no dissipation, a speed limit c, and stable localized configurations. |
| [02_energy_and_density.md](02_energy_and_density.md) | Energy and Density | Energy as deviation from uniform density. Mass as integrated density excess. Vacuum as the uniform ground state. |
| [03_propagation.md](03_propagation.md) | Propagation | Wave and knot propagation, the speed limit, why CHPT is not an aether theory (no object-medium duality), and why Lorentz covariance is a natural consequence. |

### Particle Structure

| File | Title | Summary |
|------|-------|---------|
| [04_knots_and_particles.md](04_knots_and_particles.md) | Knots and Particles | Particles as stable localized field patterns. Stability mechanisms (topological vs. dynamic). Mode structure and the mapping to known particles. |
| [05_chirality.md](05_chirality.md) | Chirality | Handedness of knots as the origin of electric charge. Interaction rules: same-chirality repulsion (Pauli exclusion), opposite-chirality nesting (binding), achiral transparency. |
| [06_null_rotors.md](06_null_rotors.md) | Null-Rotors | Massless propagating field patterns as photon analogs. Dual-state oscillation (harmonic/wave), emission scaling, PGA connection, force mediation. |

### Forces and Interactions

| File | Title | Summary |
|------|-------|---------|
| [07_gravity.md](07_gravity.md) | Gravity | Back-pressure from density depletion zones. 1/r^2 from geometry. Comparison with GR. Le Sage problem resolution. Dark matter discussion. |
| [08_electromagnetism.md](08_electromagnetism.md) | Electromagnetism | Coulomb's law from null-rotor exchange. Magnetic fields from moving charges. Perpendicular field structure. Path to deriving Maxwell's equations. |
| [09_nuclear_interactions.md](09_nuclear_interactions.md) | Nuclear Interactions | Strong force from direct knot overlap. Color charge problem (SU(3) from binary chirality?). Weak force options. Fusion and fission. |

### Emergent Physics

| File | Title | Summary |
|------|-------|---------|
| [10_conservation_laws.md](10_conservation_laws.md) | Conservation Laws | How energy, momentum, angular momentum, and charge conservation emerge from field symmetries and topological invariants. |
| [11_relativity.md](11_relativity.md) | Relativity | Special relativity automatic from Lorentz-invariant field equation. General relativity via effective metric from density gradients. Post-Newtonian tests. |
| [12_quantum_phenomena.md](12_quantum_phenomena.md) | Quantum Phenomena | The Bell's theorem crisis. Nonlocal hidden variables as the recommended escape. Wave-particle duality, uncertainty, entanglement. |

### Phenomenology

| File | Title | Summary |
|------|-------|---------|
| [13_particle_spectra.md](13_particle_spectra.md) | Particle Spectra | The Standard Model particle zoo vs. CHPT predictions. Three generations, mass hierarchy, stability and lifetimes. Current status: zero computed values. |
| [14_cosmology.md](14_cosmology.md) | Cosmology | Big Bang, expansion, CMB, dark matter, dark energy, black holes, baryogenesis. Mostly re-narration of standard cosmology. |

### Assessment

| File | Title | Summary |
|------|-------|---------|
| [15_open_problems.md](15_open_problems.md) | Open Problems | Categorized list of all unresolved issues (fatal, major, minor). Decision register. Prioritized research program. Honest final assessment. |

---

## Key Decisions Recommended in This Spec

These are positions **recommended** (not yet confirmed) throughout the documents where the original proposal was ambiguous or underspecified. Each represents a choice that the mathematical formalization must validate or revise:

1. **Lorentz invariance**: The field equation must be Lorentz-covariant. This is not an imposed constraint but a physical consequence of CHPT's ontology: because knots ARE the field (not objects moving through a medium), there is no object-medium duality, no preferred frame, and no aether problem. The math must formalize this. (See [03_propagation.md](03_propagation.md))

2. **Quantum mechanics**: CHPT must adopt nonlocal hidden variables to evade Bell's theorem. This is a **significant departure from the original proposal's deterministic local-realism claim**. Axiom 5 would need to distinguish signal locality (preserved) from ontological nonlocality (allowed). The implications of this change propagate through the entire theory and are not yet fully explored. (See [12_quantum_phenomena.md](12_quantum_phenomena.md))

3. **Black holes**: The original proposal claims "no event horizons." This spec recommends **not committing** to that claim until the effective metric is computed. The density gradient mechanism may well produce horizons. If it does, CHPT's unique contribution would be horizons without singularities (finite maximum density replacing the infinite-density point). This depends on whether the field has a maximum density — an open question. (See [14_cosmology.md](14_cosmology.md))

4. **Dark matter**: Modified gravity alone is insufficient (the Bullet Cluster and CMB power spectrum argue strongly for a particle dark matter component). CHPT should predict dark knots (achiral, non-luminous stable configurations) but **must acknowledge this is a reframing, not a solution**, until specific knot properties are derived. (See [07_gravity.md](07_gravity.md))

5. **Field type**: Recommend geometric algebra Cl(3,0,1) to naturally support null-rotor structure, chirality, and tensor gravitational waves. **This is speculative** — the choice is motivated by post-hoc fit to desired features, not derived from first principles. Other algebras (Cl(4,1), Cl(1,3), etc.) have not been evaluated. (See [01_field_axioms.md](01_field_axioms.md), [06_null_rotors.md](06_null_rotors.md))

### Warning: Unresolved Contradictions

This spec identifies several internal tensions that cannot be resolved at the conceptual level:

- **Determinism vs. Bell's theorem**: The original proposal is deterministic and local. Bell's theorem forbids this. The recommended nonlocal hidden variables approach is viable but fundamentally changes the theory's character.
- **Topological conservation vs. baryogenesis**: If baryon number is topologically exact, the matter-antimatter asymmetry cannot be generated via Sakharov conditions. Either the initial conditions contained the asymmetry, or topological conservation is approximate.
- **Gauge symmetry**: The Standard Model's gauge structure (U(1), SU(2), SU(3)) is completely absent from CHPT. How a single density field produces non-Abelian gauge symmetries is unknown and may be impossible.
- **Neutrino interactions**: If neutrinos are achiral (neutral), they should not interact via the weak force. But they do. The achiral-equals-neutral mapping may be wrong or incomplete.

These are not aesthetic complaints. They are points where CHPT's concepts, as currently formulated, conflict with established physics. Resolution requires either modifying the concepts or demonstrating mathematically that the conflicts are only apparent.

---

## Critical Path to a Physical Theory

The spec identifies a phased research program (detailed in [15_open_problems.md](15_open_problems.md)):

1. **Write the field equation** — Choose field type, write Lagrangian, derive EOM.
2. **Validate the linear sector** — Derive Maxwell's equations from linearized dynamics.
3. **Find stable knots** — Solve for soliton/topological defect solutions, compute masses.
4. **Validate gravity** — Compute effective metric, PPN parameters.
5. **Validate quantum sector** — Formulate nonlocal guidance, derive Born rule.
6. **Address nuclear physics** — SU(3) color, weak force, nucleosynthesis.

The theory lives or dies at Step 1. Without an equation, nothing can be computed or tested. With the right equation, everything follows.

---

## Notation and Conventions

Throughout this spec:
- **rho** (ρ): field density
- **rho_0** (ρ₀): background/vacuum density
- **c**: maximum propagation speed (speed of light)
- **Knot**: stable localized field configuration (= particle)
- **Null-rotor**: propagating field oscillation (= photon/radiation)
- **Chirality**: handedness of knot pattern (= electric charge sign)
- **Depletion zone**: region of below-background density surrounding a knot (= gravitational field)
- **PGA**: Projective Geometric Algebra, specifically Cl(3,0,1)
- Units: unspecified until the field equation determines natural units.
