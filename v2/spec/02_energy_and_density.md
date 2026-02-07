# 02 — Energy and Density

This chapter defines how energy relates to the field density introduced in the axioms, establishes the vocabulary for describing field configurations, and identifies the conceptual role of density gradients.

Depends on: [01_field_axioms.md](01_field_axioms.md)

---

## Density as the Primitive

The field has one primitive quantity: local density rho(x,t). Every physical observable must reduce to some functional of rho and its derivatives.

### Background Density

In the absence of any localized structure (no knots, no waves), the field is uniform:

    rho(x,t) = rho_0 = const (everywhere, always)

This uniform state is the **vacuum** or **ground state** of the field. rho_0 is the background density. It represents the "empty space" of CHPT — not truly empty, but uniformly filled.

### Knot Density and Depletion

When a knot forms, it concentrates density locally:

    rho_knot(x) > rho_0    (inside the knot)
    rho_surround(x) < rho_0  (in the vicinity of the knot)

Because total density is conserved (Axiom 3), the excess inside the knot must be exactly compensated by a deficit in the surrounding field. The knot "borrows" density from its environment.

This depletion zone is critical. It is the origin of gravity (see [07_gravity.md](07_gravity.md)) and of the long-range influence of massive objects.

---

## Energy as Non-Uniformity

CHPT defines energy operationally as **deviation from the uniform ground state**. A perfectly uniform field has zero energy in the dynamical sense — there is nothing happening, no structure, no motion.

Energy exists whenever the field is non-uniform. Specifically:

- **Gradient energy**: Energy stored in spatial density gradients. Where rho changes rapidly across space, there is energy localized at that gradient. This is analogous to elastic potential energy in a stretched medium.
- **Kinetic energy**: Energy associated with temporal changes in density (drho/dt != 0). Propagating disturbances carry kinetic energy.
- **Topological energy**: Energy locked into the structure of knots, which cannot be released without destroying the knot. This is the mass-energy of particles.

### The Energy Functional (Proposed)

The Lagrangian ([math/03_dynamics.md](math/03_dynamics.md)) gives the static Hamiltonian density with three terms:

$$ \mathcal{H} = \frac{1}{2}\langle\nabla\Psi\,\widetilde{\nabla\Psi}\rangle_0 + \frac{1}{4e^2}\sum_{i<j}\langle[R_i,R_j]^2\rangle_0 + \frac{\lambda}{4}(\langle\Psi\tilde{\Psi}\rangle_0 - \rho_0^2)^2 $$

The total static energy is:

$$ E = \int \mathcal{H} \, d^3x = E_2 + E_4 + E_V $$

This satisfies the required properties:

1. $E[\Psi = \rho_0] = 0$ (uniform field has zero dynamical energy).
2. $E \geq 0$ for all configurations (all three terms are non-negative).
3. $E$ is conserved via Noether's theorem (time-translation symmetry of the Lagrangian).

The gradient term ($E_2$) encodes twist/winding energy. The Skyrme term ($E_4$) provides soliton stability. The potential ($E_V$) penalizes deviation of $\langle\Psi\tilde{\Psi}\rangle_0$ from $\rho_0^2$.

### Update — Full Hamiltonian (A7/A8 Resolved)

The issues noted above have been resolved. The full Hamiltonian density ([math/03_dynamics.md](math/03_dynamics.md)) now includes **three terms**:

$$ \mathcal{H} = \underbrace{\frac{1}{2c^2}\langle\partial_t\Psi\,\widetilde{\partial_t\Psi}\rangle_0 + \frac{1}{2}\langle\nabla\Psi\,\widetilde{\nabla\Psi}\rangle_0}_{E_2} + \underbrace{\frac{1}{4e^2}\sum_{i<j}\langle[R_i,R_j]^2\rangle_0}_{E_4} + \underbrace{\frac{\lambda}{4}(\langle\Psi\tilde{\Psi}\rangle_0 - \rho_0^2)^2}_{E_V} $$

- **Time kinetic energy**: Now included via the explicit $\partial_t\Psi$ term (A7 resolved).
- **Skyrme term** ($E_4$): Provides soliton stability against collapse in 3D (A8 resolved). See [math/05_mass_mechanism.md](math/05_mass_mechanism.md) for the Derrick analysis.
- **Degenerate mass term** ($E_D$): Gives mass $\mu$ to the pseudoscalar and flux modes (B5 resolved). See [math/03_dynamics.md](math/03_dynamics.md), §5.

---

## Mass as Field Energy

In CHPT, mass is **Energy**. Specifically, the inertial mass of a particle is the total integrated energy of its knot configuration.

$$ m = \frac{1}{c^2} \int T_{00} \, d^3x $$

This energy comes from three sources (as seen in the Lagrangian):
1.  **Gradient Energy** ($E_2$): The energy stored in the twisting and winding of the field (spatial gradients).
2.  **Skyrme Energy** ($E_4$): Higher-order stabilization energy from interacting field currents. Essential for preventing soliton collapse.
3.  **Potential Energy** ($E_V$): The energy cost of displacing the field norm from the vacuum value $\rho_0$.

The virial theorem constrains these: $E_2 = E_4 - 3E_V$, giving the mass formula $Mc^2 = 2E_4 - 2E_V$ (see [math/05_mass_mechanism.md](math/05_mass_mechanism.md)).

Heavier particles (like protons) correspond to more complex knot topologies (higher winding numbers, more twist) which require more total field energy to sustain.

### Consequences
-   **$E=mc^2$**: This is definitional. The particle *is* a localized bundle of field energy. The conversion factor $c^2$ is set by the field's propagation speed.
-   **Inertia**: The resistance to acceleration comes from the fact that moving the knot requires propagating its internal process. The "heavier" the knot (more internal energy), the more field activity must be updated to move it.
-   **Equivalence Principle**: Gravitational mass equals inertial mass because both are the same quantity — the total field energy of the knot. The density excess that creates the surrounding depletion zone (gravitational effect) is the same density excess that resists acceleration (inertial effect). This is the equivalence principle, arising as an identity rather than a coincidence.
-   **No negative mass**: Density is non-negative, and a knot below background would be a deficit (depletion zone), not a particle.

### Unknowns

- **Quantitative mass spectrum**: Why is the electron mass exactly 0.511 MeV and the proton 938 MeV? CHPT predicts that these follow from the stable knot harmonics, but has not yet identified which harmonics correspond to which particles, nor computed any masses.
- **Mass vs. binding energy**: In standard physics, the proton mass is mostly binding energy (QCD), not quark masses. Does CHPT reproduce this? The knot model of the proton would need internal structure whose binding contributes most of the mass.

---

## Density Gradients and Forces

Wherever the density field is non-uniform, there are gradients. These gradients are the fundamental source of all forces in CHPT:

- **Gravity**: Long-range gradients from knot depletion zones. See [07_gravity.md](07_gravity.md).
- **Electromagnetic force**: Mediated by null-rotor patterns in the density field. See [08_electromagnetism.md](08_electromagnetism.md).
- **Nuclear forces**: Short-range intense gradients from overlapping knot structures. See [09_nuclear_interactions.md](09_nuclear_interactions.md).

The unifying principle: **there is only one kind of "stuff" (density), and all forces are density gradients acting on density structures.** The apparent diversity of forces comes from the different geometric arrangements of density at different scales.

### Critical Question

This is an elegant unification claim, but it must be demonstrated, not just asserted. The four known forces have qualitatively different properties:

| Property | Gravity | EM | Strong | Weak |
|----------|---------|-----|--------|------|
| Range | Infinite | Infinite | ~1 fm | ~0.001 fm |
| Strength (relative) | 1 | 10^36 | 10^38 | 10^25 |
| Charges | Mass | Electric | Color | Weak isospin |
| Mediators | Graviton? | Photon | Gluon | W/Z |

CHPT must explain how a single density field produces force behaviors spanning 38 orders of magnitude in strength, with both infinite-range and confined behaviors, mediated by different carrier structures. This is the central quantitative challenge.

---

## Vacuum Energy and the Zero-Point Problem

In the uniform field state (rho = rho_0 everywhere), the dynamical energy is zero. But the field still exists — it has density rho_0 everywhere. Is this a form of energy?

In standard physics, vacuum energy density is a notorious problem: quantum field theory predicts an enormous value (~10^113 J/m^3) while cosmological observations suggest a tiny value (~10^-9 J/m^3). This is the "cosmological constant problem," a 122 orders-of-magnitude discrepancy.

CHPT sidesteps this: the uniform background is the zero of energy by definition. rho_0 is not "energy" in the dynamical sense because it produces no gradients, no forces, no observable effects by itself. Only deviations from rho_0 are physically meaningful.

### Caveat

This works conceptually, but it raises the question: is rho_0 observable? If two regions of space had different background densities, would that be detectable? If not, then rho_0 is a gauge freedom and can be set to any value. If yes, it becomes a physical parameter whose value needs explaining.

---

## Summary

| Concept | CHPT Definition | Status |
|---------|----------------|--------|
| Energy | Deviation from uniform density | Functional proposed ([math/03_dynamics.md](math/03_dynamics.md)) |
| Mass | $Mc^2 = \int \mathcal{H} \, d^3x$ | Proposed ([math/05_mass_mechanism.md](math/05_mass_mechanism.md)); spectrum TBD |
| Vacuum | Uniform density rho_0 | Defined; observability TBD |
| Force | Density gradient acting on density structure | Conceptual; quantitative derivation TBD |
| Vacuum energy problem | Dissolved (rho_0 not dynamical energy) | Elegant if consistent |
