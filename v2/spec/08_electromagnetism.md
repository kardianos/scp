# 08 — Electromagnetism

This chapter describes how electromagnetic phenomena emerge from the interaction of chiral knots and null-rotors. It covers Coulomb's law, radiation, spectra, polarization, and Maxwell's equations.

Depends on: [05_chirality.md](05_chirality.md), [06_null_rotors.md](06_null_rotors.md), [07_gravity.md](07_gravity.md)

---

## Electric Charge as Topology

As established in [02_topology.md](math/02_topology.md) and refined in [05_chirality.md](05_chirality.md), electric charge maps to the **Topological Index (Winding Number)** of a knot, denoted $Q$:

-   $Q = +1$: Positive charge (e.g., Positron).
-   $Q = -1$: Negative charge (e.g., Electron).
-   $Q = 0$: Neutral (e.g., Neutrino, Photon).

### Distinction from Chirality
Properties like "handedness" (Left/Right) are distinct from Charge. A neutral particle ($Q=0$) can still be chiral (Left-handed Neutrino).

The electromagnetic force arises from the interaction between **Topologically Charged** knots, mediated by null-rotor exchange.


---

## Coulomb's Law

Two charged knots at rest, separated by distance r, experience a force:

    F = k * q1 * q2 / r^2

where q1, q2 are the charges (chiralities) and k is Coulomb's constant.

### CHPT Mechanism

Each chiral knot continuously emits null-rotors due to its internal oscillation (see [06_null_rotors.md](06_null_rotors.md)). The null-rotor emission pattern is anisotropic for chiral knots — it has the perpendicular structure described in the null-rotor chapter.

When two chiral knots are near each other:

- **Same chirality (like charges)**: The null-rotor patterns between them constructively interfere, creating a density enhancement in the gap. This enhancement pushes both knots outward. Repulsion.
- **Opposite chirality (unlike charges)**: The null-rotor patterns destructively interfere in the gap, creating a density reduction. The surrounding field pushes both knots inward. Attraction.

### Why 1/r^2

The null-rotor emission spreads spherically, diluting as 1/r^2 in 3D. The force, which depends on the null-rotor intensity at the location of the other knot, therefore also scales as 1/r^2.

This is the same geometric origin as gravity's 1/r^2. The difference:
- Gravity involves the static depletion zone (always attractive).
- Coulomb involves the dynamic null-rotor emission pattern (attractive or repulsive depending on chirality).

### Why EM Is Stronger Than Gravity

The null-rotor emission intensity from a knot's chiral oscillation is much larger than the density depletion from its mere presence. In magnitude: the EM force between two protons is ~10^36 times stronger than the gravitational force. In CHPT, this means the chiral oscillation produces null-rotor fluxes vastly exceeding the static depletion effect. Why?

Because the null-rotor emission is a direct, first-order effect of the knot's internal dynamics, while the depletion zone is a second-order, cumulative background effect. The ratio between them depends on the field equation's parameters.

### Unknowns

- **Derivation of Coulomb's constant**: k = 1/(4 pi epsilon_0). What combination of CHPT parameters (rho_0, c, nonlinear coupling, etc.) equals k?
- **Charge quantization**: Why is the electron charge the fundamental unit? In CHPT, this must follow from the minimum chiral winding number of stable knots.
- **Screening and running coupling**: In QED, the effective electric charge changes with distance (vacuum polarization screens the bare charge). Does CHPT produce this? The field between two knots contains null-rotor fluctuations that could screen/enhance the effective chirality.

---

## Magnetic Fields and Moving Charges

When a chiral knot moves, its null-rotor emission pattern is modified. The perpendicular emission structure (described in [06_null_rotors.md](06_null_rotors.md)) becomes dynamic:

- The "electric" component of the null-rotor emission is radial (points away from the charge).
- The "magnetic" component wraps around the direction of motion (azimuthal).

This matches the Biot-Savart law: a moving charge produces a magnetic field that circles around the current direction.

### Lorentz Force

A chiral knot moving through a region with existing null-rotor patterns (another knot's field) experiences a force that depends on:
- Its charge (chirality).
- Its velocity.
- The null-rotor pattern at its location.

The force has two components:
- Along the electric component of the null-rotor: F = qE (independent of velocity).
- Perpendicular to the magnetic component: F = qv x B (proportional to velocity, perpendicular to both v and B).

This is the Lorentz force law: F = q(E + v x B).

### Unknown

- **Derivation**: The Lorentz force law must be derived from the field dynamics, not assumed. Specifically, it must be shown that a chiral knot moving through a null-rotor pattern experiences exactly qE + qv x B.

---

## Maxwell's Equations

If CHPT is correct, the dynamics of null-rotor patterns must reduce to Maxwell's equations in the appropriate limit. Maxwell's four equations:

1. div E = rho/epsilon_0 (Gauss's law — charges source electric fields)
2. div B = 0 (no magnetic monopoles)
3. curl E = -dB/dt (Faraday's law — changing B induces E)
4. curl B = mu_0 J + mu_0 epsilon_0 dE/dt (Ampere-Maxwell law)

### CHPT Mapping

- **Gauss's law**: Chiral knots emit null-rotors radially. The total flux through a closed surface equals the enclosed chirality (charge). This is a geometric consequence of the emission being spherically symmetric for a static knot.
- **No magnetic monopoles**: Magnetic fields arise from motion (not from a separate charge type). Since magnetic components are always secondary effects of moving chiral knots, there is no source term for div B. Note: CHPT does not obviously PROHIBIT magnetic monopoles. If a knot configuration existed whose null-rotor emission had radial magnetic (rather than electric) character, it would be a magnetic monopole. Whether such configurations are stable is a question about the knot spectrum.
- **Faraday's law**: A time-varying null-rotor "wave phase" (B analog) naturally generates the "harmonic phase" (E analog) component, because the two are coupled by the oscillation mechanism. This is the null-rotor's defining dual-state oscillation.
- **Ampere-Maxwell**: Moving chirals (currents) and time-varying E fields generate B fields. This follows from the same oscillation coupling.

### Critical Unknown

**Can Maxwell's equations be rigorously derived from the CHPT field equation?** This is the litmus test for the electromagnetic sector of the theory. If the null-rotor dynamics do not reduce to Maxwell's equations in the linear (weak-field) limit, the theory fails at the most basic level.

This derivation is one of the highest-priority tasks for mathematical formalization. A successful derivation would:
1. Identify E and B with specific components of the null-rotor field.
2. Derive the four Maxwell equations from the field dynamics.
3. Extract epsilon_0 and mu_0 (and therefore c = 1/sqrt(epsilon_0 mu_0)) from the field parameters.

---

## Atomic Spectra

Atoms are composite knots: an achiral or chiral nuclear knot surrounded by lighter chiral knots (electrons). The electrons occupy specific resonant configurations around the nucleus (orbitals).

When an electron knot transitions from one orbital to another, the change in internal oscillation frequency causes it to emit or absorb a null-rotor at the corresponding frequency difference:

    f_photon = (E_upper - E_lower) / h

This produces the discrete spectral lines observed for each element.

### CHPT Contribution

The orbital structure of electrons around nuclei is already well-described by quantum mechanics. CHPT does not obviously improve on this unless it can:

1. Derive the orbital quantization from field resonance conditions.
2. Predict transition rates from first principles.
3. Explain the fine structure and hyperfine structure from knot geometry.

Without these, CHPT's account of atomic spectra is no more than a re-narration of quantum mechanics in new language.

---

## Polarization

Light (null-rotors) is polarized: the oscillation direction of the electric component can be aligned in specific orientations.

In CHPT, polarization arises naturally from the perpendicular emission structure:

- A chiral knot oscillating along a specific axis emits null-rotors polarized along that axis.
- The chiral selectivity of knots means they preferentially absorb null-rotors whose polarization matches their internal oscillation axis.

### Circular Polarization and Chirality

Circularly polarized light (where the E vector rotates as the wave propagates) maps naturally to the chiral structure of null-rotors. Left-circularly-polarized light may correspond to null-rotors with left-chiral internal oscillation, and similarly for right. This is a natural prediction: photon helicity maps to null-rotor chirality.

---

## Summary

| EM Phenomenon | CHPT Mechanism | Derived? |
|---------------|---------------|----------|
| Coulomb 1/r^2 | Null-rotor flux dilution | Conceptual only |
| Like-charge repulsion | Constructive null-rotor interference | Conceptual only |
| Unlike-charge attraction | Destructive null-rotor interference | Conceptual only |
| Biot-Savart / magnetism | Modified emission pattern from moving chiral | Conceptual only |
| Lorentz force | Knot response to null-rotor field | Not derived |
| Maxwell's equations | Null-rotor dynamics linearized | Not derived |
| Atomic spectra | Orbital resonance transitions | No improvement over QM |
| Polarization | Oscillation axis of null-rotor | Natural fit |

### Priority for Formalization

The single most important derivation in all of CHPT is: **derive Maxwell's equations from the field dynamics.** If this succeeds, the electromagnetic sector is validated. If it fails, the theory must be fundamentally revised.
