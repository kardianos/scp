# 10 — Conservation Laws

This chapter examines how the fundamental conservation laws of physics emerge from CHPT's axioms and dynamics. Conservation laws are among the most rigorously tested facts in physics; any viable theory must reproduce all of them exactly.

Depends on: [01_field_axioms.md](01_field_axioms.md), [02_energy_and_density.md](02_energy_and_density.md), [04_knots_and_particles.md](04_knots_and_particles.md), [05_chirality.md](05_chirality.md)

---

## Noether's Theorem and Its Requirements

Emmy Noether's theorem (1918) states: every continuous symmetry of a physical system's action corresponds to a conserved quantity. The canonical examples:

| Symmetry | Conserved Quantity |
|----------|-------------------|
| Time translation invariance | Energy |
| Space translation invariance | Momentum |
| Rotational invariance | Angular momentum |
| U(1) phase symmetry | Electric charge |
| SU(3) color symmetry | Color charge |

For CHPT to derive conservation laws via Noether's theorem, two things are needed:

1. **A Lagrangian/action formulation**: The field dynamics must be expressible as an action principle (extremize the integral of a Lagrangian density). This is the standard framework for all fundamental physics and is required for Noether's theorem to apply.
2. **The correct symmetries**: The Lagrangian must have the symmetries listed above.

### Status in CHPT

CHPT does not yet have a Lagrangian. It has axioms and qualitative dynamics. Translating these into a Lagrangian is a prerequisite for rigorous conservation law derivation. However, we can assess whether the axioms are COMPATIBLE with the required symmetries.

---

## Energy Conservation

### Axiom-Level Guarantee

Axiom 3 (total density conserved) directly guarantees that total field density is conserved. If energy maps to field non-uniformity (as proposed in [02_energy_and_density.md](02_energy_and_density.md)), then energy conservation follows IF the background density rho_0 is constant.

More precisely: in a closed system with no boundary, total energy is conserved if the field dynamics are time-translation invariant (the laws don't change with time). This is expected — the field equation should have fixed coefficients.

### Subtlety: Expanding Universe

In standard cosmology, energy conservation is subtle in an expanding universe. Photons redshift (lose energy) as space expands. Is this energy lost? In GR, energy conservation applies locally but not globally in a curved/expanding spacetime.

In CHPT (flat spacetime + field): if the field is expanding/diluting (see [14_cosmology.md](14_cosmology.md)), the total density is still conserved. Redshifted null-rotors have lower frequency but the energy goes into the expansion itself (kinetic energy of the dilution). Total energy is exactly conserved globally. This is arguably cleaner than the GR treatment.

---

## Momentum Conservation

### From Spatial Translation Invariance

If the field equation is the same at all points in space (no preferred location), then Noether's theorem gives conservation of momentum.

CHPT's axioms are consistent with this: the field is homogeneous in its properties (same equation everywhere). Knots and null-rotors can be at any location. No point is special.

### How Knots Carry Momentum

A stationary knot has zero momentum. A moving knot has momentum proportional to its mass (density excess) and velocity:

    p = m * v * gamma(v)    (where gamma is the Lorentz factor, if Lorentz invariance holds)

The momentum is stored in the density current — the flow of density excess through the field. When two knots collide, the total density current is conserved.

Null-rotors carry momentum p = E/c (standard photon momentum). This follows from the field dynamics if they are formulated correctly.

---

## Angular Momentum Conservation

### From Rotational Invariance

If the field equation is the same for all orientations (no preferred direction), Noether's theorem gives conservation of angular momentum.

CHPT is consistent with this: the axioms make no reference to preferred spatial directions.

### Orbital and Intrinsic (Spin) Angular Momentum

- **Orbital**: A knot orbiting another knot carries orbital angular momentum L = r x p. Standard.
- **Spin**: A knot's internal oscillation pattern carries intrinsic angular momentum. This is the knot's spin (see [05_chirality.md](05_chirality.md)). Spin is quantized because only specific internal harmonic patterns are stable.

The total angular momentum (orbital + spin) is conserved.

### Unknown

- **Spin values**: Why are fermions spin-1/2 and bosons spin-0 or spin-1? The stable knot harmonics must naturally select for these specific spin values. This has not been demonstrated.

---

## Electric Charge Conservation

### From Chirality Conservation

If chirality (the handedness of a knot) is a topological invariant — it cannot be changed by continuous deformation of the field — then chirality is exactly conserved.

In any process:
- The number of left-chiral knots minus the number of right-chiral knots is constant.
- This is equivalent to: total electric charge is conserved.

### Topological vs. Dynamic Conservation

If chirality is topologically conserved (knots cannot change handedness without being destroyed and recreated), then charge conservation is exact — as exact as the impossibility of untying a mathematical knot by continuous deformation.

If chirality is only dynamically conserved (a conservation law from symmetry, which could in principle be violated by quantum tunneling or other processes), then charge conservation is approximate. In standard physics, charge conservation is believed to be exact (no proton decay to neutral particles has ever been observed, lifetime > 10^34 years). CHPT should aim for topological conservation.

### Unknown

- **Is chirality topological or dynamic?** This depends on the field type and equation. For scalar fields, there are no topological chiral invariants in 3D. For fields valued in suitable groups (SU(2), SO(3), etc.), chirality CAN be topological. This constrains the field type.

---

## Baryon Number Conservation

In standard physics, baryon number (number of baryons minus number of antibaryons) is conserved in all known processes. Proton decay has never been observed.

In CHPT: if baryons (proton, neutron) are topological knots with a conserved topological charge (winding number), then baryon number conservation is exact.

### Connection to Charge Conservation

In the Skyrme model (a close relative of CHPT's knot concept), the topological charge IS baryon number. Electric charge is a derived quantity related to the orientation of the Skyrme field. If CHPT follows this structure, baryon number is the fundamental topological invariant, and electric charge is a consequence.

---

## Lepton Number Conservation

Lepton number (number of leptons minus number of antileptons) is conserved in the Standard Model (at the classical level; anomalies can violate it). Neutrino oscillations conserve total lepton number but violate individual flavor lepton numbers (L_e, L_mu, L_tau).

In CHPT: if leptons (electron, muon, tau, neutrinos) have their own topological invariant distinct from baryon number, then lepton number is conserved. Whether this exists depends on the knot taxonomy.

### Unknown

- **Baryon-lepton distinction**: What topological feature distinguishes baryonic knots from leptonic knots? In the Standard Model, this distinction is built in (quarks carry baryon number, leptons carry lepton number). In CHPT, it must emerge from the knot structure. A natural candidate: baryons are composite knots (three sub-knots), leptons are elementary knots. The composite/elementary distinction could define a topological invariant.

---

## CPT Symmetry

The CPT theorem states that any Lorentz-invariant local quantum field theory is invariant under the combined operation of Charge conjugation (C), Parity (P), and Time reversal (T).

In CHPT:
- **C (charge conjugation)**: Swap left-chiral and right-chiral knots. If the field equation is chiral-symmetric, C is a symmetry.
- **P (parity)**: Mirror reflection. If the field equation is parity-symmetric, P is a symmetry.
- **T (time reversal)**: Run the dynamics backward. If the field equation is time-reversal symmetric (no dissipation — Axiom 4), T is a symmetry.

Individual symmetries (C, P, T) may be broken (the weak force breaks P and CP), but the combination CPT should be preserved if the field equation is Lorentz-invariant (Option A from [03_propagation.md](03_propagation.md)).

### CP Violation

The Standard Model violates CP symmetry (observed in kaon and B-meson systems). This is essential for baryogenesis (explaining why there is more matter than antimatter). CHPT must accommodate CP violation, either:

- In the field equation itself (explicit CP violation).
- In the vacuum/ground state (spontaneous CP violation).
- In the knot dynamics (complex phases in transition amplitudes).

---

## Summary

| Conservation Law | CHPT Mechanism | Status |
|-----------------|---------------|--------|
| Energy | Total density conserved (Axiom 3) | Built in |
| Momentum | Spatial translation invariance | Consistent with axioms |
| Angular momentum | Rotational invariance | Consistent with axioms |
| Electric charge | Chirality invariant (topological?) | Depends on field type |
| Baryon number | Topological winding number | Plausible if Skyrme-like |
| Lepton number | Separate topological invariant? | Unknown |
| CPT | Lorentz invariance | Consistent if Option A |
| CP violation | TBD (needed for baryogenesis) | Not addressed |

### Key Dependency

Almost all conservation laws reduce to the question: **does CHPT have a Lagrangian formulation with the right symmetries?** Writing down this Lagrangian is the critical task. The axioms are COMPATIBLE with all required symmetries, but compatibility is not a derivation.
