# 05 — Chirality

This chapter defines chirality in CHPT, explains how it gives rise to particle/antiparticle distinctions, governs interactions between knots, and maps to charge and the Pauli exclusion principle.

Depends on: [01_field_axioms.md](01_field_axioms.md), [04_knots_and_particles.md](04_knots_and_particles.md)

---

## Definition of Chirality

Chirality ("handedness") is a geometric property: an object is chiral if it cannot be superimposed on its mirror image by any rotation. Your left and right hands are chiral — mirror images of each other, but no rotation will make one match the other.

In CHPT, knots can be chiral: the 3D+time density pattern of a knot can have a handedness. Specifically, the internal oscillation/rotation of the density within the knot defines an orientation, and in 3D, there are exactly two inequivalent orientations (left and right).

### Three Classes of Knots

1. **Left-chiral knots**: Internal pattern has left-handed orientation.
2. **Right-chiral knots**: Internal pattern has right-handed orientation (mirror image of left).
3. **Achiral (symmetric) knots**: Internal pattern is its own mirror image. No handedness.

This three-way classification is fundamental to the interaction rules.

---

## Chirality vs. Topology (The Distinction)

We distinguish between two unrelated geometric properties of a knot:

1.  **Topological Index ($Q$)**: The winding number of the field map onto the target space (e.g., Hopf invariant).
    *   **Physics Map**: **Electric Charge**.
    *   $Q = +1$ (Positive), $Q = -1$ (Negative), $Q = 0$ (Neutral).

2.  **Chirality ($\chi$)**: The "handedness" or orientation of the internal flow/process (Left vs. Right).
    *   **Physics Map**: **Particle Type / Weak Interaction**.
    *   Distinguishes Matter (Left-dominant?) from Antimatter (Right-dominant?) or determines weak interaction couplings.

### Why this is robust
*   **Charge Quantization**: Topological indices are integers by definition. This explains $Q = n \cdot e$.
*   **Orthogonality**: A neutral particle ($Q=0$) can still be chiral ($\chi = L$ or $R$). The **neutrino** is the physical proof: it is electrically neutral but maximally left-chiral (coupling to the weak force). This independence of charge and chirality is essential for the theory.

### Fractional Charges (Quarks)
If Charge is a winding number, quarks ($Q = \pm 1/3, \pm 2/3$) imply that the fundamental winding can be fractionally partitioned, or that the "elementary" winding is 3 (making quarks integers in a different base). This suggests the knot topology correlates with $SU(3)$ or comparable groups.

---

## Interaction Rules

The chirality of a knot determines how it interacts with other knots. These rules are central to the theory:

### Rule 1 — Same-Chirality Repulsion

Two knots of the same chirality (both left or both right) repel each other. They cannot overlap.

**Mechanism**: When two same-chirality knots approach, their internal oscillation patterns interfere destructively at the boundary, creating a density barrier between them. The closer they get, the stronger the barrier. This is not electromagnetic repulsion — it is a direct geometric incompatibility of their density patterns.

**Maps to**: The Pauli exclusion principle. Two identical fermions cannot occupy the same quantum state. In CHPT, this is because two identical chiral knots literally cannot occupy the same space — their patterns are geometrically incompatible.

**Critical note**: The Pauli exclusion principle applies to ALL fermions, not just same-charge ones. Two electrons repel via Pauli exclusion, but so do two neutrons (which are neutral). If chirality equals charge, then neutrons are achiral, and the Pauli exclusion mechanism (same-chirality repulsion) would not apply to them. This is a problem. Possible resolution: neutrons are composites of chiral sub-knots, and the exclusion operates at the sub-knot (quark) level. Two neutrons cannot overlap because their internal chiral sub-knots conflict. This must be verified.

### Rule 2 — Opposite-Chirality Attraction/Nesting

Two knots of opposite chirality (one left, one right) CAN overlap and share the same spatial region. When they do, they can form a stable composite state.

**Mechanism**: Opposite-chirality patterns interfere constructively at the boundary, reducing the total gradient energy. The composite has lower energy than the separated pair.

**Maps to**: Matter-antimatter combinations (which annihilate rather than bind — so this needs nuance), and more importantly, the binding of quarks with different color/flavor quantum numbers inside hadrons.

**Nuance — Annihilation vs. Binding**: If a left-chiral knot meets its exact mirror (same knot type, opposite chirality), they annihilate — the density patterns cancel and disperse as radiation. But if they are different knot types with opposite chirality, they can bind stably. Example: in a proton, the up-quark sub-knots and down-quark sub-knots have different chiralities but different structures, so they bind rather than annihilate.

### Rule 3 — Achiral Transparency

Achiral (symmetric) knots interact only weakly with other knots. They pass through chiral knots with only small perturbations from density gradient interactions.

**Maps to**: Photons (spin-1, massless, no electric charge) — though photons are described as null-rotors rather than knots in CHPT (see [06_null_rotors.md](06_null_rotors.md)).

**Does NOT map to neutrinos.** Neutrinos are chiral (left-handed only), not achiral. They pass through matter because they lack electric charge and the weak force is short-range (~10^-18 m), not because they lack chirality. See the neutrino discussion below.

### The Neutrino Problem — Resolved by Multi-Axis Chirality

The earlier draft mapped neutrinos to achiral knots (Rule 3: achiral transparency). This is **internally contradictory**:

- Rule 3 says achiral knots interact only through density gradients (gravity).
- But real neutrinos interact via the weak force — ~10^25 times stronger than gravity.
- Neutrino cross-sections (~10^-38 cm²) are tiny compared to EM (~10^-24 cm²), but enormous compared to gravity (~10^-63 cm²).

**Resolution**: Neutrinos are **NOT achiral**. They are chiral — but in a way that is *independent of electric charge*.

In the Standard Model, this is well established: neutrinos are **maximally chiral** (only left-handed neutrinos exist; only right-handed antineutrinos). They couple to the weak force *because* they are chiral, not despite it. The reason they pass through most matter is that they lack **electric charge**, not that they lack chirality. The weak force has very short range (~10^-18 m, set by the W/Z boson mass), so the neutrino must pass very close to another particle's core to interact.

**Implication for CHPT**: A single binary chirality axis (L/R) is insufficient. The theory needs **at least two independent geometric properties** that map to independent physical charges:

| Geometric Property | Physical Observable | Neutrino Value |
|---|---|---|
| Winding number $Q$ | Electric charge (or baryon number — see [mapping framework](../../proposal/mapping_framework/geometric_to_physical_mapping.md)) | $Q = 0$ (neutral) |
| Weak chirality $\chi_W$ | Weak isospin / weak force coupling | $\chi_W = L$ (left-handed only) |

The neutrino is the **proof** that these two properties are independent: it has zero electric charge but maximal weak chirality. Conversely, right-handed electrons have electric charge but zero weak isospin (they don't couple to W bosons).

**Where does weak chirality live in the algebra?** The Cl⁺(3,0,1) field has 8 components: 4 bulk (s, f₁, f₂, f₃) and 4 degenerate (j₁, j₂, j₃, p). The bulk quaternion provides the winding number ($\pi_3(S^3) = \mathbb{Z}$). Weak chirality likely lives in the **degenerate sector** — the (J, P) components that are currently non-dynamical (see [15_open_problems.md](15_open_problems.md), B6). Making this sector dynamical could:
- Give the degenerate modes their own topological structure
- Provide a geometric origin for weak chirality independent of electric charge
- Explain why neutrinos interact weakly but not electromagnetically

**Status**: Conceptually resolved (neutrinos are chiral, not achiral). Mathematically open — the degenerate sector dynamics must be specified (Roadmap item 7) before the weak chirality can be formalized. See also [proposal/mapping_framework/geometric_to_physical_mapping.md](../../proposal/mapping_framework/geometric_to_physical_mapping.md), §8.

---

## Chirality and Spin

In the Standard Model, chirality and spin are related but distinct concepts. Left-chiral and right-chiral fermions transform differently under the weak force. Spin is intrinsic angular momentum.

In CHPT, both chirality and spin arise from the internal rotational structure of the knot:

- **Chirality**: The handedness (sense of rotation) of the internal pattern.
- **Spin**: The magnitude of the intrinsic angular momentum of the internal pattern.

For spin-1/2 particles: the knot's internal pattern requires a 720-degree rotation (4-pi) to return to its initial state. This is a defining property of spinors and is consistent with the knot being a configuration in a geometric algebra.

For spin-1 particles: the pattern repeats every 360 degrees (2-pi). These would be achiral or vector-type configurations.

### Unknown

- **Spin-statistics connection**: The spin-statistics theorem (spin-1/2 particles obey Fermi-Dirac statistics; spin-1 particles obey Bose-Einstein statistics) is one of the deepest results in quantum field theory. In standard physics, it follows from Lorentz invariance + quantum mechanics + locality. CHPT must derive this from its field dynamics. Specifically, the claim that same-chirality knots repel (Rule 1) must be shown to apply to ALL spin-1/2 knots and NONE of the spin-1 configurations. This is not guaranteed and must be proven.

---

## Chirality Breaking and the Weak Force

In the Standard Model, the weak force violates parity (P) — it treats left-chiral and right-chiral particles differently. This is one of the most surprising facts in physics (Wu experiment, 1957).

CHPT can potentially accommodate this: if the field equation or the knot dynamics are not perfectly mirror-symmetric, left and right chirals would interact differently. But this raises a deep question:

### Is the Field Equation Chiral?

**Option A — Chiral Field Equation**: The field equation itself distinguishes left from right. Parity violation is built into the fundamental dynamics.

**Option B — Achiral Field Equation, Chiral Vacuum**: The field equation is parity-symmetric, but the ground state (vacuum) spontaneously breaks parity. Left-right asymmetry is a property of our universe's vacuum, not of the laws themselves.

**Option C — Achiral Field Equation, Chiral Knots**: The field equation is parity-symmetric, but specific knot solutions break parity. The weak force arises from interactions between chiral knots that are not parity-conjugates of each other.

### Recommended Path

Option B is most consistent with standard physics (the Higgs mechanism breaks electroweak symmetry, which includes parity). Option C is more novel. Both need mathematical formulation to evaluate.

---

## Summary of Chirality Mappings

| CHPT Concept | Standard Physics Analog | Confidence |
|-------------|------------------------|------------|
| Left-chiral knot | Left-handed particle (weak-interacting) | Medium |
| Right-chiral knot | Right-handed particle (weak-singlet) | Medium |
| Achiral knot | Photon / truly non-interacting mode | Medium |
| Neutral + chiral | Neutrino ($Q=0$, $\chi_W = L$) | Medium |
| Same-chirality repulsion | Pauli exclusion principle | Medium |
| Opposite-chirality nesting | Nuclear/hadronic binding | Low |
| Mirror-chiral pair | Particle-antiparticle pair | High |
| Chiral annihilation | Matter-antimatter annihilation | High |
| Chirality conservation | Weak isospin conservation | High |
| Winding number conservation | Charge / baryon number conservation | High |

### Multi-Axis Geometric Structure

The theory requires **at least two independent geometric properties** to accommodate the known particle spectrum:

1. **Winding number** ($Q \in \mathbb{Z}$): From $\pi_3(S^3)$. Maps to electric charge or baryon number (see [mapping framework](../../proposal/mapping_framework/geometric_to_physical_mapping.md), §2).
2. **Weak chirality** ($\chi_W \in \{L, R, 0\}$): Independent of $Q$. Maps to weak isospin. Likely originates from the degenerate sector of $Cl^+(3,0,1)$.

The neutrino ($Q = 0$, $\chi_W = L$) and the right-handed electron ($Q = -1$, $\chi_W = 0$) prove these are independent. A single chirality axis cannot accommodate both.

### Critical Open Questions

1. What geometric structure in $Cl^+(3,0,1)$ provides weak chirality independently of winding number?
2. How do fractional charges (quarks) emerge from integer winding numbers?
3. Does same-chirality repulsion correctly reproduce the Pauli exclusion principle for neutral fermions (via composite sub-knot chirality)?
4. What mechanism produces parity violation (left-right asymmetry) in the weak sector?
5. Is the spin-statistics connection derivable from the field dynamics?
6. Does the degenerate sector $(J, P)$, once made dynamical, provide the topological structure for weak chirality?
