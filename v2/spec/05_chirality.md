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
*   **Orthogonality**: A neutral particle ($Q=0$) can still be chiral ($\chi = L$ or $R$). This creates a natural place for the **Neutrino** (Neutral but Chiral), resolving the contradiction in the previous draft.

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

**Maps to**: Photons (spin-1, massless, no electric charge) — though photons are described as null-rotors rather than knots in CHPT (see [06_null_rotors.md](06_null_rotors.md)). Neutrinos are tentatively mapped here but this creates a serious problem (see below).

### Known Problem — Neutrino Weak Interactions

The achiral-neutrino mapping is **internally contradictory**. If neutrinos are achiral:
- They should interact only through density gradients (gravity). Rule 3 says achiral knots pass through others with "only small perturbations."
- But real neutrinos interact via the weak force — not as strongly as charged particles, but far more strongly than gravity alone. Neutrino cross-sections are ~10^25 times larger than gravitational cross-sections.

This means either:
1. **Neutrinos are NOT achiral** — they have some form of chirality that is distinct from electric charge. This requires a richer chirality structure than binary left/right.
2. **The weak force is NOT a chirality interaction** — it operates through a different mechanism entirely (perhaps knot transitions as in [09_nuclear_interactions.md](09_nuclear_interactions.md), Option B).
3. **The achiral-equals-neutral mapping is wrong** — neutral particles can still have internal chiral structure that doesn't manifest as electric charge but does enable weak interactions.

Option 3 is the most promising: in the Standard Model, neutrinos are left-chiral (they couple to the weak force) but electrically neutral. This suggests CHPT needs **multiple independent chirality types** — one for electric charge, at least one more for weak charge. See [09_nuclear_interactions.md](09_nuclear_interactions.md).

**Status**: Unresolved. This is a significant conceptual gap.

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
| Left-chiral knot | Negatively charged particle | Medium |
| Right-chiral knot | Positively charged particle | Medium |
| Achiral knot | Neutral particle (neutrino?) | Low |
| Same-chirality repulsion | Pauli exclusion principle | Medium |
| Opposite-chirality nesting | Nuclear/hadronic binding | Low |
| Mirror-chiral pair | Particle-antiparticle pair | High |
| Chiral annihilation | Matter-antimatter annihilation | High |
| Chirality conservation | Charge conservation | High |

### Critical Open Questions

1. How does binary chirality produce the multiple charge types (electric, color, weak) seen in nature?
2. How do fractional charges (quarks) emerge from a binary chiral winding?
3. Does same-chirality repulsion correctly reproduce the Pauli exclusion principle for neutral fermions?
4. What mechanism produces parity violation in the weak sector?
5. Is the spin-statistics connection derivable from the field dynamics?
