# Pre-Geometric Exploration — Initial Framing and Rationale

**Date**: 2026-05-18  
**Context**: Route 3 from LAYER_CRITIQUE.md. This sub-folder holds the early conceptual work on starting points that do not presuppose the 3+1 metric, causal structure, or continuum limit that the theory must ultimately reproduce. It is deliberately placed after the first-principles requirements document so that any concrete proposal can be evaluated against the necessary conditions in `../first_principles/EXPECTED_BEHAVIOR.md`.

---

## 1. Why Pre-Geometric?

The fixed-background problem diagnosed in LAYER_CRITIQUE.md is fundamental: once a grid (Cartesian, Voronoi, or otherwise) with a pre-defined notion of distance and time is introduced, any "gravity" that appears is a force or drift inside that fixed geometry. The medium never gets to *define* the causal structure.

A pre-geometric (combinatorial, algebraic, relational) starting point avoids the circularity by construction:

- There is no manifold, no dx^μ, no pre-assumed light cones, no embedding space.
- The fundamental objects are discrete elements (points, simplices, spin-network nodes, group elements, basis multivectors, abstract relations) together with incidence or composition rules.
- Dynamics act on the elements and on the *relations* among them.
- Dimension, signature, effective metric, and light cones are *outputs* of the collective state of the structure, not inputs.
- "Particles" appear as stable, localized, bounded patterns of relations or excitations that persist under the dynamics.
- "Density" is defined relationally (local degree, number of causal predecessors within a light-cone-like cut, algebraic multiplicity, connectivity density, etc.) and can vary spatially in the emergent sense.
- Effective propagation for small fluctuations is then defined by the same relational structure that supports the particles, automatically yielding universality and geometric (non-dispersive) behavior if the emergence is successful.

This matches the long-term preference stated in LAYER_CRITIQUE §2.3 and §3: "the only approach that does not presuppose the geometry it hopes to derive."

## 2. The User's Core Intuition in Pre-Geometric Language

"Everything is field; particles are localized field excitations; gravity is an emergent effect of field density gradients modulating effective propagation / geometry."

In a pre-geometric substrate this becomes:

- The "field" is the entire relational or algebraic structure and its evolution rules.
- A "localized field excitation" is a persistent, self-reinforcing sub-pattern (a defect, a bound cluster of relations, a coherent phase or winding in the algebraic data) that carries energy, can move relative to other patterns, and has internal degrees of freedom (hence rest mass).
- A "density gradient" is a spatial variation (in the emergent coarse-graining) in the local density of elements or in the local density of connections / causal links.
- "Effective propagation" between two distant patterns is determined by the shortest or dominant causal chains connecting them; a region with higher element or link density per emergent coordinate length requires more steps or more phase accumulation, producing the tick-rate and light-delay effects described in EXPECTED_BEHAVIOR.md.

Because the same structure generates both the particles and the causal ordering, the equivalence principle and the non-dispersive character are natural rather than engineered.

### 2.1 Gravitational Waves as Collective Medium-State Modes

A direct consequence of the hypothesis in EXPECTED_BEHAVIOR.md §3.1 is that gravitational waves are not vibrations of particles inside a geometry; they are propagating perturbations of the medium state \(S\) (density, connectivity, causal valence, or equivalent relational quantity) itself. This is the pre-geometric translation of the sound-wave analogy: the wave is a collective compression/rarefaction (or shear) of the relational structure. Because the effective causal ordering is defined by the same \(S\), the perturbation automatically affects every other excitation (light, matter, other GWs) in a universal, geometric way.

This reading is developed in detail in `MEDIUM_STATE_WAVES.md`. It imposes two especially sharp filters on any candidate substrate:
- The collective mode must be (or project onto) transverse-traceless tensor character at long range; purely scalar density waves are ruled out by observed GW polarization.
- The perturbation of \(S\) must change the emergent light cones identically for all small disturbances (the non-dispersive universality required by GW170817 and the broader §3.1 constraints).

Early exploration should therefore ask, for each relational rule set: "Can it support a propagating, tensorial modulation of its own connectivity or density that looks like a GR gravitational wave in the coarse-grained limit, while still allowing stable localized patterns?"

The direct multivector expression of this idea (density deviation \(\rho\) of \( M \), induced bivector \(\omega = \nabla M \tilde{M}/|M|^2\), force via geometric product) is developed in `MULTIVECTOR_DENSITY_GRAVITY.md`, including the explicit reduction to the Newtonian \( 1/r^2 \) law and a detailed consistency check against observed gravitational-wave results.

## 3. Minimal Requirements for a Useful Pre-Geometric Starting Point

A candidate substrate must, at minimum, be capable of addressing the open questions listed in LAYER_CRITIQUE §4 and the concrete requirements in EXPECTED_BEHAVIOR §4. Early exploration should therefore focus on:

1. **Generation of bounded, mappable structure**: the dynamics or incidence rules must produce finite, distinguishable regions rather than a featureless infinite soup. The "universe expands from nothing, but no matter the expansion, the set of what is is bounded and mapped to what we perceive" intuition must be realizable.
2. **Emergence of causal structure**: local rules must induce an effective partial order or light-cone-like relation among elements so that "before / after / spacelike" becomes meaningful in the large-scale limit.
3. **Stable localized patterns**: there must exist persistent sub-structures that behave as particles (can be tracked, carry conserved or quasi-conserved quantities, possess internal oscillation or rotation modes that define rest mass).
4. **Variable local density**: the structure must support regions of higher and lower element / connectivity density whose effect on propagation of small disturbances is geometric (same for all small modes) and universal. On galactic and larger scales this extended medium density can source the extra gravitational effects currently attributed to dark matter without new particles (see `COSMOLOGICAL_DENSITY_AND_DARK_MATTER.md`).
5. **Dimensionality and signature**: in some controlled limit the emergent geometry should be 3+1 with Lorentz signature; this is the hardest and most important filter. Early work may begin in lower or toy dimensions but must have a plausible path to 3+1.
6. **No hidden continuum assumption**: the rules must remain well-defined and discrete at every scale; the continuum / metric must arise only as a derived, coarse-grained description.

## 4. Example Classes of Substrate (Exploratory, Not Commitments)

The following are mentioned only to illustrate the space; none is adopted yet. Each will be scored against the first-principles requirements before any implementation investment.

- Causal sets or causal networks with growth dynamics (elements added according to local rules; density variations arise from growth rate).
- Spin networks or group field theory elements with evolution moves that change connectivity.
- Abstract Clifford or multivector algebras without an a-priori metric (relations and products generate effective inner products and causal structure).
- Incidence geometries or simplicial complexes with purely combinatorial dynamics.
- Quantum information or entanglement networks in which "density" is entanglement density and effective geometry emerges from entanglement structure (ER=EPR style but kept classical/algebraic at first).

Deliberately set aside (per LAYER_CRITIQUE preference): any model that begins by assuming a mechanical lattice whose discreteness is the starting point, or any soliton model whose background manifold is fixed.

## 5. Relation to the First-Principles Requirements

Every proposal in this folder must eventually map onto the list in `../first_principles/EXPECTED_BEHAVIOR.md`:

- How does the substrate guarantee non-dispersive, universal propagation in the emergent limit?
- How does a density variation (relational or algebraic) translate into the tick-rate / Shapiro / light-deflection phenomenology without introducing species dependence or frequency dependence?
- What is the precise definition of "density" (or its tensorial generalization) and of the effective metric it induces?
- Can the substrate support a propagating, transverse-traceless tensorial collective mode of its own state (see `MEDIUM_STATE_WAVES.md`) while still allowing stable localized patterns?
- Can stable composites with the required mechanical pressure profiles appear as bound states of the fundamental relational excitations?

Initial exploration is allowed to be qualitative and low-dimensional. The first decisive results will be small proofs-of-concept or analytic arguments showing that a given relational rule set can produce at least two of: (a) stable localized patterns, (b) an emergent causal ordering that varies with local density, (c) absence of built-in dispersion for collective modes.

## 6. What This Folder Will Not Do

- It will not import hedgehog, Skyrme, Cosserat, or other continuum soliton literature as the null model.
- It will not assume a background manifold "on which" the pre-geometry lives.
- It will not jump to large-scale numerical simulation before the minimal relational rules have been shown, at least conceptually, to be compatible with the requirements in EXPECTED_BEHAVIOR.md.

## 7. Immediate Next Artifacts (Planned)

- A short "candidate substrate scorecard" mapping each explored algebra/relation set to the §4 requirements of EXPECTED_BEHAVIOR.md (now explicitly including the tensor-wave filter from MEDIUM_STATE_WAVES.md).
- `MEDIUM_STATE_WAVES.md` (already created) as the living conceptual development of gravitational waves as collective medium-state excitations.
- `MULTIVECTOR_DENSITY_GRAVITY.md` (just added) giving the direct multivector-algebra expression for gravity as the geometric effect of density gradients in \( M \), with explicit Newtonian limit and GW consistency analysis.
- A one-page "density definition" note once the first promising relational notion of density (or its tensorial generalization) is identified.
- If a minimal toy model (e.g., 1+1 or 2+1 relational dynamics) yields an emergent metric with the correct qualitative behavior — including support for a transverse tensor collective mode — a small prototype implementation (language and size to be chosen for clarity, not performance) may be opened in a follow-on version directory.

---

*Pre-geometric work is high-risk, high-reward. A negative result (inability to generate 3+1 Lorentzian structure with stable particles from pure relations) would be as valuable as a positive one. The first-principles requirements give the yardstick.*