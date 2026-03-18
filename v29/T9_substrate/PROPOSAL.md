# T9: The Model Behind the Model — Substrate Derivation

## Question
What underlying physical structure, when coarse-grained or rapidly expanded
from a dense state, naturally produces L = ½(∂φ)² - V(φ₀φ₁φ₂)?

## Key Insight
The triple product P = φ₀φ₁φ₂ is the determinant of the deformation
gradient (when field index a = spatial direction a). The potential V(det F)
penalizes local volume change. This is the Lagrangian of a 3D nonlinear
elastic medium with a specific volumetric strain energy function.

The substrate must have IRREDUCIBLE 3-BODY interactions — pairwise
potentials cannot produce a triple product. The energy must depend on
the VOLUME of simplices (tetrahedra), not just edge lengths.

## Track A: Lattice Reverse-Engineering (computational)
See T9A_lattice/. Test which discrete lattice + 3-body interaction
reproduces V(P) under coarse-graining via Cauchy-Born rule.

## Track C: Expansion Cosmology (physical)
See T9C_expansion/. Start dense, expand, cool. Do braided defects
form spontaneously via the Kibble-Zurek mechanism?

## Track B: Symmetry Classification (analytical, paper work)
The Lagrangian has symmetry S₃ × (Z₂)³. Classify lattice types whose
continuum limit has this symmetry. The Cauchy-Born rule connects
lattice pair/triplet potentials to continuum strain energy:

    W(F) = (1/V₀) Σ_{bonds} φ₂(|F·d|) + (1/V₀) Σ_{triangles} φ₃(V_tri(F))

where F is the deformation gradient, d are bond vectors, V_tri is the
deformed triangle volume. If φ₃(V) = (μ/2)V²/(1+κV²), we recover V(P).

## Track D: Minimum Information Structure (theoretical)
What is the simplest cellular automaton / lattice gas that:
(i) Has 3 equivalent DOF per site
(ii) Conserves a winding number
(iii) Supports localized propagating structures
(iv) Has both symmetric and antisymmetric gradient modes

This is a search through rule spaces — potentially very large.
Defer unless Tracks A-C produce clear results.
