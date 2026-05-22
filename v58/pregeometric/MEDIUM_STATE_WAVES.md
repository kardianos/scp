# Gravitational Waves as Collective Excitations of the Medium State

**Date**: 2026-05-18  
**Status**: Conceptual theory note. Develops the hypothesis, first raised in the conversation around EXPECTED_BEHAVIOR.md §3.1, that gravitational waves are propagating perturbations in the fundamental medium state itself (density, connectivity, causal valence, etc.) rather than vibrations of particles within a pre-existing geometry. This note maps the idea onto the first-principles requirements and extracts constraints for pre-geometric substrate design.

---

## 1. The Distinction That Matters

In the post-layer-critique ontology there are two very different pictures of a "gravitational wave":

**Picture A (ruled out by the layer critique)**:  
Particles are localized field excitations living inside a passive or slowly varying background medium. A GW is a separate dynamical field (or a collective mode of the particles) that exerts a force or modulates the trajectories of those particles. The background causal structure remains essentially fixed; the wave is something that happens *to* the particles *on top of* the geometry.

This is the picture realized (approximately) by the v43/v51 depletion drift and the v39 dispersive wavefront effects: forces and speed modifications acting on solitons inside an unchanging coordinate grid.

**Picture B (the one required by the hypothesis)**:  
The sole fundamental entity is the field medium. Its local state \(S\) (node density, connectivity density, local causal valence, algebraic norm, frame field, or equivalent relational quantity) *defines* both the existence of particles and the effective causal structure experienced by every fluctuation. A gravitational wave is a small-amplitude, propagating collective excitation of \(S\) itself:

\[
S(\mathbf{x}, t) = S_0 + \delta S(\mathbf{x}, t)
\]

The particles (stable localized patterns or defects) are carried by the changing \(S\), but they are not the primary oscillating degrees of freedom of the wave. The wave *is* the geometry changing.

This is the direct analog of a sound wave in a material: the wave is a compression/rarefaction of the medium's own state. The medium elements oscillate, but the propagating entity is the disturbance in density/pressure/connectivity.

## 2. Why This Picture Is Structurally Forced by the Hypothesis

If gravity is an emergent effect of gradients in the medium state modulating effective propagation, then any *dynamic* gravitational phenomenon must be a dynamic change in that same state. There is no other place for the dynamics to live.

- A static central mass creates a static gradient \( \nabla S \). The orbital/tick-rate effects of §3.3 in EXPECTED_BEHAVIOR.md follow immediately.
- A propagating disturbance \(\delta S\) is a gravitational wave. It changes the local node density (or connectivity) as it passes, which changes the effective number of microscopic steps or the accumulated phase per macroscopic distance for *every* process defined in the medium: photon phase accumulation, internal oscillator periods of composite excitations, and the causal ordering itself.

Because the effective metric \(g_{\mu\nu}\) (or the local light-cone structure) is a function of the instantaneous \(S\), a wave in \(S\) is automatically a wave in the metric. Universality and geometric (non-dispersive) behavior are then natural consequences rather than extra constraints, provided the emergence of the continuum limit works.

The explicit multivector-algebra realization (density deviation \(\rho\) of the field multivector \(M\), induced bivector \(\omega \propto \nabla M \tilde{M}/|M|^2\)) and the full consistency check against observed GW results are given in the companion note `MULTIVECTOR_DENSITY_GRAVITY.md`.

## 3. Consistency with GR Gravitational-Wave Phenomenology

### 3.1 Speed, No Dispersion, and the GW170817 Bound

The wave \(\delta S\) propagates at the microscopic speed of the medium (or the emergent \(c\)). Photons (coherence or phase waves along causal links of the same medium) and matter particles (bound states whose internal clocks are built from the same links) experience the identical background \(S\). Therefore the effective null geodesics for the GW mode and for all other modes are the same by construction. The \(10^{-15}\) multimessenger bound is satisfied structurally, not by tuning.

Long-wavelength dispersion must be absent (or Planck-suppressed) for the same reason high-energy photons and LIGO events show no extra frequency-dependent delay beyond known effects. This constrains the response function of the medium: at macroscopic scales the collective mode must behave like a massless wave on an effective Lorentzian geometry.

### 3.2 Tensor Character, Two Polarizations, Transverse-Traceless

A purely scalar density perturbation would produce longitudinal or breathing (scalar) waves. Current observations rule these out for the detected events; only the two tensor modes (\(+\) and \(\times\)) are present, and they are transverse.

Therefore the medium state variable \(S\) that sources geometry cannot be a pure scalar. It must possess (or the dynamics must project onto) tensor degrees of freedom at long wavelength. In a pre-geometric substrate this can arise from:
- local directional connectivity (a frame or tetrad whose perturbations carry spin-2),
- the symmetric traceless part of a valence or adjacency perturbation,
- or the way isotropic changes in element density are converted into anisotropic strain by the underlying incidence rules during coarse-graining.

The restriction to the transverse-traceless sector is then a consequence of the emergent Lorentz symmetry, exactly as in linearized GR. This is a strong filter on candidate substrates: a model whose only collective density mode is scalar will fail the polarization test and must be rejected or extended.

### 3.3 Tidal Effects, Energy, and Back-Reaction

A passing \(\delta S\) wave stretches and compresses the local node (or link) density in a quadrupole pattern. Any two nearby test excitations therefore experience a time-varying proper separation even though no "force" is acting on them in the old sense. Their internal periodic processes (clocks) accumulate phase at a periodically modulated rate. This reproduces geodesic deviation and the strain \(h_{+\times}\) measured by interferometers.

Because \(\delta S\) is a real perturbation of the state that defines distances and causal ordering, it carries energy and momentum in the emergent sense. It can be absorbed, can source further excitations, and must itself be sourced by the stress-energy of other excitations (the back-reaction required by EXPECTED_BEHAVIOR §4 item 5). Nonlinearities in the medium dynamics automatically permit GW–GW interactions, memory effects, and the usual GR nonlinearities in the strong-field regime.

### 3.4 The User's Orbital / Tick-Rate Intuition, Now Dynamic

The static case (planet) creates a static \(\nabla S\). A GW is the dynamic, propagating version of the same gradient. As the wave passes a pair of satellites, it periodically alters the local density of elements (or causal links) between them. The effective number of microscopic steps per macroscopic separation therefore oscillates, producing the observed relative strain and the periodic modulation of clock rates. The wave does not primarily shake the satellites' internal constituents; it changes the medium *through which* their internal processes and the signals between them are defined.

## 4. Translation into Pre-Geometric Language

In a purely relational or algebraic substrate (no background manifold):

- "Density" is a local relational quantity: number of elements per emergent volume, average valence (degree), density of causal links crossing a cut, local algebraic multiplicity, entanglement density, etc.
- A GW is a propagating modulation of that quantity (or of the tensorial generalization of it) such that the emergent causal ordering experienced by all other elements is periodically sheared in a quadrupole pattern.
- The same incidence or composition rules that allow stable localized patterns (particles) must also support this collective mode.
- Because the causal structure itself is being modulated, every small disturbance (whether it will later be interpreted as a photon, a graviton, or a matter wave) sees the identical change in its light cones.

The emergence of 3+1 Lorentzian geometry with only two propagating tensor modes at long range then becomes a necessary property of any viable substrate that passes the polarization filter.

## 5. Concrete Constraints This Picture Imposes on Substrate Design

Any candidate pre-geometric starting point must eventually demonstrate (analytically or in controlled low-dimensional prototypes) that it can support:

1. A collective mode \(\delta S\) that is (a) sourced by localized patterns, (b) propagates at the same speed as other small disturbances, and (c) is non-dispersive at long wavelength.
2. The mode must be (or project onto) transverse-traceless tensor character; purely scalar or longitudinal density waves are insufficient.
3. The perturbation \(\delta S\) must change the effective causal ordering or light-cone structure identically for all other excitations (universality).
4. The back-reaction must be two-way: excitations source \(\delta S\), and \(\delta S\) in turn affects the propagation and internal dynamics of excitations.
5. In the continuum limit the mode must reproduce the linearized GR wave equation (or an equivalent geometric description) on the emergent metric it helps define.

Failure on item 2 (polarization) is as decisive a rejection criterion as failure to produce stable particles or an emergent 3+1 signature.

## 6. What This Picture Is Not

- It is *not* particles oscillating inside a gravitational field that exists independently of them.
- It is *not* a separate spin-2 field \(h_{\mu\nu}\) living on a fixed Minkowski background whose stress-energy then sources curvature.
- It is *not* automatically satisfied by any model that merely has "some density variable." The density variable must be the one that *generates* the causal structure, and its collective excitations must be tensorial and universal.

## 7. Relation to the First-Principles Requirements (EXPECTED_BEHAVIOR §4)

This interpretation directly supports or is required by every item:

- Item 1 (non-dispersive geometric propagation): satisfied by construction if \(\delta S\) is a mode of the same structure that defines the light cones.
- Item 2 (universality): automatic — there is only one medium state being perturbed.
- Item 3 (local Lorentz invariance): requires the microscopic substrate to be free of preferred frames in the uniform-state limit.
- Item 5 (self-consistent back-reaction): the two-way coupling between excitations and \(\delta S\) is mandatory.
- Items 6–8 (mechanical structure, rest mass, clock behavior): unchanged from the static case; the wave merely makes the geometry dynamic.

Items 4 (1/r² limit) and the strong-field requirements remain as before.

## 8. Open Questions Specific to This Picture

- What is the microscopic "stress-energy" coupling that sources \(\delta S\) from localized patterns? (Must be derived from the same relational rules, not inserted by hand.)
- How does an isotropic substrate produce a propagating mode whose long-range projection is purely transverse-traceless spin-2? (Selection rule or dynamical projection?)
- What is the nonlinear completion that allows strong-field GWs (ringdown, mergers) to remain consistent with the observed quasinormal modes?
- In a discrete substrate, what prevents the GW from "feeling" the discreteness (dispersion or birefringence) at wavelengths many orders of magnitude above the Planck scale?
- Can the same substrate support both the stable composites of §4 item 6 *and* clean tensor waves without the composites radiating or destabilizing the wave sector?

---

*This note is deliberately conceptual. Its purpose is to make the "GW = medium-state wave" reading of the hypothesis precise enough that substrate proposals can be scored against it before any implementation investment.*