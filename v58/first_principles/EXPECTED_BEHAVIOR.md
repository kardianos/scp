# Expected Behavior from Physics: Field Excitations and Emergent Geometry from Density Gradients

**Date**: 2026-05-18  
**Status**: Foundational requirements document. Derives, from established experimental phenomenology and the hypothesis that the sole fundamental entity is a dynamical field medium whose local state (density, connectivity, correlation structure) determines both particles and causal structure, the necessary properties any candidate model must possess.

This document is written in textbook style: it states what *must* be true, contrasts with what prior simulations actually produced, and lists concrete, testable constraints with quantitative anchors. It does not narrate discovery history.

---

## 1. The Physical Hypothesis (Post-Layer-Critique Posture)

The project hypothesis, sharpened by LAYER_CRITIQUE.md, is:

- The only fundamental object is a **field medium**.
- **Particles** (protons, electrons, photons, gravitons as collective modes) are **localized, stable excitations** or persistent patterns within that medium.
- **Gravity** is not a fundamental force or an added field. It is an **emergent effect** of spatial variations in the medium's own state — specifically, gradients in "node density," correlation length, stiffness, or causal connectivity — that modulate the *effective* propagation speed and causal structure experienced by *all* excitations.
- Microscopically, the propagation speed *between adjacent elements* ("local c") is constant. Macroscopically, the effective speed or clock rate varies because a coordinate distance corresponds to a different number of elements or a different phase accumulation depending on local density. Higher density → more elements per macroscopic length → signals or internal oscillators experience a different effective rate when compared across regions.

Under this hypothesis, the Einstein equivalence principle, the universality of free fall, the absence of dispersion in gravitational and electromagnetic propagation, and the precise match of observed geodesics to the Schwarzschild / Kerr / FLRW metrics are not inputs; they are *outputs* that any viable medium dynamics plus excitation spectrum must reproduce.

## 2. What the Prior Framework (v28–v57) Actually Produced

The 6-field Cosserat (and later GA multivector) simulations on a fixed Cartesian or Voronoi grid produced:

- Long-lived localized excitations ("oscillons," z-aligned braids, UUD/UDD composites) that persist under periodic boundaries and radiate under absorbing boundaries (v34 and successors).
- **Depletion response** (v43, v51): isotropic UUD composites drift toward regions of self-generated lower |P| or lower background amplitude via the asymmetric binding potential V(P). The drift is a classical force law (approximately 1/r^{1.8} in some measurements) acting on solitons *inside* an unchanging coordinate grid whose metric is Minkowski by construction.
- **Dispersive propagation modifications** (v39 wavefront and BLV analysis): when local amplitude or density-dependent κ altered wavefront speed, the effect was frequency-dependent; different Fourier components traveled at different group velocities. The BLV effective metric derived from the linearized dispersion relation around a braid background was therefore not a true geometric metric (identical null geodesics for all frequencies, all polarizations, all field species).
- Mechanical structure of the braid (v57 assessment): an extended helical flux-tube soliton with theta halo, not a compact object whose stress-energy exhibits a positive-pressure core, negative-pressure tension shell, and satisfaction of the von Laue integral condition observed in lattice-QCD proton gravitational form factors.

**Conclusion from contrast**: These are valuable analogs of binding, charge-dependent interactions, and gradient response *within* a passive fixed-background medium. They are *not* gravity. The causal structure never becomes dynamical; the medium state never defines the metric for test signals of arbitrary type; dispersion appears where universality is required.

## 3. Established Phenomenology That Must Be Reproduced

Any candidate model (pre-geometric or otherwise) must, in a suitable continuum or coarse-grained limit, recover the following to within current experimental precision. These are not optional targets; they are necessary conditions.

### 3.1 Propagation Universality and Absence of Dispersion

- Gravitational waves and electromagnetic waves must travel at identical speed to |c_GW − c_EM| / c ≲ 10^{-15} (GW170817 + GRB170817A; arrival-time difference ~1.7 s over ~40 Mpc baseline; arXiv:1710.05834).
- The effective metric experienced by light, gravitational waves, and matter must be the *same* (no birefringence, no species-dependent refraction).
- High-energy photons and gravitational waves from distant sources must show no measurable frequency-dependent arrival-time dispersion beyond standard plasma or cosmological effects (GRB constraints and LIGO/Virgo high-frequency bands).
- Local Lorentz invariance must be recovered to ≳ 10^{-17} (no preferred frame at that level or tighter, across multiple SME sectors for photons, electrons, and hadrons).

**Implication for density-gradient models**: the mapping from local medium state (density ρ or equivalent) to effective metric g_μν(ρ) must be geometric (conformal or metric rescaling that is frequency-independent) and universal across all fluctuation types.

### 3.2 Equivalence Principle and Universality of Free Fall

- All test bodies, regardless of internal composition, charge, or spin, must fall identically in an external "gravitational" field (Eötvös-type bounds ≲ 10^{-13}, improving in some sectors).
- Light deflection, Shapiro time delay, and perihelion precession must match the predictions of the Schwarzschild metric in the weak-field limit (PPN γ − 1 ≲ 2 × 10^{-5} from Cassini; tighter bounds from VLBI and pulsar timing).
- Clocks (any internal periodic process built from field excitations) at different densities must exhibit the gravitational redshift and time-dilation factors of GR when compared by light signals.

### 3.3 The Orbital / Tick-Rate Intuition Made Quantitative

Consider two satellites in circular orbits at different radii in a central density gradient (higher ρ closer to center). Local microscopic c is constant between neighboring field elements. However:

- In the higher-density region the number of elements (or the integrated phase length) per unit coordinate radius is larger.
- An internal oscillator (the "clock" defined by the particle's own field excitations, e.g., the breathing or precession period of a composite) completes its cycle after accumulating a fixed microscopic phase. Because each macroscopic time interval corresponds to more elements at higher ρ, the proper time dτ per coordinate dt is smaller: the clock "ticks slower" when viewed from afar via light signals.
- Light signals themselves traversing the same region accumulate extra phase or path length, reproducing Shapiro delay.
- Gravitational waves (tensor or collective modes of the same medium) must experience the identical effective null geodesics.

**Consequence**: any model in which "density" merely modulates a force on pre-existing particles fails. The particles *are* the excitations whose internal dynamics and whose propagation through the medium are both governed by the same local state. The effective metric must emerge for the medium's own fluctuations.

### 3.4 Long-Range Inverse-Square Behavior and Weak-Field Limit

- In the low-density, weak-gradient regime the effective force or geodesic deviation must reproduce Newtonian 1/r² gravity (or the corresponding curvature of the emergent metric) with G_eff set by the microscopic parameters of the medium.
- The sourcing of the density perturbation by energy (rest mass of excitations) must be consistent with the equivalence of inertial and gravitational mass.

### 3.5 Particle Properties and Internal Structure

- Stable localized excitations with positive rest mass must exist. Their internal field configuration must support a non-zero proper mass (energy in the particle's rest frame).
- Composite objects (baryon analogs) must exhibit internal pressure profiles compatible with stability: typically a positive (repulsive) core pressure balanced by negative (inward tension) pressure at larger radius, satisfying the von Laue condition ∫ p dV = 0 (as required by lattice QCD proton gravitational form factors and by general relativity for any self-gravitating body in equilibrium).
- The spectrum must include both fermionic and bosonic excitations (or composites) and must allow for electromagnetic-like and strong-like interaction hierarchies without fine-tuning beyond the medium dynamics.

### 3.6 Strong-Field and Cosmological Consistency (Aspirational but Necessary Long-Term)

- The same medium dynamics must support regimes in which effective curvature becomes strong (horizons, ergospheres) without introducing dispersion or violating the universality already required at weak field.
- The global structure must permit a bounded, mappable universe whose expansion or large-scale density variation is consistent with observed cosmology.
- On galactic and cosmological scales the medium state itself (rather than only localized baryonic excitations) can contribute extended density that sources additional gravitational effects, offering a possible account of the phenomenology currently attributed to dark matter (see `../pregeometric/COSMOLOGICAL_DENSITY_AND_DARK_MATTER.md`). Any such contribution must still satisfy the universality and non-dispersion requirements of §3.1 and §4.
- A deeper possibility is that the medium has an intrinsic tendency to favor locally higher density. Stable particles and composite structures can then be understood as the field’s most effective mechanisms for achieving and sustaining significantly higher local density than the free field can maintain on its own. Gravity emerges as the geometric consequence of these density-achieving processes (see `../pregeometric/PARTICLES_AS_DENSITY_ACHIEVERS.md`).

## 4. Concrete, Testable Requirements for Any Candidate Model

A model under consideration for v59+ implementation must demonstrate, analytically or numerically, that it can satisfy **all** of the following (in the appropriate limit):

1. **Non-dispersive geometric propagation** (see §3.1): the effective light cones / null geodesics for small fluctuations are independent of frequency and of the species of the fluctuation (scalar, vector, tensor) to within 10^{-15} or better.
2. **Universality**: the effective metric g_μν is the same for all excitations; any "charge" or internal quantum numbers affect only the sourcing of the medium state, not the propagation on it.
3. **Local Lorentz invariance**: in a region of locally uniform medium state, the emergent geometry is Minkowski to experimental precision (no residual preferred frame at 10^{-17}+).
4. **1/r² (or GR-equivalent) long-range limit**: the weak-field, low-density limit reproduces the Newtonian potential or the corresponding curvature with the correct scaling.
5. **Self-consistent back-reaction**: the presence and motion of excitations source changes in the medium state (density, connectivity, etc.) that in turn affect propagation for everything, including the excitations themselves.
6. **Stable composites with realistic mechanical structure** (see §3.5): the theory admits long-lived, self-bound multi-excitation states whose radial stress-energy profiles include the core positive / shell negative pressure pattern required by real hadrons.
7. **Rest mass from internal structure**: localized excitations possess a positive rest energy arising from their internal field configuration, not inserted by hand.
8. **Clock behavior under density gradients** (see §3.3): internal periodic processes of composite excitations exhibit the correct gravitational time dilation when compared across density regions via null signals.

**What has NOT been established at this stage** (explicit nulls for clarity):

- The microscopic definition of "density" or "node number" in the medium (number density of discrete elements, local algebraic norm, entanglement density, causal-link density, etc.).
- The precise functional form medium_state → g_μν or c_eff(ρ).
- The dynamical equation governing the evolution of the medium state (continuity, diffusion, wave equation on the substrate, etc.).
- Whether a pre-geometric substrate can generate an effective 3+1 Lorentzian geometry with the correct signature and dimension in the continuum limit.
- The particle spectrum and quantum numbers that emerge from the chosen substrate.
- Whether the Newtonian G and the Planck scale arise naturally or require additional parameters.

## 5. Implications for Pre-Geometric and Other Substrates

Because the requirements above are geometric and universal at the effective level, a pre-geometric starting point (no a-priori metric, no fixed grid) has a structural advantage: if an effective metric *emerges* from the collective state of the algebraic or relational substrate, then universality and non-dispersion are automatic consequences of the emergence mechanism rather than imposed constraints. The same relations that support stable localized patterns (particles) define the causal ordering experienced by fluctuations.

Two-layer explicit-medium models can serve as diagnostic prototypes: one writes the medium state ρ(x) as an independent dynamical field sourced by excitations, then defines propagation operators whose characteristics depend on ρ. Such models are still grid-based and therefore provisional, but they make the geometric-propagation question measurable and can falsify or refine the functional requirements above.

## 6. Next Steps After This Document

- Use `EXPECTED_BEHAVIOR.md` as the acceptance criteria for any proposed pre-geometric prototype or two-layer realization.
- In `pregeometric/`, articulate the minimal relational rules that could generate bounded mappable structure, emergent light cones, and stable excitations consistent with the list in §4.
- Maintain a living "requirements traceability" note that maps each experimental anchor (GW170817, Eötvös, von Laue, etc.) to the corresponding modeling constraint.
- Only after the above is stable, open a subsequent version directory containing the first concrete substrate implementations and small-scale numerical or analytic tests.

---

*This document is deliberately substrate-agnostic. Its purpose is to make the physics target unmistakable before any new code is written.*