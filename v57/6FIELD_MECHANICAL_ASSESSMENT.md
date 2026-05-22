# 6-Field Cosserat Mechanical Structure — Assessment (May 2026)

**Date**: 2026-05-15  
**Context**: Fresh recreation of the classic single z-aligned 6-field braid (v34/v55 parameters) on the current main kernel, followed by radial mechanical profiling and direct comparison to the lattice QCD proton gravitational form factors (Hackett, Pefkou, Shanahan 2024).

---

## 1. Executive Summary

A clean 6-field Cosserat simulation (N=128, L=20, T=300, exact v34-era parameters with α=β=0) was completed using the current unified kernel. Late-time analysis (t=300) of the resulting SFA shows:

- A persistent main cluster with substantial binding energy.
- Significant theta radiation (theta_rms ≈ 0.08 in the core region).
- **Multiple clusters** by late time — the original braid has fragmented or spawned secondary structures.
- The main cluster is spatially extended (R_half ≈ 3.8 on L=20).
- Radial profiles inside the main cluster show |P|(r) peaked in the core but falling over a large radius, with high theta_rms extending outward.

**Critical mechanical mismatch with the real proton**:

The lattice proton (derived from gravitational form factors A(t), J(t), D(t)) is a compact, roughly spherical object whose stress-energy tensor exhibits:
- High energy density in a small core under **positive (repulsive) pressure**.
- A clear **negative-pressure (inward tension) shell** at larger radius.
- A well-defined mechanical radius smaller than the charge radius.
- Overall satisfaction of the von Laue stability condition (∫p dV = 0).

The 6-field braid produces none of these features. It is an excellent, long-lived, radiating **current-loop / helical flux tube** soliton. It is not a compact isotropic mechanical object with the internal force balance observed in the real proton.

This is high-quality negative evidence. After more than a decade of work, the pure 6-field Cosserat (even with the best current discretization and parameters) has not spontaneously produced the mechanical structure of the proton.

---

## 2. Background — The Long Arc of the Project

### 2.1 Early Phase (oscillons and gradient response)
The project began with nonlinear wave equations supporting long-lived localized breathing structures ("oscillons"). Early work demonstrated that these objects could drift in density gradients, suggesting a possible classical mechanism for gravity. The force law appeared roughly 1/r^1.8 in some measurements, but objects were never truly stable under absorbing boundaries.

### 2.2 Introduction of the Theta Sector (V34 onward)
The addition of three angle fields θ with curl coupling (η curl(φ) ↔ curl(θ)) transformed the theory. This 6-field Cosserat system:
- Produced charge-dependent forces (same-winding objects attracted more strongly).
- Generated an extended theta halo around braids, interpreted as an electromagnetic analog (θ → A_μ).
- Allowed the definition of H_cross (a chirality measure).

The z-aligned braid became the fundamental "particle-like" object: a helical twist in the three φ fields that sources a persistent theta field via the coupling term.

### 2.3 Phase Confinement and Baryon Analogs (V41–V43)
A major conceptual step was the discovery that carrier phase offsets Δ = {0, 2π/3, 4π/3} create a confinement mechanism. At the triple overlap, P → 0 due to destructive interference. This led to the construction of UUD ("proton") and UDD ("neutron") 3-braid composites.

Key results from that era:
- UUD composites were more stable than UDD (matching real physics).
- A deuterium analog (UUD + UDD pair) showed inter-baryon attraction, force equilibration to 1:1 strong/EM ratio, and overall compaction.
- Clean gravitational response was measured on isotropic UUD protons in density gradients (V43). Gravity was shown to be a pure φ-depletion effect, while the z-braid responded oppositely (EM/anisotropy artifact).

These were the strongest "particle from field" results in the project's history.

### 2.4 The Stability Crisis (V44–V54)
Extensive parameter searches, evolutionary algorithms, and equation refinements (Cosserat strain α, curl²-hardening β, variable theta mass, etc.) failed to produce objects that remained compact and stable under absorbing boundary conditions for arbitrarily long times. All candidates eventually dispersed or fragmented when radiation was allowed to leave the system.

V54 concluded that the 6-field Cosserat equations, under both voxel and (later) Voronoi foam discretizations, do not support truly stable particles across hundreds of configurations.

### 2.5 Discretization Pivot (V55)
The Voronoi foam kernel was developed to test whether cubic grid anisotropy was the hidden culprit. After significant engineering (9.4× CPU speedup, cell-native SFA format with temporal delta compression, volview support), the foam discretization was shown to reproduce the same physics as voxels for the V34 braid (once a discrete curl bug was fixed). The stability problem persisted.

### 2.6 Algebraic Pivot (V56)
Faced with the failure of the 6-field system to produce stable compact objects, the project adopted Lengyel's projective geometric algebra ℝ(3,1,1). The new framework supplies:
- Built-in Lorentz covariance.
- Spinor structure via the relativistic quaternion (8-component even subalgebra).
- A Higgs-like vacuum manifold whose homotopy admits topological solitons.

Work moved to an 8-component multivector field with Higgs potential, later augmented by genuine Skyrme L4 terms and hard S³ projection to preserve winding. While topology survival improved dramatically with projection, the "lattice collapse" problem (soliton shrinks, gradients explode, L4 energy diverges) remained. Pion mass terms worsened the collapse.

---

## 3. The Lattice Proton as Quantitative Target

The 2024 Hackett–Pefkou–Shanahan lattice QCD calculation (arXiv:2310.08484) changed the game. For the first time, we have first-principles, flavor-decomposed gravitational form factors of the proton at near-physical pion mass. From these we obtain:

- Energy density ε(r) — positive, centrally peaked, gluons more extended than quarks.
- Pressure p(r) — **positive in the core, negative in the outer region** (the sign change is required for mechanical stability).
- D-term D(0) ≈ −3.87(97) GeV⁻² (large negative).
- Mechanical radius < charge radius.

This is no longer a vague "stable lump" requirement. It is a precise radial mechanical equation of state that any candidate field-theoretic proton must reproduce.

The v57 reference table (`proton_lattice_profiles.tsv`) was generated from the dipole parametrization of these GFFs to serve as the new gold-standard target.

---

## 4. The May 2026 Experiment — Clean 6-Field Recreation

A fresh simulation was executed using the current main unified kernel (`scp_sim.c`) with the exact parameters from the v55 "Stage 0" validation run (`recreate_v34.cfg`):

- N=128, L=20, T=300
- m²=2.25, m_θ²=0, η=0.5, μ=-41.345, κ=50
- α_cs=β_h=0 (pure V34-era physics)
- Single z-aligned braid seed (A=0.8, R_tube=3.0, ellip=0.3325, correct phase offsets)
- Periodic boundaries, f16 output

**Results at t=300 (cluster_profile analysis):**

- Main cluster: 44,167 voxels, P_peak ≈ 0.543, substantial negative potential energy.
- R_half ≈ 3.8 (significantly extended relative to the box).
- theta_rms ≈ 0.0806 in the main cluster (still strongly radiating).
- **5 clusters total** — the original coherent braid has fragmented or spawned secondary bound structures.

Radial profiles inside the dominant cluster show |P|(r) peaked near the center but declining over a large distance, with significant theta_rms persisting into the outer regions. There is no indication of a compact core under positive pressure surrounded by a tension shell.

---

## 5. What the 6-Field Actually Produces Well

The 6-field Cosserat theory has genuine strengths:

- Long-lived, high-binding, radiating braid solitons with excellent energy conservation.
- Charge-dependent (winding-dependent) forces.
- A natural electromagnetic sector (theta as vector potential, braids as current loops).
- Phase-confinement mechanism that gives proton > neutron stability in 3-braid composites.
- Clean separation of gravitational (φ-depletion) and electromagnetic (theta/ anisotropy) responses when composites are used.

These are not small achievements. The theory is rich and has produced several qualitative successes that map onto real physics.

---

## 6. What the 6-Field Has Not Produced

Despite extensive effort, the 6-field system has not produced:

- A single compact, isotropic object whose internal stress-energy tensor exhibits the positive-core / negative-shell pressure distribution measured in lattice QCD.
- True long-term stability under absorbing boundaries (all candidates eventually radiate away or fragment).
- A mechanical radius that is distinctly smaller than the overall size in the manner of the real proton.
- A D-term of the correct magnitude and sign arising naturally from the dynamics.

The objects that survive longest are extended, often cylindrical or multi-lobed, and continue to radiate theta energy outward rather than self-confining into a tight mechanical equilibrium.

---

## 7. Current Thinking (May 2026)

The lattice mechanical profile is now the strongest constraint we have. Lifetime, winding number survival, and global energy conservation are necessary but no longer sufficient criteria.

The 6-field Cosserat appears to be an excellent theory of **braids** — stable, charged, radiating helical flux tubes with interesting binding and interaction properties. It is not yet a theory of **protons** as compact mechanical objects with the observed internal pressure distribution.

The user's intuition ("multi-dimensional dynamics that create the hard shell and corresponding field forces, with a Higgs 'hole' maintained polar opposite the shell") is reasonable but requires a richer vacuum structure and topological protection than the current 6-field equations naturally provide. The v56 GA/multivector direction was an attempt to supply exactly that richer structure.

It is still possible that a clever addition *within* the 6-field framework (new potential, density-dependent coefficients, higher-derivative terms, or a fundamentally better initialization/relaxation protocol) could produce the required pressure profile. However, after the exhaustive searches of V44–V54 and the clean recreation just completed, the burden of proof has shifted.

---

## 8. What Has Been Tried (Summary)

- 3-field gravity-only models (V33 and earlier)
- 6-field Cosserat with V(P) saturation
- Phase-offset confinement for 3-braid composites (V41–V43)
- Evolutionary search for stable seeds (V40)
- Cosserat strain (α) and curl²-hardening (β) terms (V44+)
- Variable theta mass, Floquet averaging, stochastic treatments
- Voronoi foam discretization (V55) — ruled out grid anisotropy as the cause
- Full pivot to ℝ(3,1,1) multivector field + Higgs + Skyrme L4 + hard S³ projection (V56)
- Hard vs. soft constraints on |q|=1
- Pion mass terms (destabilizing)
- Proposals for resonance locking, density-dependent effective mass, frequency self-confinement

None have produced a compact object whose radial pressure profile matches the lattice data.

---

## 9. Possible Future Directions

### 9.1 Stay Inside 6-Field (High-Risk, Low-Cost)
- Systematic search for a potential or additional term that enforces the observed p(r) sign flip.
- Use the lattice ε(r) profile itself as a target for seed generation or as a soft constraint during relaxation.
- Detailed study of the 3-braid composite under the new mechanical diagnostics (does binding ever produce the required shell?).
- Absorbing boundary + controlled radiation studies to see if any 6-field object can reach true mechanical equilibrium.

### 9.2 Controlled Hybrid
- Keep the 6-field sector for the "gluonic" extended part and add a minimal additional field or constraint that supplies the compact mechanical core.
- Explore whether the theta sector can be modified (while remaining 6-field) to provide the negative pressure shell.

### 9.3 Full Commitment to Richer Algebra (v56 direction)
- Treat the lattice mechanical profile as the **primary scoring function** for all new terms and parameters in the GA framework.
- Focus on mechanisms that naturally produce a dense shell with a maintained depleted region ("Higgs hole") whose position is dynamically stabilized.
- Systematic study of the radial pressure profile as a function of Skyrme coefficient, projection strength, pion mass, and Higgs parameters.
- Development of a proper discrete stress-energy tensor and simulated GFF extraction for direct comparison.

### 9.4 Infrastructure
- Finish and harden `v57/tools/analyze_6field_mechanics.c` (and its multivector counterpart) so that every new simulation is automatically scored against the lattice table.
- Add the ability to compute simulated gravitational form factors from the discrete EMT.
- Create a library of reference lattice profiles (energy, pressure, shear, flavor decomposition) for automated comparison.
- Formalize the von Laue condition + D-term requirement in Lean 4 as a necessary condition for any viable Lagrangian.

### 9.5 Seed and Initialization
- Move away from analytical or Gaussian seeds toward data-informed radial templates derived from the lattice ε(r).
- Explore whether a slow "adiabatic" relaxation or cooling phase (varying parameters during the run) can help objects reach the mechanical equilibrium that pure dynamics do not find.

---

## 10. Open Questions

1. Can any purely 6-field configuration (single or multi-braid) ever develop a negative-pressure shell of the magnitude and radial location required by the lattice D-term?
2. Is the fragmentation seen at late time in the recent run a generic feature of 6-field braids, or can it be suppressed by better initialization or additional terms?
3. In the GA framework, which components of the 8-component multivector naturally map to the more compact "quark" contribution versus the extended "gluon" contribution seen in the lattice decomposition?
4. What is the minimal addition (if any) to the current 6-field Lagrangian that can produce the observed pressure sign flip while preserving all existing qualitative successes?
5. How should we weight the various observables (mechanical profile, topology survival, charge-dependent forces, gravitational response, deuterium binding) when deciding whether a direction is promising?

---

## 11. Recommended Immediate Next Steps

1. **Analysis of the just-completed run** — Run the full `analyze_6field_mechanics` (or enhanced cluster_profile) pipeline on the new `v55_braid_recreate.sfa` at multiple late times and produce a direct TSV comparison to the lattice reference table.
2. **Launch L=40 version** — Repeat the recreation at N=192, L=40 to check resolution dependence of the fragmentation and radial profiles.
3. **Historical 3-braid analysis** — Locate or re-create the best V41–V43 3-braid composites and subject them to the same mechanical profiling. Does explicit phase confinement ever produce the required shell?
4. **Document the mechanical target** — Make the v57 reference table and comparison code the standard first diagnostic for every new simulation in the project.
5. **Decision point** — After the above three analyses, decide whether to invest further in 6-field refinements or to treat the GA/multivector direction as the primary path forward, using the lattice pressure profile as the main optimization target.

---

This document captures the state of thinking in May 2026. The data from the clean 6-field recreation, combined with the now-quantitative lattice mechanical target, represents a significant clarification of the problem. The project has excellent tools, a clear quantitative benchmark, and two well-defined paths forward. The choice between them should be driven by the numbers, not by prior attachment to any particular formulation.