# V57 — Lattice Proton Mechanical Structure as Quantitative Target

**Date**: 2026-05
**Status**: New direction — planning + infrastructure
**Motivation**: The 6-field Cosserat (and subsequent multivector GA) models produce semi-stable localized objects (the z-aligned "braid" of V34–V55), but have not produced objects whose *internal mechanical structure* matches the real proton. Stability under periodic BC or "B winding lasts until grid collapse" is necessary but insufficient.

The 2024 lattice QCD results (Hackett, Pefkou, Shanahan, arXiv:2310.08484) now give first-principles, flavor-decomposed, spatially resolved mechanical observables for the proton:
- Energy density ε(r)
- Pressure p(r) (positive repulsive core + negative confining shell)
- Shear s(r)
- D-term (Druck term)
- Mass radius and mechanical radius (mechanical < charge radius)

These are derived from the gravitational form factors (GFFs) A(t), J(t), D(t) of the QCD energy-momentum tensor. They constitute a concrete, falsifiable target for any "particle from field" construction.

## Reference Object from Previous Work

The v55/L40 snapshots (`v55/snapshots/L40/`) show the current best 6-field Cosserat object on the Voronoi foam: the single z-aligned braid (A=0.8, R_tube=3.0, phase offsets from v28 optimization, recreated with α=β=0 to match V34).

- It is semi-stable for T≈50 under periodic BC.
- It continuously sources and propagates an extended theta halo (the "propagates theta while moving" behavior).
- Energy drift is bounded after the discrete curl fix.
- This object (and its 3-braid composites) was the basis for all charge-dependent force, gradient response, and deuterium binding results in V34–V43.

**Limitation**: This braid is a 1D-like helical current loop with an extended theta "EM" field. It does not exhibit the 3D isotropic pressure sign-flip (positive core, negative shell) required for mechanical stability in the lattice proton. The v56 Skyrme+Higgs attempts on the foam reproduced the topology problem (lattice collapse when the soliton shrinks) but still lacked a calibrated mechanical equation of state.

## Quantitative Targets (Hackett et al. 2024)

All values in MS-bar scheme at μ=2 GeV unless noted.

### Forward limits (Table 1, dipole fits)
- A(0) = 1 (momentum sum rule, satisfied)
- A_g(0) ≈ 0.501(27)   — gluon momentum fraction
- A_q(0) ≈ 0.510(25)   — total quark (u+d+s)
- J(0) = 1/2 (spin sum rule)
- J_g(0) ≈ 0.255(13)
- **D(0) ≈ −3.87(97) GeV⁻²**  (critical stability diagnostic)
  - D_g(0) ≈ −2.57(84)
  - D_q(0) ≈ −1.30(49)

### Radial structure (Breit frame, from GFF Fourier transforms — Figure 3)
- Energy density ε(r): positive, centrally peaked, falls monotonically. Gluon contribution spatially more extended than quark (larger gluonic mass radius).
- Pressure / longitudinal force density: **positive (repulsive) for small r (core)**, crosses zero, **negative (inward tension) at larger r (shell)**. The sign change satisfies the von Laue stability condition ∫ p(r) 4π r² dr = 0.
- Mechanical radius (from force density) < charge radius (~0.84 fm PDG). Gluonic contributions extend further.

These profiles are the new success metric:
A candidate soliton (braid, UUD composite, multivector hedgehog, or future object) "succeeds" when its discrete stress-energy tensor, when radially binned in the object's rest frame, reproduces ε(r), p(r), the D-term, and the two radii within uncertainties.

## Proposed Work

### Phase 1 — Infrastructure (this version)
- [x] v57/FUTURE.md (this file)
- [ ] C generator (`tools/gen_proton_profiles.c`) that produces a high-resolution reference table `data/proton_lattice_profiles.tsv` (r [fm], ε_total, p_total, ε_g, p_g, force density, etc.) using the dipole parametrization of the GFFs with parameters tuned to published radii and D(0).
- [ ] Simple C analyzer (later to move to `sfa/analysis/`) that, given an SFA (or cell-native FMSH/FCEL) frame of a relaxed object:
  - Computes the discrete EMT from the fields (φ,θ or multivector M + derivatives + potential + Skyrme terms).
  - Bins radially around the tracked centroid.
  - Outputs matching columns for direct comparison.
- [ ] Comparison metrics (χ² on profiles, extracted radii, D-term analog, location of p(r)=0 crossing, core pressure max, etc.).

### Phase 2 — Application to existing objects
- Run the analyzer on the v55 L40 braid snapshots (and earlier V34/V41/V43 proton composites).
- Quantify how far the current best objects are from the lattice mechanical structure (especially the missing negative pressure shell).
- Document the gap in a RESULTS.md.

### Phase 3 — Theory guidance
- Use the target p(r) sign-flip and ε(r) shape to design or constrain new terms in the Lagrangian (density-dependent coefficients, higher-order Skyrme-like invariants, resonance-locking or variable-mass terms from v56/PLAN_NEXT.md, etc.).
- In the ℝ(3,1,1) multivector framework: decide which components (rotor vs. boost/pseudoscalar) should map to the more compact "quark" vs. extended "gluon" contributions.
- Re-examine the role of hard |q|=1 projection vs. soft sigma constraint in light of the required mechanical stiffness.

### Phase 4 — Seed improvement
- Use the lattice ε(r) profile (or the corresponding |M| or |q| shape) as a radial template for improved initial conditions (better than pure hedgehog or analytical braid sums).
- Goal: start closer to mechanical equilibrium so that the subsequent evolution does not immediately drive the soliton into the lattice-collapse regime.

## Relation to Existing Documents

- Updates and extends `FUTURE.md` (root) items on "Proton mass decomposition analog (cf. Ji decomposition)".
- Directly addresses the stability and collapse issues documented in `v56/PION_AND_LATTICE_COLLAPSE.md`, `v56/RESULTS_SKYRME_V2.md`, and `v56/SCAN_RESULTS.md`.
- The v55 braid (this directory's reference object) is the calibration baseline.

## Success Criteria for v57

1. Reproducible C code + committed numerical table that any later analysis can reference without external data.
2. First quantitative comparison (numbers, not just "it looks stable") between a simulated object and the lattice proton mechanical structure.
3. Clear statement in a future RESULTS.md of what the current best objects get right and what they are missing mechanically.

All new code follows CLAUDE.md: C with `-lzstd`, SFA I/O where possible, no modification of `sfa/sim/scp_sim.*` or `sfa/format/sfa.h` unless explicitly authorized.

---

**Next concrete steps after the table exists**:
- Implement the SFA-side EMT radial profiler.
- Run it on the v55 braid.
- Decide whether the GA multivector direction needs a mechanical "pressure term" tuned to the lattice D-term and p(r) profile before further topology work.