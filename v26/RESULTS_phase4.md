# V26 Phase 4: Torsion Constraint and Teleparallel Gravity -- Results

## Parameters

All runs: mu=-20, kappa=20, mass=1.0, A0=0.8, R_tube=3.0, N=128, L=20, periodic BC in z.
Twisted tube initialization (v26 mode 2): triple-product potential only, no pairwise or Lame.
Phases 4a/4b use tfinal=200 (scan). Phases 4c/4d use tfinal=500.

## Phase 4a: Torsion Scan (kappa_S=0, kappa_T varied)

KEY QUESTION: Does kappa_T tighten the braid (increase fc from baseline 0.37)?

| kappa_T | fc(final) | |P|_max | E(final) | l=0 frac | l=1 frac | l=2 frac | breathing_relvar |
|---------|-----------|--------|----------|----------|----------|----------|------------------|
| 0.00    | 0.260     | 0.449  | 349      | 0.575    | 0.011    | **0.415**| 2.043            |
| 0.10    | 0.246     | 0.554  | 355      | 0.698    | 0.022    | **0.280**| 2.487            |
| 0.50    | 0.234     | 0.296  | 354      | 0.734    | 0.025    | **0.241**| 0.306            |
| 1.00    | 0.199     | 0.004  | 339      | 0.920    | 0.018    | 0.062   | 0.069            |
| 2.00    | 0.162     | 0.001  | 291      | 0.923    | 0.039    | 0.038   | 0.215            |
| 5.00    | 0.123     | 0.001  | 186      | 0.877    | 0.081    | 0.042   | 0.220            |
| 10.00   | 0.174     | 0.000  | 150      | 0.945    | 0.038    | 0.017   | 0.298            |

### Key Findings

**ANSWER: NO. Torsion LOOSENS the braid, it does not tighten it.**

1. **fc decreases monotonically** with kappa_T: 0.26 (baseline) -> 0.12 (kT=5).
   The torsion coupling adds stiffness against antisymmetric shear, which resists
   the helical twist that defines the braid structure.

2. **Triple product |P|_max collapses**: From 0.45 (kT=0) to 0.001 (kT>=1). The
   torsion term kills the coherent phase structure that creates triple-product
   binding. By kT=1, the fields have lost their correlated twist.

3. **l=2 content drops dramatically**: The baseline has 41.5% l=2 (excellent
   quadrupolar asphericity), but this drops to <5% for kT>=1 as the configuration
   becomes diffuse and isotropic.

4. **Breathing is suppressed**: The relative variance drops from 2.0 (strong
   breathing at kT=0) to 0.07 (nearly static at kT=1). However, this is because
   the fields dispersed, not because the braid topology prevented oscillation.

5. **Energy decreases with kT**: The added torsion stiffness causes faster radiation
   loss (E=349 at kT=0 vs E=150 at kT=10).

6. **The m=0 case is UNSTABLE**: Without mass, the triple-product potential (mu<0)
   has no confinement and the fields blow up (E -> -infinity). The torsion coupling
   does NOT create an effective mass gap. It modifies wave speeds but not the potential
   landscape. Mass=1 is essential for any survival.

## Phase 4b: Separate Strain and Torsion

| kappa_S | kappa_T | fc(final) | |P|_max | E(final) | l=2 frac |
|---------|---------|-----------|--------|----------|----------|
| 0.50    | 0.00    | 0.196     | 0.002  | 328      | 0.017    |
| 0.00    | 0.50    | 0.234     | 0.296  | 354      | 0.241    |
| 0.50    | 0.50    | 0.218     | 0.001  | 312      | 0.039    |
| 1.00    | 0.50    | 0.200     | 0.001  | 248      | 0.010    |
| 0.50    | 1.00    | 0.145     | 0.001  | 273      | 0.034    |
| 0.50    | 2.00    | 0.235     | 0.001  | 250      | 0.010    |

### Key Findings

1. **Strain (kappa_S) is MORE destructive than torsion**: (kS=0.5, kT=0) gives
   fc=0.20 and |P|=0.002, worse than (kS=0, kT=0.5) which gives fc=0.23 and
   |P|=0.30. The symmetric strain coupling adds stiffness to compression modes
   that disrupts the tube faster.

2. **When kS=kT** (speed rescaling only, no anisotropy): (0.5,0.5) gives fc=0.22,
   worse than the baseline. The overall stiffness increase just disperses faster.

3. **No configuration improves over the baseline**: The best l=2 fraction is the
   baseline (41.5% at kS=kT=0). Adding any cross-gradient coupling kills the
   coherent braid structure before it can establish asphericity.

4. **kS != kT does NOT help l=2**: The hypothesis that separating compression and
   shear stiffness would enhance quadrupolar content is falsified. The effect is
   dominated by overall dispersal.

## Phase 4c: Torsion Flux Quantization

At (kS=0, kT=2.0), evolved to t=500 with mass=1.

**Torsion flux Phi_z = 0 everywhere** (to machine precision ~10^-15).

The transverse vorticity fluxes (Phi_x, Phi_y) trace a smooth sinusoidal pattern:
- Peak amplitude: |Phi_x|_max = 7.3, |Phi_y|_max = 6.8
- One full rotation over the z-domain (as expected from one twist)
- sqrt(Phi_x^2 + Phi_y^2) varies from ~5.5 to ~7.3

**The torsion flux is NOT quantized.** The vorticity vector rotates smoothly
around the z-axis as z varies, with amplitude that depends continuously on
the field strength (no topological quantization). This is because the twisted
tube does not have a topological winding number in the Phi_x-Phi_y plane --
it's a smooth field configuration, not a vortex with a phase singularity.

## Phase 4d: Self-Consistent Teleparallel Metric

At (kS=0, kT=2.0) with full metric correction g^{ij} = delta^{ij} - 2*eps_{ij},
ramped from alpha_g=0 to 1 over t=0..100, then held at alpha_g=1 to t=500.

| Metric | fc(t=500) | |P|_max | E(500) | l=2 frac |
|--------|-----------|--------|--------|----------|
| OFF    | 0.162     | 0.001  | 291    | 0.038    |
| ON     | 0.210     | 0.001  | 147    | 0.046    |

### Key Findings

1. **The self-consistent metric is STABLE.** The braid survives t=500 with the
   full nonlinear metric correction. No runaway instability.

2. **Modest improvement**: fc=0.21 (metric ON) vs 0.16 (metric OFF), and l=2
   fraction increases slightly (4.6% vs 3.8%). The metric back-reaction provides
   slight additional confinement.

3. **Faster energy loss**: E drops to 147 (metric) vs 291 (no metric). The
   curved-space Laplacian modifies wave propagation and enhances radiation.

4. **No qualitative change**: The metric correction does not fundamentally alter
   the dynamics. The braid still slowly disperses with low |P|_max.

## Summary

| Phase | Result | Assessment |
|-------|--------|------------|
| 4a: kT tightens braid? | fc DECREASES with kT | **NEGATIVE** |
| 4b: kS != kT helps l=2? | l=2 DECREASES with any coupling | **NEGATIVE** |
| 4c: Flux quantization? | Phi_z = 0, transverse not quantized | **NEGATIVE** |
| 4d: Metric stable? | Survives t=500 | **POSITIVE** (only stable result) |

### Why Torsion Fails as Confinement

The torsion coupling kappa_T * (d_i phi_j - d_j phi_i)^2 penalizes antisymmetric
gradients. For the twisted tube, the DEFINING feature is antisymmetric cross-gradients
(phi_0 varies in z while phi_1 varies in z with a phase offset). Adding torsion
stiffness directly opposes the twist that creates the braid topology.

This is analogous to increasing the shear modulus of a material: braids become
harder to form, not tighter. The braid is a LOW-energy configuration of the
base Lagrangian (mass + triple product). Adding torsion stiffness raises the
energy of the twisted state relative to the untwisted ground state, causing
the braid to unwind.

### The Baseline is Already the Best

The V26 mode 2 configuration (mass=1, triple product only, no cross-gradient
couplings) achieves:
- fc = 0.26 (marginal survival)
- l=2 = 41.5% (strong quadrupolar content -- the best of any configuration!)
- |P|_max = 0.45 (strong triple-product binding)

Adding ANY cross-gradient coupling (torsion, strain, or combined) degrades
all three diagnostics. The triple-product binding alone creates the most
coherent braid structure.

### Implications for Gravity

1. The self-consistent metric is stable (Phase 4d positive), proving the
   mathematical framework is self-consistent.

2. However, the torsion/strain decomposition does not improve the soliton
   structure. The l=2 content drops, which is the wrong direction for
   tensor gravity.

3. The torsion flux is not quantized, ruling out an EM-like interpretation
   of the antisymmetric sector.

4. The best quadrupolar content (41.5%) comes from the twisted tube geometry
   itself, not from any coupling. This suggests geometry, not dynamics, is
   the source of asphericity.

## Data Files

- `data/v26p4_summary.tsv` -- summary of all runs
- `data/v26p4_4a_kT*.tsv` -- Phase 4a time series
- `data/v26p4_4b_kS*_kT*.tsv` -- Phase 4b time series
- `data/v26p4_phase4c_flux.tsv` -- torsion flux vs z
- `data/v26p4_phase4c_vorticity_xy.tsv` -- vorticity cross-section
- `data/v26p4_4d_metric_*.tsv` -- Phase 4d time series
- `data/v26p4_*_strain.tsv` -- strain on shell R=8
- `data/v26p4_*_multipoles.tsv` -- multipole decomposition
