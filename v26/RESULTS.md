# V26: Braided Solitons — Results (Phases 1-3)

## Parameters

| Parameter | Value |
|-----------|-------|
| mu | -20.0 |
| kappa | 20.0 |
| mass | 1.0 |
| A0 | 0.8 |
| lambda_pw | 0.5 (0 for mode 2) |
| eta | 0.1 (0 for mode 2) |
| lambda_L | 0.1 (0 for mode 2) |
| N | 128 |
| L | 20.0 |
| dx | 0.3150 |
| dt | 0.06299 |
| tfinal | 500 |
| Ngrid | 2.1M |

## Four Configurations

| Mode | Name | Init | Couplings | BC (z) |
|------|------|------|-----------|--------|
| 0 | twisted_tube_full | Helical twist in tube R=3 | All | Periodic |
| 1 | borromean_full | Three orthogonal rings R0=6, r_tube=2 | All | Absorbing |
| 2 | twisted_tube_tripleonly | Same as mode 0 | Triple product only | Periodic |
| 3 | oscillon_control | Spherical Gaussian sigma=3 | All | Absorbing |

---

## Phase 1: Configuration Survival

| Mode | E(0) | E(500) | fc(500) | pk_max(500) | |P|max(500) | Survived? |
|------|------|--------|---------|-------------|------------|-----------|
| 0 (twisted full) | 310.3 | 78.4 | 0.19 | 0.12 | 0.0006 | Marginal (dispersed) |
| 1 (Borromean full) | 745.5 | 235.9 | 0.92 | 0.73 | 0.387 | **YES** |
| 2 (twisted triple-only) | 575.3 | 204.8 | 0.27 | 0.20 | 0.006 | Marginal (dispersed) |
| 3 (oscillon control) | 278.6 | 226.2 | 0.98 | 0.11 | 0.001 | **YES** (oscillon) |

### Key findings

**Borromean rings (mode 1) are the clear winner.** This configuration survives
t=500 with fc=0.92 (92% of energy in core) and large triple product (|P|max=0.387).
The three orthogonal ring tubes create a stable, localized, triply-linked structure
that the dynamics preserves.

**Twisted tubes disperse.** Both mode 0 (full couplings) and mode 2 (triple
product only) lose their structure. The twisted tube topology is one-dimensional
(along z-axis) and lacks the 3D interlocking that stabilizes the Borromean
configuration. Mode 0 disperses faster due to pairwise coupling and mass term
pushing fields apart; mode 2 survives longer but still disperses by t~300.

**Oscillon control (mode 3)** survives as expected with very high fc=0.98,
confirming the code is correct and the Lagrangian supports bound states.

---

## Phase 2: Non-Breathing Verification

| Mode | mean(rho) | rel. variance | omega_rho | omega_phi | Breathing? |
|------|-----------|---------------|-----------|-----------|------------|
| 0 (twisted full) | 0.013 | 0.10 | 0.07 | 0.74 | NO (dispersed) |
| 1 (Borromean full) | 1.71 | 0.056 | **2.74** | **1.37** | **YES** |
| 2 (twisted triple-only) | 0.030 | 0.17 | 0.07 | 1.02 | NO (dispersed) |
| 3 (oscillon control) | 1.73 | 0.054 | **2.74** | **1.37** | **YES** |

### Key findings

**The Borromean soliton ALSO breathes.** It shows the same breathing frequency
as the oscillon: omega_phi = 1.37 (close to mass m=1, as expected for oscillons
where omega < m). The rho harmonic is at 2*omega = 2.74 because energy density
goes as phi^2.

**This is a negative result for the non-breathing hypothesis.** The Borromean
configuration survives but oscillates just like the standard oscillon. The
triple-product coupling creates a bound state, but the braid topology does NOT
prevent breathing. Both configurations have similar relative variance (~5%).

**Twisted tubes show no breathing** because they dispersed — the center density
is near zero, so there is nothing to oscillate.

---

## Phase 3: Strain Multipoles

Strain tensor eps_{ij} = (1/2)(d_i phi_j + d_j phi_i) sampled on a sphere at R=8.
Traceless shear |sigma|^2 decomposed into multipoles l=0, l=1, l=2.

| Mode | l=0 fraction | l=1 fraction | l=2 fraction | |sigma|^2 avg |
|------|-------------|-------------|-------------|-------------|
| 0 (twisted full) | 0.797 | 0.156 | 0.047 | 2.1e-4 |
| 1 (Borromean full) | 0.981 | 0.004 | **0.016** | 3.3e-4 |
| 2 (twisted triple-only) | 0.911 | 0.022 | **0.068** | 5.3e-4 |
| 3 (oscillon control) | **0.993** | 0.007 | **0.001** | 4.3e-5 |

### Key findings

**The oscillon is nearly perfectly spherical** as expected: 99.3% l=0, only 0.07%
l=2. This confirms the multipole analysis works correctly.

**The Borromean soliton has slightly more l=2 content** (1.6%) than the oscillon
(0.07%), a factor ~23x more aspherical. However, this is still dominated by l=0
(98%). The three orthogonal rings have cubic symmetry, which is more spherical
than a generic braid.

**The twisted tube remnants show the most asphericity** — mode 2 has 6.8% l=2
and mode 0 has 15.6% l=1. This is because the dispersed field retains some
cylindrical asymmetry from the initial tube geometry.

**Bottom line**: the braid DOES break spherical symmetry (l=2 is 23x larger than
the oscillon), but only at the few-percent level. The Borromean rings' cubic
symmetry (three mutually orthogonal planes) keeps the configuration approximately
spherical.

---

## Summary and Implications

1. **The Borromean ring configuration is a stable 3D soliton** that survives
   t=500 with 92% energy confinement and large triple product. The three
   interlocking field rings create a topologically non-trivial bound state.

2. **It breathes**, just like the oscillon. The braid topology does NOT suppress
   oscillation. This invalidates the central hypothesis of V26 that braided
   solitons would be static (non-breathing).

3. **The asphericity is modest** (l=2 ~ 1.6%), reflecting the high symmetry
   of the Borromean configuration (cubic point group). A less symmetric braid
   (e.g., trefoil knot on a torus) might show larger l=2, but the twisted tube
   already disperses before we can measure it.

4. **The twisted tube topology is unstable.** The 1D helical structure along the
   z-axis does not provide enough confinement to prevent dispersal, even with
   periodic boundary conditions.

5. **For future work**: try a toroidal braid (Method 2 in PROPOSAL.md) which
   has genuine 3D closure without relying on periodic BC. Also try larger R0_ring
   and smaller r_tube for the Borromean rings to increase the linking strength.

## Data files

- `data/mode{0,1,2,3}/v26_{name}_phase1.tsv` — time series
- `data/mode{0,1,2,3}/v26_{name}_rho_history.tsv` — center density history
- `data/mode{0,1,2,3}/v26_{name}_phase2_dft.tsv` — DFT power spectra
- `data/mode{0,1,2,3}/v26_{name}_phase3_strain.tsv` — strain on shell
- `data/mode{0,1,2,3}/v26_{name}_phase3_multipoles.tsv` — multipole decomposition
