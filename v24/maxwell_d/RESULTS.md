# V24-MD Results: Elastic Dislocation as Gauge Field

## Summary

**DEFINITIVE NEGATIVE.** The three-field oscillon has ZERO dislocation content in both 1D and 2D. All incompatibility measures (Nye tensor, Saint-Venant, Burgers vector) vanish to machine precision for smooth field configurations. This is a theorem: smooth C^2 fields cannot carry dislocations. Only singular field configurations (vortices, domain walls with discontinuous derivatives) can have nonzero dislocation density, and the v24/vortex analysis already showed vortices are unstable.

---

## Phase 1: 1D Strain Analysis

**Code**: `src/maxwell_d.c -phase 1` | Parameters: mu=-10, kappa=10, m=1.0, A=0.8, sigma=3.0

The 1D oscillon equilibrates within ~200 t.u. with peak amplitude ~0.8 and energy E~3.58.

### Strain Profile

Strain eps_{ax} = d_x phi_a for each field a = 1,2,3.

| Quantity | Value | Notes |
|----------|-------|-------|
| Peak |phi_a| | 0.576 | All three equal (symmetric) |
| Peak |eps_{ax}| | 0.177 | All three equal (symmetric) |
| Strain localization | ~6 code lengths | Same as field profile (Gaussian core) |

The strain is a smooth, bell-shaped derivative of the oscillon Gaussian. It has two lobes (positive and negative) flanking the center, with zero at x=0 (peak of phi).

### 1D Dislocation Density

**Identically zero.** In 1D, the dislocation density tensor alpha_{ij} = eps_{ikl} d_k eps_{lj} requires the Levi-Civita symbol with at least 2 spatial indices. With only one spatial direction (x), this is trivially zero for any field configuration.

---

## Phase 2: 2D Strain Topology

**Code**: `src/maxwell_d.c -phase 2` | Nx=Ny=256, L=20, mu=-20, kappa=20, m=1.0

The 2D oscillon is initialized as three identical radially-symmetric Gaussians. It equilibrates into a breathing mode with f_core > 0.99 and peak amplitude up to 0.88.

### 2D Oscillon Stability

| t | E_total | f_core | Peak amp | Status |
|---|---------|--------|----------|--------|
| 0 | 21.53 | 1.000 | 0.80 | Initial |
| 50 | 18.67 | 0.969 | 0.065 | Shedding |
| 100 | 17.08 | 0.986 | 0.88 | Breathing |
| 150 | 16.48 | 0.992 | 0.21 | Stable |
| 200 | 16.16 | 0.996 | 0.61 | Stable breathing |

Energy loss decelerates (from ~0.06/t.u. initially to ~0.006/t.u. at late times). The 2D oscillon survives and enters a stable breathing mode, consistent with v21 findings.

### Nye Dislocation Density

alpha_a = d_x(d_y phi_a) - d_y(d_x phi_a) = curl of the distortion tensor.

| t | max |alpha_1| | max |alpha_2| | max |alpha_3| | max |alpha| |
|---|----------------|----------------|----------------|--------------|
| 0 | ~1e-16 | ~1e-16 | ~1e-16 | 1.56e-16 |
| 50 | ~1e-16 | ~1e-16 | ~1e-16 | 7.21e-17 |
| 100 | ~1e-16 | ~1e-16 | ~1e-16 | 1.68e-16 |
| 150 | ~1e-16 | ~1e-16 | ~1e-16 | 8.41e-17 |
| 200 | ~1e-16 | ~1e-16 | ~1e-16 | 1.68e-16 |

**The Nye tensor is ZERO to double-precision roundoff (~1e-16).** This is the equality of mixed partial derivatives: d_x d_y phi = d_y d_x phi for smooth (C^2) fields.

### Saint-Venant Incompatibility

eta = d_yy(eps_xx) + d_xx(eps_yy) - 2 d_xy(eps_xy), where eps_xx = d_x phi_1, eps_yy = d_y phi_2, eps_xy = (d_y phi_1 + d_x phi_2)/2.

| t | max |eta| | Notes |
|---|-----------|-------|
| 0 | 4.16e-05 | Initial Gaussian — numerically exact |
| 50 | 1.28e-02 | During shedding |
| 100 | 1.89e-02 | Breathing peak |
| 150 | 1.11e-02 | |
| 200 | 9.60e-03 | |

The Saint-Venant scalar is **analytically zero** for smooth fields (proven by expanding: d_yy d_x phi_1 + d_xx d_y phi_2 - d_xy(d_y phi_1 + d_x phi_2) = 0). The O(1e-2) residual is discretization error from second-order finite differences on a dx=0.157 grid. At t=0 (analytic Gaussian), the error is 4e-5; during evolution with steeper gradients, it rises to 1e-2. This is entirely numerical.

### Burgers Vector

b_a = oint eps_{ia} dl_i around a square contour at ~2 sigma from center.

| t | b_1 | b_2 | b_3 |
|---|-----|-----|-----|
| 0 | 6.4e-3 | 6.4e-3 | 6.4e-3 |
| 100 | -7.1e-3 | -7.1e-3 | -7.1e-3 |
| 200 | 4.7e-3 | 4.7e-3 | 4.7e-3 |

All Burgers vectors are O(1e-3), consistent with discretization error on a finite-difference grid. They oscillate with the breathing mode (not growing), and are identical across all three fields (reflecting the symmetric initialization). **No quantized dislocation content.**

---

## Why the Oscillon Has Zero Dislocation Content

### Theorem (smooth fields => no dislocations)

For any C^2 scalar field phi_a(x,y):

    alpha_a = d_x(d_y phi_a) - d_y(d_x phi_a) = 0

This follows from the equality of mixed partial derivatives (Clairaut's theorem). The dislocation density is the failure of integrability of the distortion field — but smooth fields are automatically integrable.

### What WOULD have nonzero alpha?

1. **Vortex configurations**: phi = (cos(n*theta), sin(n*theta)) * f(r) has a singularity at r=0 where d_x d_y phi != d_y d_x phi. The winding number n gives the Burgers vector magnitude. However, v24/vortex showed vortices are unstable in this model (they unwind and disperse).

2. **Domain walls with kinks**: A field configuration with a line of discontinuous slope would have delta-function dislocation density along the line.

3. **Topological defects in ordered media**: In physical crystals, dislocations arise from discrete lattice symmetry breaking. The continuous fields phi_a(x,y) cannot break a discrete symmetry they don't have.

### Fundamental obstruction

The elastic gauge theory (dislocation = gauge field strength, Burgers vector = charge) requires **singular** or **multi-valued** field configurations to produce nonzero field strength. The oscillon, being a smooth, localized, oscillating field configuration, is in the "vacuum sector" of the elastic gauge theory — it has zero field strength everywhere.

This is analogous to how a smooth electromagnetic potential A_mu always gives F_munu = dA (exact), so it cannot carry magnetic monopoles (which require non-exact F). Dislocations are the elastic analog of magnetic monopoles.

---

## Conclusion

The elastic dislocation / emergent EM interpretation does not apply to oscillons:

1. **1D**: Dislocation density is trivially zero (needs >= 2 spatial dimensions)
2. **2D**: Dislocation density (Nye tensor) vanishes to machine precision (~1e-16)
3. **Burgers vector**: Zero to discretization accuracy (~1e-3)
4. **Saint-Venant incompatibility**: Zero analytically; numerical residual ~1e-2 is grid artifact
5. **Root cause**: Smooth C^2 fields have zero curl of the distortion tensor by Clairaut's theorem

The oscillon carries no elastic "charge." Emergent electromagnetism from dislocations requires topological defects (vortices, dislocations) that this model cannot stabilize.

---

## Files

| File | Description |
|------|-------------|
| `src/maxwell_d.c` | Combined Phase 1 + Phase 2 solver |
| `data/phase1_ts.tsv` | 1D time series (phi, peak, energy) |
| `data/strain_1d_t{NNNN}.tsv` | 1D strain profiles at snapshots |
| `data/phase2_ts.tsv` | 2D time series (phi, energy, Nye, SV, Burgers) |
| `data/strain_2d_t{NNNN}.tsv` | 2D strain/incompatibility snapshots |
