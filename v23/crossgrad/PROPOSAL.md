# V23-B: Cross-Gradient Oscillon Dynamics

## Thesis

Adding a cross-gradient coupling η(∂_i φ_j)(∂_j φ_i) to the three-field
oscillon Lagrangian distinguishes compression from shear in the gradient
sector. In vacuum this only changes wave speeds (c_L² = 1+η, c_T² = 1).
But around the oscillon, the cross-gradient AMPLIFIES the Hessian splitting
between symmetric (compression) and antisymmetric (shear) perturbation
sectors. At some critical η_c, the antisymmetric effective gap may close,
creating a near-massless shear resonance.

This investigation tests whether oscillons survive the addition of cross-
gradient coupling, and scans for the critical η_c where shear softening
becomes maximal.

## Mathematical Setup

### Modified Lagrangian

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)(∂_i φ_a)]
        - ½η(∂_i φ_j)(∂_j φ_i)
        - ½m²Σ_a φ_a²
        - V(φ₁φ₂φ₃)

where V(P) = (μ/2)P²/(1+κP²) as before.

### Modified Equation of Motion

    ∂²φ_a/∂t² = ∇²φ_a + η·∂_a(∇·φ) - m²φ_a - ∂V/∂φ_a

The new term η·∂_a(∇·φ) = η·∂_a(Σ_b ∂_b φ_b) couples field a to the
divergence of all fields. Note: this requires identifying field index a
with spatial direction a (the elastic interpretation).

In component form for the 1D case (x-direction only, all fields depend
on x):

    ∂²φ_a/∂t² = ∂²φ_a/∂x² + η·∂²φ_a/∂x²·δ_{a,x-component}
                 - m²φ_a - ∂V/∂φ_a

Wait — in 1D with x as the only spatial dimension, ∂_a(∇·φ) = ∂_a(∂_x φ_x).
So:
- For a=1 (x-direction): extra term = η·∂²φ₁/∂x²  (compression enhanced)
- For a=2,3 (y,z-directions): extra term = 0  (shear unchanged)

This means in 1D:
    ∂²φ₁/∂t² = (1+η)∂²φ₁/∂x² - m²φ₁ - ∂V/∂φ₁
    ∂²φ₂/∂t² = ∂²φ₂/∂x²     - m²φ₂ - ∂V/∂φ₂
    ∂²φ₃/∂t² = ∂²φ₃/∂x²     - m²φ₃ - ∂V/∂φ₃

In 1D, the cross-gradient only affects the field aligned with the spatial
axis. For the SYMMETRIC oscillon (φ₁=φ₂=φ₃), the equation becomes:

    ∂²f/∂t² = (1+η/3)∂²f/∂x² - m²f - ∂V_sym/∂f

(averaging over the three fields' roles). This just rescales the effective
wave speed for the symmetric mode.

### 3D Equations

In 3D, the full cross-gradient term:

    η·∂_a(∂_b φ_b) at grid point (i,j,k):

For a=1 (x):  η·∂_x(∂_x φ₁ + ∂_y φ₂ + ∂_z φ₃)
For a=2 (y):  η·∂_y(∂_x φ₁ + ∂_y φ₂ + ∂_z φ₃)
For a=3 (z):  η·∂_z(∂_x φ₁ + ∂_y φ₂ + ∂_z φ₃)

This requires second-order mixed partial derivatives (∂_x∂_y, etc.).
Discretize with centered differences:

    ∂_x∂_y φ₂ ≈ [φ₂(i+1,j+1,k) - φ₂(i+1,j-1,k)
                 - φ₂(i-1,j+1,k) + φ₂(i-1,j-1,k)] / (4·dx²)

### What to Test

**Phase 1 (1D, fast)**:
1. Start with the v21 1D triad code (triad1d.c).
2. Add the cross-gradient term (only affects φ₁ in 1D x-direction).
3. Initialize symmetric oscillon (φ₁=φ₂=φ₃) with Gaussian.
4. Scan η from 0 to 2.0 in steps of 0.1.
5. For each η: measure oscillon lifetime, frequency ω, energy, stability.
6. Find the maximum η for which the oscillon survives.

**Phase 2 (1D perturbation spectrum)**:
7. At each stable η: compute the linearized perturbation spectrum.
8. Decompose into symmetric (all fields in phase) and antisymmetric modes.
9. Measure the effective gap for each sector.
10. Plot gap vs η. Find η_c where antisymmetric gap closes.

**Phase 3 (3D, if Phase 1-2 are promising)**:
11. Add cross-gradient to the v21 3D triad code (triad3d.c).
12. Run oscillon at the most promising η value.
13. Measure TT strain field in the far zone.
14. Check if strain falls as 1/r or e^{-mr}/r.

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c` (base for modification)
- v21 3D solver: `/home/d/code/scp/v21/src/triad3d.c` (for Phase 3)
- v22 two-oscillon: `/home/d/code/scp/v22/src/two_oscillon.c` (for force measurement)

## Output Structure

- `src/crossgrad1d.c` — 1D oscillon with cross-gradient (Phase 1-2)
- `src/crossgrad3d.c` — 3D version (Phase 3, if needed)
- `data/` — output data files (time series, profiles, spectra per η)
- `RESULTS.md` — results and analysis

## Parameters

Base: μ=-20, κ=20, m=1.0 (v21 production)
Scan: η ∈ [0.0, 2.0] in steps of 0.1
Grid: Nx=2000, xmax=60, tfinal=5000 (1D)

Compile: `gcc -O3 -Wall -o crossgrad1d src/crossgrad1d.c -lm`
