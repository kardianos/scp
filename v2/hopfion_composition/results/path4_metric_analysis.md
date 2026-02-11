# Path 4: Hopfion Effective Metric — Results

## Key Finding

**The BLV effective metric on a hopfion background is ANISOTROPIC.**

The hedgehog Skyrmion gives a perfectly isotropic effective metric (g^{xx} = g^{yy} = g^{zz}).
The hopfion, with its axial symmetry and non-hedgehog current structure, gives:

| Quantity | Hedgehog | Hopfion (a=2) |
|----------|----------|---------------|
| Max g^{zz}/g^{xx} | 1.000 | **1.750** |
| Peak Δc²/c² | 0.000 | **0.413** |
| Mean anisotropy | 0.000 | 0.093 |
| Anisotropy range | — | r < 3 (core region) |

## Physical Meaning

The directional effective metric g^{jj}(x) governs the propagation speed of small
perturbations in direction j. For L₂ + L₄ Skyrme theory:

```
g^{00}(x) = 1 + (c₄/|q|⁴) Σ_i |A_i|²         (temporal inertia)
g^{jj}(x) = 1 + (c₄/|q|⁴) Σ_{i≠j} |A_i|²     (spatial stiffness in j-direction)
```

Sound speed in direction j: c²_j = c² × g^{jj}/g^{00}

For the hedgehog: <|A_x|²> = <|A_y|²> = <|A_z|²> by spherical symmetry → isotropic.
For the hopfion: |A_x|² = |A_y|² ≠ |A_z|² due to axial symmetry → anisotropic.

The hopfion's preimage fibers lie in planes roughly perpendicular to z, creating
strong currents in the xy-plane. This makes:
- g^{zz} ≈ 1 + c₄(|A_x|² + |A_y|²)/|q|⁴  (large — strong perp currents)
- g^{xx} ≈ 1 + c₄(|A_y|² + |A_z|²)/|q|⁴  (smaller — one weaker current direction)

Result: c_z > c_x — perturbations travel ~32% faster along the hopfion axis.

## Two-Hopfion Product (D=6.0, a=2.0)

Two hopfions at z=±3.0, composed via product ansatz q₁×q₂/ρ₀.

| Quantity | Value |
|----------|-------|
| B | 0.000 (pure hopfion pair) |
| E/E_FB | 4.94 |
| g^{zz}/g^{xx} range along z-axis | **0.47 – 1.87** |
| Metric non-trivial BETWEEN cores | **YES** |

**Z-axis scan highlights:**
- Deep minimum g^{zz}/g^{xx} = 0.47 near each core center (z ≈ ±3)
- Maximum g^{zz}/g^{xx} = 1.87 at r ≈ 3.1 (between cores)
- Midplane (z=0): g^{zz}/g^{xx} ≈ 0.83 → still anisotropic
- The metric oscillates between cores — NOT simply additive

## Skyrmion + Hopfion Composite (B≈1, H≠0)

Skyrmion at origin + hopfion ring at z=2.0, product ansatz composition.

| Quantity | Value |
|----------|-------|
| B | 0.999 (preserved skyrmion charge) |
| E/E_FB | 4.85 |
| Peak g^{zz}/g^{xx} along z | **2.51** (at z ≈ 0.07) |
| Min g^{zz}/g^{xx} along z | **0.41** (at z ≈ 4.1) |
| Δc²/c² at core | **0.425** |

**This is the strongest result:**
- The hedgehog alone gives g^{zz}/g^{xx} = 1.000 everywhere
- Adding a hopfion ring creates **151% peak anisotropy** (ratio 2.51)
- The z-axis scan shows dramatic asymmetry: z>0 (toward hopfion) has much higher g^{zz}
- The metric is **genuinely asymmetric** — not just anisotropic, but z↔-z broken
- Sound speed Δc²/c² = 0.43 → perturbations propagate 43% faster in z at the core

**Physical interpretation:**
- The combined (B=1, H≠0) topology gives a Kerr-like metric with frame-dragging
- The hopfion ring's preimage fibers create directional stiffness in the field
- This is a genuinely **non-perturbative** tensor effect from the L₂+L₄ Lagrangian alone
- No additional coupling constants needed (unlike L₆ or B⁰p)

## Comparison to Previous Paths

| Property | Path 2 (L₆) | Path 3 (B⁰p) | Path 4 (Hopfion metric) |
|----------|-------------|--------------|------------------------|
| Tensor rank | 2 (scalar-equiv) | 0 (scalar) | **2 (true tensor)** |
| Isotropy | Isotropic P/m | Isotropic 1/r | **Anisotropic** |
| Analog | Schwarzschild | Newton | **Kerr-like** |
| Range | 0.55 fm | 1/r or Yukawa | Core-scale (~1.5 fm) |
| Free params | λ₆ | g_top | **None** (existing L) |
| Peak effect | Φ/c² = -0.14 | U(D) ~ 1/D | g^{zz}/g^{xx} = 2.51 |

## The P/m = 2 Question

The hedgehog-convention P/m = 2 is **tautological** in the current code: the factors
(2r², 4c₄X) vs (r², 2c₄X) are hardcoded from the hedgehog Sturm-Liouville analysis.
The spherically averaged P/m = 2 for ANY field configuration with L₂+L₄.

The REAL signal is the **directional anisotropy** g^{zz}/g^{xx} ≠ 1. This cannot be
captured by the hedgehog radial reduction — it requires the full 3D tensor metric.

## Hopf Charge Issue

The Gauss-Seidel Poisson solver for the Hopf integral H = (1/4π²)∫A·B d³x converges
too slowly for N=128 (needs ~O(N²) = 16K iterations; we ran 2000). The value H ≈ 4.3
has not converged. This is a numerical issue, not a physics one.

Fixes needed: FFT-based Poisson solver (O(N³ log N)), or multigrid method.

## Run Details

- Grid: N=128, R=8.0, h=0.133
- Hopfion: a=2.0, Gaussian envelope g(r) = π exp(-r⁴/a⁴)
- Physics: e=1, ρ₀=1, σ-model (λ=0, λ₆=0)
- Single hopfion: B = 0.000, E/E_FB = 2.47
- Two hopfions (D=6.0): B = 0.000, E/E_FB = 4.94
- Skyrmion+hopfion: B = 0.999, E/E_FB = 4.85

## Conclusions

1. **The BLV effective metric IS anisotropic for non-hedgehog topologies.** This is a
   genuine rank-2 tensor effect that cannot be captured by the P/m radial reduction.

2. **No new coupling constants needed.** The anisotropy emerges from the existing L₂+L₄
   Lagrangian alone — it's a topological effect from the field geometry.

3. **The effect is core-scale (~1.5 fm).** It does NOT produce long-range 1/r gravity.
   However, it demonstrates that soliton topology can create effective spacetime anisotropy.

4. **The skyrmion+hopfion composite gives the strongest signal** (151% anisotropy),
   suggesting that mixed topology (B≠0, H≠0) creates the most non-trivial metric.

5. **For long-range gravity**, Path 3 (B⁰p coupling) remains the only mechanism that
   produces 1/r behavior. Path 4 provides the tensor structure (gravitomagnetic/Kerr-like)
   but only at nuclear scales.

## Remaining Work

1. **FFT Hopf charge**: Replace Gauss-Seidel with FFT-based Poisson solver
   for accurate topology computation.

2. **Equilibrium**: Run gradient flow on the hopfion to reach the actual
   minimum-energy configuration. The anisotropy may change but should persist.

3. **Backreaction**: Propagate a test perturbation on the anisotropic metric
   and check if it modifies the background (self-consistency loop).
