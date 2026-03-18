# T10B Results: Effective Metric Extraction

## Method
Linearized the EOM around the frozen (equilibrated) braid background to extract
the 3x3 Hessian matrix M²_ab(x) = m²δ_ab + ∂²V/∂φ_a∂φ_b at each spatial point.

The eigenvalues of M² determine the local effective mass for each mode.
The refractive index n = sqrt(k²+M²)/sqrt(k²+m²_vac) tells us whether
waves are slower (n>1, attractive lensing = gravity-like) or faster (n<1, repulsive).

## Key Results

### Instantaneous snapshot (frozen at T=200)

| r     | M²_min | M²_mid | M²_max | n_min  | n_max  | Anisotropy |
|-------|--------|--------|--------|--------|--------|------------|
| 0.25  | 2.74   | 2.78   | 3.88   | 1.013  | 1.044  | 0.41       |
| 2.25  | 2.83   | 2.95   | 3.45   | 1.016  | 1.032  | 0.21       |
| 3.25  | 1.86   | 3.00   | 3.19   | 0.989  | 1.026  | 0.44       |
| 4.25  | 1.38   | 2.69   | 2.75   | 0.976  | 1.014  | 0.51       |
| 5.25  | 1.49   | 2.44   | 2.47   | 0.979  | 1.006  | 0.40       |
| 6.25  | 2.15   | 2.27   | 2.27   | 0.997  | 1.001  | 0.05       |
| 8+    | 2.25   | 2.25   | 2.25   | 1.000  | 1.000  | 0.00       |

### Time-averaged (22 samples over T=20)

| r     | <M²_min> | <M²_mid> | <M²_max> | <n_min> | <n_max> |
|-------|----------|----------|----------|---------|---------|
| 0.33  | 1.90     | 2.59     | 3.02     | 0.990   | 1.021   |
| 2.33  | 1.89     | 2.63     | 3.10     | 0.990   | 1.023   |
| 4.33  | 1.45     | 2.52     | 2.60     | 0.978   | 1.010   |
| 5.67  | 2.08     | 2.28     | 2.29     | 0.995   | 1.001   |
| 7+    | 2.25     | 2.25     | 2.25     | 1.000   | 1.000   |

### At braid center (x=y=z=0)
- Background: (φ₀,φ₁,φ₂) = (-0.586, 0.578, 0.612), ρ=1.05, P=-0.207
- M² eigenvalues: (2.74, 2.78, 3.88)
- ΔM² = (+0.49, +0.53, +1.63)
- All three modes SLOWER at center (n = 1.013, 1.014, 1.044)
- Anisotropy: 41% (genuine tensor structure, not scalar conformal)

## Interpretation

### Is this gravity? MIXED SIGNAL.

**Positive signs:**
- Near the braid center (r < 2.5), ALL modes are slower → attractive
- The max eigenvalue gives n ≈ 1.04 → 4% speed reduction → gravitational lensing
- Tensor structure (anisotropy) → not just a scalar effect
- Far field returns to vacuum → localized effect

**Negative signs:**
- At r = 3-5, the MINIMUM eigenvalue drops BELOW vacuum (M² ≈ 1.4-1.9)
  → one mode is FASTER there → repulsive for that polarization
- Net effect: attractive for 2/3 of modes, repulsive for 1/3 at intermediate r
- This is BIREFRINGENCE — different polarizations feel different metrics
- Real gravity is universal (all modes bend the same way)
- Range is ~6 code units, decaying exponentially → Yukawa, NOT 1/r
- The "attractive" and "repulsive" zones may partially cancel for unpolarized waves

### The fundamental issue
The triple-product potential V(P) = (μ/2)P²/(1+κP²) creates a 3x3 mass matrix
that has different eigenvalues for different field-space directions. This is
inherent to having 3 distinct scalar fields rather than a single field with
internal symmetry. The Hessian mixes "gravity" (all modes slow) with "birefringence"
(some modes fast, some slow) in a way that's un-gravity-like.

### What would fix this?
For genuine gravity-like behavior, we need M²_ab ∝ δ_ab + (something positive),
i.e., ALL eigenvalues above vacuum EVERYWHERE. This requires:
1. A potential where ∂²V/∂φ_a∂φ_b is positive semi-definite on all backgrounds
2. OR: a different coupling (e.g., derivative coupling) that modifies the
   kinetic term rather than the mass term
3. OR: a fundamentally different mechanism (gauge field, metric backreaction)

## Files
- `data/t10b_v3_radial.tsv` — radial profile of eigenvalues
- `data/t10b_v3_xy_slice.tsv` — 2D xy slice
- `data/t10b_v3_zaxis.tsv` — z-axis profile
- `data/t10b_v3_time_avg.tsv` — time-averaged radial profile
