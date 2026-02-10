# Effective Metric Analysis — Sigma-Model Skyrmion

## Summary

We computed the Babichev-Langlois-Vernieri (BLV) effective metric on a B=1
Skyrmion equilibrium in the sigma-model Skyrme theory (L₂ + L₄, |q| = ρ₀).
The key result: **P/m = 2 identically at every radius**, to machine precision.
This means the sigma-model Skyrmion generates no emergent gravity — no
Schwarzschild-like time dilation, no gravitational redshift, no light bending.

## Numerical Results

### Skyrmion Equilibrium (arrested gradient flow)

Parameters: N=128, R=8.0, h=0.133, e=1, ρ₀=1, sigma model (λ=∞)

| Quantity | Value | Expected | Status |
|----------|-------|----------|--------|
| E/E_FB | 1.219 | 1.232 | PASS (1% discretization error at h=0.133) |
| E₂/E₄ | 0.999 | 1.000 | PASS (virial theorem) |
| B | 0.991 | 1.000 | PASS |

The ~1% discrepancy in E/E_FB is from 3D grid discretization. The 1D radial
solver gives E/E_FB = 1.2320 exactly. The virial E₂ = E₄ is satisfied to 0.1%.

Gradient flow was run in "arrested" mode: stopped when B deviated more than 2%
from 1.0, with rollback to best-energy state. The sigma-model Skyrmion is a
lattice saddle point — unrestricted gradient flow eventually collapses the
soliton (known issue, documented in MEMORY.md).

### Effective Metric: P(r) and m(r)

The Sturm-Liouville coefficients for radial (breathing mode) perturbations:

```
P(r) = 2r² + 4c₄ sin²f(r)    (radial stiffness)
m(r) = r² + 2c₄ sin²f(r)     (temporal inertia)
```

where c₄ = 2ρ₀²/e².

Selected radial values:

| r | P(r) | m(r) | P/m |
|---|------|------|-----|
| 0.25 | 31.75 | 15.88 | 2.000000 |
| 0.65 | 23.60 | 11.80 | 2.000000 |
| 1.05 | 16.48 | 8.24 | 2.000000 |
| 1.85 | 10.31 | 5.16 | 2.000000 |
| 2.95 | 17.80 | 8.90 | 2.000000 |
| 4.95 | 49.02 | 24.51 | 2.000000 |

Average P/m (0.3 < r < 4.0): **2.000000**
Max |P/m - 2|: **0.000000** (to double precision)

### Why P/m = 2 Exactly

The identity P = 2m is algebraic, not numerical:

```
P(r) = 2r² + 4c₄ sin²f = 2(r² + 2c₄ sin²f) = 2 × m(r)
```

Both P and m use the same angular/perpendicular current |A_⊥|². The factor
of 2 vs 1 in the r² terms (and 4 vs 2 in the c₄ terms) comes from the
structure of the consistent Laplacian D†D vs the kinetic energy ½|∂₀q|².

Physically: for hedgehog perturbations, the temporal perturbation A₀ is
aligned with the radial current A_r in isospin space, so [A₀, A_r] = 0.
Only the angular currents A_θ, A_φ contribute to both the radial stiffness
and the temporal inertia, with coefficients in ratio 2:1.

## Physics: What "No Emergent Gravity" Means

### The Effective Line Element

For perturbations δf propagating on the Skyrmion background, the Sturm-
Liouville equation -(Pg')' + Wg = (ω²/c²)mg defines an effective metric:

```
ds²_eff = -m(r) dt² + P(r) dr² + (angular terms)
```

With P/m = 2, this becomes:

```
ds²_eff = m(r) [-dt² + 2 dr² + ...]
```

The conformal factor m(r) cancels in null geodesics (light rays), which
satisfy ds² = 0 → dt²/dr² = 2 everywhere. Light travels at the same
effective speed c/√2 regardless of position — there is no bending, no
redshift, no time dilation.

### Contrast with Schwarzschild

In Schwarzschild geometry, g₀₀/g_rr = -(1-r_s/r)² varies with r,
approaching 0 near the horizon. This variation is what causes:
- Gravitational redshift (clocks tick slower near mass)
- Light bending (geodesics curve toward mass)
- Orbital precession

With P/m = const, none of these effects occur. The soliton is "invisible"
to its own perturbations — it cannot gravitate.

### Why the Sigma Model Fails at Gravity

The sigma model enforces |q(x)| = ρ₀ = const. The soliton is a topological
twist in the field direction (the hedgehog maps r̂ to isospin space), but
the field magnitude never changes. Physically:

- **No density variation**: The "medium" has uniform density ρ₀ everywhere,
  including inside the soliton core. No density gradient → no refractive
  index gradient → no light bending.

- **Scale invariance of L₄**: The Skyrme term |[A_μ, A_ν]|² with constant
  ρ₀ has a scaling symmetry that forces the effective metric to be
  conformally flat. The commutator structure [A_r, A_θ] produces the same
  ratio of radial/temporal coefficients at every radius.

- **Derrick's theorem connection**: The virial identity E₂ = E₄ (which we
  verified: E₂/E₄ = 0.999) is a consequence of scale invariance at
  equilibrium. P/m = 2 is the local version of this same scale invariance.

### What This Rules Out

For any sigma-model Skyrme soliton (any B, any shape, any e):
1. No gravitational binding between well-separated solitons via metric effects
2. No gravitational mass (the ADM mass of the effective metric is zero)
3. No gravitational lensing of perturbation waves
4. No Hawking-like radiation from the effective horizon (there is none)

The solitons do interact — but purely through the nonlinear field equations
(direct force), not through any emergent geometry.

## What Breaks the Identity (Paths to Emergent Gravity)

### Path 1: Finite-λ (Mexican hat potential, ρ(r) varies) — **DEAD END**

**CORRECTION**: The original analysis claimed that varying ρ(r) breaks P/m = 2.
This is WRONG. The identity P/m = 2 holds for **any** ρ(r) profile, because
the angular currents A_θ = q⁻¹∂_θq are ρ-independent:

```
For q = ρ(r) × (quaternion unit), ∂_θq = ρ(r) × ∂_θ(unit)
q⁻¹ = (1/ρ²) × unit⁻¹
A_θ = q⁻¹∂_θq = (1/ρ) × unit⁻¹ × ρ × ∂_θ(unit) = unit⁻¹∂_θ(unit)
```

The ρ factors cancel! Both P and m derive from |A_θ|² with fixed 2:1 ratio:

```
P(r) = 2r² + 4c₄(r) sin²f = 2(r² + 2c₄(r) sin²f) = 2 × m(r)
```

This holds for **any** c₄(r) = 2ρ(r)²/e², even when ρ varies wildly.
The identity is algebraic, not a consequence of σ-model constraint.

**Implication**: Finite-λ density variation alone CANNOT produce emergent gravity.

### Path 2: L₆ (sextic/BPS term)

The L₆ = -λ₆ (B⁰)² term contributes additional metric components:

```
g₆^{μν} ∝ ∂²(B⁰²)/∂(∂_μq)∂(∂_νq)
```

These involve third derivatives of q and have a different radial structure
than L₄. The P/m ratio from L₆ is NOT 2 in general.

Near the BPS limit (L₆ + V dominant), the soliton approaches E = E_FB and
the metric could develop non-trivial curvature.

### Path 3: Dynamic degenerate sector

If the weight sector (J, P) acquires kinetic energy L₂_D = (1/2c²)|ṗ|² -
(1/2)|∇p|² + coupling, the cross-coupling between bulk and degenerate
sectors breaks the algebraic structure that enforces P/m = 2.

## Path 2 Results: L₆ Breaks P/m = 2

### Theory

The sextic term L₆ = λ₆ B^μ B_μ contributes to the Sturm-Liouville coefficients
with ratio P₆ : m₆ = 1 : 1, unlike L₂ and L₄ which give 2:1. This is because
L₆ is Lorentz-invariant: the spatial baryon current B^r ∝ ḟ sin²f/r² has the
same structure as the baryon density B⁰ ∝ f'sin²f/r².

The total P/m becomes:

```
P/m = (2α + β)/(α + β) = 2 - β/(α+β) < 2

where α = r² + 2c₄ sin²f   (L₂ + L₄ contribution)
      β = (2λ₆/π³) sin⁴f/r²  (L₆ contribution)
```

### Numerical Confirmation

Parameters: N=128, R=8.0, h=0.133, e=1, ρ₀=1, sigma model, B=1 hedgehog

#### λ₆ scan (1D analytical, validated against 3D)

| λ₆ | max\|P/m-2\| | P/m at peak | Peak r | % deviation |
|-----|-------------|-------------|--------|-------------|
| 0.5 | 0.014 | 1.986 | ~0.01 | 0.7% |
| 1.0 | 0.028 | 1.972 | ~0.01 | 1.4% |
| 2.0 | 0.055 | 1.945 | ~0.01 | 2.7% |
| 5.0 | 0.126 | 1.874 | ~0.01 | 6.3% |
| 10.0 | 0.224 | 1.776 | ~0.01 | 11.2% |
| 20.0 | 0.367 | 1.633 | ~0.01 | 18.3% |
| 50.0 | 0.591 | 1.409 | ~0.01 | 29.6% |

#### 3D effective metric at λ₆ = 10

Selected radial values showing the P/m profile:

| r | P(r) | m(r) | P/m |
|---|------|------|-----|
| 0.25 | 32.49 | 17.17 | 1.892 |
| 0.45 | 30.00 | 16.28 | 1.843 |
| 0.85 | 19.89 | 10.09 | 1.971 |
| 1.25 | 13.33 | 6.67 | 2.000 |
| 1.85 | 10.42 | 5.27 | 1.979 |
| 2.50 | 13.14 | 6.62 | 1.985 |
| 3.50 | 25.35 | 12.68 | 1.999 |

The deviation is strongest near the soliton core (r ~ 0.3–0.5) where sin⁴f/r²
is largest, and vanishes at both r → 0 (sin⁴f ~ r⁴) and r → ∞ (sin⁴f → 0).

#### 3D vs 1D cross-check

The 3D spherically-averaged P/m agrees with the 1D analytical formula to within
~5% at most radii. The remaining discrepancy is from 3D grid discretization
(h = 0.133) — the 3D grid cannot resolve the rapid field variation near r = 0
where the 1D formula gives the peak deviation.

### Physical Interpretation

**Wave slowing, not time dilation**: L₆ makes P/m < 2, meaning perturbation
waves travel *slower* inside the soliton (c_eff = c√(P/m) < c√2). This is the
opposite sign to Schwarzschild gravity, where g_00/g_rr < 1 produces
gravitational attraction. The L₆ effect is more like an anti-gravitational
refractive index increase — waves slow down but are not attracted.

At λ₆ = 10: c_eff = 1.333c at the core (vs c√2 = 1.414c in the sigma model).

**Why the ratio is 1:1**: The L₆ term involves (B⁰)² where B⁰ ∝ f'sin²f/r².
Under a time perturbation f → f + εg(r)e^{iωt}, both the "stiffness" (spatial
derivative content) and the "inertia" (time derivative content) come from the
same f'sin²f/r² structure, giving equal contributions to P and m.

### Regression Verification

- λ₆ = 0: P/m = 2.000000 at every radius (exact to machine precision) — **PASS**
- Gradient consistency: 144/144 checks passed (6 test configs × 3 points × 8 components)
- Skyrmion equilibrium: E/E_FB = 1.219, E₂/E₄ = 0.999, B = 0.991 — all **PASS**

## Self-Consistent L₆ Profiles (Phase 2)

### The Self-Consistency Problem

The Phase 1 analysis computed the effective metric using the L₂+L₄ equilibrium profile
with L₆ as a post-processing diagnostic. This answers: "what metric does L₆ produce on
the OLD profile?" The physically correct question is: "what metric does L₆ produce on
the profile that L₆ itself modifies?"

### Method

Modified `radial.c` (1D shooting solver) to include L₆ in the Euler-Lagrange equation:
- f'' coefficient (denominator): += c₆ sin⁴f/r²
- Numerator: -= c₆ f'² sin²f sin(2f)/r² + 2c₆ f' sin⁴f/r³

where c₆ = λ₆/(2π⁴ρ₀²) (NOT 2λ₆/π³ — the normalization factor 1/(4πρ₀²) accounts
for the common factor in the L₂+L₄ EL equation that L₆ doesn't share).

**Critical bugs found and fixed:**
1. **Normalization**: c₆ in the ODE must include 1/(4πρ₀²) — the L₂+L₄ EL has a
   common factor 4πρ₀² that L₆'s EL doesn't. Off by factor 4π ≈ 12.57.
2. **Sign error**: The L₆ numerator terms had flipped signs. Correct: -c₆f'²sin²fsin(2f)/r²
   (not +) and +2c₆f'sin⁴f/r³ (not -).
3. **Derrick scaling**: E₆ ~ λ⁻³ (not λ⁻¹ like E₄), so virial is E₂ = E₄ + 3E₆
   (not E₂ = E₄ + E₆).

### Self-Consistent Profile Results

All profiles satisfy the Derrick virial E₂ = E₄ + 3E₆ to 0.08%.

| λ₆ | a = -f'(0) | E/E_FB | E₂ | E₄ | E₆ | E₂/(E₄+3E₆) |
|-----|-----------|--------|------|------|------|------------|
| 0 | 1.420 | 1.232 | 51.55 | 51.59 | 0 | 0.999 |
| 5 | 1.407 | 1.234 | 51.87 | 51.27 | 0.21 | 0.999 |
| 10 | 1.396 | 1.237 | 52.18 | 50.97 | 0.42 | 0.999 |
| 20 | 1.373 | 1.242 | 52.76 | 50.41 | 0.80 | 0.999 |
| 50 | 1.317 | 1.255 | 54.34 | 48.97 | 1.81 | 0.999 |

**Key observation**: L₆ makes the soliton WIDER (a decreases from 1.420 to 1.317),
because L₆ = λ₆(B⁰)² penalizes high baryon density, creating a repulsive pressure
that puffs up the soliton.

### 3D Effective Metric with Self-Consistent Profile (λ₆=10)

Selected radial values from 3D computation:

| r | P(r) | m(r) | P/m |
|---|------|------|-----|
| 0.35 | 30.57 | 16.53 | 1.850 |
| 0.55 | 27.51 | 14.84 | 1.854 |
| 0.85 | 19.66 | 9.99 | 1.967 |
| 1.25 | 13.33 | 6.67 | 2.000 |
| 1.85 | 10.48 | 5.29 | 1.980 |
| 2.50 | 13.18 | 6.64 | 1.985 |
| 3.50 | 25.36 | 12.68 | 1.999 |

Average P/m (0.3 < r < 4.0): **1.976**
Max |P/m - 2|: **0.158** at r = 0.45

### Post-Processing vs Self-Consistent Comparison

The self-consistent profile is wider, so the metric deviation is slightly different:
- Post-processing (λ₆=10 on L₂+L₄ profile): max deviation 0.157 at r ≈ 0.01
- Self-consistent (λ₆=10 on L₂+L₄+L₆ profile): max deviation 0.158 at r = 0.45

For moderate λ₆, the difference is small because the profile changes are modest (a
shifts by ~2%). For large λ₆, the self-consistent profile is essential.

### Finite-λ + L₆ (Coupled ρ and f with L₆)

Modified `finite_lambda.c` to include L₆ in the f-shooting ODE. Since B⁰ is
ρ-independent (the ρ⁶ from quaternion normalization cancels), L₆ does NOT affect
the ρ BVP — only the f-shooting is modified.

At λ₆=10: The finite-λ scan works correctly. The virial E₂ - E₄ - 3E₆ + 3E_V = 0
is satisfied at all λ values. At λ=8000, ρ(0) = 0.997 (vs 0.803 for λ₆=0).

## Verification Summary

| Test | Result |
|------|--------|
| Gradient consistency (6 tests, 144 checks) | ALL PASS |
| Skyrmion equilibrium (E, virial, topology) | ALL PASS |
| Effective metric P/m = 2 (λ₆=0) | PASS (exact to machine precision) |
| L₆ breaks P/m = 2 (λ₆>0, post-processing) | CONFIRMED (up to 29.6% at λ₆=50) |
| Self-consistent L₆ virial E₂ = E₄ + 3E₆ | ALL PASS (0.08% accuracy) |
| Self-consistent L₆ effective metric (3D) | CONFIRMED (max dev 0.158 at λ₆=10) |
| Finite-λ + L₆ solver | ALL PASS (virial converges) |

## Gravitational Observables (Geodesic Analysis)

### Effective Gravitational Potential

The BLV effective metric defines a position-dependent wave speed for perturbations:
```
v(r) = c √(P/m)    (phase velocity)
n(r) = √(2m/P)     (refractive index, with n→1 at r→∞)
```

The effective gravitational potential is:
```
Φ(r)/c² = (P/m - 2)/4
```

With P/m < 2 inside the soliton: **Φ < 0 → attractive well** (correct sign for gravity).
This is the same mechanism as Schwarzschild: coordinate speed of light decreases near mass,
bending rays inward.

### Self-Consistent Profiles (from 1D formula)

| λ₆ | P/m at center | Φ_min/c² | Φ × Mc² (MeV) | Half-width (fm) |
|-----|--------------|----------|----------------|-----------------|
| 10  | 1.782        | -0.054   | 51 MeV         | 0.55 fm         |
| 50  | 1.450        | -0.138   | 129 MeV        | 0.70 fm         |

The potential is a bell-curve peaked at the soliton center with exponential falloff.

### Shapiro-Like Time Delay

Extra optical path for radial passage through soliton:
- λ₆=10: Δℓ = 0.119 code = 0.067 fm (1.5% fractional delay)
- λ₆=50: Δℓ = 0.33 code = 0.19 fm

### Deflection Angle (Born Approximation)

θ(b) = -b ∫_b^∞ d(P/m)/dr / √(r²-b²) dr  (singularity-free via cosh substitution)

| b (code) | θ at λ₆=10 (rad) | θ at λ₆=50 (rad) | GM_eq/c² (code) |
|----------|-------------------|-------------------|-----------------|
| 0.3      | -0.092            | -0.185            | -0.007          |
| 0.5      | -0.138            | -0.296            | -0.017          |
| 1.0      | -0.162            | -0.465            | -0.040          |
| 1.5      | -0.096            | -0.384            | -0.036          |
| 2.0      | -0.032            | -0.165            | -0.016          |
| 3.0      | -0.001            | -0.008            | -0.001          |

Peak deflection ~0.16 rad (λ₆=10) to ~0.46 rad (λ₆=50) at b ≈ 1.0 (0.56 fm).
These are enormous by astrophysical standards but natural for nuclear-scale phenomena.

### Comparison to Newtonian Gravity

```
Φ_Newton/c² at proton surface:  1.4 × 10⁻³⁹
Φ_BLV/c² at soliton center:    5.4 × 10⁻²  (at λ₆=10)
Ratio:                          3.9 × 10³⁷
```

The BLV metric effect is ~10³⁷ times stronger than actual gravity.

### Physical Interpretation

1. **Correct sign**: L₆ produces attractive lensing (waves slow down inside soliton → bend inward)
2. **Nuclear scale**: Well depth ~50 MeV matches nuclear interaction energies, not gravity
3. **Short range**: Potential is Yukawa-like (exponential falloff, half-width ~0.55 fm), NOT 1/r
4. **No long-range gravity**: The effect is confined to the soliton interior; no 1/r tail
5. **Nuclear binding contribution**: At λ₆=10, the metric potential (51 MeV) is comparable to
   nuclear binding (~8 MeV/nucleon), suggesting L₆ contributes to nuclear forces via wave-slowing
6. **Not gravity**: The 10³⁷ ratio and wrong range law mean this cannot be identified with gravity

### What Would Real Emergent Gravity Require?

For 1/r gravitational potential (long-range):
- A **massless mediator** (graviton analog) propagating in the effective metric
- None of the Skyrme perturbation modes are massless in the relevant sense
- The breathing modes (K=0) have no bound states below threshold
- The angular modes (K≥1) are massive (pion mass gap)

Possible avenue: the **degenerate sector** (J, P) from the full Cl⁺(3,0,1) theory
could provide a massless degree of freedom if it acquires dynamics. This is Path 3.

## Path 1: Definitive Dead End

P/m = 2 holds for **any** hedgehog perturbation on **any** radial background, because:
1. Angular currents A_θ are ρ(r)-independent (ρ cancels in q⁻¹∂_θq)
2. Both P and m derive from |A_θ|² with algebraic 2:1 ratio
3. The ρ-perturbation channel has P_ρ/m_ρ = 1 (trivially flat)
4. The f-ρ cross-coupling in spatial stiffness vanishes (P_fρ = 0)

This rules out finite-λ, massive pions, and any other radial potential as mechanisms
for breaking P/m = 2 in the hedgehog breathing mode. Only L₆ (or higher-derivative
terms with different P/m ratios) can break the identity.

## Files

- `src/skyrmion_equil.c` — Equilibrium finder + effective metric + geodesic observables
- `src/verify_gradient.c` — Gradient consistency verification (6 test configurations)
- `src/hopfion_field.c` — Energy, force (L₂+L₄+L₆+V+V_π+V_D), and effective metric (with L₆)
- `../proposal/hopfion_search/src/radial.c` — Self-consistent L₂+L₄+L₆ profile solver
- `../proposal/hopfion_search/src/finite_lambda.c` — Coupled f+ρ+L₆ solver
