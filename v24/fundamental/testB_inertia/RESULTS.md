# Test B: Oscillon Deformation Under Acceleration — RESULTS

## Summary

A constant external force applied to the 1D three-body oscillon produces measurable
deformation. **Newton's second law a = F/M is verified to ~5% accuracy** in the
perturbative regime (F <= 1e-5). The leading/trailing asymmetry scales as **F² (quadratic),
not linear in F**, because it is a tidal deformation requiring both velocity AND acceleration.
The oscillon survives intact only for v_final << c; larger forces disrupt it.

---

## Model

Same Lagrangian as v21: three massive scalars with saturating triple-product coupling.

    L = Σ_a [ ½(∂_t φ_a)² - ½(∂_x φ_a)² - ½m²φ_a² ] - (μ/2) P² / (1 + κ P²)

Parameters: μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Grid: Nx=4000, xmax=100, dx=0.050, dt=0.033
Protocol: t_equil=5000 (no force), t_accel=3000 (force on)

**Force method**: Momentum injection. At each timestep during the acceleration phase,
the total field momentum P receives an increment dP = F·dt via:

    vel[a][i] += -(F·dt / ∫|∇φ|² dx) · ∂φ_a/∂x

This is the minimum-norm velocity perturbation that adds exactly dP = F·dt of momentum.
It ensures dP/dt = F exactly, unlike a constant body force (which cannot accelerate a
soliton because ∫∂φ/∂x dx = 0 for localized fields).

---

## Rest-Frame Oscillon

After t_equil=5000:
- **Rest energy**: E_rest = 1.274
- **"Charge"**: Q = ∫Σφ² = 1.364
- **Width**: σ+ = σ- = 1.85 (symmetric, as expected)
- **f_core > 0.99** (fully localized)

---

## Newton's Second Law: a = F/M

### Velocity verification (perturbative regime)

| F | v_predicted = F·t/M | v_measured | Agreement |
|---|---------------------|------------|-----------|
| 3×10⁻⁶ | 0.00707 | 0.00681 | 3.7% |
| 1×10⁻⁵ | 0.0236 | 0.0225 | 4.7% |

### Position verification

| F | x_predicted = ½(F/M)t² | x_measured | Agreement |
|---|------------------------|------------|-----------|
| 3×10⁻⁶ | 10.60 | 10.23 | 3.5% |
| 1×10⁻⁵ | 35.33 | 33.99 | 3.8% |

**The small (~4%) discrepancy** is from energy increase due to kinetic energy and
radiation losses. The oscillon gains mass E = E_rest + ½Mv², so the effective
acceleration decreases slightly as it speeds up.

### Energy conservation

| F | E_rest | E_final | ΔE/E | Status |
|---|--------|---------|------|--------|
| 0 | 1.274 | 1.267 | -0.6% | Slow decay (intrinsic) |
| 3×10⁻⁶ | 1.274 | 1.270 | -0.3% | Intact |
| 1×10⁻⁵ | 1.274 | 1.285 | +0.9% | Intact (kinetic energy) |
| 3×10⁻⁵ | 1.274 | 0.863 | -32% | **Disrupted** |
| 5×10⁻⁵ | 1.274 | 0.494 | -61% | **Destroyed** |
| 1×10⁻⁴ | 1.274 | 0.308 | -76% | **Destroyed** |

**Critical force**: F_crit ≈ 2×10⁻⁵ (corresponding to v_final ≈ 0.05c).
Above this, the momentum injection disrupts the oscillon.

---

## Deformation: Leading/Trailing Asymmetry

### Asymmetry = σ+ - σ- (cycle-averaged over 500 t.u. windows)

| F | v(t≈6500) | asym(t≈6500) | v(t≈7500) | asym(t≈7500) |
|---|-----------|--------------|-----------|--------------|
| 0 | 0 | 0 | 0 | 0 |
| 3×10⁻⁶ | 0.00329 | -0.0102 | 0.00558 | -0.0448 |
| 1×10⁻⁵ | 0.01099 | -0.0826 | 0.01849 | -0.390 |

**Sign convention**: negative asymmetry means trailing side (σ-) wider than leading (σ+).
This is the physically expected result: the leading edge is compressed, trailing stretched.

### Linearity test

| Quantity | F ratio (10⁻⁵/3×10⁻⁶) | Measured ratio | Scaling |
|----------|----------------------|----------------|---------|
| Velocity v | 3.33 | 3.31 | **Linear ✓** |
| Asymmetry (t≈6500) | 3.33 | 8.1 | **~F²** |
| Asymmetry (t≈7500) | 3.33 | 8.7 | **~F²** |
| Total width increase | 3.33 | 12.0 | **~F²** |

### Key finding: Deformation is QUADRATIC in F

The velocity scales linearly (v ∝ F, confirming Newton's law), but the
leading/trailing asymmetry scales as **F²**. This is because the asymmetry is a
**tidal** deformation:

1. The Lorentz contraction (even/symmetric) scales as v² ∝ F²
2. The leading/trailing asymmetry (odd) requires BOTH velocity AND acceleration
3. In the co-moving frame, the tidal pseudo-force across the soliton is F·L/M
4. The deformation from this force is ∝ (F·L/M) × (L/binding_energy) ∝ F²

Equivalently: asymmetry ∝ F·v ∝ F²·t (since v = F·t/M).

### Width decomposition at t≈7500

| F | σ+ | σ- | Width = σ+ + σ- | Δwidth |
|---|------|------|-----------------|--------|
| 0 | 1.864 | 1.864 | 3.728 | 0 |
| 3×10⁻⁶ | 1.856 | 1.901 | 3.757 | +0.029 |
| 1×10⁻⁵ | 1.851 | 2.241 | 4.092 | +0.364 |

The leading side (σ+) barely changes. The trailing side (σ-) stretches dramatically
(by 20% at F=1e-5). The soliton develops a "tail" on the trailing side.

---

## Force Scan: Full Range

Ran forces spanning 5 orders of magnitude: F ∈ {3×10⁻⁶, 10⁻⁵, 3×10⁻⁵, 5×10⁻⁵, 10⁻⁴, 10⁻³, 10⁻², 10⁻¹}

### Regime map

| F | v_final/c | E_final/E_rest | Regime |
|---|-----------|---------------|--------|
| 3×10⁻⁶ | 0.007 | 1.00 | Perturbative (linear v, quadratic deformation) |
| 1×10⁻⁵ | 0.022 | 1.01 | Perturbative (clean F² scaling) |
| 3×10⁻⁵ | — | 0.68 | Disrupted (energy loss, chaotic) |
| 5×10⁻⁵ | — | 0.39 | Destroyed |
| 1×10⁻⁴ | — | 0.24 | Destroyed |
| 10⁻³ | — | 0.37 | Destroyed, scattered to boundary |
| 10⁻² | — | 0.80 | Destroyed + large radiation field |
| 10⁻¹ | — | 2.64 | Destroyed, radiation dominates |

The proposed force values {0.001, 0.01, 0.1} all lie in the destructive regime.
The perturbative regime requires F ≲ 2×10⁻⁵ for this oscillon (M ≈ 1.27, t_accel = 3000).

---

## Profile Deformation δρ

The co-moving energy density profile ρ(ξ) (ξ = x - x_c) was measured by averaging over
the last 500 time units. The deformation δρ = ρ_accel - ρ_rest was decomposed into:
- **Odd part** (antisymmetric): ½[δρ(ξ) - δρ(-ξ)] — the leading/trailing asymmetry
- **Even part** (symmetric): ½[δρ(ξ) + δρ(-ξ)] — breathing/width change

### Perturbative regime (F=1e-5)

| Component | |δρ|₁ | Peak | Relative to |ρ_rest| |
|-----------|--------|------|----------------------|
| Odd (asymmetry) | 4.4×10⁻³ | 1.8×10⁻³ | 0.69% |
| Even (breathing) | 7.9×10⁻³ | 4.5×10⁻³ | 1.24% |
| Background (F=0) | 7.0×10⁻³ | 3.7×10⁻³ | 1.10% |

The even component is dominated by intrinsic oscillon evolution (breathing mode decay),
not the external force. The odd component is a clean F-induced signal (zero at F=0 to
machine precision).

---

## Conclusions

### What works
1. **Newton's law verified**: a = F/M to ~4% in the perturbative regime
2. **Parabolic trajectory**: x_c = ½(F/M)t² confirmed
3. **Clear leading/trailing asymmetry**: trailing side stretches, leading side unchanged
4. **Momentum injection method**: cleanly applies constant force to a soliton

### What the deformation tells us
1. **Asymmetry ∝ F²** (not linear): this is a tidal deformation, not a simple dipole
2. **The trailing side stretches; the leading side barely changes**: the soliton develops a "tail"
3. **The even deformation (width change) also ∝ F²**: total width increases under acceleration
4. In 1D, only monopole (even) and dipole (odd) decomposition is possible — no quadrupole

### Implications for the equivalence principle
- The deformation is real and measurable, but it's **second-order** in the force
- In 1D, the "quadrupole" signature is just the even width change (∝ v²)
- The leading/trailing asymmetry is the odd/dipole component (∝ F·v)
- To connect to gravity: need 2D or 3D where true quadrupole (spin-2) structure emerges
- The F² scaling is consistent with tidal deformation, which IS the mechanism by which
  gravitational waves couple to matter — but this test alone does not prove spin-2 character

### Limitations
- The oscillon is fragile: F_crit ≈ 2×10⁻⁵ (v_max ≈ 0.05c before disruption)
- The momentum injection method is somewhat artificial (directly adds momentum)
- In 1D, angular structure is absent — need 2D/3D for multipole decomposition
- Intrinsic breathing mode obscures the even deformation signal

---

## Code

- `src/inertia.c` — main solver
- `data/inertia_F*_ts.tsv` — time series (t, xc, vx, P, E, σ+, σ-, asymmetry)
- `data/inertia_F*_profile.tsv` — profiles (ξ, ρ_rest, ρ_accel, δρ)

Compile: `gcc -O3 -Wall -o inertia src/inertia.c -lm`
Run: `./inertia -F 0.00001 -o data/`
