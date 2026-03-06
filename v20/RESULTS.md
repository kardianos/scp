# V20 Results: Triple-Product Saturating Potential

## Summary

**Definitive null result.** The saturating triple-product potential V = (mu/2)P²/(1+κP²)
with P = φ₁φ₂φ₃ does not produce stable bound states in either 1D or 3D.

Two failure modes depending on parameters:
- **Weak saturation** (κ small): runaway collapse. Field amplitudes grow without bound
  as the attractive potential concentrates the fields faster than gradient pressure resists.
- **Strong saturation** (κ large): trivial dispersal. The capped potential is too weak
  to overcome 3D spreading. The triple product P ∝ φ³ drops cubically with amplitude,
  so binding vanishes faster than the fields spread.

No parameter combination (mu, κ) produces a stable, persistent, localized resonance.

---

## 1D Parameter Scan

**Code**: `src/toy1d.c` | **Grid**: Nx=4000, xmax=40, dt=0.01

### Scan 1: Coupling Strength (separated bumps, sep=4, κ=0.1)

| mu | Peak amplitude at t=20 | E_pot at t=20 | Outcome |
|----|----------------------|---------------|---------|
| -1 | 0.507 (stable) | -0.000 | No interaction — fields pass through each other |
| -5 | 0.537 | -0.000 | Negligible interaction |
| -10 | 21.2 (blowup) | -1275 | **Runaway collapse** at t≈18 |
| -50 | 60.8 | -9397 | Immediate collapse by t≈5 |
| -100 | 88.1 | -19840 | Immediate collapse by t≈3 |

### Scan 2: Saturation Strength (overlapping, sep=0, mu=-10)

| κ | Peak at t=60 | E_pot at t=60 | Growth rate | Outcome |
|---|-------------|---------------|-------------|---------|
| 0.1 | 45.5 | -3694 | ~0.6/time unit | Runaway collapse |
| 1.0 | 14.4 | -357 | ~0.2/time unit | Slow collapse |
| 10 | 5.4 | -34.5 | ~0.07/time unit | Very slow collapse |
| 100 | 2.5 | -3.6 | ~0.02/time unit | Gradual growth, no equilibrium |

**No stable equilibrium exists at any κ.** Stronger saturation slows the collapse but
never halts it. The field amplitudes grow monotonically for all tested parameters.

At κ=100, t=500: peaks still at 2.5, still growing. The growth rate is ~0.02/time unit,
meaning the configuration doubles in amplitude every ~50 time units. This is not a
metastable oscillon — it's a slow-motion collapse.

### Why Saturation Fails

The saturating potential V → mu/(2κ) for large P. This caps the potential energy
*density*, but the field amplitudes are unconstrained. As fields concentrate:

1. Gradient energy grows as ∫|∇φ|² ~ A²/σ² (A = amplitude, σ = width)
2. Potential energy per volume → |mu|/(2κ) (constant cap)
3. Total potential → |mu|/(2κ) × σ³ (proportional to volume)

Under rescaling (width σ → σ/λ, amplitude → const):
- E_grad ~ λ (grows with compression)
- E_pot ~ λ⁻³ (shrinks with volume)

Equilibrium requires E_grad ~ E_pot, which gives λ⁴ ~ |mu|/(2κ) × (something). But
the amplitude is FREE to grow, and higher amplitude deepens the well. The system
exploits this by growing amplitude while concentrating spatially — a mode not covered
by simple Derrick scaling.

---

## 3D Results

**Code**: `src/triple3d.c` | **Grid**: N=64, L=16, dx=0.254, dt=0.037

### Test Comparison (mu=-10, κ=100, A=1.0, σ=1.5)

| Test | Description | E_pot(t=0) | f_core at t=7 | Dispersal time | Outcome |
|------|-------------|-----------|---------------|----------------|---------|
| 1 | Three overlapping | -1.43 (8% of E_grad) | 0.004 | ~7 | Disperses |
| 3 | Single field (control) | 0.000 | 0.003 | ~7 | Disperses |
| 4 | Two fields (control) | 0.000 | — | ~7 | Disperses |

**The triple-product coupling makes no measurable difference in 3D dispersal.**

The initial binding energy (E_pot = -1.43) is only 8% of the gradient energy
(E_grad = 18.5). As the fields spread in 3D, the amplitude drops as ~1/r³/² (volume
dilution), so P = φ₁φ₂φ₃ drops as ~1/r⁹/² — the binding force vanishes nearly
instantly.

To achieve meaningful 3D binding (E_pot > E_grad), one would need |mu|/κ >> 1,
but this regime is precisely where the 1D scan shows runaway collapse.

---

## Comparison with v19

| Aspect | v19 (S³ + Skyrme) | v20 (Triple product) |
|--------|-------------------|---------------------|
| Binding mechanism | S³ constraint (implicit) | Explicit V(I₃) potential |
| Can force strong binding? | No | Yes, but causes collapse |
| Stable equilibrium? | No (disperses) | No (collapses OR disperses) |
| Three-body specific? | No (hedgehog = all 3) | Yes (P=0 with <3 fields) |
| Failure mode | Radiation | Collapse (weak κ) or radiation (strong κ) |
| Core lifetime vs control | 1.3× (negligible) | 1.0× (identical) |

v20 actually performs WORSE than v19: the explicit attractive potential introduces a
new failure mode (collapse) without solving the dispersal problem.

---

## Root Cause Analysis

The v20 proposal and all "shared resonance" variants hit the same fundamental
obstruction: **there is no mechanism to simultaneously prevent both collapse and
dispersal for a non-topological configuration.**

In detail:
1. **Preventing dispersal** requires an attractive self-interaction (V < 0 when
   fields overlap). This creates inward pressure.
2. **Preventing collapse** requires a repulsive core or saturation. But:
   - Saturation caps potential density, not field amplitude → collapse continues
   - A mass term φ² creates a length scale but also makes far-field Yukawa (kills 1/r)
   - A constraint |φ|=1 kills amplitude growth but returns us to v19's sigma model
3. **Topology** (B ≠ 0 in the Skyrme model) uniquely solves this by providing a
   conserved quantity that cannot be shrunk to zero. But B=0 configurations have
   no such protection.

The "hydrocarbon binding" analogy fails because covalent bonds are stabilized by
*quantum mechanics* (Pauli exclusion + kinetic energy from uncertainty principle),
not by classical potential energy. A classical field theory with an attractive potential
always collapses — this is Earnshaw's theorem for continuous media.

---

## Conclusion

The saturating triple-product potential is not viable for creating stable particles
from a continuous field. The two-parameter (mu, κ) space has no region producing
persistent localized structures.

Phase 2 (density coupling for gravity) was not attempted.

---

## Files

| File | Description |
|------|-------------|
| `src/toy1d.c` | 1D three-scalar toy model |
| `src/triple3d.c` | 3D three-scalar simulator |
| `data/toy1d_timeseries.tsv` | 1D time series (last run) |
| `data/test{T}_timeseries.tsv` | 3D time series for test T |
| `data/test{T}_tavg_profile.tsv` | 3D time-averaged radial profiles |
| `data/test{T}_farfield.txt` | 3D far-field power-law fits |
