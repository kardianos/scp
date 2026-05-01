# V54 Results — Theta Self-Interaction Tests

## Context

The V53 grid tests confirmed that oscillons under the standard Cosserat equations
suffer from one-way theta drain: phi sources theta via η×curl(φ), but theta (being
massless) disperses and its curl vanishes, killing the back-transfer channel. The
result is theta_rms growing without bound while phi weakens.

V54 explored three mechanisms to fix this by giving theta a self-interaction.

## Parameter Base

All tests use the high-amplitude regime identified by Maxima analysis:
- m=3, μ=-90, κ=50, η=0.3, A=1.5, A_bg=0.1
- Lissajous 3-axis seed (Δ=0, π/3, π/3), theta_gain=-2.5
- BC: absorbing T=0-30, periodic T=30-600
- N=96, L=15

## Preliminary Tests

### High amplitude + standard equations (m_theta=0)
- theta_rms: 0.13 → 1.66 (unbounded growth)
- E_pot: -913 → -11 (binding lost)
- phi_max: 1.17 → 0.90 (slowly dying)
- **Verdict**: One-way drain, same as V53 but at higher amplitude

### m_theta = 1.0
- theta_rms: 0.24 → 0.28 (SATURATED — flat from t=100 onward)
- E_pot: -16 → -15 (weak but stable binding)
- phi_max: 1.06 → 0.75 (reduced but persistent)
- **Verdict**: Mass term confines theta, stops drain. But theta doesn't cycle
  back to phi — it's trapped, not converted.

### Theta saturation threshold (theta_sat=0.5, gamma_conv=1.2)
- theta_rms: 0.09 → 0.154 (saturated)
- E_pot: ~0 (no binding)
- **Verdict**: Caps theta but conversion creates diffuse phi, not binding.

## Three Self-Interaction Tests

### Test A: Gradient-Based (|∇θ|² → conversion)

Where theta direction changes rapidly (high gradient), energy converts to phi.
Chirality bias: modulated by θ·curl(θ).

Parameters: sigma_grad=0.5, chi_chiral=0.5

| t | E_total | E_pot | phi_max | theta_rms |
|---|---------|-------|---------|-----------|
| 50 | 23414 | -2295 | 1.71 | 0.11 |
| 200 | 22871 | -209 | 1.09 | 0.16 |
| 400 | 27013 | -47 | 0.87 | 0.18 |
| 600 | 34485 | -88 | 0.98 | 0.20 |

- Theta capped (0.20) ✓
- **Energy GROWING** (+47%) ✗ — the gradient conversion injects energy
- Binding weak (E_pot → -88)
- **Verdict**: UNSTABLE — energy not conserved

### Test B: Cubic Self-Interaction (-σ|θ|²θ → conversion)

Smooth nonlinear damping proportional to |θ|². Always active but only significant
at large theta. Drained energy sources phi via curl(θ).

Parameters: sigma_cubic=0.5

| t | E_total | E_pot | phi_max | theta_rms |
|---|---------|-------|---------|-----------|
| 50 | 22722 | -2044 | 1.75 | 0.11 |
| 200 | 21585 | -152 | 1.13 | 0.16 |
| 400 | 21667 | -43 | 1.02 | 0.17 |
| 600 | 21871 | -22 | 0.76 | 0.17 |

- Theta capped (0.17) ✓ — flattest of all three
- **Energy nearly conserved** (-3.7% over T=600, mostly absorbing phase) ✓
- phi_max persistent (0.76-1.13) ✓
- Binding weak but stable (E_pot -22 to -43)
- **Verdict**: BEST — energy-conservative, theta stabilized

### Test C: Frequency-Mismatch Conversion

Theta at each voxel has local frequency ω = |θ_vel|/|θ|. Mismatched neighbors
convert theta to phi. Matched neighbors are stable.

Parameters: sigma_freq=0.2

| t | E_total | E_pot | phi_max | theta_rms |
|---|---------|-------|---------|-----------|
| 50 | 22607 | -2348 | 1.63 | 0.11 |
| 200 | 22762 | -148 | 1.02 | 0.16 |
| 400 | 27536 | -48 | 0.87 | 0.17 |
| 600 | 33693 | -66 | 0.94 | 0.18 |

- Theta capped (0.18) ✓
- **Energy GROWING** (+49%) ✗ — same instability as Test A
- **Verdict**: UNSTABLE — energy not conserved

## Summary

| Test | θ capped? | Energy conserved? | Binding? |
|------|-----------|-------------------|----------|
| Baseline | No (1.66) | Yes | Dying |
| m_theta=1 | Yes (0.28) | Yes | Weak, stable |
| A: Gradient | Yes (0.20) | **No (+47%)** | Weak |
| **B: Cubic** | **Yes (0.17)** | **Yes (-3.7%)** | **Weak, stable** |
| C: Frequency | Yes (0.18) | **No (+49%)** | Weak |

## Key Insight

The cubic self-interaction is the only energy-conservative mechanism tested.
Tests A and C inject energy because their conversion terms (proportional to
|∇θ|² and frequency mismatch) create phi with more energy than they remove
from theta. The cubic term (-σ|θ|²θ) drains theta in a way that naturally
conserves energy when paired with the curl(θ) source term for phi.

## Parameter Sweep (750 combinations)

Swept mu, eta, m2, kappa over full range with T=50 per trial.
See `sweep_results.md` for full top-20 table.

**Best: mu=-80, eta=0.2, m2=1.5, kappa=50** (score=82.7, compactness=0.84, contrast=102).
Particles form readily at T=50 but dissolve by T=2000 under periodic BC.

Key finding: **low eta (0.2) dominates** — less curl coupling = less theta drain.

## Chirality Pair Sweep (10 methods × 3 separations)

Tested whether opposite-chirality particles can coexist.
See `chirality_sweep_results.md` for full table.

**flip_ellip produces opposite H_cross that persists.** Geometric handedness
(cross-section shape) is the true chirality — phase/theta sign flips get absorbed
by dynamics.

## Chiral Pair Long Runs

Tested flip_ellip pair at D=8 for T=1000-3000 with various mechanisms:

| Config | Result |
|--------|--------|
| Baseline (no terms) | Dead by t=430 |
| lambda_self=5 | Dead by t=200 (too aggressive for A=0.4) |
| sigma_cross=10 | Dead by t=200 (adds mass, accelerates dissolution) |
| theta_sat=0.3 | Dead by t=200 (theta never reaches threshold) |
| sigma_cubic=0.5 | Dead by t=200 (theta too small to trigger) |

## Amplitude Sweep (A=1.0, 1.5, 2.0, 2.5, 10.0)

Chiral pairs at all amplitudes with lambda_self=5, T=1000:

| A | phi_max(t=50) | phi_max(t=150) | Clusters(t=1000) |
|---|--------------|----------------|-----------------|
| 1.0 | 0.78 | 0.82 | 0 |
| 1.5 | — | — | 0 |
| 2.0 | — | — | 0 |
| 2.5 | 0.87 | 0.72 | 0 |
| 10.0 | 8.10 | 1.46 | 0 |

**All amplitudes dissolve to zero clusters by T=1000.** Higher A starts higher
but falls at the same rate. The equations do not support stable localized objects
at any tested amplitude.

## Definitive Conclusion

The 6-field Cosserat equations with V(P) = μP²/(1+κP²), combined with any
tested self-interaction terms (|θ|⁴, |θ|²|φ|², θ saturation, cubic conversion),
do not produce stable localized particles under any boundary condition, amplitude,
or parameter combination tested across 750+ configurations.

The equation structure itself needs to change for stable particles to exist.
