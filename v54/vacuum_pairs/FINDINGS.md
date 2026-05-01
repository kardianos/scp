# V54 Vacuum Pair Extension — Findings

## What Another AI Built (CUDA kernel, scp_sim.cu lines 816-832)

Confining Lagrangian potentials implemented with **minus signs**:

```
V_cross = +(σ/2)|θ|²|φ|²     →  acc_phi   -= σ|θ|²φ     (theta makes phi heavier)
V_self  = +(λ/4)|θ|⁴          →  acc_theta -= λ|θ|²θ     (theta self-limits)
                                  acc_theta -= σ|φ|²θ     (phi makes theta heavier)
```

Parameters tested: sigma_cross=20, lambda_self=10 (testG).
Result: theta_rms stable at 0.322 ± 0.0004 for 400 TU. Energy ±0.17%.
But E_pot ≈ 0 — no binding energy detected.

## Sign Convention Error in Our V54 Postulates

Our vacuum pair postulates used **plus signs** (pair creation), the CUDA kernel
uses **minus signs** (confinement). These are opposite physics:
- Plus: phi reduces theta mass → instability above threshold (Schwinger)
- Minus: phi increases theta mass → always stable (mutual trapping)

The CPU kernel (scp_sim_v54.c) now matches the CUDA kernel's minus signs.

## Higgs VEV Attempt — FAILED

Added theta_vev parameter for Mexican hat potential V=(λ/4)(|θ|²-v²)².
Theory predicted "vacuum binding through theta holes."

**Particle tracking showed this is WRONG:**

| Config | Particles at t=200 | Largest mass | Best E_pot |
|--------|-------------------|-------------|-----------|
| sigma_cross=2, lambda_self=1, **theta_vev=0.5** | **0** | dead by t=70 | — |
| sigma_cross=2, lambda_self=1, **theta_vev=0** | **11** | 103 | -129 |

The Higgs VEV **kills oscillons faster** than the control. The theta condensate
grows everywhere and smothers the localized structure instead of creating a
binding hole around it.

Lesson: global energy diagnostics (E_pot < 0) prove nothing. Only particle
tracking with cluster detection shows whether localized bound structures exist.

## What Actually Works: Moderate Confining Cross-Potential

The control run (sigma_cross=2, lambda_self=1, theta_vev=0) maintains 4-12
particles through T=200 with mass 60-200 and E_pot -68 to -270. This is the
confining mechanism without any VEV.

Particle tracking data (control, late time t≥150):
- Average particle count: 8.2
- Average largest particle mass: 108.8
- Average E_pot of largest: -127.8
- phi_max: 0.74-0.96 (sustained)

## Three-Way Comparison (Particle Tracking Verified)

Ran identical simulations (N=64, L=15, T=200, periodic BC, braid init A=1.5)
with sfa_particle_track for ground truth:

| Config | Particles t=200 | Largest mass | Best E_pot |
|--------|-----------------|-------------|-----------|
| **Baseline** (no new terms) | **10** | **210** | **-253** |
| Confining (σ=2, λ=1, v=0) | 11 | 103 | -129 |
| Higgs VEV (σ=2, λ=1, v=0.5) | 0 | DEAD | — |

Late-time averages (t≥150):
| Config | Avg particles | Avg largest mass | Avg E_pot |
|--------|-------------|-----------------|----------|
| **Baseline** | 7.2 | **222** | **-255** |
| Confining | 8.2 | 109 | -128 |
| Higgs VEV | 0 | DEAD | — |

**The baseline wins.** The confining cross-potential fragments structures into
more but smaller pieces. The Higgs VEV destroys them entirely.

## Honest Conclusions

1. The confining cross-potential (other AI's testG) creates a stable theta
   equilibrium but does NOT improve particle binding vs baseline.

2. The Higgs VEV mechanism kills oscillons — the theta condensate smothers
   localized structures rather than creating binding holes.

3. At these parameters (m=3, μ=-90, κ=50, η=0.3), the unmodified Cosserat
   equations produce the best particle survival.

4. The v54 RESULTS.md conclusion stands: "The equation structure itself
   needs to change for stable particles to exist." The cross-potential
   and self-potential do not provide the needed change.

## What Was Learned

- Global diagnostics (E_pot, theta_rms, energy conservation) can be
  misleading. Only particle tracking with cluster detection shows the truth.
- The confining mechanism stabilizes THETA but not PARTICLES. These are
  different problems.
- The Higgs VEV creates a uniform condensate that competes with localization.
  The binding mechanism (theta hole) doesn't work because the condensate
  grows everywhere uniformly rather than forming holes at specific locations.
- The path to stable bound particles likely requires a mechanism that
  specifically concentrates energy rather than spreading it uniformly.
