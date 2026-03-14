# V22: Two-Oscillon Gravitational Interaction

## Goal

Test whether two v21 oscillons attract via the scalar gravity mediator Phi.
Measure separation vs time and compare with/without gravity.

## Setup

Two oscillons placed on z-axis at z = +D/2 and z = -D/2 (D=20).
Same v21 parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0.
Gravity: alpha=0.1, beta=0.1 (full-rho source, Option A from v21).

## Grid Parameters

| Run | N | L | dx | pts/sigma | Memory | Est. time |
|-----|---|---|----|---------  |--------|-----------|
| Production (gravity) | 400 | 60 | 0.30 | 10.0 | 7.2 GB | ~6 hours |
| Control (no gravity) | 350 | 60 | 0.34 | 8.7 | 4.8 GB | ~4 hours |

- CFL = 0.25, dt = 0.075 (N=400) / 0.086 (N=350)
- Absorbing boundary: 70-95% of L (R_abs = 42-57)
- Profile output every 50 t.u. along z-axis
- tfinal = 500

## Expected Signal

Gravitational force: F = alpha * E / (4*pi*r^2) where r = D/2 = 10.
With E ≈ 85 per oscillon (post-shedding), alpha=0.1:
- F ≈ 0.1 * 85 / (4*pi*100) ≈ 0.0068
- Acceleration: a = F/E ≈ 8e-5
- Displacement after 350 t.u. (post-shedding): dz ≈ 0.5*a*350^2 ≈ 5

This 5-unit displacement should be clearly visible relative to the control.

## Diagnostics

- separation(t): z_upper - z_lower (energy-weighted centroids)
- E_upper, E_lower: energy in z>0 and z<0 halves
- fc_upper, fc_lower: core fraction (relative to sigma-centered core)
- Phi_mid: gravity potential at midpoint (z=0)
- Q_total: monopole moment (should be ~ constant with full-rho source)
- Z-axis profiles at t = 0, 50, 100, ..., 500

## Status

- [x] Code: `src/two_oscillon.c` (compiled, validated at N=100)
- [x] Production run complete (N=400, alpha=0.1, beta=0.1, t=500, 371 min)
- [x] Control run complete (N=350, alpha=0, beta=0, t=500, 234 min)
- [x] Results analyzed — see RESULTS.md
- [x] Invalid control (N=200) killed — dx=0.60 too coarse (5 pts/sigma)

## Notes

- Initial N=200 control (dx=0.60) showed oscillons dying from insufficient
  resolution (only 5 grid points per sigma, need >= 10). Killed and restarted
  at N=350 (dx=0.34, 8.7 pts/sigma).
- Shedding phase (t < ~50): centroids move outward by ~0.5-1.0 each due to
  radiation pressure. Shorter than predicted ~150 t.u.
- Post-shedding: direct field-tail coupling DOMINATES over gravitational 1/r^2.
  Control shows merger at t~409; gravity run shows monotonic expansion.
- Key finding: gravity enhances m_eff, compactifying oscillons and suppressing
  the tail overlap that drives direct attraction. Net effect is repulsive.
- See RESULTS.md for full analysis.
