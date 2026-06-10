# v66 Smoke Test — First 3D Q-ball Persistence Runs

**Date**: 2026-06-09. **Status**: [estimate] — single resolution (N=96, L=20, f32 snapshots),
T=100, CPU kernel. Production GPU runs to follow.

## Setup

Seed: `init=qball` from `v66/results/profile_w1.39.txt` (radial solver, ω=1.39,
continuum Q=482.2, E=691.9, r_half=4.41). Grid N=96, L=20 (dx=0.421), dt_factor=0.025,
bc_type=0 (absorbing sphere, damp_width=3, damp_rate=0.01), T=100, diag_dt=0.5.
Standard parameters m²=2.25, m_θ²=0, μ=−41.345, κ=50. Configs preserved in
`/tmp/v66_smoke/*.cfg` (copies of the SPEC §9d runs).

## RUN 1 — η = 0 (exact ansatz): **PASS — first persistent particle in program history**

| metric | t=0 | t=99.5 | verdict |
|---|---|---|---|
| Q_total | 482.1846 | 482.1811 | **retention 99.99927%** (gate: >99%) |
| E_total | 691.325 | 691.299 | −0.0037% |
| s_max | 0.0687 | 0.0684 (range 0.0654–0.0705) | stationary, ±4% (gate: core within 2×) |
| r_core | 3.705 | 3.703 | frozen |
| θ_rms | 0 | 0 (exact) | θ stays identically zero, as THEORY §3 predicts |

- Late-time drain dQ/dt = −5.1e−6 → extrapolated **Q half-life ≈ 5×10⁷ t.u.**
  Baseline: every v65 real oscillon lost its core within a few t.u. (FINDINGS_STABILITY).
- s_max breathing is the lattice-discretization mode at ~0.94 rad/t (verifier FFT);
  **no power at the 3ω harmonic (4.17 rad/t)** — the v65 radiation channel is absent.
  This is the smoking-gun signature that the U(1) mechanism works as derived.
- CONCEPT.md §6 item 1 ("no localized structure persists under absorbing BC") is
  **falsified by this run** for the complexified theory. Item 2 (conserved charge) is
  supplied by the Noether Q (machine-precision conservation, drift 7e−14 over T=20
  periodic; here 7e−6/100 t.u. through an absorbing sponge).

## RUN 2 — η = 0.5 (full Cosserat coupling): coherent drain along the stable branch

| metric | t=0 | t=99.5 |
|---|---|---|
| Q_total | 482.185 | 345.89 (71.7% retained) |
| Q_phi / Q_theta | 482.2 / 0 | 319.7 / 26.2 |
| E_total | 691.3 | 502.7 (−27.3%) |
| s_max | 0.0687 | 0.0648 (transient spike 0.119 at t=8) |
| r_core | 3.705 | 3.58 |

Windowed drain rates:

| window | dQ/dt | dE/dt | E/Q at end |
|---|---|---|---|
| 0–10 | −0.000 | −0.05 | 1.4327 |
| 10–25 | −0.57 | −0.69 | 1.4368 |
| 25–50 | −1.84 | −2.54 | 1.4428 |
| 50–75 | −1.76 | −2.46 | 1.4481 |
| 75–99.5 | −1.56 | −2.18 | 1.4534 |

Interpretation [estimate]:
- ~10 t.u. quiet period while the θ dressing builds (light-crossing of the core),
  with a single s_max compression spike at t≈8; then steady θ-mediated radiation,
  exactly the channel THEORY §5 predicts (massless Θ driven at ω is radiative;
  the U(1) does not protect it).
- **The ball stays ON the equilibrium branch while draining**: measured E/Q at
  Q=346 is 1.4534; the radial-solver scan gives E/Q = 1.450 at the ω whose
  Q ≈ 352 (ω≈1.40). The drain is adiabatic evolution up the stable branch
  (ω_eff rising as Q falls), not disintegration.
- Projection: at dQ/dt ≈ −1.6 the ball reaches the evaporation threshold
  (ω≈1.438, Q≈145) after ~120–150 more t.u.; naive half-life ≈ 105 t.u.
  Whether the drain shuts off (drain ∝ source amplitude, which falls along the
  branch) or proceeds to evaporation requires a longer run — open question for
  the GPU production runs.
- Global charge bookkeeping: Q_theta = 26.2 in-box at T=100; the remaining
  deficit (~110) exited through the absorbing sponge as θ radiation — measured
  outflow, not violation (THEORY §2 boundary caveat).

## Verdicts

1. **Build PASS**: η=0 persistence gates met with 3+ orders of margin.
2. **Stability crisis (v65) resolved as designed**: charge cures Derrick; rotating
   phase kills the 3ω channel; both confirmed in 3D, not just radially.
3. **New science target identified**: the η-drain. Either the ball settles to a
   finite-η attractor (a TRUE stable particle of the full 6-field theory) or it
   evaporates at the threshold. This is the next measurement (longer T, larger box,
   GPU; vary η to test the drain ∝ η² perturbative prediction).

## Files

- diag: `/tmp/v66_smoke/qball_eta0_diag.tsv`, `qball_eta05_diag.tsv` (copy to
  v66/results/ before /tmp cleanup), SFA snapshots `/tmp/v66_smoke/*.sfa`.
- Kernel: `complex_phi=1` mode in `sfa/sim/scp_sim.{c}` + `scp_config.h` + `scp_init.h`
  (verified: real-limit 2e−15, Q-drift 7e−14, 3-lens adversarial review, 1 fix round).
- Solver: `v66/radial_qball.c`, scan `v66/results/scan.tsv` (window ω∈(1.309,1.438) abs-stable).
