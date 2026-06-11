# v70 FINDINGS — Skeptical verification of v66–v69 claims (visual + lever-arm)

**Date**: 2026-06-11. **Status**: [measured] — 9-run V100 campaign (all rc=0,
gauss_max ≤ 1.2e-13 flat, energy drift ≤ 0.05%, 10.9 min/run) + full-field visual
analysis of banked v66–v69 SFAs. Plan: `v70/PLAN.md`. New tools:
`sfa/analysis/sfa_slice.c` + `render_slices.py` (complex/gauge-aware slice/lineout
visual pipeline, in Makefile), `v70/analysis/phase_tilt.py`,
`v70/analysis/forcelaw_fit.py`. Data: `v70/results/`, SFAs `/space/scp/v70/`.

## 0. Executive summary — what survived, what changed

| v69 claim | v70 verdict |
|---|---|
| Force fields artifact-free | **CONFIRMED** (visual, full-field) |
| Fission saddle (fiss1 spherical) | **CONFIRMED** (visual) |
| Positronium "no annihilation / no open decay channel" | **CONTRADICTED** — slow annihilation measured (charge segregation −4× over T=600) |
| ℏ_eff = Q (DEBROGLIE) | **CONFIRMED at 3–5%** (measured, was an open analysis item) |
| Mediator massless / 1/r² | **CONFIRMED, strengthened**: halo flux flat ±2% over r∈[6,32]; m_A = (−0.4±1.0)e-3 |
| Force = gQE, both signs | **CONFIRMED ≤6%** at D=20 (boundary-clean regime) |
| "Absolute Coulomb calibration 10–20%" (fl runs, D=14–26) | **PARTLY ARTIFACTUAL** — naive-vacuum comparison ignores the jellium correction (−9% at D=20, −34% at r=32) and the settling transient; the 0.92→1.17 "drift" was the transient |
| Force-law exponent measurable from pair dynamics | **NOT IN THIS BOX GEOMETRY** — wall-image systematics dominate beyond D ≈ L/2, up to full sign inversion at D=40 |

## 1. Visual confirmation of banked v69 runs [measured, full-field]

Slices (rho2, s, rhoQ, th2, |E|) of `/space/scp/v69/*.sfa` (renders in /tmp/v70vis/):

- **fl14/fl20/fl26**: two intact, spherical, same-charge balls; textbook radial
  Coulomb |E| halos; θ sector exactly empty; no radiation wake or boundary artifact
  at the colormap floor. The v69 tracker inputs were clean fields. fl14's halos
  visibly bridge — consistent with its watershed classification.
- **fiss1** (Q=1745 ≈ 1.9 Q_max): single spherical ball for T=300, steady Coulomb
  ring; settling-shell radiation at t≈75–125 exits through the sponge. Supports the
  symmetric-saddle interpretation.
- **pos1**: confirmed a single heavily-overlapped pulsing composite, NOT an orbit
  (v69's own caveat). **NEW**: the ± charge lobes fade — z=0-plane segregation
  envelope decays ~4× over T=600 (plane Q⁺ 108 → 15–25, breathing-modulated;
  x-dipole moment 1016 → ~200). The gauged pair annihilates SLOWLY (~10× slower
  than ungauged tb3); "no open decay channel" is wrong as physics — annihilation
  to A-waves + free quanta is open, merely slow. [plane proxy, f16; volume
  per-lobe Q(t) is the tightening follow-up]

## 2. ℏ_eff = Q fingerprint — measured [bs03/05/07, was v68 open item #4]

In-core phase tilt k from lineouts through the moving core; v measured from the
same SFA:

| run | v_meas | k (t=0..60) | γωv(v_meas) | k/γωv | ℏ_eff = p/k | /Q |
|---|---|---|---|---|---|---|
| bs03 | 0.2977 | 0.4373 ± 0.0001 | 0.4335 | 1.009 | 493.5 | 1.023 |
| bs05 | 0.4864 | 0.8023 ± 0.0011 | 0.7738 | 1.037 | 480.1 | 0.996 |
| bs07 | 0.6678 | 1.3558 ± 0.0063 | 1.2472 | 1.087 | 457.9 | 0.950 |

The tilt is dynamically frozen at the seeded γ(v_nom)·ω·v_nom to 4 decimals over
60 t.u. — real dynamics, not a seed echo. k/γωv excess grows at the
lattice-dispersion scale ((k·dx)²/6 ≈ 3% at bs07); group velocity lags 1–5%
(lattice drag). **ℏ_eff = E₀/ω = Q(1+ε) holds to 3–5%** — the generic U(1)
soliton identity, now measured. (k = p, the "true quantum composite" alternative,
is excluded trivially: it would be ~500× larger, sub-voxel.)

## 3. Single-ball halo: the clean force-law measurement [ctrl1]

ctrl1 (N=192, L=36, ball at origin, T=150): E_r·4πr²/g vs the jellium-corrected
enclosed charge Q(1 − 4πr³/3V) (the wrapped-stencil neutralizing background,
v69 SPEC O2):

- Ratio = 1.000 ± 0.02 (point-wise ripple, NO trend) over r ∈ [6,32] — five ball
  radii. The deviations match the documented jellium correction exactly
  (e.g. r=32: measured 204, jellium 199, naive 315).
- **Yukawa mass bound: m_A = (−0.4 ± 1.0)×10⁻³** (slope of log-flux-residual,
  r ∈ [8,30]) — the mediator is massless to ~1500× below m=1.5. This replaces
  pair-dynamics exponent fits as the program's direct 1/r² statement.
- Drift control: single-ball centroid a_sys = −3.5×10⁻⁸ — three orders below the
  weakest pair signal. Lattice drift is a non-issue.

## 4. Force-law lever-arm campaign — force confirmed at D=20; wall images beyond

9 runs, N=192, L=36 (dx=0.3770 = v69 fl resolution), ω=1.42 g=0.05 pairs
(per-ball Q=315.4, M=461.2), η=0, m_θ=1.6, T=150, 31 frames. Fit: quadratic
with v0 term (absorbs the O(g²) settling kick), t ≥ 20, residual bootstrap.

| D | a_same | a_opp | estimator (s−o)/2 | naive Coulomb | est/naive |
|---|---|---|---|---|---|
| 20 | +1.788e-4 ± 0.02 | −1.829e-4 ± 0.03 | +1.808e-4 | 2.14e-4 | 0.84 |
| 26 | +1.170e-4 ± 0.02 | −0.630e-4 ± 0.02 | +0.900e-4 | 1.26e-4 | 0.71 |
| 32 | +0.899e-4 ± 0.02 | **+0.025e-4 ± 0.02** | +0.437e-4 | 0.83e-4 | 0.53 |
| 40 | +1.035e-4 ± 0.02 | **+0.867e-4 ± 0.02** | +0.084e-4 | 0.53e-4 | 0.16 |

- **D=20 (boundary-clean): force = gQE confirmed.** Same/opposite magnitudes
  symmetric to 2.3%; both match the prediction computed from the MEASURED
  jellium-corrected field (a = 2gQ·E_ctrl(D)/M = 1.91e-4) to ≤6%.
- **D ≥ 26: configuration-dependent boundary electrostatics take over.** The
  absorbing sponge clamps E in the shell (≈ conductor walls → opposite-sign
  images). Opposite-charge attraction: 0.96 → 0.61 → 0.03 → **sign-inverted
  +1.64** of |naive| across D = 20→40. Same-charge: wall images pull outward,
  inflating "repulsion" to 1.97× naive at D=40. The sign inversion at D=40 was
  predicted in-flight from the image model before flo_d40 completed — confirmed.
- A naive exponent fit over the estimators returns n = 2.9 ± 0.05 — this is the
  boundary systematic, NOT the force law. **No pair-dynamics exponent is
  measurable in a standard box beyond D ≈ L/2.** Requirement for any future
  absolute force work: wall distance ≳ D (i.e. L ≳ 2D + sponge), N=384 at D=40.
- Retroactive: the v69 fl14/20/26 "absolute calibration" compared against naive
  vacuum Coulomb; the in-box prediction at D=20 carries a −9% jellium correction
  (and the v0-less fits carried transient bias: refit fl14 a = +0.7e-5 ± 0.6e-5,
  a TRUE watershed null; fl20 0.86, fl26 0.98 of naive). The v69 conclusion
  "long-range force is Coulomb" survives, but its quoted precision was partly
  coincidental.

## 5. Verification-debt items spawned

1. Volume-integrated per-lobe Q(t) for pos1-class runs (annihilation rate law).
2. Big-box (L ≳ 2D) pair runs if an exponent-from-dynamics number is ever needed
   — or accept the halo measurement (§3) as the definitive 1/r² statement.
3. The sponge≈conductor image model is qualitative; if boundary-corrected force
   work matters, derive the actual sponge Green function or implement
   transverse-only damping (GAUGE_DESIGN risk #2 / SPEC O3 v2 option).
4. `sfa_slice`/`render_slices.py` should be the standard first look at every new
   object class — pos1 demonstrated diag-only conclusions can mislead.

## Infrastructure

- Campaign: single V100-32GB Vast.ai instance, 9 config-mode runs with eager
  auto_download, zero failures, teardown after download verification (~12 GB).
- Seeds: `gen_qball_pair`/`gen_qball_boost` built ON the instance (CPU) from the
  4-column gauged profile (generators read first 2 columns; E built by the
  kernel's init Gauss projection; O(g²) settling absorbed by fit v0 term).
- GPU kernel rejects `init=qball` (CPU-only init path) — GPU qball seeding goes
  via `init=sfa` + projection; documented here because the v69 SPEC §5.2 reads
  as if init=qball were universal.
