# v70 PLAN — Skeptical verification of the v66–v69 headline claims

**Date**: 2026-06-11. **Motivation**: the v69 flagship claim ("the long-range force
is Coulomb") rests on two clean data points (fl20, fl26) whose magnitude ratio
drifts monotonically (0.92 → 1.17) — a possible systematic, and no exponent is
measurable from two points. Several other claims (positronium non-annihilation,
fission saddle, ℏ_eff=Q) rest on scalar diagnostics only. v70 (a) re-measures the
force law with a lever arm and controls, (b) adds full-field VISUAL confirmation
of every spatial claim, (c) runs the deferred ℏ_eff=Q phase-tilt fingerprint.

## A. Visual confirmation pipeline (no new GPU time)

New tool `sfa/analysis/sfa_slice.c` (+ `render_slices.py`): complex/gauge-aware
slice + lineout extraction (rho2, s, rhoQ, th2, |E|) from 12/24/30-column SFAs.
Applied to the banked v69 SFAs: fl14/fl20/fl26 (force-law runs), pos1
("positronium"), fiss1 (fission saddle). Checks: two intact symmetric balls (no
radiation-wake centroid bias), Coulomb halo structure, dipole vs orbit identity,
spherical-saddle persistence, boundary artifacts.

## B. Phase-tilt fingerprint (banked bs03/05/07, v68 relativity ladder)

`v70/analysis/phase_tilt.py`: measure the in-core phase tilt k = dθ/dx and the
actual translation velocity v from the same SFA; test k = γωv and
ℏ_eff = p/k = E₀/ω = Q(1+ε) (Pohozaev ε = 0.032). This converts the v67
"the action quantum is the charge" interpretation into a measured number.

## C. Force-law lever-arm campaign (GPU, 9 runs)

Geometry: N=192, L=36 (dx = 0.3770, matches v69 fl runs), V100, f16 snaps,
snap_dt=5, diag_dt=0.5, T=150. Physics: ω=1.42 g=0.05 profile
(`v69/theory/gprofile_w142_g005.txt`, Q=311.5, E=456.3), η=0, m_θ=1.6,
absorbing BC. Seeds: `gen_qball_pair`/`gen_qball_boost` (24-col matter seed,
E built by the kernel's Gauss projection; O(g²) settling absorbed by a v0 term
in the fit, t<20 dropped).

| run | seed | purpose |
|---|---|---|
| ctrl1 | single ball at x=+10 | lattice/sponge drift floor (a_sys) |
| fl2_d20/26/32/40 | same-charge pairs at ±D/2 | repulsion a(D) |
| flo_d20/26/32/40 | opposite-charge pairs (qball2 ω=−1.42) | attraction a(D); net-neutral (no jellium offset) |

Estimator: a_C(D) = (a_same − a_opp)/2 — cancels additive systematics (drift,
sponge asymmetry); contact-force residual ≤7% at D=20, negligible ≥26
(F_t ∝ e^{−0.564D}). Fit: D(t) quadratic with v0 term, residual bootstrap;
exponent n from weighted log-log over D ∈ {20,26,32,40}.

**Pass gates**: n = 2.00 ± 0.15; coefficient C within ~15% of the
parameter-free 2g²Q²/4πM; |a_same| ≈ |a_opp| at matched D; ctrl1 drift
≪ a_C(40) ≈ 5×10⁻⁵.

**Falsifiers**: n far from 2 ⟹ the "Coulomb" headline is wrong (image charges /
finite-box / non-1/r mediator); a_same ≉ −a_opp asymmetry beyond contact-force
scale ⟹ unmodeled systematic; ctrl1 drift ≳ signal ⟹ all v69 single-frame-count
force numbers suspect.

## Files

- `sfa/analysis/sfa_slice.c`, `sfa/analysis/render_slices.py` — visual pipeline
- `v70/analysis/phase_tilt.py` — fingerprint test B
- `v70/analysis/forcelaw_fit.py` — campaign fit C
- `v70/cfg/common.cfg.tmpl` — run template
- `v70/results/` — diag.tsv + tracker TSVs + fit outputs; SFAs in /space/scp/v70/
