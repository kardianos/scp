# X10-pathfinder results — shell ringdown / existence test (GPU)
**Date:** 2026-07-22 · **Status:** run truncated at t≈1700/2000 (85%) by
instance failure; verdict frames secured. Data: `/space/scp/v85/topo1/gpu/`
(run3_*: frames 0/500/1000/1500, diag → t=1703.6, logs).

## Design (per ANALYSIS §2C, adapted)
Ball ω=1.4325 (ungauged profile → seeded Q=160.6; α_f = g²Q/4π = 0.032,
**a₀ ≈ 21**, ε₁ ≈ 7.6e-4, 1/ε₁ ≈ 1300 t.u.) + opposite-charge perturbation
(amplitude ×0.15 of the profile, Q ≈ −3.7, seeded at r = 21, ω = −1.4994).
N=384, L=55, g=0.05, T=2000, absorbing BC. V100-32GB (Vast).

## Campaign notes
- **Kernel fix required and applied (authorized):** `init_gauss_project` in
  scp_sim.cu was a *serial* host CG — 1362 iterations ≈ 12 s each ≈ 4.8 h at
  N=384. Parallelized with 8 OpenMP pragmas (host-side only; no physics
  change) → ~48 min at ~7× (memory-bandwidth-bound on 72 cores).
- **Two Vast instance failures**, both on the ssh7.vast.ai proxy host: gpu1
  died pre-physics; gpu2 died at t≈1700 (85%). Hourly local pulls (user
  directive) preserved everything of value both times.

## Measured: the four-frame negative-charge story

| t | total Q₋ in box | profile |
|---|---|---|
| 0 | −3.7 (seeded) | perturbation at r=21 |
| 500 | −0.967 | **shell peaked at r ≈ 20–24 (= a₀)** + outbound tail r>38 |
| 1000 | −0.504 | transient outward slosh (broad max r ≈ 28–32; < 1 binding time — settling, not evaporation) |
| 1500 | −0.331 | **peak returned to r ≈ 20–24 (= a₀), amplitude ≈ unchanged from t=1000**; losses came from the outer structure (r=28–36 bins 0.037→0.013) |

Decay of the total: τ = 765 (500→1000) → **1190 (1000→1500)** — strongly
decelerating. Core: **zero negative charge inside r = 12 at every frame**;
ball intact (+160.5, static profile). Gauss at 1.8e-13 floor throughout;
E drift −2.28% (the escaping unbound fraction), θ silent.

## Verdicts

| Question | Verdict |
|---|---|
| **F11 — does a discrete bound response structure exist?** | **PASS (existence).** A negative-charge shell formed at the *predicted* Bohr radius a₀ ≈ 21, survived ≥1500 t.u. (≈1.15 binding times), and re-concentrated onto a₀ after a sub-binding-time transient — mode behavior, not dispersal. |
| **F10 — does the cloud drain into the core?** | **Did not fire.** No negative charge within r=12 at any time; annihilation channel closed on this timescale (detuning protection sufficient; flavor protection wasn't even needed). |
| Shell lifetime | **Remanded — BC-confounded.** The mode's evanescent tail length 1/κ ≈ 21 vs sponge clearance ≈ 20: the absorbing boundary continuously eats the outer mode structure (ARTIFACT(L) per the D7 classifier). The observed drain is dominated by box loss of loosely-bound components; the surviving a₀-peaked core's intrinsic lifetime is not measurable in this box. |

## SHELL thread status after X10
The atom architecture's central bet — **orbitals as bound response modes of
the same fabric** — now has direct dynamical evidence: right radius, right
protection phenomenology (core untouched), mode-like settling. What it does
not yet have: a BC-clean lifetime, mode frequencies (needs T ≈ 1.5×10⁴ per
the ANALYSIS error budget), or occupancy structure (D13 untouched).

## Follow-up options (in value order)
1. **X10b — deeper shell, BC-clean (recommended):** Q_N ≈ 450 at g=0.05
   (α ≈ 0.089, ε₁ ≈ 6e-3, tail ≈ 9, a₀ ≈ 7.5, ρ = a₀/r_half ≈ 1.8):
   marginal scale separation but ~4.5 tail-lengths of sponge clearance in
   the same N=384/L=55 box — lifetime and settling become physical
   measurements. Same seed recipe, ~20 h V100.
2. **L-scan (D7 certification)** on the pathfinder design — 2–3 boxes,
   expensive at N≥384.
3. **Tail completion** of this run from the t=1500 frame (T+500) — low
   value; the curve's shape is established.
