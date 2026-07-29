# PRE-REGISTRATION — Phase-2 candidates vs the fitted leak/death laws

*2026-07-28, registered BEFORE any Phase-2 run. Source: REGRESS.md fits over
the 23-run mass corpus (9 uncensored deaths). Frozen here; score against
this file, do not refit after seeing the data.*

## The laws being scored

**Law A (death vs settled load; primary, CV R² = 0.97):**

    t_death = 274 · (x50 / 0.0617)^1.066          [t.u.]

**Law B (equivalent universal-current form):**

    t_death = 50 + 2.5·(x50 − 0.0617)/4.25e-4
    dM/dt(whole-life) ≈ −c0·N,  c0 = 4.25e-4 Em/t.u./voice  (±13% spread)

x50 = census M_sum(t=50)/(N·cap). **The prediction is the curve, not the
point**: seeding moves x50 (foam inflates loops/cubes), so score each sim by
its *measured* x50 at t=50. The x50 values below are best guesses from the
closest measured seeds. Late-window |dM/dt| prediction: 0.5–0.75·c0·N
(deceleration factor measured 0.45–0.75 across the corpus).

**Uncertainty band: ×/1.5** (fit σ_ln 0.077 would say ×/1.08, but foam chaos
is ±30% and the a2 ensemble contributed only one death — quoting the fit σ
alone would be dishonest).

Scoring curve (t.u.):

| x50  | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.50 | 0.60 | 0.70 |
|------|------|------|------|------|------|------|------|------|------|------|
| Law A|  458 |  706 |  959 | 1217 | 1477 | 1741 | 2008 | 2546 | 3092 | 3645 |
| Law B|  275 |  570 |  864 | 1159 | 1453 | 1747 | 2042 | 2630 | 3219 | 3808 |

The registered exception class: **wound mutual exact rings** (comp12-like)
beat Law A by ≥2.4× (comp12 bound: c_eff ≤ 1.68e-4 = 0.40·c0, alive at
5000 vs 1981 predicted).

## Candidates

### (a) Cube shell re-tuned to gates min ≥ 0.95 (vs the measured mean-0.6 cube, died 1811)

Assume the heavy-cube geometry (abar ≈ 1.586 → **x50 ≈ 0.375**).

* **Predicted: t_death ≈ 1875 (band 1250–2810); late |dM/dt| ≈ 1.7–2.5e-3;
  c_eff = 3.5–4.6e-4.**
* This is the sharpest test of the headline null: across shells v1/v2/v3 and
  all rings, seeded gate quality bought NOTHING (c0 flat from gmin 1.5e-5 to
  1.0). The consonant-skin hypothesis says perfect gates stop the leak.
* Falsification lines: t_death within the band → gate-retune buys <1.5× and
  the consonant-skin route (as a rate-level fix) is dead. t_death ≥ 4700
  (≥2.5×) or a plateau (|dM/dt| ≤ 1e-4 sustained) → the skin mechanism is
  real and the load-line null breaks. **Registered expectation: ON the load
  line (the null holds).** Confidence: this is the prediction I am most
  confident in.

### (b) Hexagonal prism N=12, struts-only

Assume ring_m-exact seeding at comp12-class load, **x50 ≈ 0.40** (score at
measured).

* **Predicted: ON the load line — t_death ≈ 2010 (1340–3010); late |dM/dt|
  ≈ 2.5–3.8e-3 (N=12); c_eff = 3.5–4.6e-4.**
* Rationale: mutual unwound closure (unwound12) sat exactly on the line;
  adding a second ring plane + struts without winding adds closure but no
  chirality. 3D-ness alone is predicted worthless at rate level.

### (c) Wound tube — two stacked co-rotating rings (comp12 links + exact m=1 axial rungs), N ≈ 24

**x50 ≈ 0.40** assumed.

* **Predicted: the exception class — BEATS Law A by ≥2.3×: t_death ≥ 4600,
  plausibly censored at T=5000; per-voice current ≤ 1.7e-4 (i.e. total late
  |dM/dt| ≤ ~0.5·c0·24 ≈ 5e-3, decreasing); late phase = parked arcs with
  no secular decay, as comp12.**
* This is the strongest *structural* bet: it doubles down on the single
  measured exception (wound + mutual). If the tube lands ON the load line
  (≤3000), the comp12 excess was a foam accident and winding-mutuality is
  not a protection mechanism — that outcome would kill route R2 at rate
  level. Registered odds: exception real, ~2:1 (single-seed evidence).

### (d) ring12 + 2 consonant chords

**x50 ≈ 0.40** assumed.

* **Predicted: between the line and the comp12 excess — t_death 2000–3400
  (point estimate ≈ 2600); c_eff 2.5–4.3e-4.**
* Rationale: chords add mutual closure paths (redundant cycles) without
  winding. If closure redundancy alone (not chirality) is what protects
  comp12, this lands high; if chirality is required, it lands on the line.
  Either outcome cleanly splits the comp12 mechanism — that split is the
  point of this arm.

### (e) Torus net, N ≈ 24

**x50 ≈ 0.35** assumed.

* **Predicted: ON the load line — t_death ≈ 1740 (1160–2610) — and the
  first N-scaling test of the current: total late |dM/dt| ≈ 0.5–0.75·c0·24
  ≈ 5–7.7e-3, i.e. per-voice current unchanged at N=24.**
* If the closed 2-surface suppresses the per-voice current (surface-tension
  clause), it beats the line — same discriminator as (a); a torus that
  merely *dies slower in % terms* while c_eff stays 3.5–4.6e-4 is NOT a
  success, it is the load line at higher mass.

## Summary of what each outcome decides

| outcome | verdict |
|---|---|
| (a) in band | consonant-skin (rate-level) dead; load line rules shells |
| (a) ≥2.5× or plateau | skin mechanism real; biggest possible upset |
| (c) ≥2.3× line | winding+mutuality confirmed as THE protection lever |
| (c) on line | comp12 excess = foam accident; no structural lever known |
| (d) high vs on-line | splits closure-redundancy vs chirality as the comp12 mechanism |
| (e) c_eff ≈ c0 at N=24 | current is per-voice, extensive in N (law hardened) |

*Single most confident statement:* a Phase-2 object that is not a wound
mutual ring will die at **t_death = 274·(x50/0.0617)^1.066 ×/1.5** with
whole-life per-voice current 3.5–4.6e-4 — regardless of its seeded gate
quality, N, box closure, or pocket pressure.
