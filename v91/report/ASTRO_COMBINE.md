---
title: "The Substrate, In Company"
subtitle: "v91 measurement report — the ASTRO decision arms and the COMBINE campaign: which effects survive in the presence of others"
author: "SCP programme, v91 (free-cell substrate, laws_V2g verbatim)"
date: 2026-08-06
lang: en
documentclass: article
classoption: [11pt]
geometry: [margin=1.05in]
mainfont: "TeX Gyre Pagella"
monofont: "DejaVu Sans Mono"
monofontoptions: [Scale=0.82]
numbersections: true
colorlinks: true
linkcolor: teal
urlcolor: teal
toccolor: black
toc: true
abstract: |
  Two config-only campaigns on the standing v91 binary (kernel and law table
  untouched; battery ALL GREEN, 93 bars, re-verified twice). The ASTRO §6
  decision arms close the (scale × medium-temperature) confound: the frozen
  UUD chord survives both medium classes (ret 0.4114 at t=5000, cool) while
  the embedded nv=48 ring folds and starves in both (warm ret 0.107 at
  t=1500, with a fed-ember plateau), and identity assertion anti-correlates
  with survival on wound rings (frozen 364 < maintained 436 < control 600).
  The COMBINE campaign — combined experiments with an ablation-on-failure
  protocol — then attributes every standing effect to its carrier: chord
  stability belongs to the coupled lock+tuner (both required); metabolism
  and luminosity belong to the radiance drive (without it the conversion
  door collapses ~200× and the object's books stop at evaporation exactly
  zero, while the chord lives on as a static hold); the parameter-free
  doublet law belongs to laws_V2g itself. First two-body measurements:
  isolated coexistence preserves every effect with pooled luminosity
  exactly 2.00× one chord and no mutual force at range 7; contact
  coexistence adds a real, probabilistic interaction — competitive
  starvation, killed by feeding-connection loss and read live by the
  spectrometer. One standing-record error was found and corrected, dated
  (the UUD brightness baseline). Nothing was adopted; both pending
  decisions remain the user's.
header-includes:
  - \usepackage{booktabs}
---

<!-- Build: pandoc ASTRO_COMBINE.md -o ASTRO_COMBINE.pdf --pdf-engine=xelatex -->

# Provenance and method

All numbers in this report are read from committed run logs of the standing
binary `v91/freecell` (2026-08-05 build, battery-green). Both campaigns are
config-only: no kernel, law-table, or format change anywhere. Registration
preceded execution in both cases — ASTRO §6.2 committed before the first
decision arm (`3c0009a`), COMBINE §1–§4 committed before the first combined
arm (`41cffb1`). Results were committed at `9ac32b0` (DECIDE), `6b93d82`
and `f9b08b7` (COMBINE). The 93-bar battery read ALL GREEN before each
results commit (`runs/BATTERY_decide.log`, `runs/BATTERY_combine.log`).
Energy drift in every run: $\le 2\times10^{-15}$ relative.

Campaign documents: `v91/ASTRO.md` (§6 arms, §6.5 results, §7 decision
brief) and `v91/COMBINE.md` (effect inventory, arms, bars, per-experiment
results §5.1–§5.7, synthesis §6). Instruments: the QATOM spectrometer with
`spec.awk`/`spec2.awk` (per-event parameter-free line residuals,
two-triangle aware), `door.awk`/`door2.awk` (per-voice door ledgers, with
the ILAG inter-object lag meter), `prof.awk` (radial $\pi$ profiles over
snapshots). An animated companion page renders the same data: *The
Substrate, In Company* (claude.ai artifact); the visual exhibit *The
Substrate, Filmed* holds the volview plates.

# The ASTRO decision arms (§6, user-directed: "more data to decide")

The two standing stability results sat at opposite corners of a 2×2 in
(scale, medium temperature) with apparatus confounded. Four single-knob
arms filled the missing cells; a fifth added the second species; an
analysis-only arm measured the identity lane's design envelope.

## Stability: scale kills, temperature does not

| body | medium | apparatus | ret\@1500 | ret\@5000 | outcome |
|:--|:--|:--|--:|--:|:--|
| chord nv=3 | warm 0.5 | frozen lock | 0.5101 | 0.4975 | lives (reference) |
| chord nv=3 | cool 0.15 | frozen lock | 0.4258 | **0.4114** | lives — P-D2 PASS |
| ring nv=48 | cool 0.15 | none | 0.0376 | — | folds, starves |
| ring nv=48 | warm 0.5 | none | **0.1069** | — | folds, fed-ember plateau |
| ring nv=48 | cool 0.15 | frozen lock | 0.0627 | — | null, 1.67× (CO-RL) |

The warm and cool rings fold on the same schedule (rms $\approx 5.3$–$5.5$
by $t=200$–$400$) and track each other to $t\approx1000$; the cool arm then
drains monotonically to 0.035 while the warm arm arrests at a **fed-ember
plateau** of ret $\approx 0.12$–$0.15$ (books: evap/cond 1.46 warm vs 3.83
cool). Feeding sets the remnant's asymptote; it does not touch the fold or
the ~75–89% loss. Scored exactly as registered: **fold and starve → the
limitation is general → identity-lane urgency confirmed.**

On the wound nv=6 ring, escalating identity assertion *shortens* life:
$t_{\mathrm{ret}<0.25}$ = 600 (control) > 436 (maintained) > **364**
(frozen). This winding wall does **not** generalize to the embedded nv=48
(CO-RL: 456 vs 460 with the lock — the fold regime is a different failure
mode, and assertion is irrelevant to it).

## The second species and the spectroscopy edges

The UDD chord is config-reachable (`tri_kind=1`) and **lives** (ret 0.3675
at t=5000) as an open chain — exactly as the kernel's own seed comment
predicted (the D–D edge $d^{*} = \pi/w_D$ exceeds contact). Its D emission
line carries $n=61$ events at residual rms 1.5% under the same
parameter-free law with no new constants.

Two registered expectations broke, informatively:

1. **Luminosity is intake-limited at per-emitter grain.** With two
   emitters but one intake, UDD's per-emitter brightness drops to
   0.55–0.77× the UUD voice. *(The original clause "dimmer in absolute
   than UUD" rested on an erroneous baseline and is retracted — see §5.)*
2. **The healthy-D locus is not a species tag.** UUD and UDD healthy D
   voices emit at 1.6710 ± 0.0052 and 1.6774 ± 0.0071 — indistinguishable.
   Species identification runs on the **line-ratio pattern**
   (absorption:emission 2.1–2.6× apart between species), the way real
   spectroscopy reads abundances. Drained voices slide blueward; loci are
   operating points.

Pooled U-emission statistics over four seeds: $n=55$ (bar $\ge 30$, pass),
mean $-0.0035$ (pass), rms 0.063 vs the 0.05 bar (**fail as registered**):
the excess decomposes cleanly into partial-death redward broadening — the
intact-topology seed alone reads rms 0.018. Seed-robustness reframed:
matter alive 4/4 seeds at t=5000; full triangle topology 2/4.

## The identity lane's design envelope (analysis-only)

From the standing full-rate raws: a door-carried parcel gid is
load-bearing if it survives ~30–500 t.u. (the DF→FD lag quartiles) at a
traffic of ~0.04 tags/t.u. per chord — computationally trivial. The
identity-bearing observable is the door's standing *direction*
(U+, U+, D−), which runs at bath loudness: living matter is distinguished
by the direction of its flow, not its rate.

# The COMBINE campaign (user-directed: combined experiments, ablate on failure)

Protocol: put more measured effects in one box than ever before
(coexistence, the lock on the flagship), and cut the standing chord's
three co-active additions (radiance; lock+tuner; registry meter) out one
at a time. Every effect signature was scored against its registered value;
any failure triggered a back-out or replicate. Eight experiments ran
(two smokes, six arms, one protocol-triggered replicate).

## The ablation ladder

| arm | cut | first edge death | ret\@5000 | verdict |
|:--|:--|--:|--:|:--|
| CO-AL | − lock − tuner | t = 36 | 0.19 (ember, t=1500) | organism dies; **the doublet law survives** (rms 1.6%/0.9%) |
| CO-AT | − tuner (lock frozen) | t = 48 | 0.1647 | frozen gauges die with their bonds (nlock 0 by t=720) |
| CO-AR | − radiance | none | **0.4366** | chord lives; economy stops |

CO-AR is the sharpest attribution of the campaign. Without the radiance
drive the chord holds (conn 1.0, all three locks, t=5000) while:

- bath door events collapse ~200× (≈115,000 → **504**);
- the object's door fires **2 events in 4000 t.u.** (vs ~190);
- the tagged books stop: cond 0.67, **evap 0.000 exactly** (vs 29.2 ~ 31.1);
- the bath's glow band at $w\approx2.35$ disappears (bath $\bar{x}$ 0.32 → 0.20).

**Stability belongs to the lock+tuner (both required — the frozen gauge
field without its tuner is barely better than no lock at all). Metabolism
and light belong to radiance.** Door traffic behaves as one currency
$\approx$ noise × radiance drive: cool-with-radiance ≈ warm-without
(≈700 vs 504 events), and the frozen chord holds across the full ~200×
traffic range measured. The registry meter carries nothing (physics-silent,
as gated).

## Two bodies, first measurements (`exp=tri2`)

**Isolated pair (sep 8.0 = L/2).** Every effect survives company:

- joint ret\@5000 = 0.4529 (band 0.35–0.65); all six lock slots held;
- both flux machines sign-correct (U+, U+, D−; turnovers 189/162 t.u.);
- pooled D-line rate **0.0220 ev/t.u. = 2.00×** the seed-matched single —
  two full metabolisms double the light, each chord inside the single band;
- **no mutual force**: after a one-time settle contraction 8.00 → 7.07,
  the true separation (from snapshots) fluctuates 6.5–7.1 for 4500 t.u.
  with slope $-6.8\times10^{-5} \pm 4.2\times10^{-5}$/t.u. — consistent
  with zero, as the π-flat far-field closure demands;
- between-zone π flat to $\le 5\times10^{-4}$ (this window's resolution);
- ILAG (inter-object door lag): 1.4× vs the memoryless null — identical to
  the *within*-object factor; a common bath rhythm, nothing pairwise.

**Contact pair (sep 2.6, born inside contact range).** Shape coexists —
sep bounded 2.39–2.82 for the full 5000, no fusion, no runaway — but the
books do not: **competitive starvation**. In seed 20260802 the ledgers run
mirror-asymmetric (T0 net +0.92, T1 net −1.96 over the window); the
loser's apex received **1** feeding event in 4000 t.u. (the winner's U
voices: 31–33) while its D voice kept emitting at the highest rate of all
six voices, load-blind — the starving chord *shines itself to death*
(x$_D$ 0.48 → 0.10 in the last 500 t.u.; the first D-voice collapse of the
frozen-chord programme). The spectrometer reads it live: the winner stays
on the balance line (locus 1.6879 ± 0.0089), the loser's line walks
blueward (1.7815 ± 0.0276) through its starvation episodes; residual rms
still 3.7% — the law holds to the end.

The protocol-triggered replicate (seed 314159) discriminates: **both**
partners live the full 5000 (joint ret 0.5268 — the best two-chord
retention measured; sep pinned at 2.606). Contact coexistence is
survivable; the starvation risk is real and probabilistic (1 of 2 contact
pairs, 0 of 2 isolated). The kill variable is feeding-connection loss —
door-net drain alone only warns (the replicate's T1 runs net −0.00079/t.u.
and holds, apex intake 18 events).

## The attribution map

| effect | carrier (measured) | survives |
|:--|:--|:--|
| chord stability | lock + tuner, both | −radiance; cool medium; coexistence (7/8 chords) |
| doublet law | laws_V2g itself | −lock (1.6%); coexistence; the dying loser (3.7%) |
| π-flat / no far force | the law (books close at the door) | two bodies: flat $\le 5\times10^{-4}$; sep slope ≈ 0 |
| flux-machine current | radiance (drive) + lock+tuner (pattern) | coexistence: all radiance-on chords sign-correct |
| luminosity | radiance; magnitude set by intake | two objects = 2.00×; per-emitter sub-linear (UDD) |
| balance line | the medium's operating point, shared | coexistence (three chords ±0.001) |
| bath glow | radiance by construction | removed with it |
| winding wall | nv6 wound-ring-specific | does not generalize to the fold regime |

Interference census: **one record error** (§5), **zero
apparatus-on-apparatus interference**, **one genuine physical
interaction** (competitive starvation). Instrument findings: the kernel's
diagnostic `sep=` meter is wrap-broken for boundary-straddling objects
(display-only; snapshot-based true separation is the standing method);
door-net drain warns, intake disconnection kills.

# Corrections to the standing record (dated)

**ASTRO §6.5.4, corrected 2026-08-06.** The UUD windowed D-line baseline
"0.0230 ev/t.u." is not reproducible under the uniform metric (D-cell DF
events per t.u., window [1000,2000], per-run cell indices). The recompute
reads 0.0110 (seed 20260802) and 0.0150 (seed 314159). Consequences:

- the clause "UDD is dimmer with twice the emitters" is **retracted** —
  UDD (0.0170) is 1.1–1.5× *brighter* in absolute, sub-linear per emitter;
- the species line-ratio discrimination survives with updated numbers
  (UUD absorption:emission 1.13–1.36 vs UDD 0.53 — 2.1–2.6× apart);
- the corrected baseline is independently confirmed by the two-chord arm:
  pooled 0.0220 = 2.00 × 0.0110.

The error was caught by the COMBINE brightness bar failing and tracing the
failure to its source — the campaign's ablation protocol working exactly
as directed, on the records themselves.

# Standing before the user

Three closures point the same way. The config-level assertion path is
exhausted: it keeps exactly one object class alive (chords), it is null on
the embedded flagship, and it preserves skeletons — 26 of 48 gauges
rock-stable through an entire starvation — never living contents. The
contents' currency (metabolism, light, feeding) is radiance-driven flow,
and flows carry no identity. A door-carried parcel gid is computationally
trivial but presumes the radiance-driven door.

The two pending decisions, unchanged and not made here:

(a) open the **parcel-carried ontological-identity lane** (REGISTRY §5.5 —
kernel work, requires explicit authorization);

(b) re-grade audit row **B4** (ABSENT → STRUCTURAL), stated as
ratios-and-pattern on a shared locus, with the corrected brightness
numbers.

Nothing was adopted in either campaign. The ratchet governs.

# Appendix: exact arms

All commands run from `v91/` on the standing binary; seed 20260802 unless
stated. The frozen-chord reference recipe:

```
./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000
  diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0
  k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18
  reg_tau=10 qatom_every=1
```

Arms (knob deltas only): D-E1 ring48 warm (`noise_amp=0.5`); D-E2 chord
cool (`noise_amp=0.15`); D-E3 frozen nv6; D-S1 UDD (`tri_kind=1`); D-S2
seeds 271828/141421; CO-T2F/T2N two chords (`exp=tri2`,
`tri2_sep=8.0`/default 2.6); CO-AR (`k_rad=0`); CO-AT (`k_tune=0`); CO-AL
(`k_cant=0 k_tune=0`); CO-RL ring48 cool + full lock stack; replicate
`seed=314159`. Smokes: tri2 T=100; ring48+lock L=36 T=30.
