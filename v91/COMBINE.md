# COMBINE — do the effects survive each other? (user-directed, 2026-08-06)

User directive (2026-08-06, verbatim): *"Please do combined experiments
to see which effects survive in the presence of others. If one effect
fails to replicate, cut out experiment additions until it can be
determine which experiments are interfearing. After each experiment,
quantify and write down in the project the results."*

Charter, per house discipline: config-only arms on the standing binary
(`v91/freecell`, battery-green; kernel `v91/kernel/freecell.c` and
laws_V2g untouched). Every §2 command and every §3 bar in this file was
written and committed BEFORE the first §2 arm ran. Results land in §5,
one subsection per experiment, quantified, **written immediately after
each harvest** (the directive's protocol). Nothing is adopted; all
decisions are the user's.

The campaign's question, sharpened: the standing v91 discoveries were
each measured with a specific apparatus stack active. Which of the
measured effects are properties of the *matter*, and which are
artifacts of — or dependencies on — a co-active experiment addition?
Two directions of attack:

- **Combination** (§2 T-arms): put more measured objects/effects in one
  box than ever before (two chords coexisting; the lock apparatus on
  the embedded flagship) and check every standing signature against its
  registered value.
- **Ablation** (§2 A-arms): the standing frozen-chord state is itself a
  stack of three inert-by-default additions running together — radiance
  (`k_rad=0.05 p_rad=4`), the harmonic lock + tuner (`k_cant=1
  k_tune=0.2`, frozen), and the registry meter (`reg_tau=10`). Remove
  one addition at a time and see which signatures each addition was
  actually carrying. This is the directive's "cut out experiment
  additions" applied *proactively*, before any failure forces it.

## 1. The effect inventory (registered signatures)

Each row is a measured standing effect with the quantitative signature
it must reproduce to count as "survived" in a combined arm. Bars come
from the standing measurements (ASTRO §4/§6, REGISTRY §4–5); they are
signatures to CHECK, never to soften.

| id | effect | signature (pass condition) | standing value |
|---|---|---|---|
| F1 | frozen-chord stability | ret@1500 ∈ [0.36,0.66]; ret@5000 ∈ [0.35,0.65]; matter alive (no voice dead unless noted); nlock=3, a_tag=1.0 | 0.5101 / 0.4975 (seed 20260802); matter 4/4 seeds |
| F2 | parameter-free doublet law | per-event w_pred=(2.9\|1.65)/(1+1.2·x_cell); per line n≥30 where statistics allow, resid \|mean\|≤0.02, rms≤5% | warm rms 1.7% (DF); cool 0.34/0.58% |
| F3 | π-flat far field | far-zone π ripple ≤ 4e-4 around every stable body | chord ±2e-5, blob ±2e-4 |
| F5 | door flux-machine | per-voice net sign pattern U+,U+,D− (U net-absorbs, D net-emits at the door); object turnover Σ Emfl/Σ gross ∈ [100,300] t.u. | +0.0012/+0.0012/−0.0020; 167–179 t.u. |
| F6 | species by light | D emission line ≥10× bath rate in its band [1.6,1.8] | 30–60× |
| F8 | universal balance line | healthy-D emission locus ∈ [1.665,1.685] | UUD 1.6710±0.0052, UDD 1.6774±0.0071 |
| F9 | radiance glow (medium) | warm bath emission band near w≈2.35, bath xbar elevated | ASTRO §4.3 |
| F10 | winding wall | identity assertion does not rescue wound unison rings (frozen ≤ maintained ≤ ctl horizons) | 364 < 436 < 600 (nv6) |

Global CONVTAG books are deliberately NOT a scored cross-arm signature:
they are bath-dominated, and the −radiance arm changes the bath's
legitimate physics — a shifted global ledger there is not interference.
Door-grain books (F5) carry the object's books instead.

## 2. Arms (exact commands, standing binary, from `v91/`; seed 20260802 everywhere)

Reference recipe = the frozen UUD chord (ASTRO §6.2 D-E2's warm
original). Every arm below is that command ± the stated knobs. RNG
note, on record before running: `exp=tri2` draws one extra seed-time
phase (T1's apex), so its runtime noise stream diverges from the
single-chord runs — comparisons are statistical, never byte-level; the
bath carve and T0's seed state are identical.

```
# S-T2   SMOKE tri2 (exp=tri2 has never been run): T=100, formats only
#        verify: SEED ntri=2, NC=2706 (i0=2700; T0=2700-02, T1=2703-05),
#        diag "| T0 ... | T1 ... | sep=", QATOM from both triangles
./freecell exp=tri2 tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=100 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 > runs/combine/smoke_t2.log 2>&1

# S-RL   SMOKE ring48+lock (cant_seed at nv=48 never run): L=36 T=30
./freecell exp=ring ring_n=48 ring_x=0.47 bath=1 noise_amp=0.15 convtag=1 L=36 T=30 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 > runs/combine/smoke_rl.log 2>&1

# CO-T2F two frozen UUD chords, FAR (sep=8.0 = L/2, footprints are
#        contact-local so this is isolated coexistence), T=5000
./freecell exp=tri2 tri2_sep=8.0 tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/combine/t2far.fcs > runs/combine/t2far.log 2>&1

# CO-T2N two frozen UUD chords, NEAR (tri2_sep default 2.6: born inside
#        contact range — nearest inter-triangle vertex gap ≈1.15 vs
#        contacts 1.53–1.65), T=5000
./freecell exp=tri2 tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/combine/t2near.fcs > runs/combine/t2near.log 2>&1

# CO-AR  chord − RADIANCE (k_rad 0.05→0; the only knob changed)
./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/combine/ar.fcs > runs/combine/ar.log 2>&1

# CO-AT  chord − TUNER (k_tune 0.2→0; frozen lock stays)
./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/combine/at.fcs > runs/combine/at.log 2>&1

# CO-AL  chord − LOCK − TUNER (k_cant=0 k_tune=0; the positive control
#        for the ablation ladder — CANTUS says this dies), T=1500
./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=1500 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=0 k_tune=0 reg_tau=10 qatom_every=1 > runs/combine/al.log 2>&1

# CO-RL  ring48 embedded COOL + the frozen-lock apparatus
#        (= ASTRO §4.4 A-F2 + ONLY the lock stack; single-knob vs A-F2)
./freecell exp=ring ring_n=48 ring_x=0.47 bath=1 noise_amp=0.15 convtag=1 L=36 T=1500 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=200 snap_every=500 snap_file=runs/combine/rl.fcs > runs/combine/rl.log 2>&1
```

Registered as NOT config-reachable (kernel work — not opened here): a
**mixed-species pair** (UUD + UDD in one box). `tri_kind` is global to
both tri2 triangles. Goes on the identity-lane design list.

Instruments: `runs/combine/spec2.awk` and `runs/combine/door2.awk` —
tri2-aware versions of the standing spec/door analyzers (parse both
per-triangle xUDD tuples; door2 adds ILAG, the inter-triangle door lag:
FD on one triangle vs the most recent DF on the other — the direct
event-grain meter for cross-object exchange).

## 3. Registered predictions and bars (BEFORE any arm runs)

**CO-T2F (coexistence, isolated).** Prediction: every effect survives —
the far field's measured absence (ASTRO §4/B6) predicts **no mutual
force**: two π-flat bodies have nothing to pull on each other with.
Bars: (a) both chords matter-alive at 5000, joint ret@5000 ∈
[0.35,0.65]; (b) F2 per chord: D-line resid rms ≤5%; (c) F8 both
healthy-D loci ∈ [1.665,1.685] and |Δ| ≤ 0.02 between chords; (d) F5
sign pattern per chord; (e) F6 pooled D-line ev-rate = 2× single
(window [1000,2000]: 0.036–0.052 vs single 0.0230) — the §6.5.4
intake-limited-brightness claim's clean test: two FULL metabolisms must
double the light that UDD's two shared-intake emitters could not; (f)
**no mutual force**: |sep(t) − 8.0| < 0.5 for all t (drift beyond that
= a medium-carried interaction the far-field result says cannot exist);
(g) F3: π ripple ≤ 4e-4 in the zone between/around the pair at t=5000.

**CO-T2N (coexistence, in contact).** Genuinely open — registered as a
classification, not a prediction: (a) **fusion** — inter-triangle bonds
form (U–U contact at d*_UU=1.447 is an available unison rung; report
the edge census and whether a 6-voice molecule holds); (b) **mutual
disruption** — either triangle dies before 5000 (this seed lives single
— death here = interaction-caused); (c) **contact coexistence** — sep
locks, both alive, spectra intact. Weak prior on record: fusion is
plausible via U–U unison; a closed 6-cycle would newly expose the
winding wall (F10) at an object built entirely from chord parts —
watch for it. Scored either way: matter-alive count, F2 on survivors,
sep(t), ILAG (does the pair exchange at the door more than chance?).

**CO-AR (− radiance).** Prediction: the chord's signatures do NOT
belong to radiance — F1 in band, F2 rms ≤5% with per-line resid means
within ±0.02 of the standing arm (the law has no k_rad term; if loci
move, radiance was hiding in the "parameter-free" constants), F5
pattern holds. The bath changes legitimately: F9 glow is REMOVED by
construction — bar: bath 2.35-band rate drops ≥3×; F6 contrast should
RISE (darker bath under the same D line) — report the new ratio.
FAILURE branch (chord dies or law breaks): radiance is load-bearing
for the goal object — a hidden dependency the isolation-era
measurements could not see; ablation cascade per §4.

**CO-AT (− tuner).** Open, two live hypotheses on record: CANTUS
measured lock-without-tune dying for *growing pairs*; whether a
FROZEN lock still needs its tuner at hold-time is unmeasured. Bars: F1
band = tuner not load-bearing for a frozen chord (then the standing
recipe carries a removable part); death/leaving band = lock–tuner
co-dependency extends to the frozen state (tune=… line stays 0.000 —
verify the tuner really is off). F2/F5 scored on whatever survives.

**CO-AL (− lock − tuner).** Expected DEATH (positive control): ret <
0.15 by t=1500 with voices dying. If it LIVES, the lock was never
load-bearing for chord stability and the CANTUS attribution is wrong —
that would be the largest finding of the campaign; immediate re-check
against REGISTRY ctl arms before believing it.

**CO-RL (lock apparatus × embedded flagship).** The two standing
claims COMPETE here, registered as a discriminator: the winding wall
(F10: assertion hurts wound rings — frozen died earliest at nv6) vs
bond-holding (the frozen apparatus keeps chord bonds alive through
5000 embedded). Bars against the single-knob reference A-F2 (cool,
no lock: fold rms 5.32@200, ret@1500 0.0376, t_{ret<0.25} ≈ 600):
fold expected to persist (rms ≤6 by t=400 — the fold is spatial, no
lock opinion); then ret@1500 **> 0.075** (2×) = the lock materially
helps embedded large matter → a no-kernel-work path the §7 decision
must weigh; ret@1500 **< 0.019** (½×) = the lock actively kills at
nv48 too → winding wall generalizes; between = null (assertion
irrelevant at this scale — the identity-lane case stands unchanged).

## 4. The ablation protocol (the directive's, made concrete)

If any effect signature FAILS in a combined arm while its standing
single-effect reference passes: re-run the combined arm minus one
addition at a time (order: registry meter → tuner → radiance → lock;
nearest-to-standing first), T=1500 probes, same seed, until the
signature returns. Report the minimal interfering set. Bars are never
softened to make a combination pass; an interference, once isolated,
is a RESULT (written in §5 with its pair identified), not a nuisance.

## 5. RESULTS (one subsection per experiment, written at harvest)

(registered empty; filled in run order)
