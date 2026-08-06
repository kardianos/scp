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

**5.1 CO-AL (chord − lock − tuner) — the positive control kills the
ORGANISM as expected, and the ablation cleanly attributes F2 to the
law, not the apparatus.** (harvested 2026-08-06, first arm home)

Death of the organized object, quantified: first channel death at
**t=36** (the unlocked triangle cannot hold a gate from the start);
the D voice — the metabolic engine — drains 0.837→0.017 by t≈300;
ret crosses the 0.15 death bar at t=560; at close 2/3 edges are
CHANNEL DEAD, ggm=0, x_D=0.001. The registered death bar (ret<0.15
by 1500) PASSES (crossed at 560). What remains is not nothing: in a
WARM medium the corpse is a bath-coupled ember whose ret fluctuates
0.09–0.26 (close 0.193) as tagged matter re-exchanges with the bath
— organization dies, stuff keeps trading. (Bath activity unchanged:
bath DF rate 11.4/t.u. vs frozen-arm 11.0 — the apparatus knobs
never touched the medium.)

**Effect-survival scoring:** F1 dead (as required for the control);
**F2 SURVIVES the lock's removal** — the parameter-free doublet law
holds on the dying remnant at resid rms **1.6% (DF, n=14) / 0.9%
(FD, n=20)**, indistinguishable from the frozen reference (1.7%/
1.4%). The law rides each cell's own load wherever a conversion
fires; the lock's real contribution is keeping the OBJECT alive
(gates, D-loading, directed flux) — not making the spectrum lawful.
(ASTRO's dying control showed DF rms 4.5% — that broader figure
accumulates over its full 5000 t.u. death tail; in AL's 1500 t.u.
window the law is tight. Event rate is the other casualty: 14 DF in
1300 t.u. vs the frozen arm's 80 in the same window — unlocked
matter is ~6× dimmer with no standing direction.) F5 dead with the
organism (no directed pattern to score).

**5.2 CO-T2F (two frozen UUD chords, isolated coexistence) — every
effect survives company; one registered bar failed and the failure
traced to a WRONG BASELINE in the standing record, not to physics.**
(harvested 2026-08-06)

State at 5000: joint ret = **0.4529** (bar [0.35,0.65] ✓), both D
voices strong (x_D 0.557 / 0.575), all six lock slots held
(nlock=6, a_tag=1.0). Degradation at the single-chord ensemble's own
rate: T0's U pair drains to ~0 by close and T1 leans (conn 0.500) —
2 degraded triangle-instances across the two t2 arms' 4, matching
the 2/4 topology fragility the SINGLE chord shows across seeds
(ASTRO 6.5.6). No interference claim is supportable from stability.

Effect-survival scoring, per registered bar:
  - **(d) F5 door flux-machine: PASS both chords.** Sign patterns
    (+0.00116,+0.00094,−0.00251) and (+0.00119,+0.00139,−0.00251);
    turnovers 189 / 162 t.u. (band [100,300]). Two flux machines run
    their standing directed currents side by side.
  - **(b) F2 doublet law: PASS with the known degradation caveat.**
    FD rms 1.5%/1.9% both chords; DF rms T0 2.7% ✓, T1 5.3% ✗ —
    decomposed: T1 early [1000,2500) rms **0.52%** (n=15), late
    **6.1%** with redward mean −0.012 (n=47) = ASTRO 6.5.5's
    partial-death broadening as T1's U voices lean out late. The law
    is tight wherever topology is intact — in company as alone.
  - **(c) F8 balance line: universality PASS, and sharpened.**
    Between the two coexisting chords |Δlocus| = 0.0011 (full
    window) / 0.014±0.019 (healthy window) — indistinguishable. But
    BOTH loci sit 0.01–0.02 above the isolated-single references,
    moving TOGETHER: the healthy-D locus is a shared OPERATING POINT
    set by the box/medium state (the law's w(x_D) at the common
    x_D), not a universal constant — cohabitants share it; the
    law-referenced residual is the true invariant (≈0 everywhere).
  - **(e) F6 brightness: PREDICTION CONFIRMED — against a CORRECTED
    baseline.** Pooled D-line rate [1000,2000] = **0.0220** ev/t.u.
    = **2.00×** the seed-matched isolated single (0.0110); each
    chord emits inside the single-chord band (0.0130 / 0.0090 vs
    singles 0.0110–0.0150). Two full metabolisms DO double the
    light. The registered numeric band (0.036–0.052) is void — it
    was anchored on ASTRO 6.5.4's "UUD 0.0230", which a uniform
    recompute cannot reproduce (see the correction recorded in
    ASTRO 6.5.4, dated 2026-08-06: UUD singles measure 0.0110 /
    0.0150 D-line, so UDD's 0.0170 was never "dimmer than UUD" —
    sub-linear per emitter, yes; dimmer in absolute, no). This arm's
    "failed replication" is the campaign's first isolated
    interferer: an erroneous standing number, caught exactly as the
    protocol intends.
  - **(f) mutual force: bar FAIL as registered; the protected CLAIM
    survives.** True separation (from snapshots — the kernel's
    diag `sep=` meter is wrap-broken when a triangle straddles the
    boundary and read 2.67 for a true 8.00; display-only, kernel
    untouched; fcs-based true-sep is now the standing method): one
    one-time settle-era contraction 8.00→7.07 by t≈500 (axial),
    then 4500 t.u. of bounded fluctuation 6.5–7.1 with slope
    **−6.8e-5 ± 4.2e-5 /t.u. (1.6σ, consistent with zero)** — no
    persistent attraction between two stable bodies at range ~7,
    as the π-flat far-field result demands. The |sep−8|<0.5 bar
    failed only because it did not budget the settle displacement.
  - **(g) F3 π between the pair: PASS at this window's resolution.**
    Between-zone profile (midpoint fixed center, [4000,5000], 101
    frames): π deviations ≤ ±5e-4 ≈ 1 sem bin-to-bin, level 1.0418
    = the global bath π; no resolvable structure between the
    bodies.
  - **ILAG (new meter): no EXTRA cross-object door coupling.**
    Cross-triangle FD-after-DF medians 30/35 t.u. vs ~45–48
    memoryless nulls = 1.4× — the SAME factor as the within-object
    lag structure (ASTRO 6.5.1's 1.4–2×); a common bath rhythm
    explains both; nothing pairwise resolved at range 7.

**5.3 CO-AR (chord − radiance) — the sharpest attribution of the
campaign: STABILITY belongs to the lock, METABOLISM and LIGHT belong
to radiance.** (harvested 2026-08-06)

**F1 PASS exactly as predicted**: ret@1500 = 0.4399 ∈ [0.36,0.66],
ret@5000 = **0.4366** ∈ [0.35,0.65], conn = 1.000, nlock=3,
a_tag=1.0 at close. The chord does not need the radiance field to
live.

Everything else went DARK — and that is the attribution, not a
failure of the chord:
  - **The conversion door is radiance-driven, ~200×.** Bath door
    events collapse from ~115,000 (radiance on, warm) to **504**
    over the same window; the object's own door fires **2 events in
    4000 t.u.** (vs ~190). Global conversion ledger: cond 688 with
    rad=0.000, vs cond 21,243 + rad 23,914 radiance-on.
  - **The object's books stop: tagged cond 0.67, evap 0.000
    (exactly), net +0.35** — against the standing balanced books
    (29.2 ~ 31.1). Without radiance the frozen chord is an
    economically STATIC hold — bonds and matter kept, no standing
    exchange current. F5 (the flux-machine sign pattern) is not
    violated; its currency is absent. **The metabolizing-matter
    picture (eat-at-U/shine-at-D, turnover 167–179 t.u.) is a
    radiance-era phenomenon at the door grain.**
  - **F2/F6/F8 unscoreable — the ablation removed the observable.**
    n=1 events per line class (residuals unquotable at n=1). No
    lines without the drive; the "loci shift?" question is moot. F9
    confirmed removed: the bath's 2.35 glow band is gone (residual
    sparse bath DF at 2.17–2.27), bath xbar 0.32→0.20.
  - **Operating mode: starved-lean even in a WARM bath** (U voices
    0.016/0.025, D 0.569) — U-side feeding is itself
    radiance-borne. Cross-arm unification: cool-with-radiance (D-E2:
    ~700 bath events) ≈ warm-without-radiance (504) — door traffic
    behaves as ONE currency ≈ noise_amp × radiance drive, and the
    frozen chord holds its band across the full ~200× traffic
    range measured.

No ablation cascade is needed: the single-knob difference IS the
isolated interferer. Design note for the identity lane (ASTRO §7a,
REGISTRY §5.5): a door-carried parcel gid presumes door traffic,
which this arm measures as radiance-conditional — at V2g defaults
the door is near-silent, so the lane's identity currency must
either ride the radiance-driven door or be carried below the atom
grain.

**5.4 CO-T2N (two frozen UUD chords, in contact) — contact
coexistence for shape; COMPETITIVE STARVATION for the books: the
first D-voice collapse of the frozen-chord programme, with
mirror-asymmetric ledgers.** (harvested 2026-08-06; replicate
launched per §4 protocol)

Registered classification lands on **(c) contact coexistence** for
geometry — sep oscillates 2.39–2.82 around ~2.55 for 5000 t.u., no
fusion into a six-voice molecule, no runaway — with an outcome the
classification didn't name for the books: **T1 dies of starvation
at close** (xUDD 0.005/0.007/**0.101** at 5000 — every prior frozen
degradation kept D ≥ 0.5; this is the programme's first D
collapse), terminal in the last ~500 t.u. (0.48@4500 → 0.10@5000;
one lock slot lost at t=4376; T0 healthy-lean 0.066/0.084/0.587;
joint ret 0.3049).

The door ledger [1000,5000] identifies the mechanism:
  - **Mirror-asymmetric books:** T0 net **+0.00023**/t.u. (+0.92
    over the window) vs T1 net **−0.00049**/t.u. (−1.96 ≈ its whole
    D store). The isolated single runs +0.0004.
  - **The loser is cut off from food, not robbed at the door:**
    T1's apex U voice received **1 feeding event in 4000 t.u.**
    (T0's U voices: 31–33); T1's D kept emitting at the full
    load-driven rate (0.0132/t.u., the highest of all six voices)
    until the store emptied — D emission is load-blind to intake,
    so the starving chord SHINES ITSELF TO DEATH.
  - **No door-to-door coupling even at contact:** ILAG 31.5/41.2
    t.u. vs ~46 nulls = the same 1.4× common-mode factor as the
    far pair and the within-object structure — the competition is
    in the shared feeding field, not in event-correlated exchange.
  - **F2 holds through the death:** T0 DF rms 2.1%, T1 DF rms 3.7%
    (n=60, mean −0.007 mild-redward) — parameter-free to the end.
  - **F8 as operating point, again:** the three healthy coexisting
    chords across both t2 arms sit on ONE box locus — 1.6877,
    1.6888, 1.6879 (±0.001) — while the starving T1's locus walks
    blueward (1.7815±0.0276), the drained-voice signature tracking
    its falling x_D. The spectrogram reads the competition in real
    time: winner on the balance line, loser walking blueward.

Interference verdict, stated with its statistics: contact-range
coexistence degrades ONE partner's books into net drain while the
other holds net gain — consistent with competition for the shared
contact-zone deliverable flux (the intake-limited-luminosity
economy at two-body grain). n=1; the seed-314159 replicate
(t2near_s314159, launched at harvest per the §4 protocol) is the
discriminator between "unlucky ensemble instance" and "contact
starves a partner"; its result lands in §5.6.

**5.5 CO-AT (chord − tuner, frozen lock kept) — the tuner is NOT a
removable part: frozen gauges without their tuner are barely better
than no lock at all.** (harvested 2026-08-06)

Registered question: does a FROZEN lock still need its tuner at
hold-time? Answer — yes, absolutely:
  - **Organization dies on the no-lock schedule:** first edge
    CHANNEL DEAD at **t=48** (the no-lock control AL: t=36); x_D
    crashes below 0.3 by **t=72**. The gauges then die WITH their
    bonds — nlock 3→2 by t=300, **0 by t=720**, a_max=0.0
    thereafter (a gauge lives on its pair slot; edge death takes
    the gauge with it — "frozen" protects against tau-decay, not
    against substrate loss). tune=0.000000 verifies the tuner
    really was off.
  - **Stability out of band:** ret 0.2211@1500, 0.1647@5000 (bar
    [0.36,0.66]/[0.35,0.65]), conn 0.333 — a slower slide than AL
    (some bond persists longer) but the same destination: the
    lean-ember mode, not a living chord.
  - **F5 erased:** the U+,U+,D− sign pattern never establishes
    (measured nets −0.00045/+0.00025/−0.00051; the D voice
    window-averages x=0.083 — drained from t≈72). No standing
    directed current without the tuner.
  - **F2 wears the death signature:** DF rms 6.6% with mean −0.026
    redward over the long decline (the 6.5.5 broadening); FD stays
    tight at 2.4%.

Attribution: what CANTUS measured for growing pairs — coupled
lock+tune holds, each half alone dies — now holds at the FROZEN
chord limit. The tuner maintains the phase alignment the gates
need; the frozen gauge field alone neither locks nor holds. The
standing "frozen chord" recipe is therefore irreducibly
lock+tuner+medium-drive (and §5.3: its metabolism additionally
needs radiance). Nothing to cascade: single-knob, cleanly
attributed.

**5.6 CO-T2N replicate (seed 314159) — the discriminator lands on
"unlucky instance": BOTH contact partners live the full 5000.**
(harvested 2026-08-06)

Joint ret@5000 = **0.5268** (the best two-chord retention measured,
above the isolated pair's 0.4529), conn = 1.000, all six lock slots
held, sep pinned at 2.606, both chords healthy-lean (x_D 0.604 /
0.577; U voices 0.05–0.10). Both flux machines sign-correct
(T0 +0.00118/+0.00171/−0.00221; T1 +0.00079/+0.00059/−0.00217).
ILAG 46/39 t.u. vs ~55–62 nulls = the same 1.2–1.6× common-mode
band. Contact-pair verdict across the two seeds, stated honestly:
contact coexistence is fully survivable (this seed) and carries a
real but probabilistic starvation risk (s20260802's partner death —
1 of 2 contact pairs, 0 of 2 isolated pairs, and the only D-voice
collapse anywhere in the frozen ensemble). The mechanism variable
is the FEEDING CONNECTION, not the door-net sign: this seed's T1
runs door-net −0.00079/t.u. yet holds (apex intake 18 events),
while the dead T1's apex had been cut to 1 feeding event in 4000
t.u. Door-net drain is a warning light; intake disconnection is the
kill.

**5.7 CO-RL (frozen-lock apparatus × embedded nv48) — the
registered NULL branch: assertion preserves STRUCTURE and is
irrelevant to the STARVATION; the winding wall does not generalize
to the fold regime.** (harvested 2026-08-06)

Scored against the single-knob reference A-F2 (cool, no lock):
  - Fold: untouched by the lock — rms 12.95→5.38 by t=200
    (reference 5.32), 5.14@400 (bar ≤6 ✓), 5.29 at close.
  - Retention: ret@1500 = **0.0627** = 1.67× the reference 0.0376 —
    inside the registered null band [0.019, 0.075], at its upper
    edge; starvation onset IDENTICAL (t_{ret<0.25} = 456 vs 460);
    books evap/cond 4.41 vs 3.83 (no improvement).
  - The structure/contents split, sharpest yet: **26 of 48 gauges
    survive the fold and hold rock-stable from t≈400 to the end**
    (nlock 48→26 by 400, then flat; a_max=1.0; 28 edges alive at
    close, 26 of them locked) while the matter drains through the
    held structure on the no-lock schedule. Asserted bond identity
    keeps the SKELETON; it cannot keep the CONTENTS — the flows
    that carry the contents carry no identity (REGISTRY §5), and
    the lock has no handle on them.
  - **F10 does NOT generalize:** the nv6 winding wall (frozen 364 <
    ctl 600 — assertion kills wound rings) does not appear at
    embedded nv48 (456 ≈ 460): the embedded ring fails by
    fold-and-starve, not by winding, and assertion is simply
    irrelevant to that mode. The winding pathology needs the
    winding.

Decision impact (ASTRO §7a): the no-kernel-work path — "save the
flagship with the lock apparatus" — is now measured CLOSED (1.67×,
null band). The identity-lane case is unchanged-to-strengthened:
config-level assertion preserves skeletons, not living contents, at
every scale tested.

---

## 6. SYNTHESIS — which effects survive in whose presence (the campaign's answer)

The attribution map, measured by combination and ablation
(campaign verdict; decisions the user's):

| effect | carrier (measured) | survives |
|---|---|---|
| F1 chord stability | the coupled lock+tuner — BOTH (AL: dead t=36; AT: dead t=48; frozen gauges die with their bonds) | radiance removal (AR 0.4366@5000); cool medium (D-E2); coexistence far AND contact (7 of 8 coexisting chords lived) |
| F2 doublet law | **laws_V2g itself** — no apparatus carries it | lock removal (1.6%); coexistence (≤5% wherever topology intact); even the starving loser (3.7%). The per-event residual is THE invariant of the spectroscopy |
| F3 π-flat / no far force | the law (stable books close at the door) | two-body: between-zone flat ≤5e-4; sep slope −6.8e-5±4.2e-5 ≈ 0 at range 7 over 4500 t.u. |
| F5 flux-machine current | **radiance (the drive) + lock+tuner (the pattern)** — AR: books stop at evap=0.000 exactly; AT: pattern never forms | coexistence (all four radiance-on locked chords sign-correct per seed) |
| F6 luminosity | radiance; magnitude set by intake (metabolism), not emitter count | two objects = 2.00× pooled (per-object conserved); per-emitter sub-linear within one object (UDD) |
| F8 balance line | the medium's operating point — box-conditional, shared by cohabitants (three chords ±0.001), NOT a constant | coexistence; degradation moves it in the KNOWN directions (drained-D blueward, dying-DF redward) |
| F9 bath glow | radiance by construction | removed with it (2.35 band gone, xbar 0.32→0.20) |
| F10 winding wall | nv6 wound-ring-specific | does NOT generalize to embedded nv48 (fold regime; lock null 1.67×) |

**Interferences found: one record error, zero apparatus
interferences, one real physical interaction.** (1) The only
effect that "failed to replicate" for apparatus reasons traced to a
WRONG BASELINE in the standing record (ASTRO 6.5.4's 0.0230 —
corrected, dated). (2) No experiment addition disturbs another's
effect: the lock+tuner, radiance, and both meters carry cleanly
separable physics — every ablation attributed its differences to
the removed knob's own carried phenomena. (3) The one genuine
cross-EXPERIMENT interaction is physical and new: **competitive
starvation at contact** — probabilistic (1 of 2 contact pairs, 0 of
2 isolated pairs), mediated by feeding-connection loss (1 vs 18–33
intake events), read in real time by the spectrometer (loser's line
walks blueward), with no door-to-door event coupling (ILAG stays at
the 1.4× common-mode everywhere).

New effects only combination could show: two-body force-free
coexistence at the π-flat floor; luminosity additivity (2.00×);
the shared box locus; competitive starvation; the
structure/contents split (26 locked bonds held through a full
starvation; E1's edge-census inversion). New instrument findings:
kernel diag `sep=` wrap-broken for straddling objects (fcs true-sep
is the standing method); door-net drain is a warning light, intake
disconnection is the kill.
