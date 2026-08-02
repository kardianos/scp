# The unification battery — ROADMAP §6.5, operational

One **law table** (`laws_*.cfg`), byte-shared by every experiment; one
**apparatus file** per experiment (`apparatus/*.cfg`) that may not touch a
law key (the runner refuses to start otherwise); **physics acceptance**,
not byte comparison. Law variants compete as whole tables — switching a
law per experiment is structurally impossible here.

## THE RATCHET RULE (standing, 2026-07-28)

1. **Every modification of the kernel or of a law table runs the FULL
   battery before it is committed.** No exceptions for "obviously inert"
   changes — inert is what the green run proves. (S2 exists because a
   zero-ablation was skipped once.)
2. **Experiments that pass and defend a claim are added to the battery**
   and become part of the gate for every future modification. The suite
   only grows.
3. **Bars encode defended claims, never snapshots of current output.**
   A bar may be sharpened by a measurement; it may not be softened to
   make a change pass.
4. **A green test leaves the suite only by explicit user decision** —
   never as a side effect of making something else work.
5. **A bar measured on one foam is a claim about that foam.** Bars whose
   claim is about the *law* are gated on a quantile over the frozen
   `SEED_PANEL` (5 foams). The panel changes only when the substrate
   does, and then for every experiment at once — enlarging it because a
   bar fails on the current five is softening under another name.

## The seed panel (ratchet action, 2026-08-02)

Every bar in this suite had only ever been measured on foam seed
`20260727`. A multi-seed scan of the two tightest bars
(`GLM_REVIEW_2026-08-01.md` Step 3) found that the standing seed is a
lucky draw for tilted transport, so two protocol changes were made under
rule 3 (sharpen) and rule 4 (explicit decision):

* **`e7_tune` is now a quantile bar** — `frac ≥ 0.75` on **≥3 of the 5**
  panel foams, replacing the single-seed floor. Strictly harder: it
  cannot be met by tuning to a seed.
* **`e3b_blob_tilt` moved to recorded, NOT gated.** It passes on 1 of 5
  foams — the standing one — and fails a different clause on each of the
  other four, drifting *backward* on `314159`. A bar that holds only on
  the foam it was measured on does not defend the claim it is written
  for. It still runs and still reports on every invocation.

This is read as frozen-foam disorder, not a law failure: `q_detune=0`
makes `e3b` drift backward on *every* seed, which is the law-failure
signature, and is not what the panel shows. S1 and livefab both predict
the variance shrinks — that prediction is `e3b`'s route back into the
gate.

Panel seed 0 keeps the bare log name (`e7_tune.log`), so the standing log
and cfg stay byte-identical to every earlier run and remain the baseline
the ratchet diffs against; the other four write `e7_tune.seed<N>.log`.
ħ-linearity is scored over the primary logs only, so the invariant stays
comparable across variants whether or not they have a panel.

```
python3 battery.py --laws laws_V2z.cfg [--jobs 10] [--only e6_pairs ...] [--skip-run]
```

Outputs `runs/<variant>/<experiment>.log` + `summary.tsv`, prints the
pass/fail table, exits nonzero on any fail.

> **`summary.tsv` is the ratchet baseline and is version-controlled.**
> The repo-root `.gitignore` excludes `*.tsv`; `runs/.gitignore` negates
> it for exactly this file. Rule 1 compares against a baseline, and an
> untracked baseline cannot be checked. Note also that `--only` rewrites
> `summary.tsv` with just the experiments it ran — on 2026-08-02 the V2g
> baseline was found truncated to 3 rows by an earlier targeted run, and
> with no history there was nothing to detect it against. **After any
> `--only` run inside a standing variant's dir, restore the baseline with
> a full run**, or use `--tag` to run in a side dir.

## The experiments (20 runs + 1 cross-check)

e1 conservation (all mechanisms live) · e2 field packet ≥0.3C · e3a heavy
blob seals · e3b light tilted blob translates · e4 curvature linearity ·
e5 CHSH · e6 rung tongue in the foam · e7 tuning curve (P4) · e8 the comma
is paid, monotone · e9 the fifth (3:2) lives · d1 double-slit fringes ·
t1 Tonomura clicks · q2 eraser ± · t4 HOM ordering · qt_lo sub-threshold
(nothing condenses) · qt_hi above-threshold (condensation) · p1 momentum
of light (ballistic centroid + forward momentum current) · g1 the
gravitational footprint (a mass maintains its extended space depression)
· g3 occultation (matter is emergently opaque; beam exits away) ·
g4 space throughput (no steady monopole: far-shell space flux is
mass-rate bookkeeping, decays, stays subdominant to the radiative
channel — the 1/r far field awaits a stable particle's internal cycle)
· LIN:
every fired grain across every log sits on the eps(w) = A0·w/2π grid
(ħ-linearity as a battery invariant).

**Recorded, NOT gated** (ratchet rule: tests enter as they pass; these
encode real claims that currently fail or sit below measurement floors —
see ROADMAP §7): p2 radiation pressure (absorbed light's momentum does
not survive conversion, deficit ~100× — S2-full criterion); g2 dense
free-fall in a frozen space gradient (below the ~1e-3 chaos floor);
halo lensing (right sign, two orders below the foam bias floor); and
since 2026-08-02 **e3b** (tilted-blob transport — seed-favourable, 1 of
5 panel foams; see "The seed panel" above). So the run line reads
**19/20 gated + 1 recorded**, not 21/21.

## The law-dependency map — what each law buys

Each probe is `laws_V2g.cfg` with **one key changed** and everything else
byte-identical (`laws_P_*.cfg`; the merged cfg in each `runs/P_*/cfg/` is
the record). A law is load-bearing if breaking it breaks exactly the
phenomena whose physics depends on it — and nothing else.

First four probes from `GLM_REVIEW_2026-08-01.md` Step 0; the remaining
four were launched in that session into `/tmp` and never landed, and were
re-run on 2026-08-02.

| probe | gated | what falls | the law it confirms |
|---|---|---|---|
| `quant_A0=0` | 14/21 | qt_lo, e7, t4, e3a, g1, g3 | **atoms at boundaries** — creates the photoelectric threshold, pays the comma, separates Bose/Fermi |
| `q_detune=0` | 17/21 | e3b (drifts *backward*, cos −0.99), e7, e9, g3 | **load flattens pitch** (detune side) |
| `s_k` ∈ {0.02…2.0} | 21/21/20/19 | only e4 (and e3b at the low edge) | **space transport** — robust over ~30× |
| `field_J=0` | ~16/21 | e2, d1, p1 all v/C=**0**, t4 | **field coupling** = c for light |
| `cap=1000` | **14/21** | e3b, **e4**, e7, **e8**, e9, g3, **g4** | **load flattens pitch (capacity side)** — and it is *not* symmetric with `q_detune`: it fails a strict superset. Curvature linearity collapses (r² 0.9933→**0.519**), the comma goes unpaid (shed \|δ₀\| 0.092→**2.38**, total 7.2→130), field flux →**0**, the fifth never locks (8→**0**), and the blob freezes solid (speed 0.00134→6e-05). `cap` is where the pitch landscape itself lives. |
| `comb_limit=1` | **16/20** | e3a, e9, g1, g4 | **the consonance comb** — collapsing to the unison kills the dense seal (e3a speed 0.00134→0.00211, g1 with it), the fifth (locked 8→**4**, t_last 60→50) and space throughput. Notably it does **not** kill e7 (3/5 seeds) or e8 — the comb's job is the dense sector's seal and the fifth, *not* the tuning curve. |
| `gamma_res_m=0` | **15/20** | e6, **e7 (0/5 seeds)**, e8, e9, g4 | **the dense rim seal.** Removing it destroys the *pair* sector wholesale: tuned pairs collapse from 31 to 7 with gg 0.684→**0.126**, e7 fails on **all five** panel foams, and the comma is not paid at all (shed \|δ₀\| exactly **0.000**, total 7.2→**0.0**) — nothing is left to shed it. The seal is what lets two dense objects hold a relationship; without it there is no pair physics to tune. |
| `w2=1.65` (=w1) | **14/20** | e4, e6, **e7 (0/5)**, e8, e9, qt_hi | **the dense/field pitch separation** — the most destructive probe of the eight. Collapsing the two pitches onto each other leaves **no** tuned pairs (n=0), no curvature fit at all, the fifth locks for one step then dies (t_last 60→**10**), and above-threshold conversion stops entirely (qt_hi qatoms 404→**0**). ħ-linearity survives on only **17** grains against 498. If dense and field ring at the same pitch there is no conversion boundary, and with no boundary there are no atoms. |

**Conservation sat at the FP floor (~1e-15) in every probe, including
every failing run.** That is the recurring signature: the One Law is
unconditional; the gates are the conditional physics, and they break in
the right places when broken.

**The map is now complete — all eight knobs probed.** Ranked by damage:
`w2` (14/20) and `cap`/`quant_A0` (14/21) are the most load-bearing, then
`gamma_res_m` (15/20), `comb_limit` (16/20), `q_detune` (17/21),
`field_J` (~16/21, but uniquely surgical — it removes light and touches
nothing dense), and `s_k` (robust over ~30×). No knob is inert, and none
breaks everything: each takes down the phenomena that depend on it and
leaves the rest standing. That pattern is the law/tuning discriminator.

*(`comb_limit=1`, `gamma_res_m=0` and `w2=1.65` were scored under the
seed-panel protocol, hence /20 gated + 1 recorded rather than /21.)*

## Variant verdict (2026-07-28, after §6.10 S2 retirement)

| variant | law | result |
|---|---|---|
| V1 | source atom, per-cycle floor | 15/17 — the floor cannot pay the comma (e8) or fit e4 |
| V2 | source atom, credit, kappa_freq=0.6 | 17/17 — the §6.9 table; the bias buys margin, not passes |
| V2z | source atom, credit, NO bias (kf=0, kr=0) | 17/17 then 18/18 with p1 — the §6.10 retirement point |
| **V2g** | **V2z + space transport (s_k=0.06, s_disp=0.3)** | **20/20 — the standing law table (footprint + occultation gated; robustness s_k ∈ [0.02, ~0.15])** |
| V2d | source atom, credit, kappa_reac=1 (derived reactive term) | 16/17 — e3b decoheres: rate transport cannot host the raw current; kept as the S2 exhibit and S2-full acceptance criterion |
| V3 | destination atom, floor | 16/17 — e8 (comma needs memory) |
| V4 | destination atom, credit | 16/17 — e3b (dense-atom tax perturbs the front) |

Annealing history (each value applied to ALL variants — mechanisms are
variant-orthogonal; full record in ROADMAP §6.9/§6.10):
L1 = the §6.4 proposal (Γ_m 0.02); L2 = Γ_m → 0.10 (unified landscape
has ~5× per-link detunes); L3 = kappa_freq dispersive bias (choir's
correction, A1); L4 = flight-loads-pitch (zero-constant kernel law —
the shelter was detuning its own pair); L5 = bias windowed by the
acceptance (blob rims feel nothing); L6 = bias × mutual gate closure
+ mob_floor 0.004 (pair bleed vs front recruitment both pass);
L7 = S2 derivation (CONSONANCE Part VII): the bias is real amplitude
physics, derived in v89/s2/ — and unnecessary at rate level after L4;
kappa_freq retired (the zero ablation had never been re-run post-L4;
the battery caught the vestige).
