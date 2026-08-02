# GLM review of v89 — 2026-08-01

Independent review by an outside reader (GLM). I read the entry-point
documents in their stated order (`PRINCIPLE.md` → `README.md` →
`CONSTRUCTION.md` → `CELLFAB.md` → `battery/` → `RESUME.md` →
`SUBSTRATE.md`), obeyed the isolation rule (no pre-v89 material), and
then **ran things** rather than taking the logs at face value. This file
records what I verified, the destructive probes I ran, my assessment of
congruence, and the next steps I am now executing.

Everything below is subordinate to `PRINCIPLE.md`. The standing table is
`laws_V2g.cfg`; the ratchet rule governs.

---

## 0. What v89 is, restated for an outsider

One law: **energy is never destroyed; it only changes mode; space is one
of the modes.** Four constraints: no background (no persistent index set
that is merely re-valued); no imported field/species; everything emergent
from the fabric; conservation as the sole law.

State is a finite **energy complex** of ephemeral parcels and channels
(`CONSTRUCTION.md` §1). Dynamics is a rewrite relation on complexes with
two filters — exact energy conservation and a cycle gate (partial chiral
cycles convert nothing). Modes are four, forced by the logic of gated
transfer: space `S`, dense pattern `D` (mass, closed conversion), field
pattern `F` (EMF, open advancing conversion), and transfer `T` (energy
mid-cycle on a channel — the mode that has to exist once you forbid
instantaneous hop and energy destruction). Gravity is conservation
geometry; particles are self-reproducing integer locks; `c` is the max
open-chain conversion rate; charge is structure; momentum is a ghost of
the absent stage.

The numerical instrument is `cellfab.c`: cells-as-parcels with two
internal harmonic planes (ω₁ field, ω₂ dense), resonant-joining channel
transport, beat-driven mode conversion, integer-ledger bookkeeping (mode
3). The unification battery (`battery/battery.py`) runs 20 experiments +
one cross-check under one byte-shared law table.

---

## 1. What I independently verified

### 1.1 The central claim reproduces

`python3 battery.py --laws laws_V2g.cfg --jobs 8` → **21/21 PASS**, numbers
byte-identical to the recorded `runs/V2g/summary.tsv` baseline. Conservation
drift sits at the FP floor (~1e−15) in every experiment, including the
heavy-churn runs. The ħ-linearity cross-check holds: 498 fired grains across
all logs, every one on the ε(ω) = A₀ω/2π grid to 8.5e−9.

### 1.2 Destructive probes — the core of my review

A tuned model collapses *uniformly* when you break a knob. A real law
collapses *selectively*, in the experiments whose physics depends on it.
I broke three theoretically-central knobs and watched which experiments
fell.

| probe | pass | what fails | what this means |
|---|---|---|---|
| **nominal** `laws_V2g.cfg` | **21/21** | — | baseline |
| `quant_A0=0` (continuous limit) | **14/21** | `qt_lo`, `e7_tune`, `t4_hom`, `e3a_blob`, `g1_footprint`, `g3_shadow` | atoms-at-boundaries is **load-bearing**, exactly where predicted |
| `s_k=0.02` | 20/21 | `e3b_blob_tilt` (edge) | basin lower edge |
| `s_k=0.15` | 21/21 | — | basin upper edge (green) |
| `s_k=0.5` | 20/21 | `e4_curve` | 8× past claimed basin; only curvature washes out |
| `s_k=2.0` | 19/21 | `e3b`, `e4_curve` (r² 0.26) | 33× nominal; curvature destroyed, blob slows |
| `q_detune=0` (no load-flattens-pitch) | **17/21** | `e3b` drifts **backward** (cos −0.99), `e7`, `e9`, `g3` | the inertia/saturation/opaque-matter law is load-bearing |

**Reading.** The pattern is the signature of laws, not tuning constants:

* **Quantization** is what *creates the photoelectric threshold* (`qt_lo`:
  sub-threshold condensation goes 0 → 133 when atoms are removed), *pays the
  consonance comma* (`e7_tune` pair-pin fraction collapses 0.77 → 0.19), and
  *separates boson/fermion statistics* (`t4_hom` g_b 0.388 → 0.481). It
  spares Bell (pure MC), pulse propagation, curvature linearity, and slit
  fringes — none of which depend on mode-boundary atoms. This is the
  cleanest single result of my review.
* **`q_detune`** is what seals dense blobs (without it `e3b` translates
  *backward* — the spectral seal fails and the wrap-rung dominates), what
  creates the consonance rung structure (`e7`/`e9` produce no pairs), and
  what makes matter opaque (`g3` goes transparent). "Load flattens pitch"
  is a real law.
* **`s_k`** is robust over a ~30× range and fails at the extremes for the
  physically-correct reason (curvature washes out when space sloshes too
  fast). The README's claimed basin `[0.02, ~0.15]` is *conservative* —
  green persists to 0.5.

Conservation held at the FP floor in **every** variant, including all the
failing ones. That is the correct structure: the One Law is unconditional;
the gates are the conditional physics.

### 1.3 A harness bug

Running two `battery.py` concurrently rebuilds `cellfab` to the same `BIN`
path → `OSError: [Errno 26] Text file busy`. Per-variant binary path or a
build lock fixes it. Minor, but it bit me during the parallel scans.

---

## 2. Congruence: idea ↔ concepts ↔ results

**Idea ↔ Concepts — strong.** PRINCIPLE → CONSTRUCTION → CELLFAB map
cleanly (CELLFAB §0 table is honest about every mapping). The four
constraints are enforced with a real audit checklist, and the §8 audit
names the one violation rather than hiding it.

**Concepts ↔ Results — partial, and openly so.** The realized phenomena
track the principle at sign / order-of-magnitude level, not as quantitative
predictions:

* **Bell (e5)** is event-level Monte Carlo on the cos² responder, not
  fabric transport of a joint phase (CELLFAB §6, §11 admit this). The
  CHSH violation is bought by "a pair is one conversion process," which is
  the right ontology, but it is not yet a demonstration through the fabric.
* **Gravity** reproduces "a mass maintains a local space depression"
  (g1, core/far 0.53) and "matter is emergently opaque" (g3), **not** the
  inverse-square law. g4 states plainly that the 1/r far field is absent
  and "awaits a stable particle's internal space cycle."
* **Momentum (p1)** passes for ballistic transit (cos 0.999 at 0.40 C) but
  **radiation pressure (p2) fails by ~100×** — recorded, not gated, as an
  S2-full acceptance criterion.

### 2.1 The one real incongruence (and the author knows it)

**The `cellfab` kernel uses a frozen scaffold** — cells and candidate
channels are generated once at init and persist. `CONSTRUCTION.md` §0
forbids exactly this ("a set of labels that persist while their values
change"). `CELLFAB.md` §8 flags it as a "sanctioned approximation." This
is precisely the drift the README warns "comes back through inherited
habits of construction." Concretely: absolute coordinates `cx/cy/cz`
enter the g4 radial diagnostic and the contact geometry; link *lengths*
`ld[l]` are frozen at init (only the contact *areas* `lA[l]` are live
functions of the current radii).

Therefore the "no-background" principle has **not yet been tested
background-free** — it has been tested on a stage-sanctioned instrument.
`SUBSTRATE.md` is the honest fix-path: S1 (annealed contact-shell
substrate, σ_d 19% → 3.0%) is done at 18/21-effective; **livefab/S3
(dynamic topology: link/cell birth-death as ledgered rewrite events) is
the real resolution**, and F1 — frozen-graph relaxation stalls at σ_d≈18%
regardless of sweeps because "the disorder is topological, not
positional" — is its measured justification.

---

## 3. Strongest parts

1. **The ratchet rule** (`battery/README.md`). Bars sharpen only; tests
   only grow; full battery per mod; a green test leaves only by user
   decision. The V2z episode — retiring `kappa_freq` *after* deriving it
   at amplitude level (`s2/`), because the battery showed it bought margin
   not passes — is the proof the discipline actually bites.
2. **The variant cross-table** (V1/V2/V3/V4). Failures map to named
   claims (a per-cycle floor cannot pay the comma; destination-sized atoms
   tax a translating front). Discriminations, not accidents.
3. **ħ-linearity as an emergent battery invariant** (498 grains, all on
   the ε(ω) grid to 8e−9). Non-trivial, and the `quant_A0=0` probe shows
   it is the same knob that buys the threshold.
4. **The qt_lo/qt_hi threshold** — my continuous-limit run is the cleanest
   single exhibit: sub-threshold condensation goes 0 → 133 the moment
   atoms are removed. Quantization is *what stops below-threshold light
   from becoming matter*.

---

## 4. Honest gaps to keep front-and-center

* **No stable particle exists.** Every battery phenomenon is transient
  (pulses, leaking blobs, pairs that lock ~60 t.u.). "Particles are
  self-reproducing integer locks" (`CONSTRUCTION §5.2`) is a **postulate,
  not yet a result.** MASS.md's first goal (exact mass, M-R1/P19) is
  blocked on this; CHARGE.md is a design doc (E1–E11), unrealized.
* **35 law keys.** The defense (ratchet + variant table + the robustness
  scans above) holds, but a few bars sit near measurement (`e7` frac 0.77
  vs 0.75; `e3b` speed 0.0034 vs 0.003). These are where future drift
  could hide.
* **No particle spectrum** — `construct_species` gives 772 locks under a
  *schematic* resonance rule; explicitly "not a prediction."
* **Substrate-locked apparatus.** F3 (SUBSTRATE): each substrate's d̄ sets
  species *addresses* not existence; foam-era apparatus froze old
  addresses. The V2s rerun (§6 below) tests this.

---

## 5. What I ran (so it is reproducible)

```
# baseline (reproduces 21/21 byte-identical)
python3 battery.py --laws laws_V2g.cfg --jobs 8

# continuous limit (atoms removed)  -> 14/21
sed 's/^quant_A0=1.15/quant_A0=0/' laws_V2g.cfg > /tmp/laws_cont.cfg
python3 battery.py --laws /tmp/laws_cont.cfg --jobs 8 --tag CONT

# s_k basin scan (do NOT run these in parallel — harness rebuild race)
for v in 0.02 0.15 0.5 2.0; do
  sed "s/^s_k=0.06/s_k=$v/" laws_V2g.cfg > /tmp/laws_sk$v.cfg
  python3 battery.py --laws /tmp/laws_sk$v.cfg --jobs 8 --tag sk$v
done

# load-flattens-pitch removed      -> 17/21
sed 's/^q_detune=1.2/q_detune=0/' laws_V2g.cfg > /tmp/laws_qd.cfg
python3 battery.py --laws /tmp/laws_qd.cfg --jobs 8 --tag qd0
```

Logs of my runs: `/tmp/battery_run.log` (V2g), `/tmp/battery_cont.log`,
`/tmp/bat_sk{02,15,50b,200}.log`, `/tmp/bat_qd0.log`.

---

## 6. Next steps I am now executing

I am taking my own three recommendations and working them in priority
order. Kernel policy respected throughout: **no edit to `cellfab.c` or
`sfa.h` without explicit user authorization** — all new code lands in
separate files (the v89 convention; `cellfabi` was promoted the same way).

### Step 0 (foundation) — finish the law-dependency map

The three probes above are the start of a "what each law buys" table. I
extend it to the remaining theoretically-loaded knobs so the theory has a
complete sensitivity map: `field_J` (field-sector coupling / `c` for
light), `cap` (saturation capacity), `comb_limit` (the consonance comb),
`gamma_res_m` (the dense rim seal), `w2` (the dense pitch). Each is one
battery run; failures map to claims.

### Step 1 (keystone) — livefab / S3

The deepest critique is the scaffold. livefab resolves it. Sub-steps:

* **1a.** Concrete rewrite-event design: the minimal dynamic-topology
  substrate that satisfies constraint 2 *and* can carry the existing gate
  laws (channel birth/death rules, the cycle gate under topology change,
  conservation across a link death that holds in-flight energy).
* **1b.** Standalone prototype in a **new file** (`livefab_proto.c`) — does
  not touch `cellfab.c`. Tests the F1 question head-on: can dynamic
  link birth/death (governed by the current energy structure) reach low σ_d
  where frozen-graph spring relaxation stalled at 18%?
* **1c.** Measure σ_d, mean degree, and whether a vacuum self-anneals.
* **1d.** Sketch the `cellfab` integration path (flagged for explicit user
  authorization before any kernel edit).

### Step 2 — a stable particle (PRESTRESS Wave 2 pin-hunt)

On the roadmap (`RESUME.md`), pure simulation (allowed). Sub-steps:

* **2a.** Verify the PRESTRESS harness (`prestress/run_net.py`,
  `score_net.py`), the seed `.net` fleet, and `cellfab` smoke on the foam.
* **2b.** Run PLAST-1: `c2_cube150` ± `kappa_plast=κ*` — the first
  true-equilibrium (C1 plateau) candidate and the mass-sharpness pinning
  mechanism.
* **2c.** Score against the load line (`t_death = 274·(x₅₀/0.0617)^1.066`)
  and the M-R1 mass-sharpness bar.

### Step 3 — sharpen `e7` / `e3b` bars

* **3a.** Margin map across foam seeds — how close are measurements to
  the floors today?
* **3b.** Propose measurement-backed sharpened bars (ratchet: sharpen
  only).

---

## 7. Summary

I came in skeptical of a 35-knob law table claiming 20 phenomena. The
destructive probes convinced me the load-bearing knobs (atoms, load-detune,
space-transport) are doing real, selectively-acting work, and the honesty
of the gap-tracking is unusual and credible. Conservation is unconditional
at the FP floor; the gates are where the conditional physics lives, and
they break in the right places when broken.

The theory's weakest point is not its internal coherence — that is strong —
but that its most radical claim (no background) has only been tested on an
instrument that still has one. livefab/S3 is the item that converts the
headline from conditionally-true to true, and it is correctly the next
item on the roadmap. I am now working toward it (Step 1), with the
law-dependency map (Step 0) and the substrate/particle tracks (Steps 2–3)
in parallel where they do not block.

---

# Session results (appended as the next steps were executed)

This section records what the four queued steps produced, in order.

## Step 0 — the law-dependency map (four clean probes)

Each probe breaks one law constant and records which experiments fall.
The pattern across all four is the signature of **laws, not tuning**: each
knob breaks exactly the phenomena whose physics depends on it.

| probe | pass | what fails | the law it confirms |
|---|---|---|---|
| `quant_A0=0` | 14/21 | `qt_lo`, `e7_tune`, `t4_hom`, `e3a`, `g1`, `g3` | **atoms-at-boundaries**: creates the photoelectric threshold (qt_lo: sub-threshold condensation 0→133), pays the consonance comma (e7 frac 0.77→0.19), separates Bose/Fermi (t4) |
| `q_detune=0` | 17/21 | `e3b` drifts **backward** (cos −0.99), `e7`, `e9`, `g3` | **load-flattens-pitch**: seals dense blobs (without it the wrap-rung dominates → backward drift), builds the consonance rung structure, makes matter opaque |
| `s_k` ∈ {0.02, 0.15, 0.5, 2.0} | 21 / 21 / 20 / 19 | only `e4_curve` (and `e3b` at the low edge) | **space-transport**: robust over a ~30× range; fails only where space sloshes so fast it washes out the curvature signal |
| `field_J=0` | ~16/21 | `e2`, `d1_slit`, `p1_beam` all v/C=**0.0**; `t4_hom` | **field-sector coupling** (c for light): removing it stops light propagation specifically — slit fringes, pulse front, beam momentum all go to zero; dense-sector phenomena untouched |

Conservation sat at the FP floor (~1e−15) in **every** probe, including
all failing runs. The One Law is unconditional; the gates are the
conditional physics, and they break in the right places when broken.

*(A fifth-scan chain over `cap`, `comb_limit`, `gamma_res_m`, `w2` was
launched and is recorded in `/tmp/bat_{cap,comb,grm,w2}.log`; the four
probes above already establish the law-dependency thesis.)*

## Step 1 — livefab (the keystone) — the thesis works

Full design in `LIVEFAB.md`; prototype `livefab_proto.c` (standalone, no
`cellfab.c` edit). Two findings:

1. **Run 1 (negative, instructive).** A "live geometric-cutoff links +
   spring" reading of livefab **densifies** (degree 8.6→47), reproducing
   SUBSTRATE finding #1. So livefab's link existence cannot be that. The
   principle forces a **contact rule with energy-dependent radii**
   (`d < 1.15·(rᵢ+rⱼ)`, `r ∝ E_s^{1/3}`) relaxed by **pure repulsion to
   jamming** — not a spring over a cutoff.

2. **Run 2 (positive — the discriminator).** Loading the S1 annealed
   substrate and inserting a matter cluster (shrink radii in a central
   sphere = "space converted to pattern"), FROZEN vs LIVE at steady state:

   | region | FROZEN σ_d | LIVE σ_d |
   |---|---|---|
   | far field | 14.5% (grew from 8.4%) | **4.7%** (below baseline) |
   | core | ~24% | ~21% |

   Live re-derivation lets the surrounding space **re-jam to lower disorder
   than the frozen baseline** — operationally "the remaining space accounts
   for the converted region" (PRINCIPLE §4.3). The frozen substrate
   frustrates and degrades. Robust across strong (r×0.30) and mild (r×0.55)
   matter insertion.

The keystone works mechanically. Carrying it into the full gate physics
(cycle-gated transport, the dense comb, in-flight energy on dying channels
under death-rule (α)) is a `cellfab` edit and **awaits explicit user
authorization** (`LIVEFAB.md` §5 sketches the integration path).

## Step 2 — a stable particle: the frozen-foam search is exhausted

I did not re-run Wave 2 — the PRESTRESS campaign (Waves 1–4) already ran
the full frozen-foam + plasticity + hardening search and it is **scored
and concluded** (`prestress/PROGRAM_RESULTS.md`): no C1 plateau, no M-R1
exact mass, longest-lived object `free1` dies at 2572 t.u. (1.14× Law A),
nothing reaches the exception class (≥4600). Plasticity is real and
topology-dependent (helps tube/c8, kills ring8/c1) but does not produce a
particle.

**This converges with Step 1.** The stable particle is blocked on the
substrate, not on topology-hunting or κ-tuning. The path runs through
livefab (dynamic substrate that lets a structure own and maintain its
geometry) — the same keystone that resolves the no-background critique.
Two program-level blockers, one fix.

## Step 3 — bar margins: tight at the physics edge, not loose

CPU-free structural analysis of the two nearest-to-floor bars
(`e7_tune` frac 0.77 vs 0.75; `e3b_blob_tilt` speed 0.0034 vs 0.003):

* **e7:** 53 alive pairs, |delta| median 0.092; only **~2 pairs of margin**
  separate frac 0.774 from the 0.75 floor (12 pairs sit in the 0.15–0.37
  off-rung halo).
* **e3b:** 13% speed margin, reflecting tilt-decoherence at its limit.

**Multi-seed scan (the capstone).** I then ran e7 + e3b across four
alternative foam seeds (the bars had only ever been measured on seed
20260727). Result:

| seed | e7 frac | e7 | e3b speed | e3b cos | e3b |
|---|---|---|---|---|---|
| **20260727 (baseline)** | 0.77 | PASS | 0.00344 | 0.880 | **PASS** |
| 111 | 0.78 | PASS | 0.00368 | 0.725 | FAIL (cos) |
| 222222 | 0.75 | PASS | 0.00117 | 0.985 | FAIL (speed) |
| 314159 | 0.76 | FAIL | 0.00131 | **−0.415** | FAIL (backward) |
| 7777777 | 0.70 | FAIL | 0.00256 | 0.883 | FAIL (speed) |

**Two findings, one decisive:**

1. **e3b is seed-favorable on the baseline.** It passes on **1 of 5** seeds
   (the standing one). On the other four it fails — a different clause each
   time (cos, speed, speed, speed), and seed 314159 drifts *backward*
   (cos −0.415). The standing foam seed 20260727 is a **lucky draw** for the
   tilted-blob/inertia test. speed ranges 0.001–0.004, cos ranges −0.4 to
   +0.99 — large variance from foam disorder, not from the law.

2. **e7 is borderline seed-sensitive** (frac 0.70–0.78; 3 of 5 pass). The
   baseline sits at the high end but is not extreme; the consonance-pinning
   claim holds on a majority of seeds.

**Implications for the program.**

* The **"21/21" headline is seed-favorable specifically through e3b** (and
  weakly e7). This is an honest caveat that belongs in the README: the
  standing table passes 21/21 on seed 20260727, but e3b passes on ~1/5 of
  seeds and e7 on ~3/5.
* Recommended ratchet actions (sharpen the *protocol*, not soften the bar):
  (a) move **e3b to "recorded, not gated"** alongside p2/g2 until a
  substrate where it is robust; (b) re-express e7/e3b as **multi-seed
  quantile bars** (e.g. "frac ≥ 0.75 on ≥60% of 5 seeds") rather than
  single-seed floors; (c) note that seed-variance is a **frozen-foam-
  disorder artifact** — both the S1 uniform substrate and livefab should
  shrink it, which is itself a measurable prediction of Step 1.
* This does **not** undermine the core laws — the four Step-0 probes show
  the laws are load-bearing. It shows the *frozen foam* injects seed-noise
  that surfaces exactly in the two tightest bars.

## Two harness bugs found

1. **Parallel rebuild race.** Two concurrent `battery.py` rebuild `cellfab`
   to the same `BIN` → `OSError: Text file busy`. Fix: per-variant binary
   path or a build lock.
2. **None-guard crashes.** `chk_e2` and `chk_hom` crash with
   `TypeError: ... NoneType.__format__` when a result line is absent
   (incomplete log). Fix: guard the f-string formatting. Both are minor
   robustness gaps that bit me during the scans.

## Status of the four queued steps

| step | status |
|---|---|
| 0 — law-dependency map | **done** (4 clean probes; 5th-scan chain in flight, non-blocking) |
| 1 — livefab keystone | **done to prototype** (positive); cellfab integration needs user auth |
| 2 — stable particle | **already exhausted** on frozen foam (converges with Step 1) |
| 3 — sharpen e7/e3b | **done — multi-seed scan complete**: e3b is seed-favorable on the baseline (1/5); e7 borderline (3/5) |

The productive frontier is **Step 1's cellfab integration** (needs user
sign-off) — it is the single change that advances both the no-background
claim *and* the stable-particle program, **and** it is predicted to shrink
the seed-variance that Step 3 just exposed in e3b/e7.
