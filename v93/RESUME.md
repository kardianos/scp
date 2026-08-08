# v93 RESUME — read this FIRST after compaction (self-contained handoff)

You are resuming the **SCP v93** campaign mid-stream. This file is the single
briefing. Read it fully, then the files it points at. Working dir: `/home/d/code/scp`.

## 0. What v93 is (one paragraph)

SCP is a no-background physics-simulation programme ("energy is never
destroyed; it only changes mode; space is one of the modes"). **v93 = THE
UNITARY DENSE CHANNEL:** give the dense (matter) sector the field sector's
transport algebra — within-mode dense transport becomes a product of UNITARY
PAIRWISE PLANE ROTATIONS (a cousin of the field's pass F) on the dense
amplitude ψ_m = √Em·e^{i·th2}, replacing the additive magnitude "want". Each
Givens hop conserves the two-cell norm exactly ⇒ ΣEm conserved to roundoff
(conservation is a theorem of the update, not a patched ledger). The cross-
term 2·Im(ψ_i*·ψ_j) is the link current J = the dense momentum. Knob `amp_nat`
(default 0 = byte-inert). The conversion DOOR stays discrete/quantized.

## 1. Reading order (after this file)
1. `CLAUDE.md` — project conventions (the banner now points at v93 as active).
2. `v93/README.md` — the design (PART I–V). NOTE its status banner was updated:
   v93 is OPENED/implemented (the "docs-only" line is stale history).
3. `v93/ITEMS_OUTCOMES.md` — the items 1/2/4 synthesis, then
   `v93/RING_DNLS.md` — **the most current results** (the retention face:
   closed-ring route A + DNLS route B, both executed to completion).
4. `v93/L1_FINDINGS.md` — the chronological findings (has a POST-REVIEW
   CORRECTION banner at top that **supersedes** the early "CONFIRMED" claims —
   the early claims were inflated by two bugs; trust the banner + ITEMS_OUTCOMES).
5. `v93/consult/REVIEW_claude.md`, `REVIEW_grok.md`, `REVIEW_opencode.md` — the
   3-reviewer consultation that found the bugs and reframed everything.
6. `v92/consult/SUBQUARK_synthesis.md` — the convergence this design came from.

## 2. Current status (honest, post-review + items 1/2/4)

| bar | status |
|---|---|
| **L1-A translation** | **RESCUED with caveats.** Linearize (amp_logate>0, drop the gate) + measure the dense-hop current `fd` (p1_meter=1), NOT the tagged-centroid meter → coherent, **seed-robust +x** translation (fd +287/+297/+277 across 3 seeds), theory-correct direction (v_g=+2J̄d·sin(kx·d)>0), speed>2.6e-3. Open: a dt-Trotter artifact (item 1) + formal tolerance-C≡Go at amp_nat>0. |
| **L1-B conservation** | **GREEN on live matter** (drift −1.8e-15 at amp_nat=2 e3b ≪1e-13). The matterless `exp=bath` runs (Em≡0) were VACUOUS — never use them to grade the channel. The long T=300 QUENCH arm is 2.15e-13 (>1e-13; needs per-step-normalized reckoning). |
| **L1-C anti-ignition** | The byte-identical-bath result was vacuous (matterless bath). The channel IS active on matter (driven-beam births 8800→4200). The §II.9 **armed ρ_coh≈0.77 coherent bath** is the real bar — UNRUN. "Self-gating τ≈0" is quantitatively wrong (gate mean ~0.1–0.2 at p_gate=8); the real guarantee is unitarity (can't raise ΣEm). |
| **QUENCH-3 + item-4 retention** | **NEGATIVE, precisely diagnosed.** A *localized* winding cannot be retained by the linear unitary channel: it diffracts the packet and |ψ|→0 edges are phase-slip sites. Independent of the door. Even after fixing the one-way evap door + empty-cell reset, R2d≈0.02. **Needs topological support (closed ring, |ψ|>0 on cycle) + possibly DNLS nonlinear binding.** |
| **Retention face (RING_DNLS.md)** | **EXECUTED, both routes.** Route A: vacuum C6 ring HOLDS the winding (W=+2 exact ≥200 t.u. — first retained matter winding; existence proof for topological support), but the medium kills it in ~10 t.u. via contact dephasing (both laws; leak = accretion, throttled non-monotonically by the q_detune barrier); the vacuum ceiling is the kernel's own sequential hop-sweep order (Trotter family, same as item 1). Route B: the q_detune detuning IS a real self-trap nonlinearity (qd0 melts / qd≥0.6 condenses / deep corner envelope-frozen) but it binds ENERGY not PHASE (retention dead at every corner) and kills mobility (kx arms stall). Spontaneous incoherent condensation under the unitary channel = a creation-adjacent discovery (no lock/gate/door; additive can't). amp_door refill UNTESTED (dark bath — door never fired). |

**Bottom line:** pass U is a sound, conservative, battery-safe unitary
primitive that genuinely transports coherently (when linearized + measured by
the current). It does NOT bind local structure, so localized-winding retention
fails. The unitary channel is good at TRANSPORT, bad at RETENTION-of-local-
structure.

## 3. The code (v93 additions, all byte-inert at default 0)

C kernel of record: `v93/kernel/freecell.c`. Go A/B mirror: `v93/fab/`.
Search by these comment markers (line numbers shift; markers don't):

- **`amp_nat`** (default 0): unitary dense channel strength. Struct field,
  parser `amp_nat`, header print. 0 = additive want path (byte-inert).
- **pass U** (`/* pass U: v93 UNITARY DENSE CHANNEL */`): the unitary hops.
  Reconstruct (m1,m2)=√Em·(cos th2, sin th2); **in-pass precession**
  (`Local clock precession IN-pass`) rotates ψ_m by w2e·dt; per slot (canonical
  i-side) Givens-rotate by shau_[s]; writeback (Em=|ψ|²; th2=arg ψ only if
  Em>1e-12 — the empty-cell-reset fix). Mirrors pass F's hop form. Go mirror
  in `fab/step.go`.
- **τ_s branch** (`v93 UNITARY DENSE CHANNEL: fold the want`): in pass 2,
  τ_s = amp_nat·base·√(gᵢⱼgⱼᵢ), OR (amp_logate>0) amp_nat·base (linear, no gate);
  clamp 0.5.
- **`amp_logate`** (item 1, default 0): >0 drops the phase-dependent gate from
  τ_s → linear Schrödinger hop. This RECOVERS coherent transport (the gate was
  an artifact source post-empty-cell-fix).
- **`hop_order`** (face A, default 0 = byte-inert): the pass-U hop schedule.
  0=sequential canonical sweep; 1=Strang symmetric (forward τ/2 + reverse τ/2,
  time-reversal-invariant — best retention); 2=randomized (per-step Fisher-Yates
  shuffle, separate `hop_rng_s`/`hopRng` stream — best dt-invariance). Refactor:
  collect active canonical slots into `actslot_`/`actslot`, dispatch by schedule;
  `apply_uhop`/`applyUhop(sl, fac)` is the hop primitive (fac=1 full, 0.5 Strang
  half). Header line `# v93 hop schedule`.
- **`amp_door`** (default 0): the arg(ψ) door. At condensation (pass 6,
  `v93 arg(psi) door (§II.7)`): coherent merge arg(ψ_m_new)=arg(√Em_old·e^{i th2}+√(amp_door·d1)·e^{i aph}),
  |ψ_m| fixed by Em. CONTINUOUS mix weight (not boolean).
- **symmetric reverse door** (`field_inject`, `v93 symmetric reverse door`):
  when amp_door>0, matter→field evap composes √(amp_door·dE)·e^{i th2} into fa,
  |fa|=√(e+dE). Fixes the one-way door. Go: `fab/sim.go:fieldInject`.
- **`seed_mw`** (item 4 + retention face, default 0): seeds a MATTER
  azimuthal vortex th2=m·atan2(dy,dx) — on the blob (item 4) AND on the
  exp=ring cycle (`v93 ring-retention face` marker; overrides the chain
  lock). Mirrored in Go (params/blob/ring; the Go ring also gained the
  previously-unmirrored MOTION-#31 kx tilt).
- **pass 6 clock skip** (`Byte-inert (amp_nat==0 takes the original branch)`):
  when amp_nat>0, the out-of-pass `th2 += w2e·dt` is skipped (precession is
  in-pass U instead).

Knob wiring (each): struct field + default 0 + `else if (!strcmp(k,"..."))`
parser in `cfg_defaults`/`set_kv` (C) or `params.go` (Go); header print in the
config-dump. Always wire BOTH kernels (C is record, Go is the A/B mirror).

## 4. How to build / run / measure (CONCRETE)

Build: `cd v93 && make all` (freecell + fabrun + battery + volview).
- C: `./freecell k=v ...` ; Go: `./fabrun k=v ...`
- **Battery (the ratchet gate):** `./battery` → expects ALL GREEN 87 at
  amp_nat=0/amp_door=0. Run before every commit that touches kernels.
- **Detached runs** (e3b/QUENCH exceed the 2-min tool limit): use
  `setsid bash -c './script > out 2>&1; touch out.done' </dev/null >/dev/null 2>&1 & disown`
  then poll with instant `ls`/`cat`/`kill -0 $PID`. (Plain `&` makes the bash
  tool hold the pipe for 120s — always setsid+redirect.)

**Key experiments + meters:**
- **e3b (L1-A blob translation):** `exp=blob bath=1 L=16 dt=0.02 T=80
  diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=2 amp_logate=1`.
  Result line `# RESULT blob_drift speed=… cos_to_kdir=…`. **WARNING: the
  tagged-centroid `blob_drift` measures the WRONG direction** (tracks draining
  tagged cells). Use `p1_meter=1` and read `# RESULT p1x sp=… fl=… fd=… gm=…`
  — `fd` is the dense-hop current (the correct +x translation signal).
- **QUENCH-3 spin retention:** `exp=slit slit_mask=3 L=64 T=300 sigma=8
  slit_sy=8 amp=2 slit_srcx=32 kx=0 spin_m=2 seed=20260802 diag_every=2000
  snap_every=1000 bath=1 snap_file=runs/quench/NAME.fcs …`. Analyze:
  `./fcsdump -mode cells NAME.fcs | python3 report/analyze_winding.py /dev/stdin SRCx SRCy KX`
  (SRC=(32,32) for QUENCH; (8,8) for the L=16 blob). Outputs per-frame W_A
  (field winding), W_th2/R2d (matter winding), m-spectrum. **build fcsdump
  first:** `go build -o fcsdump ./cmd/fcsdump`.
- **Matter-winding retention (item 4):** add `seed_mw=2 kx=0 k_rad=0 wf_on=0`
  (no field/door); analyze R2d about (8,8).
- **Conservation (L1-B) on live matter:** run the e3b (it HAS Em), read
  `# RESULT drift_rel`. NOT `exp=bath` (vacuous, Em≡0).
- **Linearize (item 1):** add `amp_logate=1`.

## 5. CRITICAL lessons / gotchas (do not re-repeat)

1. **Go `math.Sincos` returns (sin,cos), NOT (cos,sin).** C's `sincos(x,&ss,&cc)`
   gives ss=sin,cc=cos. In Go write `ss, cc := math.Sincos(x)`. A swapped
   assignment is a ~90° rotation bug. (Found once in the precession; the
   C≡Go-at-strength abx exists to catch this — but see #2.)
2. **No C≡Go A/B has been run at amp_nat>0.** All amp_nat>0 results to date
   are C-only. At amp_nat>0 C≠Go is EXPECTED (tolerance, not byte — chaotic
   resonance amplifies 1-ulp FP differences, §IV.12). Build a relative-
   tolerance abx before trusting any amp_nat>0 number as "the" result. The user
   has said C/Go diffs are fine for now (reasonably explained), so don't block
   on it, but flag it.
3. **`exp=bath` is matterless** (`freecell.c` "exp=bath: nothing to seed" ≈3938);
   with wf_on + no beam, Em≡0 forever → the unitary channel never engages.
   NEVER grade L1-B/L1-C on `exp=bath`. Use a live-matter run (e3b blob).
4. **The tagged-centroid `blob_drift` meter is misleading** for the unitary
   channel (shows -x when the true current is +x). Use `fd` (p1_meter=1).
5. **The empty-cell `th2→0` reset was a real bug** (atan2(0,0)=0 phase-flattened
   the medium → faked the L1-A "coherence window"). It's FIXED now (preserve
   old th2 when Em<1e-12). Don't reintroduce.
6. **`amp_nat>0` BYPASSES the additive want+inflight** (passes 3–5 no-op). This
   is a bigger physics change than "swap debit for rotation" — flight load,
   rough conversion, parcel-borne phase all go dark for dense. Intentional.
7. **p_gate=8 makes the random-phase gate mean ~0.1–0.2, NOT ≈0** — so the
   channel is NOT strictly inert on incoherent matter (it scrambles
   diffusively). Anti-ignition = unitarity, not "self-gating τ≈0".
8. **Harvest-script pitfall:** `grep … drift | head -1` catches the column
   HEADER, not the RESULT value; awk `$(NF-3)` on analyze_winding output grabs
   the wrong column. Always re-derive from the fcs via analyze_winding, or use
   precise `awk '/RESULT drift_rel/{print $4}'` / `grep -oE "R2d=[0-9.-]+"`.

## 6. Conventions (non-negotiable)
- **Ratchet:** every kernel change runs the FULL battery before commit; bars
  sharpen by measurement, never soften to pass. Full battery must be ALL GREEN
  87 at amp_nat=0/amp_door=0 (the invariant surface).
- **Byte-inertness:** new knobs default 0 and must reproduce the prior
  surface byte-for-byte at 0. Verify with a conserve-diff or the full battery.
- **Both kernels:** C (`kernel/freecell.c`) is record; Go (`fab/`) is the A/B
  mirror — keep them in sync (Sincos arg order, `sFree` not `S_FREE`,
  `P.AmpDoor!=0` not boolean, lut vs libm choices).
- **Commits:** one per face/fix, detailed message (the project runs
  commit-per-step; the user has authorized commits as part of v93 work).
- **Pre-v89 ban:** never consult versions before v89.

## 7. The next face (when the user wants to continue)

Faces **A/B/C below are EXECUTED (2026-08-08)** — see `FACE_A.md`,
`FACE_B.md`, `FACE_C.md` for the data. Outcomes in brief:
- **A. Symmetrized hop schedule — DONE.** New byte-inert knob `hop_order`
  (0=seq default, 1=Strang, 2=randomized) in BOTH kernels. Battery ALL
  GREEN 87. LIVE vacuum ring holds W=+2 to t=1000 under Strang/random
  (sequential degrades to W=+1); frozen Strang holds to ~t=108 then succumbs
  to the independent q_detune condensation. e3b dt-spread 4%(seq)->1.3%(Strang)
  ->0.3%(rand). **Strang = best retention; random = best dt-invariance.**
- **B. Lit-bath refill + arg(ψ) door — DONE, prediction REFUTED.** noise_amp
  lights the medium (cond=41-52 vs 0 dark). In the bath, contact dephasing
  (A.2 ~8 t.u.) beats the door's beat cycle (~7 t.u.). Fair test (vacuum +
  Strang + noise_amp): door traffic (~28 fires) destroys W by t~14-16
  REGARDLESS of amp_door. Composing a random-phase field "coherently" is
  still random injection. amp_door is NOT the retention rescue; the ceiling
  is contact dephasing + the door's own random-phase injection.
- **C. Condensation lane — DONE.** Self-trapped hoards persist to t=300
  (stable mass formation, no lock/gate/door). Deep corner (qd3.6 amp2):
  Em_max~7 (3×cap), frozen envelope. 189 hoards in a long-tailed hierarchy
  (0.5→7.2). **Radiance SELECTS sizes**: V3a tax truncates the spectrum at
  the fixed point (Em*~1.55); deep self-traps resist. Creation-adjacent.

The retention face (§2 last row, `RING_DNLS.md`) EXECUTED both prior routes.
What it opened (A/B/C above, now all done):

Also still open (lower priority per user): formal tolerance-C≡Go abx at
amp_nat>0; armed L1-C (ρ_coh≈0.77 coherent bath); a COHERENT phase-matched
field driver (the one case where amp_door might still help — untested).

## 8. Recent git (the v93 arc)
```
(HEAD)  face C: condensation lane (mass formation; radiance selects hoards)
42a6aab face C: condensation lane — stable hierarchical mass formation
f70f238 face B: lit-bath refill + arg(psi) door — NOT a retention rescue
2f55e78 face A: symmetrized hop schedule (hop_order) — fixes Trotter artifact
(retention face: ring seed_mw + RING_DNLS.md, routes A+B executed)
d178813 items 1/2/4 writeup + battery ALL GREEN 87
1239bd4 item 4 (matter-winding retention): unitary scrambles localized winding
c7e77e2 item 2 (group velocity): tagged-centroid meter wrong; fd is seed-robust +x
6c7adba item 1 (linearize-τ): gate was artifact; linear recovers transport
871337d post-review honesty pass (L1_FINDINGS correction, README, Re→Im)
cc16db7 review FIXES: empty-cell reset + symmetric reverse door + amp_door cont.
227a832 FIX Go precession sin/cos swap + lut mirror
32138ec arg(ψ) door face + 3-reviewer consultation
6dcc0c8  STEP ZERO (v92→v93 carry, baseline ALL GREEN 87) — the byte-inert anchor
```
`git show 6dcc0c8:v93/runs/conserve_c.log` is the canonical byte-inert
reference for amp_nat=0 diffs.

## 9. User's standing preferences this session
- Wants honest measurement over positive claims; invited the external review
  and endorsed walking back overclaims.
- Authorized full v93 implementation + kernel edits + commits ("no additional
  authorization required").
- Less worried about C/Go FP differences right now (reasonably explained);
  prioritize physics over the tolerance-abx.
- Most interested in: linearize-τ, group-velocity reframe, slot-borne/retention
  (items 1/2/4 — DONE). Next likely: closed-ring retention substrate / DNLS.
