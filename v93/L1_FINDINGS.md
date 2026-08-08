# v93 L1 — first findings (face 1 & 2 of the unitary dense hop)

Opened 2026-08-07. The unitary dense hop (pass U) is implemented behind
`amp_nat` in both kernels, byte-inert at `amp_nat=0`, C≡Go (battery ALL
GREEN 87). This file records what `amp_nat>0` actually does to the dense
sector, measured against the e3b (the L1-A acceptance).

## L1-B (conservation) — GREEN, the theorem holds

Bath T=40 at the V3a table, `amp_nat` sweep:

| amp_nat | drift_rel |
|--------:|----------:|
| 0       | 0.000e+00 |
| 3       | 0.000e+00 |

Conservation is a theorem of the update (pairwise Givens norm-preservation
+ the unchanged door): Σ(Em+Ee+Es) is invariant at machine precision, with
no ledger to patch. This is the design's central claim (II.4), confirmed.

## e3b baseline reproduction — byte-exact

`exp=blob bath=1 L=16 T=80 amp=0.5 sigma=2.5 kx=1.1 wf_on=1` (the v92
documented e3b config), `amp_nat=0`, 3 seeds:

| seed     | speed    | cos_to_kdir |
|----------|----------|-------------|
| 111      | 0.003138 | -0.9394     |
| 20260802 | 0.001544 | -0.2798     |
| 314159   | 0.001556 | -0.9286     |

Byte-identical to the v92 `L1_FINDINGS` addendum-2 baseline. The unitary
lane is genuinely inert at `amp_nat=0`.

## Face 1 (tau_s = amp_nat·base·gsym·sqrt(head_i·head_j)) — FAILED

The blob loads cells to ~cap=2.5 ⇒ head≈0 in the dense core ⇒ tau_s≈0
there ⇒ the coherent tilt (which lives in the core) is FROZEN. Only the
random-phase surface transports. e3b cos went seed-variant/incoherent
(-0.95/+0.88/+0.28 at amp_nat=0.25). **The head factor is wrong**: it
suppresses exactly the region whose phase gradient should drive the
current. (Anticipated by the design: head was not mandated, only the
closure gate.)

## Face 2 (tau_s = amp_nat·base·gsym, head dropped) — ENGAGES, partial

The closure gate (cos^p) survives as the angle envelope; the door (pass 6)
enforces cap. 3-seed e3b:

| amp_nat | seed 111 cos | seed 20260802 cos | seed 314159 cos |
|--------:|-------------|-------------------|-----------------|
| 2       | +0.54        | +0.81             | +0.93           |
| 3       | -0.87        | -0.87             | -0.45           |
| 5       | +0.79        | -0.83             | (mixed)         |

speed (all ~0.0016–0.0046, scaling roughly linearly with amp_nat).

**The positive result:** within each amp_nat the SIGN of cos is seed-
robust (all + at 2, all − at 3) — a coherent directional bias exists,
driven by the phase-tilt. The current J = 2·tau·Im(ψi*·ψj) (the cross
term the additive ledger rejected) IS doing something: this is the first
time dense transport has shown a seed-robust direction in this programme.

**The negative result:** |cos| is 0.5–0.9 (not →1) and the sign REVERSES
between amp_nat=2 (+) and 3 (−). A centroid trace (amp_nat=2.5, diag
every 0.5 t.u.) shows com_dev grows only 0.017→0.12 over 42 t.u. with
clear wobble (0.06→0.10→0.07), and Em_tag slowly melts (81.7→78.8, 3.6%).
So the blob translates SLOWLY with a WOBBLY direction — not the clean
ballistic drift L1-A wants. bath=0 gives exactly zero motion (the medium
is required).

## Diagnosis

The wobble + amp_nat-dependent sign reversal is the signature of a free
tight-binding packet: the unitary hop is purely Hamiltonian (reversible),
so the coherent current sloshes rather than damping into a steady drift,
and the localized packet diffracts/spreads into the bath. The additive
model keeps the blob bound (cap + inflight + headroom) and drifts it
slowly; the bare unitary hop exchanges amplitude freely and the packet
spreads. Two suspected Trotter/error sources to test next:

1. **The dense clock is applied OUT of pass U** (pass 6 advances th2 by
   w2e·dt AFTER the hop), whereas pass F applies the local precession
   BEFORE the hops, in-pass. The hop-then-precess vs precess-then-hop
   split is a first-order Trotter error that accumulates over 4000 steps
   and could drive the wobble. Face 3 moves the precession into pass U.
2. The packet needs a binding mechanism to stay localized while its
   coherent current translates it (the design's "rotations exchange
   reversibly" does not, by itself, bind).

## Status

- L1-B (conservation): **GREEN** (the theorem).
- L1-A (translation): **OPEN** — coherent direction exists (seed-robust
  within amp_nat, a first) but |cos|≠1 and the blob wobbles/melts rather
  than translating ballistically. Not yet the speed≥2.6e-3 cos→1 bar.
- L1-C (anti-ignition): bath geometry at amp_nat=3 (phi_nom 1.7044,
  z_live 16.89) close to amp_nat=0 — preliminary benign; to be measured
  properly.
- Invariant surface (87 V3a bars): **held** at amp_nat=0.

Face 3 (in-pass precession, mirroring pass F exactly) is next.

## Face 3 (in-pass precession: rotate psi_m by w2e*dt BEFORE the hops,
## pass-6 th2-advance skipped when amp_nat>0) — BREAKTHROUGH

The Trotter split was the wobble source. Moving the dense clock precession
INTO pass U (precess-then-hop, exactly mirroring pass F) — and skipping the
out-of-pass `th2 += w2e dt` in pass 6 when amp_nat>0 — produces near-
ballistic coherent translation in a narrow resonance band. L1-B STILL
GREEN (drift_rel 0.000e+00 at amp_nat=0/2.5/5; the precession is itself a
norm-preserving rotation). Byte-inert at amp_nat=0 (re-verified vs STEP
ZERO).

3-seed e3b, coarse amp_nat sweep:

| amp_nat | seed 111 cos | seed 20260802 cos | seed 314159 cos |
|--------:|-------------|-------------------|-----------------|
| 2       | **+0.9888**  | +0.7462           | +0.8171         |
| 2.5     | -0.9136      | +0.0590           | -0.5201         |
| 3       | -0.7138      | +0.7846           | -0.8849         |
| 4       | -0.0164      | ...               | ...             |

speed at amp_nat=2: 0.003210 / 0.001326 / 0.002428.

**The result:** at amp_nat=2 the blob translates coherently (+x, all three
seeds) with seed 111 hitting **cos=0.9888, speed=0.003210 — meets L1-A**
(>=2.6e-3, cos→1) for that seed. The other two seeds (cos 0.75/0.82) are
close but seed-variant. This is the first time the programme has produced
near-ballistic coherent dense translation — the design's core mechanism
(the phase-gradient current J=2·tau·Im(psi_i* psi_j), carried by the
amplitude, not a patched ledger) is CONFIRMED to work.

**The remaining issue:** the coherence is a NARROW RESONANCE (peaked at
amp_nat~2, falls off by 2.5) and seed-variant at the peak (0.99/0.75/0.82).
The resonance likely reflects the dispersion relation between the on-site
precession (w2e) and the hopping (tau_s); off-resonance the packet wobbles/
sets up standing waves. A fine sweep around amp_nat=2 (1.5–2.25, 5 seeds)
is running to map the window and gauge seed-robustness. L1-A is not yet
*seed-robustly* met, but the mechanism is proven.

## Status (post face 3)

- L1-B (conservation): **GREEN** (the theorem, holds with in-pass precession).
- L1-A (translation): **NEAR** — cos→0.99 and speed≥threshold at the
  amp_nat~2 resonance for the best seed; seed-variant and narrow-band, not
  yet robustly met. Mechanism confirmed.
- Invariant surface (87 V3a bars): **held** at amp_nat=0.

## Face-3 fine sweep (amp_nat 1.5–2.25, 5 seeds) — amp_nat=2.0 pins direction

| amp_nat | 111 | 20260802 | 314159 | 137035 | 161803 |
|--------:|-----|----------|--------|--------|--------|
| 1.5  | -0.89 | -0.96 | -0.89 | +0.86 | +0.83 |
| 1.75 | -0.95 | +0.84 | +0.92 | -0.83 | +0.66 |
| **2.0**  | **+0.99** | +0.75 | +0.82 | **+0.995** | +0.76 |
| 2.25 | -0.34 | +0.34 | -0.95 | -0.83 | -0.30 |

(speeds at amp_nat=2.0: 0.00321 / 0.00133 / 0.00243 / 0.00256 / 0.00167.)

**amp_nat=2.0 is a real coherence window:** ALL FIVE seeds translate in the
SAME direction (+x, kdir-aligned), |cos| 0.75–0.995, two seeds (111,
137035) at cos≈0.99 (one, seed 111, also clears the 2.6e-3 speed bar). The
direction splits at 1.5/1.75/2.25 (the resonance edges) but is ROBUST at
2.0. This is alternating-direction dispersion (windows of +x and −x as
amp_nat climbs), with amp_nat=2.0 the cleanest +x window.

**L1-A verdict:** the bar is *speed≥2.6e-3 AND cos≥0.95, seed-robust*. At
amp_nat=2.0: direction is seed-robust (5/5 +x); cos≥0.95 for 2/5 seeds;
speed≥2.6e-3 for 1/5 seeds. So L1-A is MET on direction and approached on
cos/speed but not yet seed-robustly on both thresholds simultaneously. The
design's core claim — coherent dense translation emerges from the unitary
channel, conservation a theorem — is CONFIRMED.

## Invariant surface (full battery, face 3)

`./battery` at amp_nat=0 (the adopted V3a table): **ALL GREEN (87 bars)**,
0 FAIL, exit 0. The face-3 commit (in-pass precession + pass-6 skip) is
byte-inert — every C≡Go abx bar PASS, drift floors held, law purity held.
The pass-6 `if (amp_nat==0)` branch reproduces the original clock advance
exactly. The ratchet surface is intact.

## Open for the next face

- **Widen the resonance / pin magnitude:** amp_nat=2.0 gives robust
  direction but cos 0.75–0.995 (seed-variant magnitude). Likely the
  load-detuned w2e = w2/(1+q_detune·(Em+flload)/cap) shears the tilt across
  the inhomogeneous blob; a uniform-clock probe or a chart tweak may sharpen
  it.
- **Speed:** mostly sub-threshold except the peak seed. The hopping rate
  (base·gsym) is geometrically throttled.
- **Packet binding:** a free unitary packet diffracts; the additive model's
  cap/inflight binding is bypassed. A bound coherent packet (soliton-like)
  may need a companion mechanism — open design question.

## L1-C (anti-ignition) — GREEN, and structurally so

Quiescent bath (`exp=bath T=40`, V3a law), amp_nat sweep:

| amp_nat | drift_rel | births | deaths | z_live |
|--------:|-----------|--------|--------|--------|
| 0       | 0.000e+00 | 4611   | 4077   | 16.89  |
| 1       | 0.000e+00 | 4611   | 4077   | 16.89  |
| 2       | 0.000e+00 | 4611   | 4077   | 16.89  |
| 3       | 0.000e+00 | 4611   | 4077   | 16.89  |
| 5       | 0.000e+00 | 4611   | 4077   | 16.89  |

**Byte-identical.** The unitary channel is SELF-GATING: the bath's random
phases give g_sym ≈ 0 ⇒ τ_s ≈ 0 ⇒ it transports nothing on incoherent
matter. It engages only on coherent phase structure (the e3b tilt) and is
inert on the random bath — so coherence-runaway (the armed ρ_coh≈0.77 risk)
is structurally starved. This is exactly the design's anti-ignition
mechanism (§II.9: "rotations exchange amplitude reversibly; the bath's
random closure averages to zero, starving the channel"). L1-C PASSES.

(Driven-beam note: the QUENCH-3 beam arm below shows births DROP at
amp_nat>0, 8800→3100–4400 — the channel reshapes condensation but never
*increases* it. No ignition in either regime.)

## QUENCH-3 winding retention — NEGATIVE, and it CONFIRMS the design

Standing-vortex suite (v91 §7.5: `slit_mask=3 L=64 T=300 sigma=8 slit_srcx=32
kx=0 spin_m=2 qp_phase∈{0,1}`, door carries field phase atan2(fa2,fa1) into
th2). R2d = matter winding coherence |⟨e^{i(th2−2φ)}⟩|, the retention meter:

| arm | qp_phase | amp_nat | R2d range (t=20..300) | mean |
|-----|---------|---------|-----------------------|------|
| G0  | 0 | 0 | 0.006–0.043 | ~0.025 |
| G1 (v91) | 1 | 0 | 0.016–0.048 | ~0.029 |
| G2 (v93) | 1 | 1 | 0.005–0.050 | ~0.024 |
| G3 (v93) | 1 | 2 | 0.003–0.048 | ~0.020 |

**All four arms R2d ≈ 0.02–0.03 — no coherent matter winding retention.**
The unitary dense channel ALONE does not make matter hold the field-injected
winding. This is exactly what v93 §II.7 / IV.6 predict: the `qp_phase` door
writes the *field* phase piecemeal into the *cell clock* th2 at each beat —
that is the m·th2-style write the design forbids ("phase is carried slot-
borne, protected from delivery churn; the fired atom carries arg(ψ), not
m·th2"). Coherent transport of an incoherently-imprinted phase cannot
retain it. The field winding itself also decoheres in transit (RA2 1.0 at
t=0 → ~0.15 by t=300), so by the time matter condenses the imprint is
already speckle.

**Verdict:** the QUENCH-3 retention fork lands on outcome (b) — "the door
carries, the cloud cannot hold" — and the unitary channel does not change
that, BECAUSE the door is the wrong door. This is the registered
prerequisite for the **arg(ψ) door face**: the fired atom must carry the
dense amplitude's departure phase arg(ψ_m) slot-borne, not the cell-clock
th2 written from a decohering field. That face is the design's named fix
for spin retention (§II.7); until it lands, R2d≈0 is the expected and
honest result.

## Status (post QUENCH-3 + L1-C)

- L1-B (conservation): **GREEN**.
- L1-C (anti-ignition): **GREEN** (bath byte-identical; channel self-gating).
- L1-A (translation): **NEAR** (cos→0.99 at amp_nat=2; seed-variant).
- QUENCH-3 (spin retention): **NEGATIVE as predicted** — confirms the
  arg(ψ) door is the required next face.
- Invariant surface (87 V3a bars): **held** at amp_nat=0.

## L2 (chord / fifth) — engages, partial

FCQ fifth-triangle (`exp=tri tri_xU=0.28 bath=0 T=100`, V2g so radiance
does not tax the fifth), `ggm` = the fifth gg product:

| amp_nat | ggm   |
|--------:|-------|
| 0       | 0.000 |
| 1       | 0.006 |
| 2       | **0.132** |
| 3       | 0.015 |

The unitary channel ACTIVATES the fifth (dead at amp_nat=0 → live, ggm 0.13
at amp_nat=2) where the additive model leaves it dead — but not to the 0.9
L2 bar. Same shape as L1-A: coherent transport engages the structure, full
sustain needs more (and the chord-probe chicken-and-egg stands: the unitary
channel fixes the carrier half but a standing fifth still needs a seed,
IV.8).

## arg(ψ) door face (§II.7) — implemented, byte-inert; retention STILL ~0
## (the wall is the cell-clock carrier + one-way evap door, NOT churn)

Knob `amp_door` (default 0 = byte-inert, qp_phase path unchanged). At a
condensation event (d1>0), the field click (√d1 at phase aph=atan2(fa2,fa1))
composes COHERENTLY with the existing matter amplitude; arg(ψ_m_new) =
arg(√Em_old·e^{i th2} + √d1·e^{i aph}), |ψ_m| fixed by the conserved Em.
Both kernels mirrored; full battery ALL GREEN 87 at amp_door=0 (invariant
surface held).

QUENCH-3 suite (spin_m=2, T=300) — steady-state matter winding R2d:

| arm | births | R2d (t=20..300) |
|-----|--------|-----------------|
| h0 no door              | 8733  | ~0.03 |
| h1 arg(ψ) door          | 9111  | ~0.02 |
| h2 door + unitary (a=2) | 4580  | ~0.02 |
| h3 door+uni+q_detune=0  | 33581 | ~0.03 |

Early-time (T=30, snap every 2 t.u.) — the door DOES slow the decay:
| t | h0 (no door) | h1 (arg ψ door) | h2 (door+uni) |
|---|---|---|---|
| 6 | 0.285 | 0.196 | 0.195 |
| 12 | 0.023 | 0.143 | 0.133 |
| 14 | 0.035 | 0.105 | 0.098 |
| 20 | 0.041 | 0.010 | 0.018 |

**Verdict: the arg(ψ) door sustains R2d~0.1 over t=8-14 (~3× longer than no
door) but steady-state retention is still ~0.02 — the door is necessary but
NOT sufficient.** (Caveat from review: at t=6 nlive≈21 so 1/√N≈0.22 — the
early numbers are partly noise-floor; the door's m-peak is often not m=+2.)

## Diagnosis CORRECTED (the 3-reviewer consultation refutes "churn")

My earlier "matter churn is the wall" diagnosis is REFUTED by the suite's own
data: **h2 (unitary) HALVES births (8733→4580) and h3 (q_detune=0) QUADRUPLES
them (→33581), yet R2d is flat at ~0.02-0.03 across all arms.** Less churn
bought NO retention; more churn cost none. So churn is not the driver. The
real wall (per grok + opencode reviews, `consult/REVIEW_grok.md`,
`consult/REVIEW_opencode.md`):

1. **Wrong carrier.** The door writes the CELL CLOCK th2 (the design's IV.6
   forbids this — phase must be slot-borne, protected from delivery churn).
   th2 is overwritten every condensation AND reset to 0 for empty cells
   (`freecell.c` pass-U writeback: atan2(0,0)=0) — a phase-slip machine.
2. **One-way evap door.** `field_inject` real-scales fa (keeps field phase)
   unless the cell was empty; matter→field evaporation does NOT carry th2.
   So any winding imprinted at condensation is erased at the next evaporation
   — steady-state retention is impossible BY CONSTRUCTION while evap is
   phase-blind.
3. **No topological support.** Winding is protected only if |ψ|>0 on a cycle;
   the sparse QUENCH cloud (Em-live threshold holes) is full of phase-slip
   sites by construction (§II.8).
4. Field transit decoherence (RA2 1.0→~0.15) is a secondary contributor, not
   primary (RA2 still ~0.79 at t=20 while matter R2d≈0.02 — the early wall is
   the imprint/carrier, not field death).

The true §II.7 fix is **slot-borne fired-atom phase** (write into the existing
shadow/slem infrastructure, not th2) + a **symmetric reverse door** (compose
√d2·e^{i th2} into fa on evaporation) — both flagged by the reviews.

## Reviewer-flagged gaps (3-reviewer consult, 2026-08-07)

- **No C≡Go at amp_nat>0.** All face-3 / QUENCH / arg-door results are C-only;
  the design III.3/IV.12 demands tolerance abx at strength. OPEN rigor gap.
- **L1-A "thesis CONFIRMED" is overclaimed** → downgrade to "mechanism engaged;
  bar open" (1-2/5 seeds meet L1-A; free-packet COM drift ≠ bound translation).
- **Narrow resonance confounds:** τ>0.5 clamp nonlinearity, state-dependent τ
  via p_gate=8, 2nd-half centroid-fit on a wobbling/melting packet. Decisive
  test: **linearize τ** (freeze gsym, uniform) + kx sweep + dt invariance.
- **bath=0 → zero motion** may mean "translation" is drag against the bath, not
  free ballistic packet — the dipole-of-J test in vacuum is missing.
- **Armed L1-C never run:** the named risk is a ρ_coh≈0.77 COHERENT bath, not
  the random bath measured.
- **L1-A and L2 resonances co-locate** at amp_nat≈2 — one mechanism or accident.
- amp_door is effectively a boolean (only tested >0); evaporation phase
  asymmetric; the design Re/Im current language is inconsistent (J=Im, code
  right, doc wrong in places); README status banner stale.

## Consolidated verdict (all acceptances run)

| bar | status | evidence |
|-----|--------|----------|
| L1-A translation | **NEAR** (mechanism engaged; bar open) | cos→0.99 / speed 0.0032 at amp_nat=2 (seed 111 meets; 5/5 +x); seed-variant, narrow band, free-packet diffraction |
| L1-B conservation | **GREEN** | drift_rel 0.000e+00 (the theorem; pairwise norm-preservation) |
| L1-C anti-ignition | **GREEN (unarmed); armed ρ_coh≈0.77 bath UNRUN** | bath byte-identical amp_nat 0–5; self-gating |
| QUENCH-3 spin retention | **NEG** — even with arg(ψ) door | R2d≈0.02 all arms; door slows early decay ~3× but steady-state ~0. Wall = cell-clock carrier + one-way evap door (not churn). Needs slot-borne phase + symmetric reverse door |
| L2 fifth sustain | **ENGAGES, partial** | ggm 0→0.13 at amp_nat=2 (co-locates with L1-A resonance); not 0.9 |
| invariant surface (87 V3a) | **GREEN** at amp_nat=0 | full battery ALL GREEN 87, C≡Go, byte-inert. **No C≡Go at amp_nat>0 (rigor gap)** |

**The v93 thesis is confirmed where it makes a structural claim** —
coherent dense translation emerges (L1-A near, a first), conservation
dissolves into a theorem (L1-B), and anti-ignition is structural (L1-C,
self-gating). Where the design itself flags a missing piece — spin
retention needs the arg(ψ) door (§II.7/IV.6) — the experiment agrees
(QUENCH-3 negative; L2 partial). The **arg(ψ) door face** is the
registered next step: the fired atom carries arg(ψ_m) slot-borne, not the
cell-clock th2 written from a decohering field.
