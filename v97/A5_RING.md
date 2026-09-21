# A5_RING — the registered A.5 arm, unbundled (2026-08-13, PLAN step 2)

Config-only on the carried Go kernel; one knob per arm; graded with
the step-1′ exact-ledger meters (`REGRADE.md`). Base: leads COMMON
(`L=16 dt=0.02 T=240 amp_nat=2 amp_logate=1 hop_order=1
seed=20260802 snap_every=200`) + `exp=ring ring_n=6 seed_mw=2 bath=1
wf_on=1 k_rad=0.05 promo=2`. Script: `runs/a5.sh`; streams/logs in
`runs/a5/`.

## rg2a — ring + bath + V3a + promo=2, nothing else

Drift −2.1e-15 (FP floor). **cond = 0.000** — the intake door never
fired: the A.5 diagnosis (dark bath, workfn shut) reproduces under
promo. Unlike vacuum p2_ring, the outtake is now alive at a trickle:
rad = 0.841 over 240 t.u., fired in ~3 atom-grain bursts (visible as
isolated `boxPO` spikes, e.g. 0.275 at t=52), with qcnvF = 0.714
still parked in credit. The ring drains into the bath (census Em
4.88 → 2.38, n 6 → ~4 — the A.2 leak), staying `parked` only within
the 0.20 bar. H1 grade: **1/5** (parked), `-tag` identical.

**Instrument finding (measured here, fixed in the meter):** promo>0
books every door event at the oriented link's **far cell**
(`fab/promo.go` dest), so the machine's own deposits land one cell
outside its census — motif ⟨O⟩ read 0.000 while the box ledger
carried the full rad total (ΣboxPO = 0.8405 = RESULT rad exactly).
`dipole -dilate 1` adds one ring of link-neighbours to the I/O book
(default 0; direction acknowledged — it can only raise motif I/O;
it is census-scope correction, not a softer bar). With `-dilate 1`
the motif sees its outtake (⟨O⟩ = 0.0115 in the best window); the
verdict stays 1/5 — the blocker is intake, which is physics, not
scope.

**Binary answers (pre-registered):** door fires? — outtake yes
(grain-limited trickle), intake no. Motif stays C6? — no: it leaks
toward n≈4 but survives the window as a connected motif.

## rg2b — + amp_door=1 (the A.5 "untested refill", now under promo)

Drift −1.9e-15. **cond = 0.000 again** — with no intake events there
is no phase to compose; the refill stays untested-by-darkness, as
predicted. (Not byte-identical to rg2a: rad = 0.852 vs 0.841 —
amp_door's reverse-evap compose engages on the radiance trickle.)
H1 grade 0/5 in the best window (the ring's drain crossed the 0.20
parked bar there). The blocker is unchanged: **no field in the bath
⇒ the workfn gate never opens.** Next: light the field (rg2c).

## rg2c — + noise_amp=0.4 (light the field honestly; no bath_xb)

*Pre-registered before results:* cond > 0 for the first time in this
line (field present, workfn opens on lit cells); the foam will
box-fill toward its own load (Theorem 3.4 — noise is g=0), which is
why the grade is read on the **tagged ring** (`-tag -dilate 1`), not
the foam. Open questions the arm decides: does the tagged motif's
Φ/Em move toward the 0.015 bar as its x climbs (PHI_BAR Statement A
says it can at x ≳ 0.37); do I and O land fore/aft (sep) on a
door-live ring; J clause expected to FAIL as a loop (helix is rg2d).

**Results (drift −3.4e-14; cond = 120.4, rad = 4.57):**

| clause | tagged ring (-tag -dilate 1) | foam census |
|---|---|---|
| Φ > 0 | **PASS — Φ/Em = 0.0154** (⟨I⟩=0.64 ⟨O⟩=0.13, exact) | FAIL 0.0027 |
| I/O separated | **PASS** — sep 1.14, polar 0.57 (8/26 books) | FAIL (polar 0.249) |
| parked | FAIL — Em 4.25→5.97 (+29%); box ×4.0 (noise keeps pumping ambient) | — |
| J steady | FAIL — Jdip≈Circ, dir-stable 0.69 (loop, as registered) | — |
| Es at const M | FAIL — tracks the box-fill dM (downstream of park) | — |

**Both pre-registrations HIT, and the first Φ-clause PASS of the
programme is on the books** — a lit bath at g = 0 clears the Φ bar on
a small tagged motif, exactly as `theory/PHI_BAR.md` Statement A
computes (in the very stream whose *foam* census fails by census
dilution, Statement C). The intake is flood-not-balance (I/O = 5.0),
so the ring grows: parked/Es wait on the drive stabilizing (longer T
and/or milder noise — a follow-up arm), and J waits on the helix
(rg2d). 2/5 tagged.

## rg2d — helix tilt kx=1.1 — NULL BY APPARATUS

**Byte-identical to rg2c** (cond = 120.438297, drift bit-equal). Root
cause found in the seed code, not the physics: `fab/run.go` applies
the kx tilt first and then `seed_mw` **assigns** `th2 = m·φ`
("overrides the chain lock"), so with `seed_mw=2` the tilt is
silently discarded. The composed helix seed (`th2 = m·φ − kx·x`)
needs a `fab/run.go` edit — **kernel-file change, needs user
sign-off; registered in PLAN §4**, noting it is behavior-preserving
for every previously-meaningful combo (kx+seed_mw is today a silent
no-op).

## rg2dp — the tilted UNWOUND ring (config-reachable polarity probe)

`seed_mw=0 kx=1.1` + the rg2c stack. *Pre-registered:* H1's J clause
needs a dipole, not a winding (winding is H3): a ring with a pure
x-phase-gradient carries J along +x on both arcs — Jdip should rise
against Circ. Φ/sep should behave like rg2c (same door, same bath).
Open: whether the tilt survives contact long enough to matter (B.3
measured ~4 t.u. stall on blobs; the bonded ring may differ).

**Results (drift −3.4e-14; cond = 131.4, rad = 2.95):** **2/5** —
sep PASS (1.12, polar 0.61) and parked PASS (Em 4.15→4.51, the
unwound ring grows slower). Φ FAIL at 0.0066 (I/O = 6.7 — less
outtake than the wound ring: rad 2.95 vs 4.57). **J FAIL, and this
answers the open question: the seeded tilt does NOT pin J**
(dir-stable 0.69, Jdip ≈ Circ — B.3's stall reproduces on the
bonded ring in a noisy bath). Across every arm of this campaign the
J clause now fails in the same band, dir-stable 0.54–0.69 vs 0.80:
**direction-pinning is the one clause no config reaches.**

## rg2f — the composed helix (seed_mw=2 + kx=1.1, post seed-fix)

Enabled by the user-authorized `fab/run.go` seed-order fix (battery
ALL GREEN 87 after — `REGRADE.md` §4b). SEED lines confirm the
composition (`winding` then `tilt applied on top`). Drift −3.4e-14;
cond = 120.2, rad = 4.71.

**2/5 — with the campaign's strongest Φ pass: Φ/Em = 0.0196** (rate
0.0042/t.u.), sep PASS (1.46, polar 0.59). Parked FAIL (Em 3.06→5.59
— still filling toward cap under the noise drive). **J FAIL:
dir-stable 0.65 — the composed helix does not pin J either.**

**Campaign synthesis on the J clause:** seven arms — vertex control,
oriented door, oriented+space-leg, hot, noisy, tilted-unwound,
wound, composed helix — all land at dir-stable 0.54–0.69 vs the
0.80 bar. **No seeded phase structure pins the direction of J in a
lit bath; the pin has to be dynamical** (a door-locus memory: the
signed per-link promotion credit with FLOW-τ decay — the grok
consult's preferred face — or the one-bit latch it refines). That
face needs its own brief and user sign-off; it now stands on seven
config nulls plus `theory/DOOR_RULE.md` (no linear local door is
Δm-selective, so the pin cannot come from door geometry alone).

## rg2e — park probe: rg2c stack at T=1000

*Pre-registered:* the ring climbs while noise pumps ambient; parked
(dM/dt) can only pass if the ring saturates near cap (cap + evap
defend) while a ≥100 t.u. window stays within the 0.20 bar. If it
parks AND Es-at-constant-M still fails only via box-fill tracking,
the Es clause needs a stabilized ambient (drive-off pool), which is
the next design conversation, not a knob spree.

**Results (drift −3.6e-13, consistent with T-scaled roundoff;
cond = 538, rad = 140):** 2/5, best window at the END (t=900–1000).
The **box saturates** (|ΔEmTot|/EmTot = 0.014 — noise drive balanced
by foam radiance), and with ambient stabilized the **Es clause
passes for the first time on a lit arm** (coreThru 1.72 vs dM 0.57,
not tracking), sep passes, and J posts its closest approach
(dir-stable 0.72, Jdip/Circ = 1.96 — a dipole at last, still under
the 0.80 bar). But the ring itself **overshot and is draining**
(Em 4.14→3.26 late; I/O flipped from 5.0 early to 0.24 late): it
*transits* balance instead of sitting at it. That is qgraph
Thm 3.4 in vivo — at g = 0 there is no restoring root above
ambient; the prediction registered in `theory/PHI_BAR.md` §B.

## Campaign verdict

Every H1 clause now passes somewhere in this campaign — and never
together:

| clause | passes in | fails everywhere else because |
|---|---|---|
| Φ > 0 | rg2c, rg2f (lit ring at x > 0.37) | dark bath / grain / census dilution |
| I/O separated | rg2c, rg2dp, rg2e, rg2f | no door traffic |
| parked | rg2a, rg2dp | climb (early lit) or drain (late lit) |
| J steady | **never** (0.54–0.72 across all seven arms) | no seeded structure pins direction |
| Es at const M | rg2e (saturated box) | dM tracking while not parked |

**H1's emptiness is now localized to exactly two missing dynamical
mechanisms**, each a named face with theory behind it, neither
authorized:

1. **The J-pin** — a door-locus memory (signed per-link promotion
   credit with FLOW-τ decay, or the one-bit latch it refines).
   Seven config nulls say seeds cannot do it;
   `theory/DOOR_RULE.md` says no local door geometry can either.
2. **The restoring intake** — the stimulated g-term (qgraph
   Thm 3.5), whose absence rg2e exhibits as overshoot-and-drain.
   This is the H2 face (`theory/PHI_BAR.md` §B); g ≈ 0.3 for a
   parked hold 8% above ambient.

Everything else on the H1 card is apparatus, and the apparatus now
works: lit bath + tagged small motif + exact ledgers + book
dilation.

---

## Dated correction (2026-08-13 evening, grok state-check — `consult/REVIEW_grok_state.md`)

Grok re-graded all seven streams (dilate on/off, sign-runs) and
confirmed the clause table reproduces — and refuted the "exactly two
missing mechanisms" reading. Three corrections adopted:

1. **The Φ pass is dilate-carried.** With `-tag` only, rg2c reads
   ⟨Φ⟩=0 (0/26 books) and rg2f falls to 0.0049 FAIL. Under
   `-dilate 1`, 78% of I and 91% of O sit on halo cells (dilated I is
   16% of whole-box intake) — the one-ring dilation imports the
   bath's own box-fill books, the inverse of the census-dilution
   error PHI_BAR §C names. The dest-scope defect rg2a found is real;
   the fix overshot. **There is a Φ-shaped signal under a same-day
   convention; there is not yet a Φ pass on the machine's own
   books.** The honest correction is dest-only booking, not a
   neighbourhood halo.
2. **rg2e's caption was wrong.** "Thm 3.4 in vivo" requires the ring
   near xb; measured x peaked at 0.53 and sat at 0.28→0.22 late. The
   late drain is **supply death** (intake dies when the foam
   saturates), not a missing restoring root — the ring never reached
   ambient. The deciding arm (quench/milder noise after fill) was
   named in PLAN §2 and never run. rg2e is NOT an exhibit for the
   g-term.
3. **The J verdict rides the SKIP meter, and DOOR_RULE was
   mis-cited.** Snapshot Jdip sign-runs are 1–4 frames on every lit
   arm — the H0 current-SKIP shape — so dir-stable 0.54–0.72 cannot
   decide a face until `pJ`/`fd` exists and is H0'd on e3b.
   DOOR_RULE is a theorem about Δm (winding, H3), not about the
   first moment of J (H1); it does not license "the pin must be
   dynamical." The lawful core that remains: `promoPick` steers by
   the sloshing snapshot compass (a steer defect in the write path).

**Restated localization (three items, not two faces):** (1) book
scope — instrument debt, dest-only; (2) matched intake at constant
ambient — apparatus (the un-run quench park probe) before any H2
talk; (3) door steer / J-pin — face candidate only after a current
meter exists. The g-term stays H2 and must not be opened to buy H1.

**Adopted queue (grok, replication before faces):** 0) dest-only
book re-grade of rg2c/f/e with pre-registered outcome bands;
1) three seeds of rg2c and rg2f; 2) the quench park probe;
3) `pJ` + H0 against e3b `fd`; 4) optional promo=3 lit ring
(DOOR_RULE's nonlinear escape). No face briefs until 0–2 are read.
