# v97 PLAN — post-leads next steps (2026-08-13)

**Authority:** user, this session — "First ensure the plan is written
down. Then execute on it." Inputs: the fresh-eyes review (claude, this
session) as corrected by the resumed grok consult
(`consult/REQUEST_next_steps.md` → `consult/REVIEW_grok_next.md`).
Scope per user: **Go kernel only; C parity out of scope for now.**
Config-only + offline-meter work. **No law change. No bar change.
No kernel face is opened by this plan.**

## 0. The corrected findings this plan stands on

1. **The heuristic I/O meter INFLATES Φ** (grok, verified against the
   kernel's own RESULT lines): on `p2_space` the exact door ledger is
   cond = 0.5249 / evap+rad = 14.185 over the whole run, while the
   graded ⟨I⟩ = 0.36/frame is ~40× the door — hop noise counted as
   conversion, faking I/O balance and simultaneity. `hunt/convert.go`
   uses the exact `promo_i/promo_o` columns only on snaps where a
   motif cell is nonzero; zero frames (atom grain) fall back to the
   heuristic. The 10× Φ shortfall is real and **worsens** on exact
   meters; the exact books say the leads arms were draining blobs
   (I ≪ O), not near-balanced cycles.
2. **The Φ/Em = 0.015 bar is NOT the g ≳ 0.3 statement.** H1's Φ is a
   per-snap integral (Δt = 4 t.u. on the leads) and the adopted law
   carries `comp = 1 + q_detune·x`. A small uniformly-converting
   motif with I ≈ O at x ≳ 0.37 passes the frame bar at g = 0
   (`0.2·x³·(1+1.2x) > 0.015`). What needs the stimulated cross term
   g ≳ 0.3 (qgraph Thm 3.5) is a **restoring hold above ambient — an
   H2 statement**, and a law change not to be smuggled in for H1.
3. **p2_ring failed three ways, not one:** dark vacuum (no field →
   cond = 0), **atom grain** (qcnvF = 0.3896 leftover credit —
   radiance never accumulated one atom at x = 0.325), and loop
   topology (Jdip/Circ = 0.02). "Only dark" was a slogan.
4. Archive EMPTY stands; sep 2/106 is a real geometric fact; the
   leads' sep+Es co-occurrence under `promo=2` is a real lead. The
   J clause remains graded on snapshot J (H0's own SKIP); dir-stable
   0.54–0.69 is a mix of meter slosh and real per-event re-picking.

## 1. Step 1′ — honest meters, then re-grade (offline `hunt/` only)

- `hunt/convert.go`: when the FCS schema carries `promo_i/promo_o`,
  use them **always, zeros included** (schema-level flag, not the
  per-frame opportunistic scan). Heuristic remains the fallback for
  archive streams without the columns.
- Report Φ/Em **also as a rate** (÷ mean frame Δt) in the clause
  detail, so the bar's snap-dependence is visible. **The 0.015 bar
  itself does not move** (user-level decision, deferred).
- Report both ⟨min(I,O)⟩ and min(⟨I⟩,⟨O⟩) + the I/O ratio, so mean
  imbalance and asynchrony stop being conflated.
- New H0 class `exact` on `runs/leads/p2_space.fcs` (+ `p0_ctl` as
  the no-columns control): the exact-graded totals must reproduce
  the kernel's own ledger (cond = 0.5249, evap+rad = 14.185); the
  heuristic totals must show the documented inflation. Add a Go test
  with those constants so the meter cannot regress.
- Re-run the standing H0 suite (archive streams have no promo
  columns → must be unchanged).
- Re-grade all ten leads arms; write `runs/regrade/H1_exact.out`.

**Pre-registered predictions (record hits/misses in REGRADE.md):**
exact Φ/Em drops (≲1e-3 on p2_space-class arms); a screaming I/O
imbalance appears (I ≪ O); sep/polar on promo>0 arms survives (it is
door-geometry, not hop noise). If sep dies on the exact ledger, the
sep clause was hop-shaped too and LEADS.md gets a dated correction.

## 2. Step 2 — the registered A.5 arm, unbundled (config-only)

Base (leads COMMON): `L=16 dt=0.02 T=240 amp_nat=2 amp_logate=1
hop_order=1 seed=20260802 snap_every=200`.

- **2a** `rg2a`: `exp=ring ring_n=6 seed_mw=2 kx=0 bath=1 wf_on=1
  k_rad=0.05 promo=2`. Nothing else. One seed. Grade `-h1` and
  `-tag` (the painted cycle, so the foam stays out of the
  denominator). Binary questions: does the door fire (cond > 0,
  rad > 0)? does the motif stay C6?
  *Prediction:* door likely still dark or grain-blocked (A.5
  measured cond = 0 in this bath; the ring sits at x = 0.325 where
  one radiance atom never accrues) — that outcome names the blocker
  as apparatus (field supply / grain), not g.
- **2b** `rg2b`: 2a + `amp_door=1` (the A.5 "hoped-for metabolic
  refill", recorded untested where it matters). Still no rad_dual,
  no SW, no bath_xb, no helix.
- **2a-lit** `rg2c`: if 2a/2b are dark — light the field honestly:
  2b + `noise_amp=0.4` (field in the foam; workfn keeps true vacuum
  dark; **no `bath_xb`** — that is an energy hose). `-tag` grading
  keeps foam box-fill out of the motif's books.
- **2c-helix** `rg2d`: only if doors are live and Φ behaves but the
  J clause fails as a loop — add the kx tilt **alone** on top of the
  best-of-2a–2c. Geometry-pin before any latch: if a helix with live
  doors still sloshes through 180°, the compass is the bug; if it
  holds, no memory (latch face) was needed.

Each arm: check drift at the FP floor; one knob delta per arm;
results per-arm into `A5_RING.md` **before** the next arm runs.

## 3. Step 3 — theory (Go/symbolic, `v97/theory/`, no kernel)

- **3a′ (corrected Φ theorem, replaces F2-as-written):** Under V3a
  (`comp = 1+q_detune·x`), a uniformly converting motif with I = O
  and snap Δt has `Φ_frame/Em = k·x^{p−1}·Δt·(1+q_detune·x)`;
  exceeds 0.015 for x ≳ 0.37 at defaults ⇒ H1's Φ clause is
  passable at g = 0. Theorem 3.4 still forbids a *restoring* root
  above xb at g = 0 ⇒ g ≳ 0.29 is the **H2** threshold
  (u* = 0.669/0.62, Thm 3.5). Write with a Go check.
- **3b (edge-door selection rule):** M_ab for a door supported on
  one oriented edge of C_N vs the vertex door (qgraph L5 measured:
  flat over Δm, mean 0, sd 3.49 on C12). *Pre-registered
  prediction:* edge support concentrates transfer at Δm = ±1
  (dipole selection rule). This pre-registers the falsifiable bar
  for any future slot-door face.
- **3c (closure-rule graph-selection toy):** deferred — P3's claim,
  not the H1 blocker.

## 4. Step 4 — faces (NOT opened here; recorded so nobody reinvents)

Only if 2a–2c produce live I and O that will not pair in time or
will not pin a direction, the candidate faces, each needing its own
brief + explicit user sign-off:

- **Intake pairing** (grok: belongs on the charter's face list as
  its own line — the leads showed orientation alone buys sep, not
  matched intake).
- **Door-locus credit:** a per-link running integral of *signed*
  promotion with FLOW-τ decay (the vector credit v95 declined —
  the objection flips at the door locus) — preferred over a one-bit
  latch (which can lock a drain). A slow-compass EMA of `cellDir`
  is the minimal variant.
- **`pJ` integrated hop-vector accumulator** (physics-inert `fab/`
  write, needed to retire the snapshot-J SKIP): deferred until
  sign-off; must be H0'd against e3b `fd` before it grades H1.
- **Composed helix seed** — **AUTHORIZED and done** (user,
  2026-08-13: "kernel edit accepted"). `fab/run.go` seed-order fix +
  the promo-header parity repair (pre-existing 6-bar break from the
  leads commit, proven by stash-rerun). Battery **ALL GREEN 87**
  after both (`REGRADE.md` §4b). The composed-helix arm is rg2f.
- The g-term (stimulated intake, qgraph L4/Thm 3.5): an **H2/law
  face**, not an H1 patch.

## 5. Self-deception watchlist (from the consult, kept in force)

- Exact meters make Φ *smaller*; any "sharpening" that raises Φ is a
  soften in disguise.
- No `bath_xb` to inhabit H1 (energy hose). No seed-SW sold as a
  maintained organized field (`sw_*` is a t=0 seed; it dies).
- No grading J on new accumulators without an e3b H0 first.
- No raising `EmFloor` / halo-cropping until Φ/Em passes on honest
  books; report floor-sensitivity as a diagnostic only.
- One knob per arm; a PASS on a six-knob stack is unattributable.
- H1 stays H1; the g-term stays H2; 0.015 and 0.80 do not move.

## 5b. Post-consult queue (user, 2026-08-13 evening: "Execute all,
but increase the size of the simulation, decrease the timestep, in
any re-grade or re-simulation")

Grok's state-check queue (`consult/REVIEW_grok_state.md`), adopted:

0. Dest-only book re-grade of rg2c/f/e — offline `hunt/` (`-dest`:
   replayed promoPick picks, not a neighbourhood halo). Bands
   pre-registered: rg2c dest-only Φ/Em ≥ 0.015 ⇒ apparatus claim
   survives; 0.000–0.005 ⇒ there was no Φ pass.
1. Three seeds of rg2c and rg2f.
2. Quench park probe (mild steady drive `noise_amp=0.1`, long T —
   no mid-run noise gate exists in config; the steady-mild variant
   answers the same question: park at/below xb with I≈O at stable
   ambient?).
3. `pJ` integrated hop-vector (physics-inert `fab/` meter knob,
   default 0, battery-gated) + H0 against a fresh e3b-style `fd`
   run.
4. promo=3 on the lit tagged ring (DOOR_RULE's phase-gate escape),
   one seed.

**New-run convention (user directive):** `L=24 dt=0.01
snap_every=400` (frame cadence stays 4 t.u. so the Φ frame bar is
comparable). Existing-stream re-grades keep their recorded L/dt.
No face briefs until 0–2 are read.

**EXECUTED — all five items. Verdict: `REPLICATE.md`.** Item 0 hit
the kill band (no Φ pass ever, on honest books); item 1: 0/6 seeds;
item 2: mild drive starves, strong drive floods, no balance seen;
item 3: pJ H0 GREEN (snapshot J can invert sign; J clause now grades
the integrated channel; battery ALL GREEN 87 after the edit);
item 4: promo=3 closed. rf_pj gives the first trustworthy J verdict:
net directed current, direction unpinned — the steer face's honest
exhibit.

## 6. Status

| item | state |
|---|---|
| Plan written | this file |
| 1′ meters + H0 `exact` + re-grade | **done** — `REGRADE.md`, `runs/regrade/` (leads 0/10; archive re-check 0/259; H0 10× PASS; hunt+calc tests green) |
| 2a/2b/2a-lit/2c arms | **done** — `A5_RING.md`, `runs/a5/` (7 arms + T=1000; first Φ-clause PASS of the programme on the lit tagged ring; J-pin and restoring-park localized as the two missing mechanisms; composed-helix seed fix + parity repair, battery ALL GREEN 87) |
| 3a′/3b theory | **done** — `theory/PHI_BAR.md`, `theory/DOOR_RULE.md`, `theory/calc/` (3b prediction REFUTED: edge door is anti-selective; selection rule = DFT of door profile) |
| Faces | not opened |
