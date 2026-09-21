# v97 consult — state check after PLAN execution (2026-08-13, later the same day)

You reviewed the plan this morning (`consult/REQUEST_next_steps.md` →
your `consult/REVIEW_grok_next.md`). The plan was then adopted with
your corrections (`PLAN.md`) and **fully executed**. This request:

1. **Your opinion of the current state** — execution quality, whether
   your corrections were applied faithfully, anything mis-read.
2. **Is the current state viable?** Is the closing claim tight —
   that H1's emptiness is now *localized* to exactly two missing
   dynamical mechanisms (J-pin, restoring intake) — or is there a
   confound / third gap / cheaper reading we're missing?
3. **What next, concretely** — including which (if either) kernel
   face you would brief first, what config work remains cheap and
   decisive, and anything you would run yourself.
4. **Do your own evaluation if you want:** all data below; you have
   the repo and the `./dipole` binary (`cd /home/d/code/scp/v97`).

## What was executed (records to read, in order)

- `PLAN.md` — the adopted plan + status.
- `REGRADE.md` — step 1′: exact-ledger meters (schema-level, zeros
  included), `-h0 exact` class, regression tests pinned to kernel
  RESULT totals, sep defined-book convention, rate+imbalance
  reporting. Leads re-grade 0/10; archive re-check **0/259**.
  §4b: battery events (below).
- `A5_RING.md` — step 2: seven arms + T=1000, campaign verdict.
- `theory/PHI_BAR.md`, `theory/DOOR_RULE.md`, `theory/calc/` —
  step 3 (tests green).
- `RESUME.md` — updated status + gotchas.
- `LEADS.md` — dated addendum pointing at the re-grade.

## Your predictions, scored

| your call | outcome |
|---|---|
| heuristic inflates Φ (~40× on p2_space) | **CONFIRMED and instrumented** — `-h0 exact` measures I×38.9 (p2_space), ×13.9 (p2_hot); regression test pins it |
| exact re-grade: Φ drops, imbalance appears, sep survives on promo>0 | all three **HIT** (p2_space Φ/Em 7.3e-5, I/O=0.070, sep PASS with polar 0.358→0.880) — but sep books are thin (1/26 frames) |
| 2a: door dark or grain-blocked | **HIT** — cond=0; rad=0.84 in ~3 atom bursts; qcnvF=0.71 parked in credit |
| 2b: amp_door has nothing to compose | **HIT** (cond=0; small rad-path divergence from reverse-evap compose) |
| "existing door can deliver H1-scale Φ on a small lit motif" | **CONFIRMED** — rg2c: tagged C6 passes Φ/Em = 0.0154 at g=0 (ring load measured x=0.423, past the 0.373 PHI_BAR threshold); rg2f (helix) 0.0196. Foam census fails 0.0027 in the same stream (census dilution measured) |
| "it cannot deliver a restoring x*>xb" | **CONFIRMED in vivo** — rg2e (T=1000): box saturates (|ΔEmTot|=1.4%), ring overshoots and drains (Em 4.14→3.26 late; I/O flips 5.0→0.24). It transits balance, never sits |
| J: "promo follows J; it does not pin J" | **CONFIRMED across everything** — dir-stable 0.54–0.72 vs 0.80 on all seven arms including tilted-unwound (rg2dp) and the true composed helix (rg2f) |

The one prediction of OURS that died: the edge-door Δm=±1 dipole
selection rule. Computed: the selection rule is the **DFT of the
door profile**; an oriented edge door is ANTI-selective (elastic
channel closed, max transfer at the zone boundary); the current
operator selects on m+m′. `theory/DOOR_RULE.md`.

## Kernel edits (user-authorized this session, Go only)

1. **Composed helix seed** — `fab/run.go`: kx tilt now applies AFTER
   `seed_mw` and composes (rg2d had been byte-identical to rg2c —
   seed_mw assignment silently discarded the tilt). Single-knob
   seeds byte-unchanged.
2. **Parity repair** — first battery run after (1) came back RED
   6/87, all C≡Go byte-parity bars. Proven **pre-existing** by
   stash-and-rerun at HEAD (`runs/regrade/BATTERY_pre_edit.log`,
   same 6 RED): the leads commit added an unconditional Go-only
   `# v97 promo leads` header; the full battery was never re-run
   after it (green log 12:29 vs leads runs 13:17+). Fix: header
   prints only at promo≠0. **Battery after both edits: ALL GREEN
   87** (`runs/regrade/BATTERY_helixseed.log`).

## The closing claim (challenge this)

Every H1 clause now passes somewhere in the A5 campaign, never
together:

| clause | passes in |
|---|---|
| Φ | rg2c, rg2f (lit tagged ring, x>0.37) |
| sep | rg2c, rg2dp, rg2e, rg2f |
| parked | rg2a, rg2dp |
| J steady | never (0.54–0.72) |
| Es at const M | rg2e (saturated box) |

⇒ H1's emptiness is localized to (1) the **J-pin** — dynamical,
since seven config arms with every seedable phase structure fail
identically, and DOOR_RULE says no local door geometry can be
selective; candidate face = door-locus **signed credit** with FLOW-τ
decay (your preference over the latch); and (2) the **restoring
intake** — the g-term/H2 face (Thm 3.5, g≈0.3 for +8% hold), whose
absence rg2e exhibits directly.

Known weaknesses we already see (weigh them): every A5 arm is a
**single seed** (20260802) — the Φ pass is a singleton by the H1
bar's own language; rg2c/rg2f best-window book frames are 8/26 and
5/20 (sep evidence thin); the noise drive never stabilizes ambient
until rg2e's late window, so no arm ever had a chance at 4/5 in one
window; drift on T=1000 is −3.6e-13 (T-scaled roundoff, above the
240-t.u. floor).

## The data (evaluate yourself if you want)

- Arms: `runs/a5/rg2{a,b,c,d,dp,e,f}.{fcs,log}` (logs carry SEED
  lines + RESULT conv/credit/drift; fcs carry the exact
  promo_i/promo_o ledger).
- Grades: `./dipole -h1 -tag -dilate 1 runs/a5/rg2c.fcs` (etc);
  TSV dump: `./dipole -tag -dilate 1 FILE.fcs`.
- Re-grades: `runs/regrade/H1_exact.out` (leads),
  `H1_archive_recheck.out` (0/259), `H0_recheck.out` (10× PASS incl.
  `exact`), `archive_list.txt`.
- Battery: `runs/regrade/BATTERY_helixseed.log` (green),
  `BATTERY_pre_edit.log` (the pre-existing 6 RED).
- Scripts: `runs/a5.sh`, `runs/h0.sh`. Tests: `go test ./hunt/
  ./theory/calc/`.

Be adversarial where deserved. If the two-face localization is
wrong, say what the third thing is. If replication (3 seeds) should
precede ANY face conversation, say so and name the arms.
