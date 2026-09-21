# REPLICATE — PLAN §5b queue results (2026-08-13 night)

All new runs at the user-directed scale: `L=24 dt=0.01
snap_every=400` (frame cadence 4 t.u. preserved). Grading:
`-tag -dest` (dest-only book scope, `hunt/dest.go`; `-dilate`
demoted to diagnostic after item 0 hit the kill band). Streams:
`runs/rep/`, `runs/pj/`. Grades: `runs/rep/{RC,RF}_GRADES.out`,
`runs/regrade/DEST_regrade.out`.

## Item 0 — dest-only re-grade of the old arms: KILL band hit

rg2c Φ/Em = **0.0012** (pre-registered "no Φ pass" band 0.000–0.005),
rg2f 0.0057, rg2e late 0. The L=16 "first Φ-clause PASS" was carried
by the dilate halo (`REGRADE.md` §6). Retracted.

## Item 1 — replication, 3 seeds × {lit ring rc, helix rf}: 0/6 pass Φ

| arm | cond (whole run) | Φ/Em dest | I/O | book frames | clauses |
|---|---|---|---|---|---|
| rc_s20260802 | 35.7 | 0 | 12.4 | 0/26 | 1/5 (parked) |
| rc_s111 | 139.8 | 0.0028 | 10.9 | 1/23 | 0/5 |
| rc_s271828 | 80.1 | 0.0071 | 4.3 | 3/20 | 2/5 (sep, parked) |
| rf_s20260802 | (helix) | 0 | 15.8 | 0/26 | 1/5 (parked) |
| rf_s111 | (helix) | 0.0022 | 1.85 | 1/26 | 1/5 (sep) |
| rf_s271828 | (helix) | 0 | 0.54 | 0/22 | 0/5 |

**Grok's pre-registration ("if 2/3 miss, stop saying the apparatus
works") is exceeded: 6/6 miss.** Best Φ/Em = 0.0071, half the bar.
Combined with item 0: **the programme has never had a Φ-clause pass
on honest books.** Additional structure: cond is strongly
seed-sensitive (35.7–139.8 — feeding is contingent, not robust);
intake is one-sided nearly everywhere (I/O 4–16 except two arms);
the door almost never fires both ways within one snap (book frames
0–3 per window).

## Item 4 — promo=3 (phase-gated door): branch CLOSED

p3ring: cond = 91.6 — it did **not** starve at L=24 (the L=16
starvation 0.17-vs-0.53 was a feeding artifact of the small box) —
and it bought nothing: 0/5, Φ/Em 0.0016, sep 0.122, J 0.58.
DOOR_RULE's nonlinear escape hatch neither starves nor selects at
this apparatus. Closed as "no benefit."

## Item 3 — pJ meter: H0 GREEN, and it indicts every prior J number

`pj_meter=1` (Go, byte-inert, battery ALL GREEN 87 after): per-cell
integrated hop current; identity Σcells ≡ `p1fd` holds at the f32
floor (ratio 1.0001–1.0002, 3 seeds). On e3b-style arms the
integrated current is seed-robustly **+x** (316–360) while the
snapshot bilinear reads **−78 to −92 ± ~250 — the wrong sign**. All
dir-stable numbers previously graded (0.54–0.72) were on a meter
that can invert direction. The H1 J clause now grades on the
integrated channel whenever the columns exist (same bars; labeled
fallback otherwise). First honest J read: `rf_pj` (below).

## Item 2 — quench park probe: the ring STARVES at mild drive

`noise_amp=0.1, T=800` (drift −1.6e-13): cond = **0.69 over 800
t.u.** (vs 36–140 at noise 0.4) — the mild ambient sits below the
workfn/grain feeding threshold. The tagged ring compresses (n→3,
x up to 0.56 at t≈100), then drains: Em 4.9 → 1.5 by t=300, **dead
(n = 0) by ~t=400**. Grade: 0/5, "no live motif" in the best window.

So the drive window is measured at both ends: **0.4 floods**
(I/O = 4–16, ring grows, box outruns it), **0.1 starves** (doors
nearly dark, ring dies). No self-balancing point was observed at
either end — the shape Thm 3.4/3.5 predicts at g = 0, now with an
honest exhibit on each side. A rigorous "no park anywhere" claim
would need a 0.2–0.3 drive sweep; that is the natural next config
arm if anyone wants it, and it was not part of this queue.

## Campaign verdict (the full §5b queue, honest meters, both scales)

**H1 is EMPTY on every meter generation, at both scales, 0/N
everywhere.** What the queue nailed down:

1. **Φ has never passed on honest books** (item 0 kill band + 0/6
   replication). The dominant killers, in order: I/O asynchrony
   (0–3 book frames per window — the door fires one-sided per
   snap), one-sided intake (I/O = 4–16), seed-contingent feeding
   (cond 36–140 across seeds of one config), and a narrow drive
   window (starve↔flood with no observed balance).
2. **The J failure is real and now trustworthy** (rf_pj): rings
   carry net directed integrated current (pJnet/Σ|pJ| = 0.71 — not
   loops), but the direction is unpinned (dir-stable 0.58 on the
   H0-validated channel). The steer face ("the door listens to a
   bad clock") has its honest exhibit.
3. **Orientation itself is the one robust positive from promo**:
   polar 0.77–1.00 wherever the door fires, every arm, both scales.
4. **promo=3 closed** (neither starves nor selects at L=24).

Face queue standing (unchanged order, now with honest exhibits;
each needs its own brief + sign-off): (1) **door-locus signed
credit** (steer); (2) **H2 g-term** — wants the 0.2–0.3 drive sweep
first; instrument debt: **source-booked ledgers** (dest replay is a
lower bound). `CHARTEXP.md` (noted, no action): the 3+1 chart
projection is a candidate hunt/ANLZ instrument under the same
meter-only contract; its chart=4 null independently corroborates
DOOR_RULE from the composition side.

## rf_pj — first honest J-clause read (integrated fd channel)

Helix + `pj_meter=1`, graded on the H0-validated integrated channel:

- **⟨|pJnet|⟩/⟨Σ|pJ|⟩ = 0.71** — the not-a-loop sub-clause passes
  easily: the ring's integrated hop current is predominantly a net
  directed flow, not circulation. (The "wrong topology for the J
  clause" worry about rings was a snapshot artifact.)
- **dir-stable = 0.58 on the honest meter** — the direction
  genuinely wanders frame to frame. This is the campaign's first
  *trustworthy* J verdict: the J-pin failure is **real physics**,
  not meter slosh. The steer-face case ("the door listens to a bad
  clock" — `promoPick` re-picks from snapshot J every event, and
  the resulting current direction is unpinned) now has an honest
  exhibit instead of a SKIP-tainted one.

Also 0/5 overall (Φ/Em = 0, I/O = 14, book frames 0/26 — consistent
with the item-1 seeds).

## Standing instrument debt (named, needs sign-off)

Dest-only replay is a **lower bound** (picks replayed from frame-end
state; mid-interval events may book off-set) and dilate an upper
bound. The definitive fix is **source-booked ledgers**: index
`pI/pO` by the door-owning cell (spatial deposit unchanged, box sums
unchanged) — physics-inert, but it changes recorded column
semantics, so it waits for explicit sign-off.
