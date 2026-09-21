# v97 promo leads — oriented door (2026-08-12)

**Question:** can a Go-only, grainy, default-off change make the H1
set findable?

**Answer:** it changed the *structure* — the two clauses that never
co-occurred in the archive now co-occur — but **H1 is still EMPTY**.
J does not lock, and Φ/Em stays ~10× below the bar.

Knob: `promo` (0 = vertex door, byte-inert). Go only. C ignores it.
Conservation held at the FP floor on every arm (|drift| ≤ 1.7e-13).

| promo | meaning |
|---|---|
| 0 | control (vertex cond / evap / rad, as before) |
| 1 | deposit Em / inject field on the J-oriented neighbour |
| 2 | also move the Es pull/return on that neighbour (space-leg) |
| 3 | (2) + cycle-closure stall (`\|cos Δθ\|^p` < 0.05 → fall back) |

Acceptance is still the five H1 clauses, unsoftened. Scoreboard:
`runs/leads/H1.out`.

---

## Arms (T=240, L=16, `amp_nat=2 amp_logate=1`, seed 20260802)

| arm | extra | H1 | what moved |
|---|---|---|---|
| **p0_ctl** | kx=1.1 V3a | 2/5 | polar 0.14 — **no split** (vertex) |
| **p1_slot** | same, promo=1 | 2/5 | **sep+polar PASS** (0.29); Es tracks dM |
| **p2_space** | promo=2 | **3/5** | sep+polar + Es-at-constant-M (the archive never had this pair) |
| **p3_close** | promo=3 | 2/5 | sep holds; closure starved cond (0.17 vs 0.53) |
| **p2_hot** | amp=2 wf=0 promo=2 | **3/5** | polar **0.60**; Es cycle PASS; Φ/Em 8e-4; J slosh |
| **p2_noise** | noise=0.4 promo=2 | **3/5** | same pair; Em climbed to 600 (not a particle) |
| **p2_ring** | C6 vacuum | 1/5 | door never fired; still a loop |
| **p2_dual** | rad_dual=1 promo=2 | **3/5** | same pair as p2_space; dual did not lift Φ |
| **p2_jm** | j_metric=1 promo=2 | 2/5 | motif melted (Em 100→63) |
| **p2_tight** | tight+hard tax | 0/5 | drained (Em 18→5) |

Drift: all arms ≤ 1.7e-13.

---

## What the leads proved

1. **Orientation is real.** Control polar = 0.14 (FAIL). Every promo>0
   blob with door traffic passed I/O spatial split. That clause was
   2/106 in the archive.

2. **The space-leg is real.** `promo=2` is what makes Es throughput
   at constant M pass *together with* the split. `promo=1` split I/O
   but Es still tracked dM — exactly the `i1` miss. The warning was
   right: slot-smear without Es on the same atom is 3/5 forever.

3. **H1 still empty, for two named reasons.**
   - **Φ/Em ~ 0.001–0.004** vs bar 0.015. I and O rarely overlap in
     the same snap (⟨I⟩ and ⟨O⟩ can be decent while ⟨min(I,O)⟩ is
     small). Dual intake did not fix it.
   - **Jdip dir-stable 0.54–0.69** vs bar 0.80. Promo follows J; it
     does not pin J. j_metric melted rather than locked.

4. **Closure (promo=3)** was too tight for this apparatus: it cut
   cond without buying a lock.

5. **Noise/hot** raise raw throughput but add mass or slosh; they
   do not mint a species.

---

## What this does *not* authorize

H1 is still EMPTY. Do not adopt `promo` into a law table. Do not
soften Φ/Em or dir-stable to pass. The next design, if any, has to
**pin a fore/aft pair** (so J stops sloshing and I,O fire as a
matched cycle), not add more isotropic traffic.

A candidate that was *not* tried here: a slot that remembers its
last successful promotion (a one-bit “this link is the machine”)
instead of re-picking max J every beat. That is a different face.

---

**Addendum (2026-08-13, exact-ledger re-grade):** the ⟨I⟩/⟨O⟩ and
sep/polar numbers above were graded by the mode-crossing heuristic,
which the exact `promo_i/promo_o` ledger shows inflates I by ~14–39×
(hops counted as conversion). On honest books H1 is still EMPTY
(0/10), the split *narrows* to arms where the door fired both ways
(and strengthens there: p2_space polar 0.88), and the Φ failure is
re-diagnosed as near-dead one-sided intake, not pairing noise.
Record: `REGRADE.md`, `runs/regrade/H1_exact.out`.
