# v103 gate g1 — decision (can a medium with no container hold three dimensions?)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing.

- Tool: `cas/geo3d.c` (C). T0 passes; see `cas/out/t0.out`.
- Grid: `cas/run_g1.sh`, 40 runs (8 sets × L = 6, 8, 10, 12, 16). It
  ran 20:25–22:09 with no errors.
- Tables: `cas/tab.go` (Go), output in `cas/out/tab.out`.
- Post-hoc readings, labelled as such: `cas/posthoc.go`, output in
  `cas/out/posthoc.out`.

**Frame (R14).** The medium is bare quanta and their contacts. There
is no container, and the places are the quanta. A gap is a missing
quantum that the medium closes over, read from the medium's own shape.

## Ending: **G — the gate passes**, in six sets: a, ac, b, bc, ab, abc. All are crystals, and all have contact-range defects

| set | holds 3-D | gap readable and persistent | phase | defect footprint | excess readable | assembles from the melt |
|---|---|---|---|---|---|---|
| none | ✘ (an expander) | ✘ | disordered | — | ✘ | ✘ |
| c | ✘ (an expander) | ✘ | disordered | — | ✘ | ✘ |
| a, ac | ✔ | ✔ | crystal | contact | ✔ | ✘ |
| b, bc | ✔ | ✔ | crystal | contact | ✔ | ✘ |
| ab, abc | ✔ | ✔ | crystal | contact | ✔ | ✘ |

**In one sentence:** with no container, a medium of bare quanta holds
three dimensions only by being rigid. Either flat squares (a) or the
octahedral link (b) is enough, and with neither (none, c) it loses
dimension entirely.

**In three dimensions a missing quantum is a readable thing, as it
cannot be in one.** The medium closes over it and is left with one
compact defect, frozen where it was made. Its "less where" reaches
across the whole closed medium.

## 1. The medium (M, from the flat start; L = 12 | L = 16)

| set | cheapest re-forming of the flat medium | accepted flips | ⟨d⟩ | λ2 | d_H (12→16) | d_s (12→16) | ladder fit d_H, d_s |
|---|---|---|---|---|---|---|---|
| none, c | 0 (every flip free) | 99 % | 4.52 \| 5.06 | 1.54 \| 1.54 | 7.6 | 400 (none), 173 (c) | 6.5, 60–81 |
| a, ac | **72t** | **0** | 9.005 \| 12.003 | 0.268 \| 0.152 | **3.00** | **3.05** | 3.01, 3.12 |
| b, bc | **144t** | **0** | same | same | 3.00 | 3.05 | 3.01, 3.12 |
| ab, abc | **216t** | **0** | same | same | 3.00 | 3.05 | 3.01, 3.12 |

- **No re-forming is free in any set with a or b.** At L = 16 there
  are 614,400 flip choices in the flat medium, counting each flip from
  both ends. Every one costs at least 72t, at a temperature of t, and
  not one was accepted in 4000 sweeps. The
  medium holds three dimensions because it is a frozen crystal
  (ψ_T = ψ_O = 1).
- **With nothing holding it, the medium melts into an expander.**
  ⟨d⟩ grows like ln N (4.5 → 5.1 from N = 1728 to 4096) and λ2 stays
  at 1.54. This is the known result that free local re-forming makes
  expanders, now measured on this medium.
- c alone (no triangles) holds nothing.

## 2. Assembly (A: melted, then cooled to t; reported)

- **The melt is a random medium:** d_H ≈ 7.6 in every set.
- **Nothing assembles.** After cooling from 64t to t, every set is a
  glass: ψ_T 0.11–0.18, ψ_O = 0, d_H 7.6–8.1, accepting 0.3–0.4 % of
  flips at t.
- **No 2-D triangular patches formed** without c. Triangles per
  quantum stay ≤ 0.012, against my expectation (§5).
- **So holding is not forming.** These rules keep a 3-D medium that is
  already there, and do not make one.

## 3. One gap (G) and one excess (X), L = 12 | L = 16

| set | gap: cost | gap: mismatched contacts, spread | excess: cost | excess: mismatched, spread | gap δ⟨d⟩ | excess δ⟨d⟩ | square deficit, gap \| excess |
|---|---|---|---|---|---|---|---|
| a | 144t (relaxes to 120t at L = 16) | 47, 5 | 76t | 34, 4 | −0.052 \| −0.046 | −0.022 \| −0.019 | +28 \| +28 |
| ac | 144t | 27, 3 | 80t | 34, 4 | −0.046 \| −0.038 | −0.022 \| −0.019 | +36 \| +28 |
| b | 288t (→ 264t) | 39, 5 | 192t | 30, 4 | −0.051 \| −0.044 | −0.015 \| −0.013 | +28 \| +36 |
| bc | 288t | 27, 3 | 198t | 30, 4 | −0.046 \| −0.038 | −0.015 \| −0.013 | +20/+36 \| +36 |
| ab | 432t (→ 284t) | 34, 4 | 288t | 34, 4 | −0.052 \| −0.042 | −0.019 \| −0.015 | +28/+20 \| +31/+33 |
| abc | 432t | 27, 3 | 292t | 34, 4 | −0.046 \| −0.038 | −0.021 \| −0.018 | +36 \| +29/+30 |

(Contact counts and spreads are at L = 16; the spread is in contacts.)

- **The gap is readable and persistent in all six sets.**
  - In every sample at both sizes, its mismatched contacts form one
    cluster, with a spread of 3–5 contacts.
  - The quiet medium has no mismatch at all, so the cluster is the
    gap and nothing else.
  - Its centre is read from the state alone.
- **It is frozen where it was made.**
  - Where c is off, a few downhill flips relax the closure right after
    it forms (a: 144 → 120t; ab: 432 → 284t).
  - After that nothing moves, because a gap in a rigid crystal has no
    classical way to travel.
- **The excess is readable and persistent too, and cheaper than the
  gap** (76–292t against 144–432t).
  - In ab and abc it keeps flickering: about one flip in 10⁵ toggles
    back and forth beside it.
  - It stays one cluster with a spread of 4.
- **Defect footprint: contact.** The mismatched contacts reach at most
  2 from the centre (r_def 0–2) at every size.

**Signs.**
- **Squares:** both the gap and the excess break squares (square
  deficit +20 to +36), as expected.
- **Where:** both make the medium smaller. The gap's δ⟨d⟩ is −0.04 to
  −0.05, and the excess's −0.013 to −0.022. **So the medium's "where"
  answers to the gap and to the excess with the same sign.** The size
  differs: the gap counts 2–3.4 times as much.
- **This is against my registered expectation** that the excess would
  make more where.
- **My reading (an interpretation, not a measurement).**
  - The perfect crystal is the most spread-out arrangement of a
    six-contact medium (⟨d⟩ 12.0, against 5.1 for the random medium).
  - Any defect of either sign adds shortcuts.
  - So "where" is greatest in the perfect medium, and anything that
    is not the perfect medium has less of it.

## 4. The reach of a gap's "less where"

Δm(r) is how much closer, on average, the rest of the medium is to a
quantum at distance r from the gap.

- **Registered reading.** Δm(r) is negative at every distance and
  decays across the medium.
  - Over the registered window 2 ≤ r ≤ L/2, the fit Δm ∝ r^−p gives
    p = 0.85–1.19 in every passing set.
  - That window ends at r = 8 at L = 16. Further out it keeps falling
    to the bookkeeping floor, about −0.003: the shift from having one
    quantum fewer.
- **But it is not a power law** (post-hoc, `cas/out/posthoc.out`).
  - **Δm collapses on r/L.** In ac, the cleanest case (its closure
    never relaxes), Δm at r/L = 1/4 is −0.104, −0.108, −0.106 and
    −0.104 at L = 8, 10, 12 and 16. At r/L = 1/2 it is −0.042, −0.041,
    −0.039 and −0.038.
  - The other sets collapse within about 20 %. Their closures relax
    differently at each size.
- **So the gap's "less where" has no decay of its own: its reach is
  set by the size of the closed medium.** The registered p ≈ 1 is the
  shape of a finite torus, not a 1/r law.
- **Why, by the medium's geometry.** Counted in contacts, a crystal's
  distance is taxicab-like, and very many paths tie for shortest. So
  the shortcuts where the medium closed over the gap shorten a fixed
  share of all paths, however far away they start. That was the
  mechanism behind my registered expectation p ≈ 0. Over the finite
  window the number came out ≈ 1. The collapse shows the expectation
  was right about the mechanism and wrong about what the window would
  read.
- **The mean distance moves accordingly.** A single gap shortens the
  whole medium's mean distance by 0.09 at L = 6 and 0.046 at L = 16
  (set a). The share it touches thins slowly as the medium grows.
- **Field footprint ΔG(r): reported, not read.** It compares a source
  spread over the gap's centre (18–28 quanta) with a one-quantum
  reference. So it mixes source shape with any change in the medium.
  It falls from −0.08 to −0.19 at the centre to about 10⁻⁴ at the far
  side. A matched reference would be needed to read it.

## 5. What I expected (PREREG §7, written before running), against what happened

| expected | outcome |
|---|---|
| none and c melt into an expander and fail H | **held** (d_H 7.6, d_s 170–400) |
| every set with a or b is rigid ("about 110t for a straight flip under a") and holds as a crystal | **held.** The cheapest flip is 72t (a), 144t (b) and 216t (ab), with zero flips accepted. The 110t estimate was for straight flips; bent ones are cheaper |
| nothing assembles; triangle patches without c | **half held.** Nothing assembles, but no triangle patches formed (≤ 0.012 per quantum) |
| the gap is readable, persistent and frozen, contact defect, less where | **held** |
| the excess makes more where; the square deficit has the same sign for both | **half held.** The square deficit has the same sign, but **the excess also makes less where** |
| distance footprint p ≈ 0; field footprint degeneracy-free | **not as registered.** p reads ≈ 1 over the window. Post-hoc, the profile collapses on r/L: no decay of its own, as the mechanism predicted. ΔG cannot be read without a matched reference |
| ending G (crystal, contact) for every set with a or b | **held** |

## 6. Checks

- **T0 passes.**
  - Local energy changes equal full recomputation exactly (61,279
    flips).
  - The flat medium reads its exact ⟨d⟩, diameter and λ2.
  - The instruments read 2.00 on a 2-D medium and far from 3 on a
    random one.
  - A zero-energy medium of L³ − 1 quanta exists (a sheared
    circulant).
  - **One gate passed with little margin:** T0.3a, detailed balance at
    N = 8, gave TV 0.0095 against a bar of 0.01. Its N = 9 companion
    passed with ample margin (max 0.0032).
- **All 40 runs** completed, every CHECK passed (six contacts, simple,
  one piece), every Lanczos run converged, and every .err file is
  empty.
- **Green's-function NaN, disclosed.**
  - In the melted sets (none and c at L = 10, 12, 16), the whole
    medium is mismatched. The gap's centre is then every quantum, the
    source is uniform, and the solve returns NaN (12 solves).
  - This happens only where the gap is already unreadable (P1 fails),
    so it decides nothing.
- **One presentation fix in `tab.go` during the grid, before the final
  tables.** Where the gap is not readable, the footprint class now
  shows "—". No bar changed.
- **The proposal and acceptance** are those tested in T0.3. The flips
  keep the medium in one piece, and it never split.

## 7. What the gate says

1. **A medium with no container can hold three dimensions, but here
   only as a frozen crystal.**
   - One local rule is enough: flat squares, or an octahedral link.
   - The holding is rigidity: the cheapest re-forming costs 72t or
     more at temperature t.
   - Free re-forming loses dimension entirely.
   - Holding is not forming: from a melt, nothing assembles.
2. **In three dimensions, "something with gaps" works literally.**
   - With nothing but quanta and contacts, a missing quantum is a
     readable, compact thing: one defect cluster, read from the
     medium's shape, with no label and no slot.
   - One dimension could not show this (c2 needed pieces with
     content).
3. **A gap makes less where, and so does an excess.** In this medium
   the response of "where" to a defect does not carry the defect's
   sign. Only its size differs (the gap 2–3.4 × the excess).
4. **The less-where of one gap reaches the whole closed medium.** Its
   profile scales with the medium's size and does not decay on its
   own.
   - This is a property of distance counted in contacts in a crystal.
     It is not yet a 1/r law, and not a force.
   - Whether two gaps attract through it is c3's question.
5. **Every passing set is a crystal.** By Kleinert and Zaanen, a
   crystal's forces between curvature sources confine, and 1/R needs
   dislocations to move. And in a frozen crystal a gap has no
   classical way to travel, which matters for lane a.

## Limits

- **Five sizes, up to 4096 quanta.** The collapse on r/L spans
  L = 8–16 only.
- **The random process only.** The geometry is classical, and the gaps
  cannot move classically.
- **Temperature t and 2t per unit deviation** are c2's scale, not
  scanned. The rigid sets are far below any melting point at that
  scale.
- **The excess is made by a registered rule** (three contacts within
  distance 2 of a quantum). Other insertions are not tested.
- **The field footprint has no matched reference.**

## Then

- **c3 (3-D), two gaps in the frozen crystal (abc,** the strictest
  set, where c also excludes the 2-D zero).
  - Their energy and their shared less-where against separation, gap
    with gap, excess with excess, and gap with excess.
  - This is the sign test of §3, now as an interaction.
- **a2 (3-D), the gap as a quantum particle on the frozen crystal.**
  - It cannot move classically. Does coherent re-forming (tunnelling
    between closures) give it a band?
  - Then spin, exchange and binding in three dimensions.
- **g2, a medium that keeps orientation but lets dislocations move**
  (the Kleinert–Zaanen nematic). The closest candidate here is b
  without a.
  - It asks whether three dimensions can be held while defects stay
    mobile, the condition for a 1/R force.
- **A matched-reference field footprint,** if c3 needs it.
