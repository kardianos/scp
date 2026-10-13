# v103 lane c step c3 — decision (two defects in the medium that holds three dimensions)

**Written:** 2026-10-13, after the registered runs. Gates and bars:
`PREREG.md`. A measurement in models; it rules nothing.

- Tool: `cas/c3pair.c` (C). Compiled clean with `-Wall -Wextra`, no
  change from the registered text. T0 passes; see `cas/out/t0.out`.
- Grid: `cas/run_c3.sh`, 96 processes (8 sets × L = 8, 12, 16 × 4
  placements). It ran 2026-10-12 20:46:54 to 2026-10-13 00:13:12 with
  no errors: every process reached DONE, every .err file is empty, and
  all 2552 CHECK lines passed.
- Tables: `cas/tab.go` (Go), output in `cas/out/tab.out`.
- Post-hoc readings, labelled as such: `cas/posthoc.go` (Go), output in
  `cas/out/posthoc.out`.

**Frame (R14).** The medium is bare quanta and their contacts, g1's
medium unchanged. There is no container, and the places are the
quanta. A gap is a missing quantum that the medium closes over, and an
excess is a quantum the medium takes in. Both are read from the
medium's own shape, with no label. Separation and "where" are counted
in contacts through the medium.

## Ending: **U — unsigned at contact**, in all six scored sets (a, b, ab, ac, bc, abc), with contact reach

| set | GG (L = 12 \| 16) | XX | GX | ending | energy reach R_E (12 \| 16) | at R = 8 |
|---|---|---|---|---|---|---|
| a | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 4 \| 3 | contact |
| b | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 3 \| 3 | contact |
| ab | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 4 \| 3 | contact |
| ac | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 3 \| 3 | contact |
| bc | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 3 \| 3 | contact |
| abc | attract \| attract | attract \| attract | annihilate \| annihilate | **U** | 4 \| 3 | contact |

none and c melt, as in g1, and are not scored.

**In one sentence:** in the frozen three-dimensional medium every pair
of defects attracts at contact, whatever their kinds, and a gap and an
excess at contact undo each other exactly. No pair feels any energy
beyond three or four contacts.

**The second result is not energy:** the two defects' "less where"
overlaps at every separation measured, falling with separation, and it
is unsigned too. Nothing in this medium charges energy for it.

## 1. The medium and the scored sets

- **The six sets with a or b** hold the flat crystal unchanged at every
  size and placement: no contact changed, no mismatch. They are scored.
- **none and c** re-form into g1's expander (mismatch everywhere) and
  are not scored.
  - Their E_int is 0 (none, no energy) or noise within ±4t (c).
  - At L = 16, R = 8 there is no quantum at distance 8 from p1 in an
    expander, so those 8 media were skipped (the 8 SKIP lines in the
    ledger).

## 2. One defect alone (g1's, reproduced)

Medians over placements; energies in t, made → after the quench.

| set | gap, L = 12 \| 16 | excess, L = 12 \| 16 |
|---|---|---|
| a | 144 → 88 \| 144 → 104 | 76 \| 76 |
| b | 288 → 264 \| 288 → 264 | 192 \| 192 |
| ab | 432 → 280 \| 432 → 334 | 288 \| 288 |
| ac | 144 \| 144 | 80 \| 80 |
| bc | 288 \| 288 | 198 \| 198 |
| abc | 432 \| 432 | 292 \| 292 |

- Each single is one compact mismatch cluster.
- The quench lowers a gap by two flips (median) in a, b and ab, and
  never moves an excess.

## 3. The interaction energy at contact (R = 1; window median over placements; t)

E_int = E_pair − E_1 − E_2 + E_0. The bars are ≤ −2t for attract and
≥ +2t for repel.

| set | GG, L = 12 \| 16 | XX | GX |
|---|---|---|---|
| a | −32 \| −65.5 | −32 \| −36 | −164 \| −180, healed 4/4 \| 4/4 |
| b | −248 \| −204 | −124 \| −124 | −456 \| −456, healed 4/4 \| 4/4 |
| ab | −188 \| −212 | −200 \| −200 | −568 \| −622, healed 4/4 \| 4/4 |
| ac | −120 \| −120 | −44 \| −40 | −224 \| −224, healed 4/4 \| 4/4 |
| bc | −240 \| −240 | −121 \| −121 | −486 \| −486, healed 4/4 \| 4/4 |
| abc | −284 \| −360 | −48.25 \| −195.5 | −724 \| −724.5, healed 4/4 \| 4/4 |

- **Every value is far past the bar.** The smallest attraction is
  −32t, sixteen times the 2t bar.
- **The class agrees at L = 12 and 16 in every set and kind,** and at
  L = 8 too.
- **Two like defects at contact share their mismatch.** For example,
  in abc two gaps together cost 360t less than apart, out of 864t;
  two excesses cost 195t less, out of 584t.
- **A gap and an excess at contact annihilate:** in 72 of 72 scored
  media (6 sets × 3 sizes × 4 placements).
  - It happens at the making. The excess's Gibbs pick at t, among the
    triples within 2 of p2, took exactly the three contacts the gap's
    closure had made. The medium was the flat crystal at once
    (E = 0 before any quench), and stayed so through the window.
  - So E_int is minus the two singles' energies (in abc,
    432 + 292 = 724).
  - **The where comes back exactly:** S_P = 0 and B_P = 0 for every
    healed pair.
  - **No GX at R ≥ 2 healed** (0 of 288 scored media). From two
    contacts apart they bind strongly (quenched medians −20 to −304t at
    R = 2) but stay two defects, because the frozen crystal cannot
    bring them together.

## 4. The reach of the energy

- **R_E** (the largest R with any |median E_int| ≥ 2t) is 3 at L = 16
  in every set, and 3 or 4 at L = 12.
- **At R = 4 (L = 16) every quenched median is exactly 0,** and at
  R = 8 every kind is below 2t. The reach is contact.
- **Exact locality:** no pair at R ≥ 6 has E_int^q ≠ 0, in any set or
  size.
- The window medians differ from 0 there only by flicker, at most
  0.25t.
- **The energy has no sign by sublattice.** Opposite-sublattice
  separations (R = 1, 3) and same-sublattice ones (R = 2) all attract.
  - Seven medians are weakly positive (≥ 2t), all at R = 3 or 4 and
    none at L = 16:
    - at R = 3: a GG +2t and a GX +4t (L = 12), b GG +6t (L = 12),
      b GX +10t (L = 8), ab GX +10t (L = 12);
    - at R = 4, the half-period of L = 8: a GG +8t and b GG +2t.

    They sit beside attractions down to −155t at the same R and follow
    no sublattice pattern.
  - The sign that sublattice gave c1's Coulomb on its container has no
    counterpart here: a closure's contacts join quanta of the same
    sublattice, so the medium is not two-coloured near a defect.

## 5. Where

**S_P = ⟨d_P − d_0⟩** is taken over the original quanta present in
all four media.

- **Unsigned.** Every unhealed pair makes less where, for all three
  kinds: S_P < 0 in 44/44 GG, 44/44 XX and 36/36 GX pairs per set at
  L = 12 and 16 (tab §9).
- **Sub-additive.** The mean ρ_W = S_P/(S_1 + S_2) over R ≥ 4 is 0.887
  to 0.941 in every set and kind (tab §9).
  - Two shortcuts do not both shorten a path that can use only one of
    them.
  - ρ_W rises toward 1 with separation: 0.41–0.87 at R = 1, and
    0.92–0.97 at R = 8.
- **The interaction in where reaches past the energy.** I = S_P − S_1 −
  S_2 is positive at every separation measured (the pair makes less
  less-where than its singles), while the energy is exactly zero
  beyond R = 4.
- **Post-hoc (`posthoc.out` §2), the relative interaction
  I/|S_1 + S_2|** falls with separation:

  | kind | R = 2 | R = 4 | R = 8 |
  |---|---|---|---|
  | GG | 0.271 | 0.148 | 0.065 |
  | XX | 0.278 | 0.119 | 0.050 |
  | GX | 0.257 | 0.114 | 0.046 |

  - The values are at L = 16, means over the six sets.
  - That is a factor of 1.8–2.5 per doubling of R.
  - This is three points per kind, the last at the half-period of the
    torus, with no fit registered. It is not offered as a law.
- **No collapse on R/L (registered, tab §8).**
  - At R/L = 1/4 and 1/2, I falls with L (for example GG in a at 1/2:
    0.0114, 0.0072, 0.0049 at L = 8, 12, 16).
  - Post-hoc, the relative interaction at fixed R varies with L by up
    to about 25% (GG at R = 2: 0.216, 0.235, 0.271). At fixed R/L = 1/2
    it varies by about 35% (0.102, 0.081, 0.065).
  - Neither collapse is clean at these sizes.
- **Less space between (B, tab §7).** B is the change in distance
  between the quanta within 2 of p1 and those within 2 of p2.
  - **Every kind at every R ≥ 2:** B_P < 0.
  - **At R = 1:**
    - GG −0.26 to −0.38;
    - GX 0 (healed);
    - XX −0.002 to −0.091, except b at L = 12 and 16 (+0.005 and
      +0.003), essentially zero.
  - **At range the shortening is additive.** At R = 8 (L = 16), B_int
    ≤ 0.031 for every set and kind, while B_P is:
    - GG −0.48 to −0.74, 6–9% of the 8.3–8.5 contacts between;
    - GX −0.28 to −0.51;
    - XX −0.04 to −0.23.
  - **At contact the shortening is strongly sub-additive.** B_int is
    +0.17 to +0.34 for GG and up to +0.41 for GX. It is not small
    there, against expectation.

## 6. What the pairs become

- **Clusters.** At R ≤ 4 nearly every pair is one mismatch cluster,
  and a healed GX is none (tab §10). At R = 8 nearly every pair is two.
- **Persistence.** The cluster count is unchanged from the quench to
  the end of the window in all but five media, all in abc: XX at R = 4
  (L = 8 and 12), 6 and 8 (L = 16), and GX at R = 6 (L = 16).
  - The excess's flicker moves a mismatched quantum across the
    "within 2" merge rule.
  - Their E_int^q is exactly 0, so it is the reading, not a merger.
- **The centre distance D_cc does not read the separation.**
  - A closed-over gap has no centre quantum. Its centre set (the
    quanta of least eccentricity to the cluster) is broad.
  - D_cc is the least distance between two such sets.
  - Two gaps 8 apart read D_cc = 4 in every placement of ac, bc and
    abc, and 4–7 elsewhere. Excess pairs read 7–9.
  - So D_cc measures the instrument's breadth, not "less space
    between". B is the reading of that (§5).

## 7. Persistence in the window (`posthoc.out` §1)

- **Singles are frozen, as in g1:** no contact changes, except the
  excess's flicker in ab and abc (10/12 and 9/12 X1, 56/60 and 33/60
  X2), and two gap media in a (one G1, one G2).
- **Pairs flicker in every set, near contact.** In a, b, ac and bc the
  flicker is at R ≤ 4, where the two closures share mismatch. The only
  exception is three media in a at R ≥ 6 (two GG, one GX).
  - In ab and abc, XX and GX also flicker at R ≥ 6 (9–12 of 12 media
    per kind): the excess's own flicker.
  - At most 12 contacts change in any window.
- **Healed GX never moves** (0 of 72).
- **The flicker never changes a class.**
  - At R = 1 the window and quenched medians differ by at most 20t
    (b XX: −124 against −104).
  - Anywhere, they differ by at most 34.5t (ab XX, L = 8, R = 2: −76.5
    against −42).
  - Every value at contact is −32t or below, sixteen times the bar.

## 8. What I expected (PREREG §7, written before running), against what happened

| expected | outcome |
|---|---|
| none and c melt and are not scored | **held** |
| locality: E_int^q = 0 exactly beyond contact; R_E ≤ 6; contact reach at R = 8 | **held.** R_E 3–4; exactly 0 at R ≥ 6 in every pair |
| GG attract at R = 1; XX attract for the same reason | **held**, in every set and size |
| GX annihilates at R = 1 (likely, not guaranteed) | **held**, 72 of 72, at the making |
| GX attracts by partial overlap at R ≥ 2 | **held** (R = 2: quenched medians −20 to −304t) |
| ending U at contact in every holding set, with contact reach | **held** |
| every unhealed pair makes less where (S_P < 0), unsigned | **held** |
| sub-additive (ρ_W < 1) | **held** (0.89–0.94 over R ≥ 4) |
| B_P < 0, B_int small | **half held.** B_P < 0 everywhere except XX at R = 1 in b (+0.005, +0.003: zero). B_int is small at range but **not at contact** (+0.17 to +0.41) |
| a healed GX gives S_P = 0 exactly | **held** |
| frozen: no contact changes, except flickers beside excesses in ab and abc | **not as expected.** Singles held to that (bar two gap media in a). **Pairs flicker near contact in every set,** at most 12 contacts, without changing a class |

## 9. Checks

- **T0 passes** (2026-10-12, 2.5 min).
  - T0.1: local energy changes equal full recomputation (30,427
    flips).
  - T0.2: identity-ordered tables 15 / 3712, independent of the order
    the contacts were entered.
  - T0.3: restricted quench = full quench in 70/70 cases, ending at a
    local minimum.
  - **T0.4: pairs at R = 16 are exactly their singles** in all 8 sets
    for GG, XX and GX (energies add, contacts are the disjoint union).
  - T0.5: healing gives back the flat medium and every distance.
  - T0.6: the instruments read the known states.
- **96 of 96 processes** completed; 2552 CHECK lines passed (six
  contacts, simple, one piece, distinct identities); no FAIL; empty
  .err files.
- **The quench restriction is exact** (T0.3), so E_int^q = 0 at
  R ≥ 6 is an exact statement, not a sampled one.
- **One presentation fix in `tab.go` during the grid** (before the
  final tables). The R = 8 reach column said "contact" while R = 8 was
  still missing; it now says "(missing)". No bar changed.
- **The grid ran at nice 10** beside another workload of the user's,
  from 21:18. Each process is single-threaded and deterministic, so
  this changes no result. It took 3 h 26 min against the registered
  4–6 h estimate.

## 10. What the step says

1. **In three dimensions, with no container, defects attract at
   contact whatever their kinds.**
   - Gap with gap, excess with excess, and gap with excess all bind,
     by sharing the medium's mismatch.
   - The pair acts at contact only (three or four contacts), as c2's
     did in one dimension.
   - **Unlike c2, the unlike pair does not repel.**
2. **A gap and an excess are each other's undoing.**
   - At contact the medium takes the excess back into the gap's
     closure. It is left exactly flat, and the "where" both had removed
     returns exactly.
   - "Something with gaps" here has an anti-gap: medium where there
     was none, which cancels a gap to nothing.
   - In a frozen crystal this happens only where they are made
     touching. From two contacts apart they bind but cannot meet.
3. **No sign.**
   - The energy does not distinguish like from unlike, and it does not
     alternate with sublattice.
   - The where does not distinguish either: every defect and every
     pair makes less where.
   - Of the two things the user reads into the earlier results (a
     signed charge and an unsigned gap–gap interaction), this medium
     shows only the unsigned one.
4. **The where-interaction outreaches the energy.**
   - The energy is exactly local, as stated before running: in a
     frozen crystal it can only count mismatch at contact.
   - The overlap of two defects' less-where is nonzero at every
     separation measured. It falls by about 2 per doubling of R
     (post-hoc).
   - Nothing charges energy for it. Turning it into a force needs a
     medium whose defects move (g2) or a coherent medium (a2), as
     PREREG §1 said.

## Limits

- **Three sizes, up to 4096 quanta,** and R ≤ L/2 on a closed torus.
  The where falloff has three points per kind, the last at the
  half-period.
- **The random process only, at temperature t.** The crystal is frozen,
  so pairs cannot approach or separate; E_int is a static
  interaction, not a measured force.
- **The excess rule is g1's,** and GX is made gap first. Healing at
  R = 1 depends on the excess's reach (within 2 of p2) covering the
  gap's closure. It did in every case, but the rule was registered,
  not derived.
- **Four placements per size.** The crystal is uniform, so placements
  differ only in the direction of p2 from p1 and in the closures.
- **The observer's bookkeeping** (identities, p1 and p2) orders
  candidates and bins readings. The dynamics never reads it.

## Then

Not registered; for the user to choose.

- **g2, a medium that keeps three dimensions but lets defects move**
  (the Kleinert–Zaanen nematic; the closest candidate here is b
  without a).
  - This is the condition for the where-interaction to become a
    force.
  - It also lets a gap and an excess find each other from afar.
- **a2, the gap as a quantum particle on the frozen crystal.**
  Coherent re-forming (tunnelling between closures) could give it a
  band. Then two coherent gaps could interact through more than
  contact.
- **The signed sector.**
  - c3 finds no sign in a medium of bare quanta: gap and excess
    annihilate rather than repel or pair as charges.
  - If a signed charge is real in this frame, it has to come from
    something the bare medium does not carry: an orientation or a
    winding (R5 is unruled, and winding is never seeded by default).
    That question is the user's to rule.
