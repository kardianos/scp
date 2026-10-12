# v103 lane c, step c3 — pre-registration (two defects in the medium that holds three dimensions)

**Written:** 2026-10-12, before any run of this step. **Nothing has
been run:** no compile, no gate, no shakedown, no registered run.

**Authority:** the user, 2026-10-12: "Setup the next combintorial step,
and commit (no coauthor) and push, but do NOT run yet. I have another
process running, so you can read and write, but do not run anything,
not even a smoke test yet."

The step was named in g1's decision (`../../gate/g1/DECIDE.md`,
"Then"): "c3 (3-D), two gaps in the frozen crystal … Their energy and
their shared less-where against separation, gap with gap, excess with
excess, and gap with excess. This is the sign test of §3, now as an
interaction."
- The README's older c3 entry ("tweaked links in 2-D", from c2's Then)
  is overtaken. g1 has since made a 3-D medium of bare quanta, so the
  two-defect question is asked there, in three dimensions and with no
  container.
- The mobility question (fracton binding) needs defects that can move.
  It goes with g2 and a2.

**A measurement in models; it rules nothing.**

## 0. Frame check (R14, `../../FRAME.md`)

| question | answer for c3 |
|---|---|
| **What are the places?** | **The quanta themselves.** This is g1's medium, unchanged: contact lists and nothing else. |
| **Can a place exist with no medium in it?** | **No.** A removed quantum leaves no slot: its six contacts are re-paired among its neighbours at once. An added quantum takes three existing contacts. |
| **How is distance counted?** | **Through the medium:** contacts along the shortest path. The separation R of the two defects is the number of contacts between the two quanta where they are made, counted in the medium before either is made. Nothing else supplies a separation. |
| **What is a gap, and how is it read?** | **A missing quantum**, read from the medium's own shape: the mismatched contacts and their clusters (g1). An **excess** is one quantum more. **A pair that heals** is read by the medium having no mismatch left. No label (R13): the dynamics never reads where a defect was made. |
| **Container caveat** | **None** (g1's medium). The observer's bookkeeping is disclosed in §8: quantum identities order candidate lists and ties, and the two creation sites bin one reading (B). |

## 1. Why this step, and what it can and cannot show

**What g1 left.**
- The medium holds three dimensions only as a frozen crystal (a or b).
- One gap is a readable, frozen, contact-range defect.
- The gap and the excess both make less where.
- The gap's less-where reaches the whole medium (post-hoc: it
  collapses on r/L).

**The user's reading (FRAME):** two things are going on, a signed
charge and an unsigned gap–gap interaction.
- In 1-D, c2 found like-with-like at contact in the random process.
  Gaps cluster with gaps and excesses with excesses, and unlike pairs
  repel.
- c3 asks the same question in three dimensions, in a medium with no
  container: the sign test, now as an interaction.

**Three readings.**
1. **The pair's energy against separation:**
   E_int(R) = E_pair − E_1 − E_2 + E_0. Its sign at contact for GG,
   XX and GX is like against unlike.
2. **The shared less-where.**
   - Δ₂(u, v) = d_pair − d_1 − d_2 + d_0, for every pair of quanta
     present in all four media. It asks whether the medium holds less
     where with both defects than the sum of each alone.
   - "Less space between gaps" (FRAME), counted through the medium:
     the B readings of §4.
3. **What the pair becomes:** two defects, one, or a healed medium. A
   gap plus an excess restoring the medium is the dictionary's
   antiparticle line.

**What it cannot show (said before running).**
- **Energy here is local.**
  - Every term (squares, links, triangles) lives within two contacts.
  - In a frozen crystal, two defects whose mismatches never come
    within reach of each other share no term. So E_int = 0 exactly
    beyond contact, by construction.
  - **c3 cannot find a long-range force in this medium.** It measures
    where the contact range ends; T0.4 gates the exact zero at one far
    separation.
- **The less-where is long-range (g1), but nothing in a, b or c charges
  for distance.** So less-where costs no energy here.
  - Whether a geometric interaction can become a force needs defects
    that move (g2) or a coherent medium (a2).
  - c3 measures the geometric part.
- **Defects cannot move classically here.** The separation is set when
  the defects are made. Attraction is read as energy, never as motion.

**Known physics (structure and method only, R8).**
- **Point defects in crystals.**
  - A vacancy and an interstitial within a small recombination volume
    annihilate (a Frenkel pair recombines).
  - Vacancies bind into divacancies at nearest-neighbour separation.
  - In continuum elasticity, point defects interact at long range as
    1/r³ (Eshelby).
  - The g1 medium has no elastic far field in its energy, because the
    energy counts only mismatch. So no 1/r³ tail is expected.
- **Fracton–elasticity duality** (Pretko and Radzihovsky 2018, cited
  in g1). Vacancies and interstitials are the conserved charges of
  crystal elasticity, and in a frozen crystal they are immobile. c2
  found that structure in one dimension.

## 2. The medium (g1's, unchanged)

- **N quanta and their contacts:**
  - six contacts per quantum, all paired;
  - two quanta touch at most once, and no quantum touches itself.
- **Re-forming:** g1's local flips, proposal and Metropolis
  acceptance at t = 1.
- **Candidates:** a (flat squares), b (octahedral link), c (no
  triangles), at 2t per unit deviation squared, in **all eight sets**.
- **The code is carried from `../../gate/g1/cas/geo3d.c` verbatim**
  (medium, energy, flips, the gap and excess candidates, the Gibbs
  pick). It is extended by:
  - **quantum identities** (oid), carried through relabelling, so
    that the same quantum can be matched across media;
  - **candidate lists in identity order;**
  - **the quench** (§3).
- **By g1, the sets fall into two groups.**
  - none and c melt into an expander, so the defects there are
    unreadable. They are controls and are not scored.
  - a, ac, b, bc, ab and abc hold as frozen crystals.

## 3. Making two defects

**Sizes.** L = 8, 12, 16, so N = 512, 1728, 4096. The bars use 12 and
16; L = 8 is reported for the collapse readings.

**The medium M.**
- The flat start, then 500 sweeps at t, gives s0.
- In the holding sets nothing changes (g1); the tool reports the
  contacts changed against the flat start.
- s0 uses one stream, the same for every placement.

**Placements.**
- K = 4 placements per (set, L): k = 0, 1, 2, 3.
- Site p1 is uniform over the quanta.
- Separations R ∈ {1, 2, 3, 4, 6, 8} with 2R ≤ L: 1–4 at L = 8, up to
  6 at L = 12, and up to 8 at L = 16.
- Site p2 is uniform among the quanta at distance exactly R from p1 in
  s0, counted through the medium.

**One defect.**
- **A gap at p:** remove p, and re-pair its six contacts among its
  neighbours by Gibbs weight at t over the valid pairings (g1).
- **An excess at p:** a new quantum takes three pairwise-disjoint
  contacts with both ends within distance 2 of p, chosen by Gibbs
  weight at t (g1).

**Kinds.**

| kind | made |
|---|---|
| G1, X1 | a gap / an excess at p1, alone |
| G2, X2 | a gap / an excess at p2, alone |
| GG | G1, then G2 |
| XX | X1, then X2 |
| GX | G1, then X2 (the order is fixed; XG is not run) |

**Shared closures.**
- Each defect has its own random stream. Candidate lists are
  enumerated in the quanta's identity order.
- So a defect made beside another and the same defect made alone are
  the same closure, whenever the other is out of reach (T0.2, T0.4).

**The quench.** After making, apply the strictly downhill flip with the
lowest energy change, ties broken by the identities of its four quanta.
Repeat until none is left.
- **Why.** In g1, closures relaxed by whichever downhill flip came up
  first (a: 144 → 120t at one size and not at the other). Two defects
  far apart must relax exactly as they do alone, so that E_int and Δ₂
  read the interaction and nothing else. The quench makes the
  relaxation a function of the state, and steepest descent on separate
  regions decomposes into the regions' own descents.
- **Restriction, exact.** Only flips whose middle quantum lies within 3
  contacts of a mismatched quantum are tried.
  - A flip's energy change is counted over its four quanta and their
    contacts.
  - If none of those touches a mismatched contact, the flip starts
    from zero local energy and cannot go down.
  - T0.3 checks the restriction against the full enumeration.
- The quench is a relaxation toward zero temperature. The window then
  samples at t.

**Window.** 400 sweeps at t, read every 100 sweeps (4 samples).
Accepted flips are reported. Persistence is the number of contacts
changed between the quench and the end.

**Seeds.**
- Base = 2000003·L + 7·(4a + 2b + c) + 1.
- Stream(k, R, role) = base·1000003 + 100000·k + 100·R + role, expanded
  by splitmix64 into xoshiro256**.
- Defects at p1, and their windows, use R = 0.
- Roles:
  - 1: medium;
  - 2: site p1;
  - 3: site p2;
  - 4 to 7: making G1, X1, G2, X2;
  - 8 to 14: the windows of G1, X1, G2, X2, GG, XX, GX.
- The excesses take identities N (the first) and N + 1 (the second).

**Processes.** One per (set, L, k): 8 × 3 × 4 = 96, with 12 at a time,
heaviest first. Each process makes G1 and X1 once, then G2, X2, GG, XX
and GX at every R.

## 4. Readings

**Every medium, at the quench and at every sample:**
- E;
- mismatched contacts and the quanta touching them (D);
- clusters (D members within 2 joined) and their spread;
- ⟨d⟩, the diameter and the number of pieces.

**A pair at the quench.** When there are exactly two clusters, D_cc is
the medium distance between their centres. A cluster's centre is the
set of quanta whose largest distance to the cluster's members is
smallest, read from the state.

**Healed** means no mismatched contact is left. Then E = 0 in every
set.

**Pair readings.** Here E_0 = E(s0) and d_0 is the distance in s0.
- E_int^made = ΔE_made(pair) − ΔE_made(1) − ΔE_made(2). It is exact, at
  creation.
- E_int^q, after the quench (exact), and E_int^w, from window means.
- **The shared where.** Over every unordered pair of original quanta
  present in all four media (s0, 1, 2, pair):
  - S_X = ⟨d_X − d_0⟩ for X = 1, 2, P. S_P is the pair's less-where.
  - The additive prediction is S_1 + S_2, and ρ_W = S_P / (S_1 + S_2).
  - I = S_P − S_1 − S_2 = ⟨Δ₂⟩.
- **Between the two defects.** Over ordered pairs (u, v) with u within
  2 of p1 and v within 2 of p2 in s0:
  - B_0 = ⟨d_0⟩, and B_X = ⟨d_X − d_0⟩;
  - B_int = B_P − B_1 − B_2.
  - **B_P < 0 is "less space between" the two defects** in the frame's
    sense, counted through the medium.

## 5. Units

- Energy is in t. A sweep is 3N proposals.
- Means and medians are over the K = 4 placements, as each table
  states.

## 6. Gates (T0, before any registered run), in `cas/out/t0.out`

1. **The carried medium**, at L = 6, in every set, from flat and
   melted starts:
   - g1's invariants, plus distinct identities;
   - local energy changes equal the full recomputation (2000 flips per
     start);
   - identities survive a gap and an excess (flat start).
2. **Candidate tables in identity order**, at L = 6 in every set:
   - 15 valid pairings and 3712 triples;
   - local equals full for every candidate;
   - the same tables and the same closure (same stream) whatever order
     the medium's contacts were entered in.
3. **The restricted quench equals the full quench**, at L = 6 and 8, in
   the seven sets with an energy. The configurations are a gap, an
   excess, GG at R = 1, XX at R = 1 and GX at R = 2.
   - Same flips, same energy change, same contacts.
   - It ends at a local minimum: the full quench then finds nothing.
4. **Out of reach, a pair is its two singles exactly**, at L = 16 and
   R = 16, in every set, for GG, XX and GX:
   - the made energies add;
   - the quenched energies add;
   - the contacts added and removed are the disjoint union of the
     singles'.
5. **Healing is reachable and readable**, at L = 8 in every set. A
   gap, then an excess on exactly the three contacts its closure made,
   must give:
   - no mismatch, E = 0;
   - the flat ⟨d⟩ exactly;
   - every distance between the remaining original quanta unchanged.
6. **The instruments on known states**, at L = 8 in every set:
   - the flat medium reads no mismatch, its exact ⟨d⟩ and diameter
     3L/2;
   - one gap reads one cluster, and the matrix mean equals ⟨d⟩;
   - two gaps on opposite sides (R = 12) read two clusters, with a
     centre distance in [8, 12].

## 7. Bars and the ending

**Scored sets.** A set is scored if its medium s0 is the unchanged flat
crystal (no contact changed, no mismatch) at L = 12 and 16, for every
placement. Expected: a, ac, b, bc, ab, abc. none and c are controls.

**The class of each kind at contact.** For each scored set, kind κ and
L ∈ {12, 16}, let v be the median over placements of E_int^w at R = 1.
- **attract** if v ≤ −2t;
- **repel** if v ≥ +2t;
- **none** otherwise;
- **annihilate** if the pair medium is healed at the end of the window
  in at least 3 of 4 placements. This counts as attract.
- **unreadable** if fewer than 3 placements completed.
- **The class must agree at L = 12 and 16.** Otherwise it is
  size-dependent, and the set's ending is M.

2t is one unit deviation, the smallest energy step of the medium.

**The ending of a scored set:**

| ending | condition |
|---|---|
| **U** | GG, XX and GX all attract: unsigned at contact |
| **L** | GG and XX attract; GX does not: like with like (c2's 1-D result) |
| **Q** | GX attracts; neither like pair does: charge-like |
| **N** | nothing interacts at contact |
| **M** | anything else |

**The step's ending** is the letter shared by every scored set.
Otherwise it is "mixed by set", with the table.

**Reported with the ending, not scored:**
- **Energy reach.**
  - R_E is the largest R at which any kind has |median E_int^w| ≥ 2t.
  - The reach is **contact** if every kind is below 2t at R = 8
    (L = 16), and **extended** otherwise.
  - Also the number of pairs with E_int^q ≠ 0 at R ≥ 6.
- **Where.**
  - The sign of S_P for every kind (unsigned if all are negative;
    healed pairs are excluded).
  - The additivity class from the mean ρ_W over R ≥ 4: additive within
    [0.95, 1.05], sub-additive below, super-additive above.
  - B_P and B_int against R.
  - The collapse of I and B_int at R/L = 1/4 and 1/2 across
    L = 8, 12, 16.
- **What pairs become:** clusters at the quench and at the end, D_cc
  against R, and healed counts by kind and R.
- **Persistence:** the contacts changed in the window.

**What I expect** (written before running):
- **none and c melt** and are not scored.
- **Locality.** E_int^q = 0 exactly beyond contact. I expect R_E ≤ 6,
  and the contact reach at R = 8.
- **GG attract at R = 1.** Two closures side by side share mismatch, so
  together they break fewer squares than apart. **XX attract** for the
  same reason.
- **GX annihilates at R = 1.**
  - The excess is made within distance 2 of p2, which is a contact of
    p1, so the three contacts the gap's closure made are likely within
    its reach.
  - Taking exactly those restores the flat medium. Its energy change is
    minus the gap's (−144 to −432t), and every partial overlap leaves
    some mismatch. So I expect it to be the lowest candidate by a wide
    margin, and the Gibbs weight at t to pick it.
  - This is not guaranteed: if one closure contact has an end beyond
    distance 2 of p2, the excess cannot take it.
  - At R ≥ 2, GX attracts by partial overlap.
- **Expected ending: U at contact in every holding set, with contact
  reach.**
  - Like c2, it acts at contact only.
  - Unlike c2, the unlike pair attracts. c2's unlike pairs repelled in
    the random process with headroom.
- **Where.**
  - Every unhealed pair makes less where (S_P < 0) for all three kinds:
    unsigned, as in g1.
  - Sub-additive (ρ_W < 1). Two shortcuts do not both shorten a path
    that can use only one of them.
  - B_P < 0: less space between the defects, and B_int small.
  - A healed GX gives S_P = 0 exactly: the where comes back.
- **Frozen.** No contact changes in the window, except flickers beside
  excesses in ab and abc (g1).

## 8. Disclosed before running

- **Nothing has run, not even a compile.**
  - The user asked that nothing be run. The tool (`cas/c3pair.c`, C)
    and the tables (`cas/tab.go`, Go) were written and read through,
    but never compiled.
  - **The first step at go** is to compile with `-Wall -Wextra` and run
    `go vet`.
  - A compile fix that does not change behaviour is made without
    amendment. Any other change is recorded here as a dated amendment
    before T0.
- **The quench is new relative to g1,** which relaxed by sampling at t.
  It takes the steepest branch, not a random one.
- **The observer's bookkeeping.**
  - Quantum identities order candidate lists and quench ties, and
    match quanta across the four media for Δ₂.
  - The creation sites p1 and p2 bin the B reading.
  - The dynamics never reads them. Readability and healing are read
    from the state alone.
- **The window is shorter than g1's** (400 sweeps, 4 samples). g1
  measured acceptance ≤ 1.5 × 10⁻⁶ in the holding sets, so here the
  window checks persistence; it does not measure motion.
- **The excess rule is g1's.** Other insertions are not tested, and the
  order of GX is fixed (gap first).
- **Four placements per size.** In the holding sets the crystal is
  uniform, so placements differ only in the direction of p2 from p1
  and in the closures.
- **Nothing is scanned.** t, 2t, the separations, K and the sizes are
  registered here.
- **In the melted sets,** E_int and Δ₂ are noise, because the medium
  re-forms everywhere. They are reported only.
- **Estimated cost (unmeasured).** From g1's timing (about 0.37 s per
  sweep at L = 16 with b):
  - the heaviest process (L = 16, b on) takes about 1.5 h;
  - the grid takes 4–6 h at 12 processes;
  - an L = 16 process holds about 100 MB of distance matrices.

## 9. Reproduce

From `cas/`:
1. `./run_t0.sh` compiles the tool and runs the gates. It must end with
   `T0 PASS`.
2. `./run_c3.sh` runs 96 processes, 12 at a time. It refuses to start
   unless T0 passed. Progress goes to `out/run.log`, and the marker
   `out/run.done` appears at the end.
3. `go run tab.go > out/tab.out`.
4. Write `../DECIDE.md`.

## Then (after DECIDE)

- **A contact-range ending (U, L or Q)** sends the force question to:
  - g2, mobile defects in a medium that keeps orientation (the
    Kleinert–Zaanen nematic);
  - a2, a coherent gap on the crystal;
  - a matched-reference field footprint, if needed.
- **If GX heals:** this is lane a's "annihilation of a void with an
  excess". The energy released (E_G + E_X) is the conversion. Light is
  not matter, so the released energy is counted and nothing else.
