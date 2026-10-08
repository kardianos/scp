# v102 fabric — pre-registration (a firm nucleon from the fabric?)

**Written:** 2026-10-08, after `../options/DECIDE.md` §12 and before
any result in this directory. The only numbers here are table
constants, cluster coefficients already computed in `../lock/` and
`../options/`, and algebra.
**Authority:** user, this session:

> "Pursue these ideas, see where they work or if they fall down. The
> goal is to not achive a "cellular clumping" but a firm proton or
> neutron from underlying fabric."

**Scope.** A measurement campaign in models. No kernel or law table is
edited. Nothing is mapped onto a measured number: no MeV, no fm, no
1836.

---

## 1. The two ideas under test (from the user's note)

| | reading | model |
|---|---|---|
| **V** | A lock's link is underlying energy. It never holds energy itself; it only mediates. | The fiber is the effective chain: the S-law clusters' h2 + r3 (exact to fourth order in g), with no real space states. The microscopic S-law, whose link holds real quanta, is run beside it as the check. |
| **Q** | The fiber cannot be hot. Its excitations are lock quanta. | Formation is scored by counting lock quanta, not by ground-state population. |

The working point is the best bundle of `../options/` (T0, C-lock,
X1), at g/Δ = 0.1.

## 2. The target

The target is the structural content of `v98/UNFOLD.md` §7 and
`v98/PRINCIPLES.md` §3–6. Not a fit. A **candidate** is an object in
the fiber; each candidate is scored on:

| | criterion | bar |
|---|---|---|
| N1 | **From the fabric.** An excitation or defect of a uniform fiber. Not a cluster of cells built for it, and not a boundary placed for it. | by construction; stated per candidate |
| N2 | **One size.** Its spin density is localized, and the size does not grow with the fiber. | the decay length changes by < 10% between the two longest fibers run |
| N3 | **Spin ½.** | its ground multiplet is a doublet |
| N4 | **Firm.** It is the lowest state of its conserved sector, and there is a gap to the next level of that object. | gap > 0.1 J |
| N5 | **A ladder.** Discrete excited levels of the same object, below the continuum or as resonances. | reported |
| N6 | **Translation by conversion.** It moves while field does not travel. | its band (the spread of its n positional states on an n ring) has width W > 1e-6 J **at t = 0** |
| N7 | **Two occupancies of one lock (p and n).** "Complete" is n and "one open slot" is p, both spin ½ (`v98/PRINCIPLES.md` §6). | both members exist, with spin ½ |
| N8 | **Formation from the fabric.** Some process of the model makes it from field without placing it. | ≥ 1% per incoming photon or per event, under V |

## 3. Algebra fixed before any run

1. **Parity (the handshake lemma).**
   - Under R6 a cell's ceiling counts its link ends. A fiber graph at
     its ceiling therefore holds Σ caps = 2 × (number of links) field
     quanta, which is always even.
   - A conversion removes two quanta into one space quantum, so the
     field count stays even at every moment.
   - N Jones quanta in two modes have total spin ≤ N/2, with
     S ≡ N/2 (mod 1).
   - **So a complete object always has integer spin. Spin ½ needs an
     odd deficit below the ceiling.**
   - With charge read as the occupancy deficit (v98), **a neutral
     ("complete") spin-½ object cannot exist. N7 fails before any
     run.** It holds for C-lock and C-geo alike; ends and junctions do
     not escape it.
2. **The lock quantum (magnon)** is a triplet, by the fiber's SU(2).
   **N3 fails for it by symmetry.** It is still scored on the rest.
3. **A hole is pinned at t = 0.**
   - A hole is one cell at cap − 1. With t = 0 every cell conserves
     N_i + Σ_l c_l exactly (`../lock/PREREG.md` algebra 1).
   - So the hole's position is a conserved quantity, its n positional
     states are exactly degenerate, and **W = 0. N6 fails for the hole
     under T0, by algebra.** The run checks this and then measures
     what transport t buys, and what it costs.
4. **Hole couplings** (g/Δ = 0.1, from the clusters):
   - J_end/J = 1.955 on the hole's two links;
   - first conversion on a hole link costs 0.2275;
   - the three-site remainder centred on the hole is |r3(2,1,2)| = 0.224 J,
     four times the bulk's 0.052 J. Strings through the hole are cheap,
     so r3 is kept.

## 4. Experiments

**F1. V against the S-law.** Rings n = 6 and 8: the effective chain's
spectrum against the microscopic one. Report which levels V removes:
the space-carrying ones. Also the one-hole ring, n = 6.

**F2. Counting lock quanta (Q).**

- (a) One photon on the n = 4 lock (the E5 machinery, exact). Sort
  each source's final ring states into:
  - the ground state;
  - one lock quantum (the lowest triplet multiplet);
  - more;
  - space-carrying levels, which are V's mismatch.

  Report the yield of exactly one lock quantum per photon.
- (b) A sudden start: the energy above the ground state in units of
  the gap, per cell, on effective rings n = 8 and 10 (Néel and three
  Jones starts). That is the density of lock quanta made.

**F3. Candidate: the lock quantum (magnon).** Score N1–N8 from
`../lock/` (L3 walk, L4 mass) and F2. No new run beyond F2.

**F4. Candidate: the hole (one open slot).**

- (a) Effective rings with one hole, n = 6 … 14 (Lanczos, h2 + r3
  including (2,1,2)). Report:
  - the ground spin;
  - the lowest levels by total spin;
  - the magnon threshold (ground state of the hole ring plus the
    no-hole ring's gap);
  - which levels lie below that threshold (bound).
- (b) Effective open chains, DMRG, with R6 end cells (cap 1, which
  screens the edges) and the hole in the middle; L_bulk = 24, 40, 56.
  Report:
  - the spin profile of the S^z = ½ ground state;
  - its decay length (N2);
  - the S^z = 3/2 level (N5).
- (c) Microscopic S-law, one-hole ring n = 6 (Σc ≤ 3): ground spin and
  lowest levels against (a).
- (d) Translation:
  - microscopic one-hole ring n = 6, transport t (T2 rule, the hole's
    only way to move), t ∈ {0, 0.1J, 0.3J, J, 3J, 10J, 0.1Δ};
  - report the band width W of the six lowest S^z = ½ states, and the
    ground spin S;
  - N6 is gated at t = 0. The t needed for W > 0, and whether the
    hole stays spin ½ there or grows a polarized cloud, are reported.
- (e) Formation: a hole is made when a photon and a lock quantum convert
  at the boundary and the space leaves (E5, F1). Report that per-photon
  rate. Then V: the boundary link holds no energy, so the space cannot
  leave. N8 is gated under V.

## 5. Endings (exactly one is claimed in `DECIDE.md`)

| ending | condition | meaning |
|---|---|---|
| A | some candidate passes N1–N4, N6 and N8 under V and T0 | a firm spin-½ object from the fabric, which walks and forms. The neutron is absent by §3.1 |
| B | some candidate passes N1–N4, but N6 or N8 needs a cost (a transport constant, or real link occupation) | a firm spin-½ object exists, but it cannot walk or form under the readings |
| C | no candidate passes N1–N4 | the fabric under V hosts no firm spin-½ object |

N7 fails by §3.1 in every ending unless charge is re-read. That
re-reading is the user's call, not this page's.

## 6. What does not count

- Placing a hole by hand and calling it formation.
- Reading a decay length, ξ, or any length as a cell count.
- Choosing g, t or a reading because a criterion passes.
- Mapping anything onto 938 MeV, 0.84 fm or 1836, or onto quark
  species.
- No hedgehog, no radial profile seeded or scored.
- Softening a bar after a result.
