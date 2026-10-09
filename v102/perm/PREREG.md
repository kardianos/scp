# v102 perm — pre-registration (re-readings and decisions, alone and in combination)

**Written:** 2026-10-08, after `../fabric/DECIDE.md` and before any
result here except the two-ring check of W listed in §4.
**Authority:** user, this session:

> "Do commit, then on each re-readings: try each in isolation, then
> combine each together in each permemtation and then all toegher. Same
> with decisions. I'm not going to do well making these decisions. Try
> each in permeitatinos in isolation and then together. See what works
> and what doesn. Record results."

**Scope.** A measurement campaign in models. No kernel or law table is
edited, and no choice is ruled. Nothing is mapped onto a measured
number.

---

## 1. The factors

The carried readings are V (the link holds no energy) and Q (count lock
quanta). The working point is `../options/` B-label-cool (T0, C-lock,
X1), at g/Δ = 0.1. The candidate is the hole of `../fabric/`.

**The re-reading (the charge question):**

| code | reading | model |
|---|---|---|
| R0 | the fiber's SU(2) is spin; charge is the occupancy deficit | `../fabric/` as is |
| **I** | the fiber's SU(2) is isospin. A hole's charge is `I^z + ½` (which plane holds the open slot) | the same Hamiltonian, relabelled |
| **W** | v98's named widening, 3+4 (`v98/PRINCIPLES.md` §4: "only if 3+2 cannot host handedness and the n/p slot together"). Each quantum carries spin ½ × isospin ½ (four modes). Conversion pairs a quantum across a link into a singlet in both (P = ε ⊗ ε). Charge is `I^z + ½` | the M = 4 S-law (`cas/graphm.py`) |

**The decisions:**

| code | decision | registered point | also reported |
|---|---|---|---|
| **T** | field transport inside a lock, bounded, so a hole can walk (T2 rule: refused into a cell at its cap) | t = J | t = 0.1 J, 3 J |
| **S** | space may leave through a lock's boundary link (a limited exception to V) | 1-D space lead at `t_s = s_k = 0.06` | |

The third decision, charge, is the re-reading axis.

**The combinations (all run):**

| | none | T | S | T + S |
|---|---|---|---|---|
| **R0** | baseline | T | S | T+S |
| **I** | I | I+T | I+S | I+T+S |
| **W** | W | W+T | W+S | W+T+S |

**I + W is not a separate cell.** W already reads one of its two
doublets as isospin. I applied to W only swaps the two names, so I+W is
W by symmetry, and this is checked numerically, not assumed (§3 P3).

## 2. What is scored

For each cell, two things.

**K: the lock's health.** The fiber at its ceiling, without a hole, must
satisfy all three:

- a unique singlet ground state;
- a gap that survives the longest chain run (> 0.1 J);
- ⟨Σc⟩/n ≤ 0.1 (the hold).

**The candidate's card** (N1–N8, the bars of `../fabric/PREREG.md` §2):

| | criterion | bar |
|---|---|---|
| N1 | from the fabric | |
| N2 | one size | decay length changes < 10% between the two longest chains |
| N3 | **spin ½** | |
| N4 | firm | gap > 0.1 J |
| N5 | a ladder | a bound excited level below the threshold "hole + one lock quantum". Reported, not gated |
| N6 | walks | band width W > 1e-6 J |
| N7 | **p and n**: one object with two occupancies, both spin ½ | |
| N8 | forms | ≥ 1% per resonant photon |

Under I, the fiber's SU(2) is isospin. N3 then asks for a spin from
elsewhere.

**Tiers:**

| tier | condition |
|---|---|
| **Full** | K and N1–N4, N6, N7, N8 all pass |
| **Partial** | K and N1–N4 pass; at least one of N6, N7, N8 fails |
| **Fail** | K fails, or any of N1–N4 fails |

N5 is reported in every tier.

## 3. Algebra fixed before the runs

- **P1 (I is a relabelling).**
  - Every number under I equals the same number under R0.
  - Under I the hole's ground doublet is (p, n): one object, two
    occupancies, the same size. The split is exactly 0 by SU(2). v98
    asks for a small split, so that is reported.
  - The fiber has no second SU(2), so **spin ½ is absent: N3 fails
    under I in every cell,** and N7 passes.
- **P2 (parity under W).**
  - Each quantum carries (½, ½). A fiber at its ceiling holds an even
    number N of them, so S and I are both integers.
  - A hole (N odd) has S and I both half-integer. Its lowest possible
    multiplet is (½, ½): p↑, p↓, n↑, n↓.
  - With charge `I^z + ½`, p and n are both holes. The neutron wall of
    `../fabric/` §5.1 does not apply.
  - Whether the ground multiplet *is* (½, ½) is measured.
- **P3 (I + W = W).** The pair form ε ⊗ ε is symmetric under swapping
  the two doublets, so W's spectrum is symmetric under spin ↔ isospin.
  Checked on the n = 4 ring.
- **P4 (T moves the hole).**
  - At t > 0 the hole's position is no longer conserved, so W > 0 is
    expected.
  - The lock without a hole is untouched by T: the T2 rule refuses
    hops into full cells, so a full fiber has nowhere to hop to.
  - At fourth order, conversion opens holes that could hop. The 3+2
    margin for that is t ≲ 0.17 Δ = 8.8 J (`../lock/` L5b). t = J sits
    9× inside it.
- **P5 (S is formation).**
  - The 3+2 numbers exist: 12.1% (g/Δ = 0.1) and 24.9% (0.2) per
    resonant photon (`../fabric/` F2a).
  - By total S^z (and I^z) conservation, the hole made carries the
    photon's polarization. Under I and W, a photon's state therefore
    decides p or n.

## 4. Known before this page

`cas/out/w_lock.out` was run first, to see whether W has a lock at all
(g/Δ = 0.1). Effective rings, microscopic n = 4 as the check:

- **3+2 control:** n = 4 and n = 6 reproduce `../lock/`.
- **W, n = 4 and 6:** the ground state is unique and lowest in sector
  (0, 0). The first excitation is 0.0190 (n = 4) and 0.0144 (n = 6),
  and it is a 15-fold multiplet:
  (S, I) = (1,1) ⊕ (1,0) ⊕ (0,1).
- **The microscopic W ring, n = 4,** agrees: E₀ −0.127331 against
  −0.127193, and first gap 0.0191 against 0.0190.

K under W is not yet settled. That needs long chains.

## 5. Experiments

| code | what | instrument |
|---|---|---|
| **X1** | K under W: DMRG gaps of W chains with R6 cap-1 ends, L_bulk = 16, 24, 32; the ground state's (S, I) | `weff.py` + DMRG |
| **X2** | The hole under W at t = 0: rings n = 4, 6 by (S, I); chains, DMRG, L_bulk = 16, 24, 32 — the profile, the decay length (N2), E(3/2, 3/2) and the other multiplets against the threshold (N5) | |
| **X3** | T under 3+2 and W: mobile-hole effective chains (site = full or hole) with hop t. Rings (exact): band width, ground multiplet. Chains (DMRG): the decay length of the hole-centred profile, and the S = 3/2 level against the threshold. t = 0.1 J, J, 3 J | `mobile.py` |
| **X4** | S under W: the per-photon hole-making probability, by the S-matrix of `../options/cas/e5_cool.py` with the channel space truncated to the ring's low levels. The truncation is calibrated on 3+2 against the exact E5 number | |
| **X5** | Reported, not gated: two holes, energy against separation (3+2 and W, t = 0 and t = J). Do holes bind? (This is the JOINT class, `v98/UNFOLD.md` §6) | |
| **X6** | I and I+W: the relabelling checks (P1, P3) | |

Every cell of §1's table is scored from X1–X6 and `../fabric/`.

## 6. What does not count

- Picking a cell, t or a reading because a criterion passes.
- Reading a decay length as a cell count.
- Mapping onto 938 MeV, 0.84 fm, 1836, quark species, or the pion.
  Multiplets are named by (S, I) only.
- Softening a bar after a result.
