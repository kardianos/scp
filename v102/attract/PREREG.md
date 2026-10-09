# v102 attract — pre-registration (R8; decisions 2–3 by test; the search for an attraction)

**Written:** 2026-10-08, after `../perm/DECIDE.md` and before any
result in this directory. One feasibility timing (W with a mobile hole,
8 cells, χ = 128: 229 s, truncation 3e-5) was run first, to know the
runs were affordable. It is not a result.
**Authority:** user, this session:

> "Allow the law change for 1. Your judgement (or test) for 2-3. For (4)
> yes please proceed."

**Scope.**

- No kernel or law table is edited. R8 changes what the theory may
  assume; a kernel step still needs explicit authorization.
- Nothing is mapped onto a measured number. Channels are named by
  (S, I) only.

---

## 1. Ruling R8 (user)

**The 3+4 widening (W) of `../perm/` is adopted for v102.**

- Each field quantum carries spin ½ and isospin ½: four modes per
  cell.
- Conversion pairs a quantum across a link into a singlet in both
  (P = ε ⊗ ε).
- Energies are those of R1–R7: the load law in the total count per
  cell, Ω, Δ.
- Charge is read as `I^z + ½` on an open slot. The open slot ("the
  hole") is the nucleon candidate: a (½, ½) quartet.

The readings V and Q, the label C-lock, and the boundary conversion X1
are carried.

## 2. Decisions 2 and 3, settled by test (rules fixed here)

**Decision 2: transport inside a lock (T).**

| option | statement |
|---|---|
| T-none | no transport in a lock (the hole is pinned exactly) |
| T-bound | transport t > 0 inside a lock, T2 rule (refused into a cell at cap). A cell may empty |
| T-floor | as T-bound, and transport may not empty a lock cell (a lock cell keeps ≥ 1 quantum) |

Rule:

1. Among the options that let the hole walk (N6), keep those where:
   - the hole's ground multiplet stays (½, ½);
   - N2 passes for the moving hole (registered fit, two longest W
     chains);
   - two holes do not collapse into a fiber cut (B2).
2. If more than one remains, take the one with fewer added statements.
3. Its bound is stated as the measured window, t < t_max. It is not
   placed at a chosen point.

**Decision 3: space leaving through the lock's boundary (S).**

| option | statement |
|---|---|
| S-none | V everywhere; no route found to form a hole |
| S-boundary | space may leave a lock only through a link with one end outside the lock (flow `s_k`) |

Rule: adopt S-boundary if:

- it is the only formation route measured;
- it leaves a lock in cold light untouched;
- under the decision-2 outcome, the holes it makes do not cut the
  fiber (B2).

Otherwise S-none.

**Tests for the decisions:**

- **B1. The moving W hole.** TeNPy, leading order, mobile sites,
  L_bulk = 16 and 24, χ = 128, t ∈ {0.1, 1, 3} J0. Report:
  - the ground multiplet;
  - the hole-centred profile and its decay length (N2, registered fit
    over each chain's own range, as `../perm/` X3);
  - the conditional profile and the 90% radius, reported beside it.
- **B2. Two holes under each T option.**
  - T-bound is microscopic (`../perm/` X5: collapse measured).
  - T-floor is the leading-order mobile model, whose local space holds
    N = 2 or 1 only, so the floor is exact.
  - Report whether two holes stay two (A3 below gives their binding).

## 3. Item 4: the search for an attraction

An attraction is anything that holds two objects at a distance. It
would show up as a ladder (a bound excited level of the hole) or a
nucleus (a bound pair of holes).

- **A1. A ladder from transport (a spin polaron?).**
  - Mobile hole: 3+2 at t ∈ {0.1, 1, 3, 5, 7} J0, L_bulk = 24 and 40;
    W at t ∈ {0.1, 1, 3} J0, L_bulk = 16 and 24.
  - Margin m(L) = E(next multiplet) − E(ground) − gap(L), with gap(L)
    the no-hole gap at the same length. The next multiplet is S = 3/2
    for 3+2, and (3/2, 3/2) and (3/2, ½) for W.
  - Extrapolate m in 1/(L+1)².
  - **Bound** if the extrapolated margin is < −0.01 J0 and m < 0 at
    the longest length.
  - The ground multiplet must stay the hole's (½ or (½, ½)), or the
    point is reported as past the Nagaoka onset.
- **A2. A ladder from the fourth-order terms (larger g).**
  - 3+2 static hole, all-order h2 + r3 (`../fabric/cas/f4b_dmrg.py`
    and `../lock/cas/l4_mass.py`'s chain), at g/Δ ∈ {0.1, 0.2, 0.3},
    L_bulk = 24 and 40, χ = 64.
  - The same margin and bar as A1.
  - W at all orders is out of reach (r3 has about 10⁴ operator terms);
    stated, not run.
- **A3. A nucleus: two holes.** Leading order, floor exact (T-floor or
  static).
  - Binding B = E(2 holes) + E(0) − 2E(1), per channel:
    - 3+2: S^z = 0 and 1;
    - W: sectors (2S^z, 2I^z) = (0, 0), (2, 0), (0, 2), (2, 2), which
      separate the (S, I) channels (0,0), (1,0), (0,1), (1,1) by
      recurrence.
  - Static pairs (t = 0) at separations d = 2…6, chain L_bulk = 24.
  - Mobile pairs at t = 0.1 and 1 J0, L_bulk = 16 and 24 (W) or 24 and
    40 (3+2). Also the pair's size: the hole–hole distance
    distribution.
  - **Bound** if B < −0.01 J0 at both lengths, B does not shrink
    toward 0 with L, and the pair's mean separation does not grow with
    L.
  - **Reported, not gated:** which (S, I) channel binds first. v98
    lists the deuteron as the first JOINT row; nothing is mapped.
- **A4. W's other channels.** Static W hole, the (3, 1) sector at
  L_bulk = 24 and 32, to see whether (3/2, ½) lies below
  (3/2, 3/2).

## 4. Endings for item 4 (exactly one is claimed in `DECIDE.md`)

| ending | condition |
|---|---|
| L | a ladder: some A1/A2/A4 case binds, with the ground still the hole's multiplet |
| N | a nucleus: some A3 channel binds |
| LN | both |
| 0 | neither: no attraction at any tested point |

## 5. What does not count

- Choosing t, g or a channel because it binds.
- Reading a pair size or a decay length as a cell count.
- Mapping onto measured binding energies, 938 MeV, 0.84 fm, 1836,
  quark species or the pion.
- Softening a bar after a result.
