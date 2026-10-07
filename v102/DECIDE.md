# v102 — decision

**Written:** 2026-10-07, after `gates/WEDGE.md` and `gates/TWO.md`.
**Status:** this page closes v102 at its gates. It does not open
Phase 1. It does not edit a kernel.

---

## 1. Endings, as registered

| ending | condition | consequence |
|---|---|---|
| A | both gates pass | Phase 1 named, with the bars in PREREG §4 |
| B | G-wedge fails | stop; the exchange is a choice |
| C | G-wedge passes and G-two fails | stop; name the one sentence that fixes `N` |

**The decision is C.**

## 2. What v102 established

1. **The conversion operator is derived.** Under R1 (scope as
   accepted) and R2, the only conversion into a scalar space mode is
   the antisymmetric pair across a link,
   `A_ij = b_{i1}b_{j2} − b_{i2}b_{j1}`. Its on-cell version vanishes
   for bosons. That a conversion is a link process, as v89 §3 says,
   is derived here, not assumed.
2. **The exchange is derived.** Its virtual exchange is exactly
   `(2g²/Δ) S_i·S_j` for every `N`, with no biquadratic term at
   leading order. It is antiferromagnetic because the lock is field
   (`Δ > 0`).
3. **A microscopic model confirms it.** A ring with explicit space
   modes on its links reproduces the spin-1 Heisenberg spectrum
   (0.044% miss at n = 6, g/Δ = 0.1), with a singlet ground state.
4. **Transport is excluded by the ceiling.** Transport's virtual
   exchange is ferromagnetic. It vanishes only at an occupancy
   ceiling. Quanta hoard, so cells sit at that ceiling (W4 rests on
   v98 §1's licensed ceiling).
5. **`N` is set by the ceiling's value, which nothing fixes.** Odd
   values give no gapped lock.

Of the two A+B premises left open in `v99/reopen/OPTIONS.md`, one is
now derived: the antiferromagnetic, bilinear link exchange (W1–W3).
The other, the occupation per cell, reduces to a single sentence:

> **The occupancy ceiling of a lock is one quantum per link end.**

## 3. What is not true

- `N = 2`, `S = 1`, or ξ ≈ 6.03 is derived.
- `n = 6`, or any ring size, is selected. Ring size is still not
  selected: `v99/reopen/OPTIONS.md` E2 found `E₀(n)` concave on even
  rings. At fixed total quanta `Q`, the ring length is `Q/N_c`.
- Spatial spin ½ is derived. The SU(2) here is internal, acting on
  polarization.
- Phase 1 is authorized.

## 4. How to move past this page

**Rule the sentence in §2.** Then:

- W4's ceiling exists, and has a value.
- `N = 2` on the ring (degree 2).
- The fiber is the spin-1 antiferromagnetic Heisenberg ring with
  `J = 2g²/Δ`.

Phase 1 then opens with the PREREG §4 bars:

1. Exact diagonalization or DMRG of the derived chain, from the
   microscopic model, not the effective one.
2. One isolated line.
3. ξ measured, not read as a cell count.
4. No second line below the two-magnon threshold.
5. The `(g/Δ)²` dependence of ξ measured, and not scanned to place
   anything.
6. The slit recovered from the unlocked, linearized field.

**The size gate.** Under the gate as written ("force `n*`"), a ruled
v102 still fails literally, because the forced scale is ξ and not
`n`. Under `v99/reopen/OPTIONS.md` C ("one isolated line plus
coupling-independent ratios"), it would pass: J is the unit, and the
only continuous dependence enters at `(g/Δ)²`. Which reading of
"size" governs is the user's call.

**To reverse W1:** show R1 was meant for transport only. Then `δ` and
`ε` both survive, and their ratio is a choice.

**To reverse W4:** show v98 §1's occupancy ceiling is not a principle
of this programme. Then the net sign is a ratio of constants.

## 5. Premise ledger

| premise | status |
|---|---|
| R1: two amplitudes per cell, U(2)-invariant coupling | ruled (user) |
| R2: quantized transport (bosons) | ruled (user); reverses v89 law 1 |
| R5: integer `W` retired as the lock | ruled (user) |
| antiferromagnetic bilinear link exchange | **derived** (W1–W3) |
| an occupancy ceiling exists | licensed (v98 §1), measured in v90; used by W4 |
| ceiling = one quantum per link end | **open: the one sentence** |
