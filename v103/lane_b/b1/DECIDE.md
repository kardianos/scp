# v103 lane b, step b1 — decision (where classical begins)

**Written:** 2026-10-11, after the registered runs. Gates: `PREREG.md`.
A measurement in models; it rules nothing.

- Tool: `../../tools/bx.c` (C), passing T0 at 1, 4, 8 and 16 threads.
- Runs: `cas/run_b1.sh`, 84 points.
- Tables: `cas/tab.go` (Go), output in `cas/out/tab.out`.

**Frame caveat (R14, ruled after this step was registered).** b1's
medium is a standard lattice model: quanta held in a fixed ring of
cells. A fixed set of cells is a container, "something in something".
So b1 measures known quantum mechanics in the container frame. It does
not test the gap frame. See `../../FRAME.md`.

## Ending: **C** (complementary)

No point (rule set, ring, occupation) is two-regime.
- Wherever the bulk is a classical, phase-coherent field, the missing
  quantum is smeared through it.
- Wherever the missing quantum is a discrete thing, here or there, the
  bulk has no phase at all.
- Raising the occupation does not escape this.

| rule set | bulk g1(L/2) | void β | regime, at every n̄ and both rings |
|---|---|---|---|
| none, c (free; U = 0, no cap) | 1.0000 | −1.0000 | classical bulk, smeared void |
| a, ab, ac, abc (number lock: cap = n̄) | 0.0000 | +1.0000 | discrete void, no phase |
| b, bc (U = 4t) | 0.70 → 0.88, rising with n̄ | −0.14 → −0.07 | neither |

Bars: the bulk is classical if g1 ≥ 0.9, and the void is discrete if
β ≥ 0.5.

## Trends

**b and bc (stiffness U = 4t).** These are the only sets that move
with occupation.

| ring | g1, n̄ = 1 → 8 | β, n̄ = 2 → 8 |
|---|---|---|
| L = 4 | 0.782 → 0.876 | −0.137 → −0.088 |
| L = 6 | 0.703 → 0.847 | −0.109 → −0.069 |

- β never turns positive. A stiff medium smears the void less than a
  free one, but it never makes it discrete.
- g1 is lower on the longer ring at every n̄. That is the
  one-dimensional signature: phase order decays as a power law. So
  "classical bulk" at fixed stiffness is itself a finite-size statement
  here.

**c (two planes) changes nothing measured.**
- With U = 0 and with U = 4t, the c-on values of g1, σ_N and β equal
  the c-off values to every printed digit.
- That is the known theorem for SU(2)-symmetric bosons: their ground
  state is the fully polarized multiplet, so its density and coherence
  equal the one-plane model's.

**a (number lock) freezes the bulk completely.**
- σ_N = 0 and g1 = 0. The void is an exact unit: β = 1 to every digit.

## What I expected (PREREG §5, written before running), against what happened

| expected | outcome |
|---|---|
| a on: g1 = 0, β = 1 | **held** |
| a off, b off: g1 = 1, β = −1 | **held** |
| b moves between the two | **held**, but it moves only part-way, and β stays below 0 |
| c changes little | **held**, and exactly, by SU(2) |
| ending C | **held** |

## Disclosures

1. **Z at U = 0 (sets none and c) is a roundoff artifact. No bar uses
   it.**
   - For free bosons, b_k|G⟩ = 0 exactly at every k ≠ 0.
   - The computed ground state carries roundoff of about 1e-9. That
     residue at k ≠ 0 projects onto the (N−1) ground state, and the
     strict band-minimum comparison then picks its random weight.
     The printed values are 0.75, 0.65, 0.0, 0.61, 0.016 and 0.002.
   - Analytically, Z = 1 exactly for these sets. b_0|G⟩ is the (N−1)
     ground state for one plane, and also for two planes at
     2S^z = 0 → −1.
   - g1 and β use no spectral routine and are unaffected. The fix, if
     Z is ever scored at U = 0, is a relative norm threshold before
     the spectral call. Any fix reruns T0.
2. **Z for the number-locked two-plane sets (ac, abc) is small,
   0.08 → 0.**
   - That is physical. With cap = n̄ and no exchange term, the locked
     two-plane medium is spin-degenerate, and the void's weight spreads
     over the degenerate manifold.
3. **Four points ended at the restart cap rather than at the 1e-9
   residual:**
   - b, L = 4, n̄ = 8: 3.5e-7;
   - b, L = 6, n̄ = 3: 1.2e-8;
   - bc, L = 4, n̄ = 6: 4.0e-8;
   - bc, L = 4, n̄ = 8: 1.1e-8.

   The eigenvector error is about residual/gap ≤ 1e-6, far below the
   4-decimal values and the bars. Every other point is ≤ 2e-9, and
   every .err file is empty.

## What lane b now says

- **In known quantum mechanics on a lattice, occupation alone does not
  split the world into a classical bulk and quantum voids.**
  - Number and phase trade off exactly as known physics says.
  - A medium coherent enough to be a classical field dissolves a unit
    deficit into itself (β = −1).
  - A medium locked enough to keep the deficit whole has no field
    (g1 = 0).
- **Stiffness moves a medium along that trade-off. It does not leave
  it.**
- **This is the container frame's answer.** In it, a void is a missing
  amplitude in fixed cells, and an amplitude can spread over the cells.
- **In the gap frame (R14), the question changes.** A gap is medium that
  is not there, not an amplitude missing from a slot. So whether a gap
  can be smeared at all is not yet asked.

## Limits

- **One dimension and small rings** (L = 4, 6). Phase order in 1-D is
  only power-law.
- **One stiffness, U = 4t.** It was registered as the natural unit, not
  scanned.
- **A quantum medium was assumed.** Whether a medium with no ħ put in
  makes quantum voids is a later step.

## Then

**b2, in the gap frame.**
- The medium's own bonds are its geometry, and a gap is a missing
  quantum. Ask:
  - can a gap be smeared, that is, can a superposition of networks
    spread a gap's place?
  - does the network's large-scale geometry become classical (smooth,
    self-averaging) while gaps stay whole?
- It shares c2's network tool once that tool exists and is gated.
