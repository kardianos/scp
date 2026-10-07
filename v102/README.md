# v102 — SPIN-ONE FIBER

**Opened:** 2026-10-07.
**Authority:** user, this session. After `v99/reopen/OPTIONS.md`, the
user accepted three of the five A+B premises as rulings ("I'm okay
with those changes. Go a bit further."). They are R1, R2, and R5
below. The other two are gates. **This document is the plan.** Gates:
`PREREG.md`. A pass or a stop is claimed only in `DECIDE.md`.

**No kernel is carried.** `v97` is not edited. Opening v102 is not a
law-table change and does not authorize a mesh. The rulings change
what the theory may assume. They do not change the code.

**Status:** the gates were run on 2026-10-07, ending C.

- G-wedge **passes**: the antiferromagnetic, bilinear link exchange is
  derived from R1 and R2.
- G-two **fails**: a cell's occupation equals the lock's occupancy
  ceiling, whose value nothing fixes.

One sentence would fix it: "the occupancy ceiling of a lock is one
quantum per link end." Phase 1 is not opened. See `DECIDE.md`.
Reproduce:
`python3 -I v102/gates/cas/wedge.py`, then
`python3 -I v102/gates/cas/two.py`.

---

## Thesis under test

A lock is the gapped quantum of a fiber chain. Each cell of the chain
carries quanta of a two-plane (Jones-vector) field. Each link converts
field into space through the link's handedness. If that is the whole
story, the fiber is a spin-1 antiferromagnetic ring (family 268 in
`github.com/openai/math`), and its correlation length (ξ ≈ 6.03
grains) is forced, with no continuous coefficient.

Two things stand between the rulings and that claim:

1. **G-wedge.** Is the link exchange antiferromagnetic and bilinear
   because of the premises, or by choice?
2. **G-two.** Is there one spin per cell, `N = 2` quanta, because of
   the premises, or by choice?

---

## Rulings (user, 2026-10-07)

| | ruling | what it reverses or retires |
|---|---|---|
| R1 | Each cell's field is two complex amplitudes, one per transverse plane (a Jones vector), with a U(2)-invariant coupling (the wording accepted in `v99/reopen/FEASIBILITY.md` §7 A). Every coupling of the field is SU(2)-invariant, and transport is U(2)-invariant. A conversion changes the number of field quanta, so it cannot also keep the U(1). | The kernel's single complex amplitude per cell. |
| R2 | Transport is quantized. The quanta are those of the classical wave, so they are bosons. | v89 law 1, "transport is never quantized." |
| R5 | Integer `W` is retired as the lock. The lock candidate is the gapped quantum of the fiber chain. | v98 §4's "winding is an integer" as the lock, and the v99 wrap thesis. |

Not ruled at opening: the occupation per cell and the form and sign
of the link exchange. Those were the gates.

| | ruling, after the gates | closes |
|---|---|---|
| R6 | **The occupancy ceiling of a lock is one quantum per link end** (user, 2026-10-07: "That's okay, can you continue with that modification?"). | `DECIDE.md` §2's open sentence. With it, W4's ceiling exists, `N = 2` on the ring, and an end cell of an open chain (degree 1) holds one quantum. |

Phase 1 is opened under R1, R2, R5, and R6: `PHASE1.md`, with gates
in `phase1/PREREG.md`.

**Phase 1 stopped** (2026-10-07; `phase1/DECIDE.md`).

- **The lock is a phase.** Transport must stay below roughly
  t/Δ ≈ 0.16–0.26; above that, the field converts into a ferromagnet.
- **Inside the phase, the fiber is a Haldane chain** with one line and
  no bound pair.
- **ξ ≈ 6 grains only at leading order.** It drifts continuously with
  g/Δ and t/Δ: 6.06, 5.79, and 5.07 at g/Δ = 0.05, 0.1, and 0.2.
- **Two registered bars were defective** (P1 at n = 4, and P6(c)).
- **The two physics misses are real:** the transport threshold, and
  ξ's drift with g/Δ.

**Law step** (2026-10-07; `law/DECIDE.md`): **ending C.**

- The current table does not host the lock.
  - Space has no pitch, so a field pair resonates with it and decays.
  - g has no counterpart: the table's conversion is a rate law.
  - Even borrowing the retired dense rung, 3-D transport is about 10×
    too fast.
- The lock needs four new law statements (S1–S4). Each is the user's
  ruling, and a law-table change.
- The two defective Phase 1 bars were re-registered on new cases, and
  both passed.

## What is carried

- `v89/PRINCIPLE.md`, in full. Where the kernel and the principle
  disagree, the principle governs. That matters here, because the
  kernel converts on a cell while §3 converts through a link.
- `v89/README.md` law 2, "a channel is a joint process of its two
  ends."
- v98 §1, "locked field cannot travel as itself," and the unfold
  filter.
- The v99 representation split.
- v100 §6: no hedgehog and no radial profile.
- v101 rules 4 and 5: no continuous ratio to break a tie, and no cell
  count read off a result.

## What is not carried

The v97 kernel as the law. The two-atom credit as an occupation. A
scanned coupling ratio. Integer winding as the lock.

---

## Order

| step | what | opens when |
|---|---|---|
| G-wedge | link exchange: operator, form, sign, net sign | `PREREG.md` written |
| G-two | occupation per cell | `gates/WEDGE.md` decided |
| Phase 1 | derived chain: line, ξ, no second line, slit regression | both gates pass; not authorized by this page |

## What not to do

1. Do not consult any version before v89.
2. Do not edit the v97 kernel or `v97/fab`.
3. Do not read ξ as a cell count, and do not take `n = 6` from it.
4. Do not scan `g/Δ`, `t/U`, or a biquadratic angle to place a line.
5. Do not use the two-atom credit (a credit that lapses at two atoms)
   as an occupation.
6. Do not map anything onto 938 MeV, 0.84 fm, or 1836.
7. Do not seed or score a hedgehog or a radial profile.
