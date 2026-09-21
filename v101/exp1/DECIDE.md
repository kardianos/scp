# Experiment 1 — decision

**Written:** 2026-09-21, after `cas/out/succession.out`.
**Status:** this page closes the succession enumeration. It does
not edit experiment 2's gates. It does not edit a kernel.

Gates were locked in `PREREG.md` before the count.

---

## 1. What was asked

Which frequency-closed words on the locked alphabet reappear as a
mutual tail call, up to cyclic rotation, at cap 5, stably from
cap 4. A pass is one class or two. Three or more, or a set that
grows with the cap, is a zoo. Empty is a miss. A partial cycle
must not appear.

## 2. What came back

Meters are clean.

| meter | result |
|---|---|
| triangle \((3{:}1,\ 1{:}3,\ 1{:}1)\) closes | yes |
| lone \(2{:}1\) closes | no |
| known reverse pair is a fixed point | yes |
| partial-cycle failures inside \(C(5)\) | 0 |

| cap | classes |
|---:|---:|
| 4 | 476 |
| 5 | 12383 |
| growth | 11907 |

Joint-length histogram at cap 5: 4→3, 5→6, 6→41, 7→264, 8→846,
9→3108, 10→8115.

**Verdict: FAIL_ZOO.** \(C(5)\) properly contains \(C(4)\), and
\(\lvert C(5)\rvert \ge 3\). Both zoo conditions fire. Either
would have been enough.

The atom ledger was not used to rank. No ratio was introduced.

## 3. What the classes are

The first classes in `cas/out/classes.tsv` are unison symbols
wrapped around a cancelling pair such as \(1{:}2\ 2{:}1\). A
cancelling pair frequency-closes, and unison symbols keep the
product at 1, so every place the unisons can sit is another
necklace. Mutual landing is then common. The cap manufactures
new places to sit. That is the growth from 476 to 12383.

The known third-triangle pair is in the set. It is one class
among thousands, not a selector.

## 4. What is now true, and what is not

**True:**

1. On this alphabet, with rotation and swap as the only
   equivalences, mutual landing does not pick a succession.
2. The set is not stable under the pre-registered refinement.
3. Non-closing words were excluded. The structural slit control
   holds. Fringe loci were not recomputed.

**Not true:**

- A mass succession was found.
- The triangle is the unique landing.
- A ratio could be added to keep one class. The pre-registration
  forbids that, and this page does not do it.

## 5. How to disagree

To reverse the zoo: name an equivalence that is already written
in `PREREG.md` §2 and under which \(C(5)=C(4)\) with size 1 or 2.
Rotation and swap are already applied. An equivalence added after
this count is a new experiment.
