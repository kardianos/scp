# Experiment 1 — pre-registration

**Written:** 2026-09-21, before the enumeration.
**Status:** commitment. Thresholds here are not revised after the
counts exist. `DECIDE.md` applies them.

---

## 1. The objects (locked)

Alphabet, the intervals named in the corrected monodromy checks
(`v89/CHARGE.md` §7.1, `v89/KIMI_REVIEW_2026-08-02.md` §2.2). Not a
scan.

| symbol | \(p{:}q\) |
|---|---|
| unison | \(1{:}1\) |
| octave | \(2{:}1\), \(1{:}2\) |
| third | \(3{:}1\), \(1{:}3\) |
| fifth | \(3{:}2\), \(2{:}3\) |

A **word** is a finite sequence of these symbols, length from 2
through \(L\). A word of length 1 is not a cycle.

A word **frequency-closes** when \(\prod (p_i/q_i) = 1\) exactly,
as a rational. The zero branch of the corrected composition then
has monodromy 0, which is an integer multiple of \(2\pi\). That is
a completed cycle. A word that does not frequency-close is a
partial cycle: it converts nothing, books no atoms, and is deleted.
It cannot be a fixed point. Fractional branches on a word that does
close are charge, not a second closure test. They are not scored
here.

A **necklace** is a word up to cyclic rotation. The canonical
representative is the lexicographic minimum among rotations.
Symbols compare in the alphabet order above.

A necklace is **painted** when every symbol is equal. A painted
necklace is a seeded ring. It is removed before pairing. Unison
words are painted.

## 2. The move (locked)

Two necklaces \(N_1, N_2\) (the same necklace allowed twice) **land
mutually** when some rotation \(A\) of \(N_1\) and some rotation
\(B\) of \(N_2\) satisfy

\[
\mathrm{last}(A) = \mathrm{first}(B)
\quad\text{and}\quad
\mathrm{last}(B) = \mathrm{first}(A).
\]

The move deletes both words. Each books two action atoms. The fresh
pair is the pair just landed: the last act of each is the first act
of the other. A **fixed point** is an unordered pair
\(\{N_1, N_2\}\) that lands mutually. The atom count is recorded
and is not a rank. No coefficient is introduced to prefer one pair.

Equivalence is that unordered pair of canonical necklaces. Swap is
the same pair. Reversal is not identified with rotation.

## 3. Bounds (locked)

Production cap \(L = 5\). Refinement cap \(L = 4\). Every word of
length \(2..L\) is enumerated. No sampling.

## 4. Meters

Instrument: `cas/succession.py`. A dirty meter makes the verdict
`INCONCLUSIVE`, which is not a pass.

| meter | clean if |
|---|---|
| triangle | \((3{:}1,\ 1{:}3,\ 1{:}1)\) frequency-closes, and \((2{:}1)\) does not |
| known pair | the unordered pair of necklaces \((3{:}1,\ 1{:}3,\ 1{:}1)\) and \((1{:}1,\ 1{:}3,\ 3{:}1)\) is a fixed point at \(L=5\). If the filters drop it, the instrument is wrong |
| partial | every word that appears in a fixed point frequency-closes. The count of failures is 0. This is the slit control this experiment can actually run: a non-closing residue is not a mass. Fringe loci are not recomputed |

## 5. Verdict

Let \(C(L)\) be the set of fixed-point classes at cap \(L\).

| verdict | condition |
|---|---|
| **PASS** | \(C(5) = C(4)\) and \(\lvert C(5)\rvert \in \{1, 2\}\) |
| **FAIL_EMPTY** | \(C(5)\) is empty |
| **FAIL_ZOO** | \(\lvert C(5)\rvert \ge 3\), or \(C(5)\) properly contains \(C(4)\) |
| **INCONCLUSIVE** | a meter in §4 is dirty |

A pass of two classes is a short ladder. It is not a licence to
insert a ratio and keep one of them. A zoo is a miss, including a
zoo that appears only when the cap grows.

## 6. What this experiment does not claim

A pass is one stable succession, or two. It is not a nucleon, not
a cell count, and not a hedgehog. A fail closes the tail-call
enumeration at this alphabet and these caps. It does not edit
experiment 2.
