#!/usr/bin/env python3
"""v101 experiment 1 — mutual tail-call fixed points.

Thresholds live in ../PREREG.md. This script does not rank classes
and does not introduce a ratio.

Run: python3 v101/exp1/cas/succession.py
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

# Alphabet order is the comparison order for necklaces.
ALPH = (
    (1, 1),
    (2, 1),
    (1, 2),
    (3, 1),
    (1, 3),
    (3, 2),
    (2, 3),
)

TRIANGLE = ((3, 1), (1, 3), (1, 1))
TRIANGLE_REV = ((1, 1), (1, 3), (3, 1))


def closes(word: tuple[tuple[int, int], ...]) -> bool:
    prod = Fraction(1)
    for p, q in word:
        prod *= Fraction(p, q)
    return prod == 1


def necklace(word: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    rots = [word[i:] + word[:i] for i in range(len(word))]
    return min(rots)


def painted(word: tuple[tuple[int, int], ...]) -> bool:
    return all(sym == word[0] for sym in word)


def lands(a: tuple[tuple[int, int], ...], b: tuple[tuple[int, int], ...]) -> bool:
    """Some rotations meet last(A)=first(B) and last(B)=first(A)."""
    brots = [b[i:] + b[:i] for i in range(len(b))]
    for i in range(len(a)):
        ar = a[i:] + a[:i]
        for br in brots:
            if ar[-1] == br[0] and br[-1] == ar[0]:
                return True
    return False


def classes_at(cap: int) -> set[tuple]:
    closing: list[tuple[tuple[int, int], ...]] = []
    for length in range(2, cap + 1):
        for raw in product(range(len(ALPH)), repeat=length):
            word = tuple(ALPH[i] for i in raw)
            if closes(word) and not painted(word):
                closing.append(word)

    # One canonical word per necklace.
    necks: dict[tuple, tuple] = {}
    for word in closing:
        key = necklace(word)
        necks.setdefault(key, key)
    keys = list(necks)

    # Index rotations by (first, last) so mutual landing is a lookup.
    bucket: dict[tuple, list[tuple]] = defaultdict(list)
    for key in keys:
        for i in range(len(key)):
            rot = key[i:] + key[:i]
            bucket[(rot[0], rot[-1])].append(key)

    found: set[tuple] = set()
    for a in keys:
        for i in range(len(a)):
            rot = a[i:] + a[:i]
            want_first = rot[-1]
            want_last = rot[0]
            for b in bucket.get((want_first, want_last), ()):
                # bucket stores necklaces that HAVE a rotation with that
                # (first, last). Confirm the junction against a itself.
                if lands(a, b):
                    found.add(tuple(sorted((a, b))))
    return found


def fmt(word: tuple[tuple[int, int], ...]) -> str:
    return " ".join(f"{p}:{q}" for p, q in word)


def main() -> None:
    tri_ok = closes(TRIANGLE)
    oct_closed = closes(((2, 1),))
    print(f"METER triangle closes {tri_ok}")
    print(f"METER lone_octave closes {oct_closed}")

    c4 = classes_at(4)
    c5 = classes_at(5)
    known = tuple(sorted((necklace(TRIANGLE), necklace(TRIANGLE_REV))))
    known_ok = known in c5
    print(f"METER known_pair {known_ok}")
    print(f"COUNT L4 {len(c4)}")
    print(f"COUNT L5 {len(c5)}")
    print(f"GROWTH {len(c5 - c4)}")

    # Partial-cycle control: every stored necklace was built from closing words.
    partial_fail = 0
    for pair in c5:
        for word in pair:
            if not closes(word):
                partial_fail += 1
    print(f"METER partial_failures {partial_fail}")

    meters_clean = tri_ok and (not oct_closed) and known_ok and partial_fail == 0
    if not meters_clean:
        verdict = "INCONCLUSIVE"
    elif len(c5) == 0:
        verdict = "FAIL_EMPTY"
    elif c5 == c4 and len(c5) in (1, 2):
        verdict = "PASS"
    else:
        verdict = "FAIL_ZOO"
    print(f"EXP1 {verdict}")

    # Length histogram of the production set. A class contributes its
    # two necklace lengths (one twice when the pair is a self-landing).
    hist: dict[int, int] = defaultdict(int)
    for pair in c5:
        hist[sum(len(w) for w in pair)] += 1
    hist_s = " ".join(f"{k}:{hist[k]}" for k in sorted(hist))
    print(f"HIST joint_length {hist_s}")

    lines = ["class\tpair\n"]
    for i, pair in enumerate(sorted(c5)[:40]):
        shown = " || ".join(fmt(w) for w in pair)
        lines.append(f"{i}\t{shown}\n")
        if i < 15:
            print(f"CLASS {i} {shown}")
    if len(c5) > 40:
        print(f"CLASS truncated {len(c5) - 40}")
    (OUT / "classes.tsv").write_text("".join(lines))
    print(f"WROTE {OUT / 'classes.tsv'}")


if __name__ == "__main__":
    main()
