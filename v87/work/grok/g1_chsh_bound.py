#!/usr/bin/env python3
"""G1: CHSH |S|<=2 under (R)^(L)^(MI) by exhaustive check on response signs.

For any fixed lambda, A(a),A(a'),B(b),B(b') in {+1,-1}:
  s = A*B + A*B' + A'*B - A'*B'
and |s| = 2 always (exactly one of B+/-B' has abs 2).

Also confirm max of |E| form under independent average is 2.
"""
from itertools import product

def s(Aa, Aap, Bb, Bbp):
    return Aa*Bb + Aa*Bbp + Aap*Bb - Aap*Bbp

vals = []
for Aa, Aap, Bb, Bbp in product([+1, -1], repeat=4):
    vals.append(s(Aa, Aap, Bb, Bbp))

assert set(vals) == {-2, 2}, set(vals)
print("G1 OK: s(lambda) in {-2,+2} for all 16 deterministic assignments")
print("  count +2:", vals.count(2), "  count -2:", vals.count(-2))
print("  therefore |integral rho s| <= 2 for any probability rho  [D]")

# Algebraic identity check
# s = A(B+B') + A'(B-B'); |B+B'| and |B-B'| are (2,0) or (0,2)
for Bb, Bbp in product([+1, -1], repeat=2):
    sp, sm = Bb + Bbp, Bb - Bbp
    assert (abs(sp), abs(sm)) in [(2, 0), (0, 2)]
print("G1 OK: complementarity of (B+B') and (B-B')")
