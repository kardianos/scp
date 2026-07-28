#!/usr/bin/env python3
"""G2: unrestricted measurement dependence realises the PR box, |S|=4.

(R)+(L): outcomes A=A(a,λ), B=B(b,λ). No (MI): ρ(λ|a,b) free.

Construction: λ labels four 'strategy atoms' {00,01,10,11} corresponding to
the four setting pairs. For each setting pair (a,b), put all mass on the λ
that forces A*B = +1 when (a,b)≠(1,1) and A*B = −1 when (a,b)=(1,1).
Responses remain local functions of (own setting, λ).
"""
import numpy as np

# Settings encoded 0/1. λ ∈ {0,1,2,3}.
# Predetermined local tables: for each λ, specify A0,A1,B0,B1 ∈ {±1}
# We only need, for each (a,b), some λ with A_a B_b = target.

# Explicit response tables (local):
# λ=0: A0=+1,A1=+1, B0=+1,B1=+1  → products: 00:+, 01:+, 10:+, 11:+
# λ=1: A0=+1,A1=+1, B0=+1,B1=-1  → 00:+, 01:-, 10:+, 11:-
# λ=2: A0=+1,A1=-1, B0=+1,B1=+1  → 00:+, 01:+, 10:-, 11:-
# λ=3: A0=+1,A1=-1, B0=+1,B1=-1  → 00:+, 01:-, 10:-, 11:+
# We need for PR: E00=E01=E10=+1, E11=-1.
# λ=1 gives E11=-1 and E00=+1, E10=+1, but E01=-1 — not all at once for one λ.
# So MD: different λ for different setting pairs.

A = {
    0: {0: +1, 1: +1},  # λ0
    1: {0: +1, 1: +1},  # λ1
    2: {0: +1, 1: -1},  # λ2
    3: {0: +1, 1: -1},  # λ3  -- wait need a λ with A0 B1 = +1 and variety
}

# Cleaner: for each setting pair, pick a dedicated λ with correct product.
# λ_ab determines only the relevant outcomes; define full tables:

# λ labeled by (a,b) pair it is "for":
# For (0,0): want A0*B0 = +1. Set A0=+1, B0=+1, and dummy A1=B1=+1
# For (0,1): want A0*B1 = +1. A0=+1, B1=+1
# For (1,0): want A1*B0 = +1. A1=+1, B0=+1
# For (1,1): want A1*B1 = -1. A1=+1, B1=-1

tables = {
    # λ: (A0, A1, B0, B1)
    "for_00": (+1, +1, +1, +1),  # A0*B0=+1
    "for_01": (+1, +1, +1, +1),  # A0*B1=+1
    "for_10": (+1, +1, +1, +1),  # A1*B0=+1
    "for_11": (+1, +1, +1, -1),  # A1*B1=-1
}

def product(table, a, b):
    A0, A1, B0, B1 = table
    A = A0 if a == 0 else A1
    B = B0 if b == 0 else B1
    return A * B

# Conditional: P(λ|a,b) = 1 on the matching atom
cond = {
    (0, 0): "for_00",
    (0, 1): "for_01",
    (1, 0): "for_10",
    (1, 1): "for_11",
}

E = {}
for a, b in [(0, 0), (0, 1), (1, 0), (1, 1)]:
    lam = cond[(a, b)]
    E[(a, b)] = product(tables[lam], a, b)

S = E[(0, 0)] + E[(0, 1)] + E[(1, 0)] - E[(1, 1)]
print("G2 correlators:", E)
print("G2 S =", S)
assert S == 4
print("G2 OK: unrestricted MD yields PR box S=4  [D]")
print("  Note: responses are local (R)+(L); only ρ(λ|a,b) carries nonlocality of correlation.")
print("  Door A without a bound on MD is unphysical (nature does not give S=4).")
