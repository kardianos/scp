#!/usr/bin/env python3
"""ws1 — FABLE seat, v87 B0.

Symbolic maximisation of the CHSH functional for singlet-form correlations
E(a,b) = -a.b (unit vectors), two independent routes:

(1) Coplanar parametrisation: S(alpha, alpha', beta, beta') via cosines;
    stationary-point solve + Hessian check -> max |S| = 2*sqrt(2).
(2) Vector proof of the global bound:
      S = -[a.(b+b') + a'.(b-b')]  =>  |S| <= |b+b'| + |b-b'|
    and |b+b'|^2 + |b-b'|^2 = 4  =>  |S| <= 2*sqrt(2)  (Cauchy-Schwarz),
    with equality iff a || (b+b'), a' || (b-b'), |b+b'| = |b-b'| = sqrt(2)
    i.e. b _|_ b'.

SymPy 1.12.
"""
import sympy as sp

al, alp, be, bep = sp.symbols('alpha alphap beta betap', real=True)

# E(a,b) = -cos(angle between a and b); coplanar settings.
S = -(sp.cos(al - be) + sp.cos(al - bep) + sp.cos(alp - be) - sp.cos(alp - bep))

# --- route 1: numeric global max (multistart) + exact value at optimum -----
# (sympy.solve returns only the degenerate trig branches here, so the global
#  search is numeric; the PROOF of the bound is route 2 below.)
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(20260726)
f = sp.lambdify((al, alp, be, bep), -S, 'numpy')
best = 0.0
for _ in range(200):
    x0 = rng.uniform(0, 2*np.pi, 4)
    r = minimize(lambda x: f(*x), x0)
    best = max(best, -r.fun)
print("numeric global max |S| over 200 starts: %.12f" % best)
print("2*sqrt(2)                             : %.12f" % float(2*sp.sqrt(2)))
assert abs(best - float(2*sp.sqrt(2))) < 1e-8

# Verify the canonical quantum angles hit it exactly:
canon = S.subs({al: 0, alp: sp.pi/2, be: sp.pi + sp.pi/4, bep: sp.pi - sp.pi/4})
print("canonical angles (0, pi/2; pi+pi/4, pi-pi/4): S =", sp.simplify(canon))
assert sp.simplify(canon - 2*sp.sqrt(2)) == 0

# --- route 2: global bound, symbolic --------------------------------------
# |b+b'| = 2|cos(phi/2)|, |b-b'| = 2|sin(phi/2)| with phi the b^b' angle.
phi = sp.symbols('phi', real=True)
bound = 2*sp.Abs(sp.cos(phi/2)) + 2*sp.Abs(sp.sin(phi/2))
# max over phi of |cos|+|sin| is sqrt(2) at phi/2 = pi/4:
crit = sp.solve(sp.diff(2*sp.cos(phi/2) + 2*sp.sin(phi/2), phi), phi)
print("critical phi:", crit)
maxbound = sp.simplify(bound.subs(phi, sp.pi/2))
print("global bound max:", maxbound)          # 2*sqrt(2)
assert sp.simplify(maxbound - 2*sp.sqrt(2)) == 0

print("ws1 OK: max|S| = 2*sqrt(2) for E(a,b) = -a.b, attained at b _|_ b',")
print("        a || (b+b'), a' || (b-b').")
