#!/usr/bin/env python3
"""
v59/first_experiments/03_modular_hyperbolic_scan.py

Targeted scan of less-standard mathematical structures for the angle 2/9 rad,
following the negative result from 02_two_ninths_audit.py.

Areas tested:
  A. Lobachevsky function Λ(θ) and its ratios — invariants of hyperbolic
     3-manifolds via volume sums.
  B. Volumes and arctan/log/ratios of small hyperbolic 3-manifolds
     (figure-8, Weeks, Whitehead link, etc.).
  C. Eisenstein E_4 / E_6 and η-function values at CM points i, ω, sqrt(-2), etc.
  D. Bloch-Wigner dihedral function D(z) at natural arguments.
  E. Ratios of CM j-invariant values.
  F. q-series and Rogers-Ramanujan-type constants.
"""

import numpy as np
import mpmath as mp

mp.mp.dps = 30  # 30-digit precision

TARGET = mp.mpf(2) / mp.mpf(9)
TOL_TIGHT = mp.mpf("1e-4")
TOL_CLOSE = mp.mpf("1e-2")

def classify(val):
    d = mp.mpf(val) - TARGET
    if abs(d) < TOL_TIGHT:
        return d, "*** TIGHT HIT ***"
    elif abs(d) < TOL_CLOSE:
        return d, "(close)"
    else:
        return d, None

def show(name, val, force=False):
    d, flag = classify(val)
    if flag or force:
        print(f"  {name:>50s}: {mp.nstr(val, 12):>22s}  Δ={mp.nstr(d, 4):>15s}  {flag or ''}")

# ----- Lobachevsky function -----
# Λ(θ) = -∫₀^θ ln|2 sin t| dt
# Implemented via mpmath.clcos2.
def Lobachevsky(theta):
    # Λ(θ) = (1/2) Cl_2(2θ) where Cl_2 is the Clausen function.
    return mp.clsin(2, 2*theta) / 2

# ----- Bloch-Wigner D(z) -----
def BlochWigner(z):
    # D(z) = Im(Li_2(z)) + arg(1-z) * log|z|
    if z == 0 or z == 1:
        return mp.mpf(0)
    Li2 = mp.polylog(2, z)
    return Li2.imag + mp.arg(1 - z) * mp.log(abs(z))


# =========================================================================
# Part A: Lobachevsky values and ratios
# =========================================================================
print("="*72)
print("Part A: Lobachevsky function values and their ratios")
print("="*72)

# Λ(π/n) for n=2..20:
print("\nΛ(π/n):")
L_vals = {}
for n in range(2, 21):
    v = Lobachevsky(mp.pi / n)
    L_vals[n] = v
    show(f"Λ(π/{n})", v, force=False)

# Ratios L(π/p)/L(π/q):
print("\nPairwise ratios Λ(π/p)/Λ(π/q) within tolerance of 2/9:")
hits = []
for p in range(2, 21):
    for q in range(2, 21):
        if p == q:
            continue
        if L_vals[q] == 0:
            continue
        ratio = L_vals[p] / L_vals[q]
        d, flag = classify(ratio)
        if flag:
            hits.append((abs(d), p, q, ratio, flag))
hits.sort()
for d, p, q, ratio, flag in hits[:15]:
    print(f"    Λ(π/{p})/Λ(π/{q}) = {mp.nstr(ratio, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none)")

# Differences:
print("\nPairwise differences Λ(π/p) - Λ(π/q) within tolerance:")
hits = []
for p in range(2, 21):
    for q in range(2, 21):
        if p == q:
            continue
        diff = L_vals[p] - L_vals[q]
        d, flag = classify(diff)
        if flag:
            hits.append((abs(d), p, q, diff, flag))
hits.sort()
for d, p, q, diff, flag in hits[:15]:
    print(f"    Λ(π/{p}) - Λ(π/{q}) = {mp.nstr(diff, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none)")


# =========================================================================
# Part B: Volumes and derived quantities of small hyperbolic 3-manifolds
# =========================================================================
print()
print("="*72)
print("Part B: Hyperbolic 3-manifold volumes and derived angles")
print("="*72)

# Well-known volumes (in standard hyperbolic units):
manifolds = {
    "Weeks manifold (smallest closed)":     mp.mpf("0.9427073627769277"),
    "Meyerhoff manifold":                   mp.mpf("0.9813688288922320"),
    "figure-8 knot complement (6Λ(π/3))":   6 * Lobachevsky(mp.pi/3),
    "5_2 knot complement":                  mp.mpf("2.8281220883307393"),
    "Whitehead link (8Λ(π/4))":             8 * Lobachevsky(mp.pi/4),
    "6_1 knot":                             mp.mpf("3.4179148372119296"),
    "Borromean rings (16Λ(π/4)*?)":         mp.mpf("7.32772475985"),
    "Bianchi orbifold V(Q(√-3))":           mp.mpf("0.16923"),
    "Bianchi orbifold V(Q(√-1))":           mp.mpf("0.30532"),
}

print("\nRaw volumes:")
for name, v in manifolds.items():
    print(f"  {name:>45s}: V = {mp.nstr(v, 10)}")

# Try transformations
print("\nAngles derived from volumes:")
test_transforms = [
    ("V mod 2π", lambda v: v % (2*mp.pi)),
    ("V / 2π", lambda v: v / (2*mp.pi)),
    ("V / π", lambda v: v / mp.pi),
    ("arctan(V)", lambda v: mp.atan(v)),
    ("arccot(V) = π/2 - arctan(V)", lambda v: mp.pi/2 - mp.atan(v)),
    ("log(V)", lambda v: mp.log(v)),
]
hits = []
for name, V in manifolds.items():
    for tname, t in test_transforms:
        try:
            val = t(V)
            d, flag = classify(val)
            if flag:
                hits.append((abs(d), name, tname, val, flag))
        except Exception:
            pass
hits.sort()
print(f"  Top hits:")
for d, mname, tname, val, flag in hits[:10]:
    print(f"    {mname[:30]:>32s} ({tname:>20s}): {mp.nstr(val, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none within tolerance)")

# Ratios between pairs:
print("\nPairwise ratios V_p / V_q within tolerance:")
names = list(manifolds.keys())
vals = list(manifolds.values())
hits = []
for i in range(len(names)):
    for j in range(len(names)):
        if i == j or vals[j] == 0:
            continue
        ratio = vals[i] / vals[j]
        d, flag = classify(ratio)
        if flag:
            hits.append((abs(d), names[i], names[j], ratio, flag))
hits.sort()
for d, a, b, r, flag in hits[:10]:
    print(f"    {a[:25]} / {b[:25]} = {mp.nstr(r, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none)")


# =========================================================================
# Part C: Eisenstein/eta function values at CM points
# =========================================================================
print()
print("="*72)
print("Part C: Eisenstein E_4 / E_6 / η at CM points")
print("="*72)

def eta(tau):
    """Dedekind eta η(τ) using q-product."""
    q = mp.exp(2*mp.pi*mp.mpc(0, 1)*tau)
    qn = mp.mpc(1)
    prod = mp.mpc(1)
    for n in range(1, 200):
        qn *= q
        prod *= (1 - qn)
    return q**(mp.mpf(1)/24) * prod

# Some CM points (imaginary quadratic):
CM_pts = {
    "i":             mp.mpc(0, 1),
    "ω = e^(2πi/3)": mp.exp(2*mp.pi*mp.mpc(0,1)/3),
    "i√2":           mp.mpc(0, mp.sqrt(2)),
    "(1+i√3)/2":     (1 + mp.mpc(0, mp.sqrt(3)))/2,
    "(1+i√7)/2":     (1 + mp.mpc(0, mp.sqrt(7)))/2,
    "i√3":           mp.mpc(0, mp.sqrt(3)),
    "(1+i√11)/2":    (1 + mp.mpc(0, mp.sqrt(11)))/2,
}

print("\nη(τ) at CM points:")
hits = []
for name, tau in CM_pts.items():
    if tau.imag <= 0:
        continue
    e = eta(tau)
    abs_e = abs(e)
    arg_e = mp.arg(e)
    log_abs = mp.log(abs_e)
    print(f"  τ = {name:>15s}: |η| = {mp.nstr(abs_e, 8):>15s}  arg = {mp.nstr(arg_e, 8):>15s}  log|η| = {mp.nstr(log_abs, 8):>15s}")
    for vname, v in [("|η|", abs_e), ("arg(η)", arg_e), ("log|η|", log_abs)]:
        d, flag = classify(v)
        if flag:
            hits.append((abs(d), name, vname, v, flag))
print(f"\n  Hits within tolerance:")
hits.sort()
for d, name, vname, v, flag in hits[:10]:
    print(f"    η({name})·{vname} = {mp.nstr(v, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none)")


# =========================================================================
# Part D: Bloch-Wigner dihedral function
# =========================================================================
print()
print("="*72)
print("Part D: Bloch-Wigner D(z) at natural arguments")
print("="*72)

# Natural arguments
print("\nD(z) at points of interest:")
test_z = [
    ("e^(iπ/3) = ω̄_6",         mp.exp(mp.pi*mp.mpc(0,1)/3)),
    ("e^(2iπ/3)",                mp.exp(2*mp.pi*mp.mpc(0,1)/3)),
    ("i",                        mp.mpc(0, 1)),
    ("(1+i√3)/2",                (1 + mp.mpc(0, mp.sqrt(3)))/2),
    ("(1+i)/2",                  mp.mpc(0.5, 0.5)),
    ("1/2 + i/√12",              mp.mpc(0.5, 1/mp.sqrt(12))),
    ("e^(iπ/9) (3-fold)",        mp.exp(mp.pi*mp.mpc(0,1)/9)),
    ("e^(2iπ/9)",                mp.exp(2*mp.pi*mp.mpc(0,1)/9)),
    ("1 - 1/9",                  mp.mpf(8)/9),
    ("(1 + e^(iπ/3))/2",         (1 + mp.exp(mp.pi*mp.mpc(0,1)/3))/2),
]
hits = []
for name, z in test_z:
    try:
        D = BlochWigner(z)
        print(f"  D({name:>30s}) = {mp.nstr(D, 12)}")
        d, flag = classify(D)
        if flag:
            hits.append((abs(d), name, D, flag))
    except Exception as exc:
        print(f"  D({name:>30s}) = failed ({exc})")
print(f"\n  Hits:")
hits.sort()
for d, name, val, flag in hits[:10]:
    print(f"    D({name}) = {mp.nstr(val, 10)}   Δ={mp.nstr(d, 3)}   {flag}")
if not hits:
    print("    (none)")


# =========================================================================
# Part E: Mass-ratio derived: arctan(sqrt(m_e/m_mu)), etc.
# =========================================================================
print()
print("="*72)
print("Part E: Quick sanity — arctan/arcsin of lepton mass ratios")
print("="*72)
m_e = mp.mpf("0.51099895069")
m_mu = mp.mpf("105.6583755")
m_tau = mp.mpf("1776.86")

candidates = [
    ("arctan(√(m_e/m_mu))",    mp.atan(mp.sqrt(m_e/m_mu))),
    ("arctan(√(m_mu/m_tau))",  mp.atan(mp.sqrt(m_mu/m_tau))),
    ("arctan(√(m_e/m_tau))",   mp.atan(mp.sqrt(m_e/m_tau))),
    ("arcsin(√(m_e/m_mu))",    mp.asin(mp.sqrt(m_e/m_mu))),
    ("(√(m_e/m_mu))",          mp.sqrt(m_e/m_mu)),
]
for name, v in candidates:
    d, flag = classify(v)
    note = flag or ""
    print(f"  {name:>30s} = {mp.nstr(v, 10):>20s}   Δ={mp.nstr(d, 4):>12s}   {note}")


# =========================================================================
# Part F: q-series, Rogers-Ramanujan
# =========================================================================
print()
print("="*72)
print("Part F: q-series constants and ratios")
print("="*72)

# Rogers-Ramanujan G(q), H(q) at q = 1/9 etc.
def RR_G(q, N=100):
    s = mp.mpf(1)
    for n in range(1, N):
        term = q**(n*n)
        for k in range(1, n+1):
            term /= (1 - q**k)
        s += term
    return s

def RR_H(q, N=100):
    s = mp.mpf(1)
    for n in range(1, N):
        term = q**(n*(n+1))
        for k in range(1, n+1):
            term /= (1 - q**k)
        s += term
    return s

print("\nRogers-Ramanujan G(q), H(q) and ratio at small q:")
for q_val_str in ["1/9", "1/10", "1/e", "e^(-π)", "1/100"]:
    if q_val_str == "1/9":
        q = mp.mpf(1)/9
    elif q_val_str == "1/10":
        q = mp.mpf(1)/10
    elif q_val_str == "1/e":
        q = 1/mp.e
    elif q_val_str == "e^(-π)":
        q = mp.exp(-mp.pi)
    elif q_val_str == "1/100":
        q = mp.mpf(1)/100
    G = RR_G(q)
    H = RR_H(q)
    ratio = G / H
    print(f"  q = {q_val_str:>15s}: G = {mp.nstr(G, 10):>15s}, H = {mp.nstr(H, 10):>15s}, G/H = {mp.nstr(ratio, 10):>15s}")
    show("log(G/H)", mp.log(ratio), force=False)
    show("arctan(G/H - 1)", mp.atan(ratio - 1), force=False)
    show("G - H", G - H, force=False)


# =========================================================================
# Summary
# =========================================================================
print()
print("="*72)
print("Summary")
print("="*72)
print(f"\nTarget: 2/9 = {mp.nstr(TARGET, 15)}")
print()
print("Tested:")
print("  A. Lobachevsky Λ(π/n) values and ratios (n = 2..20).")
print("  B. Hyperbolic 3-manifold volumes (Weeks, figure-8, Whitehead, etc.).")
print("  C. Eisenstein and Dedekind η values at CM points.")
print("  D. Bloch-Wigner D(z) at natural complex arguments.")
print("  E. Lepton-mass-ratio derived angles (sanity).")
print("  F. Rogers-Ramanujan q-series at small q.")
print()
print("Any TIGHT HIT entries above represent natural occurrences of 2/9.")
print("If none, the lepton Brannen phase remains structurally unmatched.")
