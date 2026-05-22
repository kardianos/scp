#!/usr/bin/env python3
"""
v59/first_experiments/04_muon_electron_ratio_scan.py

Scan candidate substrates for the muon-to-electron mass ratio
TYCHO_TABLE T1.5:  m_mu/m_e = 206.7682830(46).

Strategy: ratios of two natural invariants of a substrate. Lower bar than
Koide (only two invariants needed, no triple alignment). Any substrate
that produces a ratio of ~206.77 with zero or one free parameter is a
candidate.

Tolerance: experimental precision is 5e-10. We adopt:
  TIGHT       <  1e-4  (suggestive of a genuine match)
  CLOSE       <  1e-2  (worth noting)
  COINCIDENCE >= 1e-2  (random)

Categories tested:
  A. Famous near-miss numerologies (Nambu 3*alpha^-1/2, etc.)
  B. exp(rationals), log_10(rationals), log_e(rationals)
  C. Ratios of finite group orders and irrep dimensions
  D. Combinatorial sequence values: factorials, binomials, Catalan numbers
  E. Sphere kissing numbers, n-dim ball/cube volume ratios
  F. alpha-based combinations
  G. Ratios involving Riemann zeta(s) and known constants
"""

import math
from itertools import product

TARGET = 206.7682830  # m_mu/m_e, PDG 2024
ALPHA_INV = 137.035999084  # fine structure constant^-1

TIGHT = 1e-4
CLOSE = 1e-2

def rel_diff(val):
    return (val - TARGET) / TARGET

def classify(val):
    rd = rel_diff(val)
    if abs(rd) < TIGHT:
        return rd, "*** TIGHT HIT ***"
    elif abs(rd) < CLOSE:
        return rd, "(close)"
    return rd, None

def show(name, val, force=False):
    rd, flag = classify(val)
    if flag or force:
        print(f"  {name:>50s}: {val:>18.10f}   rel.Δ = {rd:>+11.4e}   {flag or ''}")


print("="*72)
print(f"TARGET:  m_mu/m_e = {TARGET}")
print("="*72)


# =========================================================================
# Part A: Famous literature near-misses
# =========================================================================
print()
print("-"*72)
print("Part A: Famous literature near-misses")
print("-"*72)

candidates = [
    ("3 * alpha^-1 / 2  (Nambu 1952)",    3 * ALPHA_INV / 2),
    ("3 * alpha^-1 / 2 - 0.5",            3 * ALPHA_INV / 2 - 0.5),
    ("3 * alpha^-1 / 2 + 1.214",          3 * ALPHA_INV / 2 + 1.214),
    ("alpha^-1 * 3/2 * (1 + 5/750)",      ALPHA_INV * 1.5 * (1 + 5/750)),
    ("alpha^-1 + 69.7",                   ALPHA_INV + 69.7),
    ("(alpha^-1)^(5/4)",                  ALPHA_INV ** (5/4)),
    ("9 * sqrt(alpha^-1) * (small)",      9 * math.sqrt(ALPHA_INV) * 1.962),
    ("2 * pi * alpha^-1 / (2*pi-2.5)",    2 * math.pi * ALPHA_INV / (2*math.pi - 2.5)),
    ("(pi^2 + e) * 7^2 - ?",              (math.pi**2 + math.e) * 49 / 2.86),
    ("207 (= 9 * 23)",                    207),
    ("207 - 0.2317",                      207 - 0.2317),
    ("4^4 - 2^(5.66)",                    256 - 2**5.66),
    ("206.7682830 (target)",              TARGET),
]
for name, val in candidates:
    show(name, val, force=True)


# =========================================================================
# Part B: exp / log of small rationals
# =========================================================================
print()
print("-"*72)
print("Part B: exp(p/q) and log_n(integer) sweeps")
print("-"*72)

# log(206.77) = 5.3315
# So exp(p/q) = 206.77 means p/q = 5.3315
print(f"\nFor reference: log_e({TARGET}) = {math.log(TARGET):.10f}")
print(f"               log_10({TARGET}) = {math.log10(TARGET):.10f}")
print(f"               log_2({TARGET}) = {math.log2(TARGET):.10f}")

# Scan rationals
hits = []
for q in range(1, 50):
    for p in range(1, 300):
        # exp(p/q)
        if abs(p/q - 5.33) < 0.5:  # range filter for speed
            for base, bname in [(math.e, "e"), (2, "2"), (math.pi, "π")]:
                val = base ** (p/q)
                rd, flag = classify(val)
                if flag and flag == "*** TIGHT HIT ***":
                    hits.append((abs(rd), f"{bname}^({p}/{q})", val, flag))
hits.sort()
print(f"\nTight hits for natural-base^(rational):")
if hits:
    for d, name, val, flag in hits[:15]:
        print(f"  {name:>20s} = {val:>15.10f}  rel.Δ = {rd:+.4e}  {flag}")
else:
    print("  (none with p,q in scan range)")

# Closer scan around log = 5.33
print(f"\nFine scan: bases and exponents near log({TARGET}):")
# Look for log(206.77) in the form ln(p) + q or similar
# Compute log_e(206.77) and see if it equals any simple expression
L = math.log(TARGET)
print(f"  log_e({TARGET}) = {L:.10f}")
candidates = [
    ("16/3",                       16/3),
    ("ln(207)",                    math.log(207)),
    ("5 + 1/3",                    5 + 1/3),
    ("3 + e/(pi-1)",               3 + math.e/(math.pi-1)),
    ("4 + ln(4)",                  4 + math.log(4)),
    ("(7+pi)/(2 + 1/pi)",          (7+math.pi)/(2 + 1/math.pi)),
    ("pi + 2",                     math.pi + 2),
    ("ln(2)*7.69",                 math.log(2)*7.69),
    ("phi^3.5",                    ((1+math.sqrt(5))/2)**3.5),
]
print(f"  Candidates for log_e of target:")
for name, val in candidates:
    delta = val - L
    print(f"    {name:>25s} = {val:.10f}  Δ_log = {delta:+.4e}")
print(f"  Inverse check: exp(16/3) = {math.exp(16/3):.4f}  (target {TARGET:.4f})")


# =========================================================================
# Part C: Finite group orders and irrep dimension ratios
# =========================================================================
print()
print("-"*72)
print("Part C: Finite group orders and irrep ratios")
print("-"*72)

# Orders of small finite simple groups and selected representations:
groups = {
    "A_5":               60,
    "PSL(2,7)":          168,
    "A_6":               360,
    "PSL(2,11)":         660,
    "M_11":              7920,
    "A_7":               2520,
    "PSL(2,13)":         1092,
    "S_5":               120,
    "S_6":               720,
    "S_7":               5040,
    "S_8":               40320,
    "A_8":               20160,
    "SL(2,3)":           24,
    "GL(2,3)":           48,
    "PGL(2,7)":          336,
    "Sz(8)":             29120,
    "He (Held)":         4030387200,
    "Cube group":        48,
    "Icosahedral":       60,
    "Z_207":             207,
    "Z_2 x Z_3 x Z_23":  138,
    "S_3 x S_4":         144,
    "Z_9 x Z_23":        207,
}

hits = []
for name1, n1 in groups.items():
    for name2, n2 in groups.items():
        if name1 == name2:
            continue
        ratio = n1 / n2
        rd, flag = classify(ratio)
        if flag:
            hits.append((abs(rd), name1, name2, ratio, flag))
hits.sort()
print(f"\nClose hits in group-order ratios:")
if hits:
    for d, n1, n2, r, flag in hits[:10]:
        print(f"  {n1:>15s} / {n2:>15s} = {r:>15.6f}  rel.Δ = {(r-TARGET)/TARGET:+.4e}  {flag}")
else:
    print("  (no group-order ratio within 1%)")


# =========================================================================
# Part D: Sphere kissing numbers, n-ball volume ratios
# =========================================================================
print()
print("-"*72)
print("Part D: Sphere kissing numbers and ball/cube ratios")
print("-"*72)

kissing = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 272,
           10: 336, 11: 438, 12: 756, 16: 4320, 24: 196560}

print("\nKissing numbers and their pairwise ratios near target:")
hits = []
for n1, k1 in kissing.items():
    for n2, k2 in kissing.items():
        if n1 == n2:
            continue
        ratio = k1 / k2
        rd, flag = classify(ratio)
        if flag:
            hits.append((abs(rd), n1, n2, k1, k2, ratio, flag))
hits.sort()
if hits:
    for d, n1, n2, k1, k2, r, flag in hits[:10]:
        print(f"  K_{n1}/K_{n2} = {k1}/{k2} = {r:.6f}   rel.Δ = {(r-TARGET)/TARGET:+.4e}   {flag}")
else:
    print("  (no kissing-number ratio within 1%)")

# n-ball volume V_n = pi^(n/2) / Gamma(n/2 + 1)
print(f"\nUnit n-ball volumes V_n:")
for n in range(1, 25):
    V = math.pi**(n/2) / math.gamma(n/2 + 1)
    show(f"V_{n}", V, force=False)

# Ratios of n-ball volumes:
print("\nRatios V_n/V_m near target:")
hits = []
for n1 in range(1, 26):
    for n2 in range(1, 26):
        if n1 == n2:
            continue
        V1 = math.pi**(n1/2) / math.gamma(n1/2 + 1)
        V2 = math.pi**(n2/2) / math.gamma(n2/2 + 1)
        ratio = V1 / V2
        rd, flag = classify(ratio)
        if flag:
            hits.append((abs(rd), n1, n2, ratio, flag))
hits.sort()
if hits:
    for d, n1, n2, r, flag in hits[:10]:
        print(f"  V_{n1}/V_{n2} = {r:.6f}   rel.Δ = {(r-TARGET)/TARGET:+.4e}   {flag}")
else:
    print("  (no V_n/V_m ratio within 1%)")


# =========================================================================
# Part E: Combinatorial sequences
# =========================================================================
print()
print("-"*72)
print("Part E: Factorial, binomial, Catalan ratios")
print("-"*72)

print("\nFactorials n!:")
fac = [math.factorial(n) for n in range(15)]
print(f"  {[ (n, fac[n]) for n in range(15) ]}")
print()
print("Ratios n!/k! near target:")
hits = []
for n in range(15):
    for k in range(15):
        if n == k or fac[k] == 0:
            continue
        r = fac[n] / fac[k]
        rd, flag = classify(r)
        if flag:
            hits.append((abs(rd), n, k, r, flag))
hits.sort()
if hits:
    for d, n, k, r, flag in hits[:10]:
        print(f"  {n}!/{k}! = {r:.4f}   rel.Δ = {(r-TARGET)/TARGET:+.4e}  {flag}")

print()
print("Binomials C(n,k) near target:")
hits = []
for n in range(1, 30):
    for k in range(n+1):
        c = math.comb(n, k)
        rd, flag = classify(c)
        if flag:
            hits.append((abs(rd), n, k, c, flag))
hits.sort()
if hits:
    for d, n, k, c, flag in hits[:10]:
        print(f"  C({n},{k}) = {c}   rel.Δ = {(c-TARGET)/TARGET:+.4e}   {flag}")

# Catalan numbers C_n = C(2n, n)/(n+1)
print()
print("Catalan numbers C_n:")
cat = [math.comb(2*n, n)//(n+1) for n in range(15)]
print(f"  {[(n, cat[n]) for n in range(15)]}")
hits = []
for n in range(15):
    rd, flag = classify(cat[n])
    if flag:
        hits.append((abs(rd), n, cat[n], flag))
if hits:
    for d, n, c, flag in hits:
        print(f"  C_{n} = {c}   rel.Δ = {(c-TARGET)/TARGET:+.4e}   {flag}")

# (n!)^2 / k! style
print()
print("(n!)^2 / m! ratios:")
hits = []
for n in range(15):
    for m in range(20):
        mfac = math.factorial(m)
        if mfac == 0:
            continue
        r = (fac[n] ** 2) / mfac
        rd, flag = classify(r)
        if flag:
            hits.append((abs(rd), n, m, r, flag))
hits.sort()
if hits:
    for d, n, m, r, flag in hits[:10]:
        print(f"  ({n}!)^2/{m}! = {r:.6f}   rel.Δ = {(r-TARGET)/TARGET:+.4e}   {flag}")
else:
    print("  (none within 1%)")


# =========================================================================
# Part F: alpha-based combinations
# =========================================================================
print()
print("-"*72)
print("Part F: alpha-based combinations")
print("-"*72)

print(f"\nalpha^(-1) = {ALPHA_INV}")
candidates = [
    ("alpha^(-1) * 3/2",                ALPHA_INV * 3/2),
    ("alpha^(-1) * 2/sqrt(2) * c1",     ALPHA_INV * math.sqrt(2)),
    ("alpha^(-1) * pi/sqrt(pi^2-2)",    ALPHA_INV * math.pi / math.sqrt(math.pi**2 - 2)),
    ("alpha^(-1) + 70 - 1/3",           ALPHA_INV + 70 - 1/3),
    ("alpha^(-1) * (1 + 1/2 + ...)",    ALPHA_INV * 1.50885),
    ("4*pi*alpha^(-1)/8.36",            4 * math.pi * ALPHA_INV / 8.36),
    ("alpha^(-1) * 24/(4*pi)",          ALPHA_INV * 24/(4*math.pi)),
]
for name, val in candidates:
    show(name, val, force=True)


# =========================================================================
# Summary
# =========================================================================
print()
print("="*72)
print("Summary")
print("="*72)
print(f"\nTarget: {TARGET}")
print(f"Tightest experimental precision: 5e-10 relative")
print(f"This scan's TIGHT tolerance: {TIGHT} relative")
print(f"This scan's CLOSE tolerance: {CLOSE} relative")
print()
print("See output above for any TIGHT or CLOSE hits.")
print()
print("Notes:")
print("  - log_e(target) = {:.6f}, not a simple rational".format(math.log(TARGET)))
print("  - target / alpha^-1 = {:.6f}".format(TARGET / ALPHA_INV))
print("  - target / (3/2 * alpha^-1) = {:.6f}".format(TARGET / (1.5 * ALPHA_INV)))
print("  - target as 9 * 23 = 207, off by 0.112%")
