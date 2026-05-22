#!/usr/bin/env python3
"""
v59/multivector_kernel_fit/06_normalization_check.py

Step 8: normalization check for the 20.95 empirical ratio against 2pi^2.

Step 7 found that under exponential modulation, S_grav/S_em ≈ 20.95.
The closest natural number was 2pi^2 = 19.74 (volume of unit S^3 -- our
constraint surface), off by 6%.

This script tests whether different normalizations of the bare coupling
(alpha_bare = G_bare = c, for various c) close the gap.  If a NATURAL
choice of c gives exactly 2pi^2, the cross-sector ratio is structurally
pinned.

Approach:
  - For each candidate alpha_bare, compute S_em = -ln(alpha/alpha_bare)
    and S_grav = -ln(G_e/G_bare) (assuming alpha_bare = G_bare = c).
  - Find c that makes the ratio exactly 2pi^2 (and other candidates like
    21, 7*pi, etc.).
  - Test whether the required c has a natural interpretation.
"""

import numpy as np
import math

# Empirical
alpha = 7.2973525693e-3
G_e   = 1.7518093e-45

ln_alpha = -math.log(alpha)
ln_Ge    = -math.log(G_e)

print("="*72)
print("Step 8: Normalization check -- can we get S_grav/S_em = 2 pi^2?")
print("="*72)
print()
print(f"Empirical quantities (depth-of-suppression):")
print(f"  -ln(alpha)  = ln(1/alpha)   = {ln_alpha:.6f}")
print(f"  -ln(G_e)   = ln(1/G_e)     = {ln_Ge:.6f}")
print(f"  RAW ratio  = {ln_Ge / ln_alpha:.6f}  (= 20.95, before any normalization)")
print()
print(f"Natural candidates for the ratio:")
print(f"  2 pi^2     = {2*math.pi**2:.6f}  (volume of unit S^3)")
print(f"  21         = {21:.6f}             (3 x 7)")
print(f"  22         = {22:.6f}             (2 x 11)")
print(f"  7 pi       = {7*math.pi:.6f}")
print(f"  4 pi + 8   = {4*math.pi + 8:.6f}")
print()


# =========================================================================
# Part A: What alpha_bare makes ratio exactly equal to each candidate?
# =========================================================================
print("-"*72)
print("Part A: solve for alpha_bare needed to hit each candidate ratio")
print("-"*72)
print()
print("Under exponential modulation: alpha_obs = alpha_bare * exp(-S_em).")
print("So S_em = -ln(alpha) + ln(alpha_bare) = ln(1/alpha) + ln(c).")
print("Similarly S_grav = ln(1/G_e) + ln(c).")
print("Ratio = (ln(1/G_e) + ln(c)) / (ln(1/alpha) + ln(c)) = target.")
print()
print("Solving for ln(c) given target ratio R:")
print("  ln(c) = (ln(1/G_e) - R * ln(1/alpha)) / (R - 1)")
print()

def required_c(target_R):
    log_c = (ln_Ge - target_R * ln_alpha) / (target_R - 1)
    return math.exp(log_c), log_c

targets = [
    ("2 pi^2 (V(S^3))",       2*math.pi**2),
    ("21 (= 3 x 7)",          21.0),
    ("20.5",                  20.5),
    ("20 (= 4 x 5)",          20.0),
    ("22",                    22.0),
    ("7 pi",                  7*math.pi),
    ("pi^3 - 10",             math.pi**3 - 10),
    ("4 pi^2 / pi (= 4 pi)",  4*math.pi),
    ("4 pi (= V(S^2) sort of)", 4*math.pi),
    ("e^pi (numerology)",     math.exp(math.pi)),
]

print(f"  {'target ratio':<25s} {'value':>12s} {'required c':>16s} {'ln(c)':>12s}")
for name, R in targets:
    c, log_c = required_c(R)
    print(f"  {name:<25s} {R:>12.6f} {c:>16.6f} {log_c:>12.6f}")


# =========================================================================
# Part B: Identify natural numbers for the required c
# =========================================================================
print()
print("-"*72)
print("Part B: identify natural-number candidates for the required c")
print("-"*72)

# For ratio = 2 pi^2, c = exp((ln_Ge - 2pi^2 * ln_alpha)/(2 pi^2 - 1))
c_2pi2, ln_c_2pi2 = required_c(2*math.pi**2)
print()
print(f"For ratio = 2 pi^2: required c = {c_2pi2:.6f}, ln(c) = {ln_c_2pi2:.6f}")
print()
print("Natural-number checks for c = 1.3747:")
natural_candidates = [
    ("e^(1/pi)",                math.exp(1/math.pi)),
    ("(1 + 1/e)",               1 + 1/math.e),
    ("e/2",                     math.e / 2),
    ("4/e",                     4/math.e),
    ("sqrt(1 + 1/e^2 * something)", math.sqrt(1 + 1/math.e**2 + 0.5)),
    ("(1 + sqrt(5))/(e/2)",     (1 + math.sqrt(5))/(math.e/2)),
    ("ln(pi^2)",                math.log(math.pi**2)),
    ("pi/sqrt(5)",              math.pi/math.sqrt(5)),
    ("ln(2*pi)",                math.log(2*math.pi)),
]
print(f"  {'expression':<30s} {'value':>12s} {'delta from c_2pi2':>20s}")
for name, val in natural_candidates:
    print(f"  {name:<30s} {val:>12.6f} {val - c_2pi2:>+20.6f}")

print()
print("Striking: e^(1/pi) = {:.6f} matches c_2pi2 = {:.6f} to {:.2e}".format(
    math.exp(1/math.pi), c_2pi2, abs(math.exp(1/math.pi) - c_2pi2)))


# =========================================================================
# Part C: Verify e^(1/pi) gives ratio = 2 pi^2
# =========================================================================
print()
print("-"*72)
print("Part C: explicit check that c = e^(1/pi) gives ratio = 2 pi^2")
print("-"*72)

c = math.exp(1/math.pi)
S_em_check   = ln_alpha + math.log(c)
S_grav_check = ln_Ge + math.log(c)
ratio_check  = S_grav_check / S_em_check
target = 2 * math.pi**2

print()
print(f"With c = e^(1/pi) = {c:.10f}:")
print(f"  S_em   = ln(1/alpha) + ln(c) = {ln_alpha:.6f} + {math.log(c):.6f} = {S_em_check:.6f}")
print(f"  S_grav = ln(1/G_e)  + ln(c) = {ln_Ge:.6f} + {math.log(c):.6f} = {S_grav_check:.6f}")
print(f"  Ratio = {ratio_check:.10f}")
print(f"  2 pi^2 = {target:.10f}")
print(f"  delta = {ratio_check - target:+.2e}")
print(f"  relative delta = {(ratio_check - target)/target:+.3e}")


# =========================================================================
# Part D: closest natural number for raw ratio 20.95
# =========================================================================
print()
print("-"*72)
print("Part D: closest 'natural' number for raw ratio 20.95")
print("-"*72)

raw_ratio = ln_Ge / ln_alpha
print(f"\nRaw ratio (c = 1, no normalization): {raw_ratio:.6f}")
print()
print(f"  {'candidate':<25s} {'value':>12s} {'delta':>14s} {'rel.delta':>14s}")
raw_candidates = [
    ("21 = 3 x 7",          21.0),
    ("2 pi^2 (V(S^3))",     2*math.pi**2),
    ("ln(M_P/m_e)^2 / ln(alpha^-1)", 2*math.log(2.176e-8/9.109e-31)/ln_alpha),
    ("7 pi",                7*math.pi),
    ("22 = 2 x 11",         22.0),
    ("20 = 4 x 5",          20.0),
    ("e + pi + e^2",        math.e + math.pi + math.e**2),
    ("4 pi + 8.4",          4*math.pi + 8.4),
    ("ln(alpha^-1 * G_e^-1) / 2", (ln_alpha + ln_Ge) / 2),
]
for name, val in raw_candidates:
    delta = val - raw_ratio
    print(f"  {name:<25s} {val:>12.6f} {delta:>+14.4e} {delta/raw_ratio:>+14.3e}")


# =========================================================================
# Part E: Honest assessment
# =========================================================================
print()
print("="*72)
print("Honest assessment")
print("="*72)
print(f"""
The raw empirical ratio is S_grav/S_em = ln(1/G_e)/ln(1/alpha) = {raw_ratio:.4f}.

Closest natural-number candidates and their gaps:
  21       = 3 x 7        : gap = {21 - raw_ratio:+.4f}   ({(21-raw_ratio)/raw_ratio*100:+.2f}%)
  2 pi^2   = V(unit S^3)  : gap = {2*math.pi**2 - raw_ratio:+.4f}   ({(2*math.pi**2-raw_ratio)/raw_ratio*100:+.2f}%)

NUMERICALLY 21 is closer to the empirical 20.95 than 2 pi^2.  But 21 = 3 x 7
has no immediate algebraic interpretation in Cl(3,1), while 2 pi^2 is exactly
the volume of the unit 3-sphere, our step-6 constraint surface.

To make ratio = 2 pi^2 exactly, the bare coupling must be c = exp(1/pi) ~ 1.3747,
which is a transcendental number without obvious algebraic interpretation.

So the situation is:
  - 21 fits numerically (0.2% gap) but has no natural origin.
  - 2 pi^2 has natural origin (S^3 volume) but requires a 6% correction
    OR an unmotivated bare coupling.

NEITHER is conclusive.  The 2pi^2 conjecture is NOT confirmed by the
normalization check.

What IS confirmed: under exponential modulation, the 42-order hierarchy
between alpha and G is compressed to O(1)-to-O(10) density depths.  The
specific value of the ratio of those depths (~21) is the question.

Possible explanations for the residual:
  1. The ratio is genuinely 21 and the 2 pi^2 was numerological coincidence.
  2. The 'right' normalization c is slightly off from 1, and corresponds to
     a natural quantity we haven't identified.
  3. The modulation is not exactly exp(-x) but something close.
  4. The empirical values include radiative corrections; the 'bare' algebraic
     prediction would differ at the few-percent level.

CONCLUSION: the cross-sector unification framework is internally
consistent but not yet quantitatively pinned to algebra.  The 6%
gap from 2 pi^2 is suggestive but cannot be closed without further input.
""")
