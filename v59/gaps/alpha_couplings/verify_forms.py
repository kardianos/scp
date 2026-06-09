#!/usr/bin/env python3
"""
verify_forms.py  —  Gap G3/G2: verify the v59 conjectured structural forms for
the fine-structure constant against CODATA/PDG to full precision, and pin down the
exact exponents/constants the v59 fits use.

This script is purely a VERIFICATION harness: it confirms (or refutes) the
arithmetic of the conjectured identities. It makes NO claim that the forms are
derived. Overfitting analysis lives in overfit_scan.py.

Reference values (2026 best):
  alpha(0)^-1   = 137.035999084     (CODATA 2018 / 2022; g-2 consistent)
  alpha(M_Z)^-1 = 127.951           (PDG: 1/127.951 +- 0.009 in MSbar at M_Z)
  sin^2(theta_W)(MSbar, M_Z) = 0.23121
  g_W           = 0.6517            (from G_F, M_W; PDG-derived)
  M_W           = 80.3692 GeV       (PDG 2024 world avg)
  M_Z           = 91.1876 GeV
"""
import math

PI = math.pi

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------
ALPHA0_INV   = 137.035999084
ALPHA0       = 1.0 / ALPHA0_INV
ALPHA_MZ_INV = 127.951            # PDG MSbar(M_Z)
ALPHA_MZ     = 1.0 / ALPHA_MZ_INV
SIN2_W_MSBAR = 0.23121            # PDG MSbar(M_Z)
GW_PDG       = 0.6517
MW_PDG       = 80.3692            # GeV
MZ_PDG       = 91.1876            # GeV

def reldev(pred, ref):
    return abs(pred - ref) / abs(ref)

def line(label, pred, ref, unit=""):
    print(f"  {label:<46s} pred={pred:.8g}{unit:<4s} ref={ref:.8g}{unit:<4s}"
          f"  reldev={reldev(pred,ref)*100:.4f}%")

print("="*78)
print("G3 / G2 verification of v59 conjectured structural forms")
print("="*78)

# ---------------------------------------------------------------------------
# 1. EW-scale form:  alpha(M_Z) = 25/(324 pi^2) = (5/(18 pi))^2
# ---------------------------------------------------------------------------
print("\n[1] EW form: alpha(M_Z) = 25/(324 pi^2) = (5/(18 pi))^2   [conj, emp]")
a_mz_form = 25.0 / (324.0 * PI**2)
a_mz_form2 = (5.0 / (18.0 * PI))**2
assert abs(a_mz_form - a_mz_form2) < 1e-18, "arithmetic identity (5/18pi)^2 = 25/324pi^2 FAILED"
print("    identity (5/(18pi))^2 == 25/(324 pi^2): OK [thm]")
line("alpha(M_Z) value", a_mz_form, ALPHA_MZ)
line("alpha(M_Z)^-1", 1.0/a_mz_form, ALPHA_MZ_INV)
# sqrt form
print(f"    sqrt(alpha_MZ) = {math.sqrt(a_mz_form):.8f}  vs 5/(18pi) = {5/(18*PI):.8f}")

# ---------------------------------------------------------------------------
# 2. The "5 sqrt(alpha)" gauge form and its equivalence to alpha(M_Z) value
#    SM tree: 4*pi*alpha = g_W^2 * sin^2(theta_W).  With sin^2 = 2/9, g_W^2 = 5 sqrt(alpha):
#       4 pi alpha = 5 sqrt(alpha) * 2/9  =>  sqrt(alpha) = 5/(18 pi).
#    Also the SM-DEFINITIONAL linear law: g_W^2 = 4 pi alpha / sin^2 = 18 pi alpha.
#    A line (18 pi alpha) and a sqrt (5 sqrt alpha) meet at exactly ONE alpha.
# ---------------------------------------------------------------------------
print("\n[2] g_W^2 = 5 sqrt(alpha):  is it a law, or alpha(M_Z) in disguise?")
sin2 = 2.0/9.0
def gW2_linear(a):   return 18*PI*a          # SM tree (definitional), holds at every scale
def gW2_sqrt(a):     return 5*math.sqrt(a)   # the v59 conjecture
# intersection: 18 pi a = 5 sqrt(a) -> sqrt(a) = 5/(18 pi)
a_cross = (5.0/(18*PI))**2
print(f"    linear law g_W^2 = 18 pi alpha (def. e=g_W sinθ_W, sin^2=2/9): exact at all scales")
print(f"    the two forms intersect at a single alpha = {a_cross:.8f} (1/{1/a_cross:.3f})")
print(f"    -> g_W^2 = 5 sqrt(alpha) is the single VALUE alpha(M_Z)={a_cross:.6f} in disguise")
print(f"    {'alpha^-1':>10s} {'18 pi a':>10s} {'5 sqrt a':>10s}  agree?")
for ainv in (137.0, 128.0, 127.951, 120.0, 100.0):
    a = 1.0/ainv
    l, s = gW2_linear(a), gW2_sqrt(a)
    agree = "YES" if abs(l-s) < 1e-3 else "no"
    print(f"    {ainv:>10.3f} {l:>10.5f} {s:>10.5f}  {agree}")
# the actual g_W at this alpha:
gW_at_cross = math.sqrt(gW2_sqrt(a_cross))
line("g_W from 5 sqrt(alpha_MZ)", gW_at_cross, GW_PDG)

# ---------------------------------------------------------------------------
# 3. IR form: -ln(alpha) + 2*alpha = pi^2/2 = 8 pi^2 / 16
#    Decompose: pure instanton -ln alpha = pi^2/2 (1.4% off), vs the +2alpha fix.
# ---------------------------------------------------------------------------
print("\n[3] IR form: -ln(alpha) + 2 alpha = pi^2/2 = 8 pi^2 / dim Cl(3,1)  [conj]")
rhs = PI**2 / 2.0
assert abs(8*PI**2/16 - rhs) < 1e-15
print(f"    identity 8 pi^2/16 == pi^2/2: OK [thm];  pi^2/2 = {rhs:.8f}")
# (a) pure instanton form, no +2alpha
a_inst = math.exp(-rhs)
print(f"    (a) PURE instanton  alpha = exp(-pi^2/2):  alpha^-1 = {1/a_inst:.4f}"
      f"   reldev(alpha0) = {reldev(1/a_inst, ALPHA0_INV)*100:.3f}%")
# (b) full form with +2alpha, solved self-consistently for alpha
#     f(a) = -ln a + 2a - pi^2/2 = 0
def f(a): return -math.log(a) + 2*a - rhs
# Newton from a0:
a = ALPHA0
for _ in range(80):
    fp = -1.0/a + 2.0
    a = a - f(a)/fp
print(f"    (b) full form -ln a + 2a = pi^2/2 solved:  alpha^-1 = {1/a:.6f}"
      f"   reldev(alpha0) = {reldev(1/a, ALPHA0_INV)*100:.5f}%")
# residual of the full form AT the measured alpha0:
res0 = f(ALPHA0)
print(f"    residual of full form at measured alpha0: {res0:.3e}  (the '3.6e-5' headline)")
# What correction makes pure pi^2/2 exact at alpha0? (= the reverse-engineered term)
corr = rhs - (-math.log(ALPHA0))
print(f"    correction needed for PURE form at alpha0: pi^2/2 + ln(alpha0) = {corr:.6f}")
print(f"        2*alpha0 = {2*ALPHA0:.6f}   (this is the reverse-engineered '+2alpha')")
for nm, v in [("1/70",1/70),("1/64",1/64),("1/(8pi^2)",1/(8*PI**2)),
              ("alpha0",ALPHA0),("2 alpha0",2*ALPHA0),("3 alpha0",3*ALPHA0)]:
    print(f"        candidate {nm:<10s} = {v:.6f}  (off by {v-corr:+.5f})")

# ---------------------------------------------------------------------------
# 4. Downstream EW masses from the forms (sanity on the chain)
# ---------------------------------------------------------------------------
print("\n[4] Downstream check: M_Z/M_W = 3/sqrt(7) from cos^2 = 7/9")
mz_mw = 3.0/math.sqrt(7.0)
line("M_Z/M_W", mz_mw, MZ_PDG/MW_PDG)

print("\n" + "="*78)
print("SUMMARY (verification only; derivation status is in FINDINGS.md)")
print("="*78)
print(f"  EW form alpha(M_Z)=25/(324pi^2):  reldev {reldev(a_mz_form, ALPHA_MZ)*100:.4f}%  [emp tight]")
print(f"  IR pure -ln a = pi^2/2:           reldev {reldev(1/a_inst, ALPHA0_INV)*100:.3f}%  [mechanism form]")
print(f"  IR full -ln a + 2a = pi^2/2:      reldev {reldev(1/a, ALPHA0_INV)*100:.5f}%  [needs fitted +2a]")
print(f"  g_W^2=5 sqrt(a): one-point coincidence with SM linear law -> NOT a law")
