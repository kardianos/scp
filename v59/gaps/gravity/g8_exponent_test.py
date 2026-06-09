#!/usr/bin/env python3
"""
v59/gaps/gravity/g8_exponent_test.py

G8 -- gravity MAGNITUDE MECHANISM.  Is  G_e = (21/16) alpha^21  a mechanism or
numerology with a striking exponent?

This extends synthesis/gravity_magnitude_test.py with the SPECIFIC mechanism
tests the gap calls for:

  (A) Pin the exponent and prefactor against CODATA (reproduce the headline).
  (B) "PRODUCT over 21 generators = alpha^(21/2)" counting test:  does a
      coherent product / determinant / top-form over the 21 Spin(7) generators
      naturally give the power 21 (or 21/2)?  Test several explicit countings:
        - product of 21 amplitudes each sqrt(alpha)         -> alpha^(21/2)   (g_grav)
        - g_grav^2                                          -> alpha^21       (G_e)
        - determinant of a 21x21 coupling matrix ~ sqrt(alpha) I -> alpha^(21/2)
        - top-form / volume of a 21-dim coupling parallelepiped  -> alpha^(21/2)
      and contrast with the ADDITIVE gauge counting g_W^2 = 5 sqrt(alpha).
  (C) Falsifier scan:  is 21 uniquely selected, or do neighbouring exponents /
      other algebra dims (so(7)=21 vs so(8)=28 vs Cl(3,1)=16 vs G2=14) fit too?
      Report the prefactor each forces, and whether it is O(1)+low-denominator.
  (D) Alternative hierarchy origins:  dimensional transmutation exp(-2pi/(b alpha))
      and instanton exp(-S) -- can EITHER reproduce 1.75e-45 with a SM-plausible
      coefficient?  (the e^{-S} route many BSM models use.)

A "mechanism" (not numerology) would: (i) fix the exponent 21 from a counting that
is forced, AND (ii) fix the prefactor to a unique structural rational, AND (iii)
say why gravity (not gauge) takes the product.  We test how many of these hold.
"""

import math
from fractions import Fraction

mPl   = 1.220890e22
alpha = 1 / 137.035999084          # alpha(0) IR
aGe   = (0.51099895 / mPl) ** 2    # alpha_G(electron)
dimSpin7, dimCl31, dimG2, dimSpin8 = 21, 16, 14, 28

print("=" * 76)
print("G8  PART A -- pin exponent and prefactor")
print("=" * 76)
pred = (dimSpin7/dimCl31) * alpha**dimSpin7
n_fit = math.log(aGe / (dimSpin7/dimCl31)) / math.log(alpha)
c_at_21 = aGe / alpha**21
print(f" alpha_G(e) CODATA          = {aGe:.6e}")
print(f" (21/16) alpha^21           = {pred:.6e}   ratio {pred/aGe:.5f}  ({abs(pred/aGe-1)*100:.3f}%)")
print(f" exponent forcing 21/16     = {n_fit:.5f}   (= dim Spin(7) = 21 to {abs(n_fit-21)/21*100:.4f}%)")
print(f" prefactor forced AT n=21   = aGe/alpha^21 = {c_at_21:.5f}")
print(f"   nearest rationals: 21/16={21/16:.4f} (err {abs((21/16)/c_at_21-1)*100:.2f}%), "
      f"17/13={17/13:.4f} (err {abs((17/13)/c_at_21-1)*100:.2f}%), 4/3={4/3:.4f} (err {abs((4/3)/c_at_21-1)*100:.2f}%)")

print("\n" + "=" * 76)
print("G8  PART B -- does a PRODUCT/DET/TOP-FORM over 21 generators give alpha^(21/2)?")
print("=" * 76)
sa = math.sqrt(alpha)

# B1 coherent product of 21 amplitudes each sqrt(alpha)
prod21 = sa ** dimSpin7
print(f"\n B1 coherent product  prod_{{i=1}}^{{21}} sqrt(alpha) = alpha^(21/2) = {prod21:.4e}")
print(f"    (g_grav). Squared = alpha^21 = {prod21**2:.4e}  -> matches G_e order.")

# B2 determinant of a 21x21 diagonal coupling matrix g = sqrt(alpha) * I_21
import numpy as np
G = sa * np.eye(dimSpin7)
detG = np.linalg.det(G)
print(f"\n B2 det( sqrt(alpha) I_21 ) = alpha^(21/2) = {detG:.4e}")
print("    A DETERMINANT over the 21 generator directions reproduces alpha^(21/2)")
print("    EXACTLY (det of n-dim scalar matrix = scalar^n). This is the cleanest")
print("    'product over all generators' realization: the coupling is the VOLUME")
print("    (top exterior form) of the 21-dim generator coupling lattice.")

# B3 top-form / volume of parallelepiped with 21 edges each length sqrt(alpha)
vol = sa ** dimSpin7
print(f"\n B3 top-form vol(e_1 ^ ... ^ e_21), |e_i|=sqrt(alpha) = alpha^(21/2) = {vol:.4e}")
print("    Same number. The Lambda^21 (top) form of the Spin(7) Lie algebra is")
print("    1-dimensional; its 'norm' with each leg ~ sqrt(alpha) is alpha^(21/2).")

# B4 contrast with additive gauge counting
print(f"\n B4 CONTRAST gauge (ADDITIVE):  g_W^2 = 5 sqrt(alpha) = {5*sa:.4f}")
print("    Gauge = SUM over a 5-generator subset (5 = 21-16 = dual Coxeter).")
print("    Gravity = PRODUCT/DET/VOL over ALL 21 generators.")
print("    => structurally consistent: additive index (linear in count) vs")
print("       multiplicative volume (exponential in count). The det/top-form")
print("       IS a natural object that turns 'all 21 generators' into a power 21.")
print("    HONEST: this shows alpha^(21/2) is the NATURAL value of a 21-leg")
print("    determinant -- but it does NOT derive WHY gravity = that determinant.")
print("    No Lagrangian yet maps the graviton coupling to det(generator-coupling).")

# B5 where could a 21-fold product come from dynamically? a one-loop / anomaly
#    argument would give a SUM of logs (additive in count), not a product. A
#    product needs a SIMULTANEOUS (coherent) coupling to all 21 -- e.g. a
#    Wess-Zumino / theta-like top-form term int Tr(F^21) on a 21-form. Note
#    dim Spin(7)=21 is ODD, and Tr(F^k) top-forms live in even degree; a 21-form
#    on the 21-dim group manifold is the volume form (Haar), which is the
#    natural 'all generators at once' object. Flag this as the most promising
#    Lagrangian candidate (a topological volume term), still UNBUILT.
print("\n B5 dynamical origin:  a SUM-over-loops gives additive (log) counting;")
print("    a PRODUCT needs a COHERENT all-21-generator object. The natural one is")
print("    the Haar/volume top-form on the 21-dim Spin(7) group manifold")
print("    (a Wess-Zumino-like int over the group). This is the most plausible")
print("    Lagrangian seed for the product, but it is NOT built/derived here.")

print("\n" + "=" * 76)
print("G8  PART C -- is 21 uniquely selected? (falsifier scan over algebra dims)")
print("=" * 76)
print(f" {'exponent n':>22s} {'prefactor c=aGe/alpha^n':>24s}  {'nearest rat (den<=20)':>22s}")
for label, n in [("dim G2 = 14", 14), ("dim Cl(3,1) = 16", 16),
                 ("dim Spin(7) = 21", 21), ("dim Spin(8) = 28", 28),
                 ("18", 18), ("19", 19), ("20", 20), ("22", 22), ("23", 23)]:
    c = aGe / alpha**n
    fr = Fraction(c).limit_denominator(20)
    flag = "  <== O(1), structural" if 0.3 < c < 3 else ""
    print(f" {label:>22s} {c:24.4e}  {str(fr)+' = '+format(float(fr),'.4f'):>22s}{flag}")
print(" => ONLY n=21 lands the prefactor in the O(1) window. n=14,16,28 are off by")
print("    huge factors. So the EXPONENT 21 is strongly selected; the data DEMANDS")
print("    an exponent of exactly 21 (= dim Spin(7)) to get an O(1) coefficient.")
print("    This is the robust part. The prefactor's exact rational (21/16 vs 17/13")
print("    vs 4/3) is NOT uniquely pinned -- weaker, as flagged.")

print("\n" + "=" * 76)
print("G8  PART D -- alternative hierarchy origins (transmutation / instanton)")
print("=" * 76)
target = aGe
# D1 dimensional transmutation: ratio ~ exp(-2pi/(b alpha)). Solve for b.
b_eff = -2*math.pi / (math.log(target) * alpha)
print(f"\n D1 dim transmutation  (m/M_Pl)^2 ~ exp(-2pi/(b alpha)):")
print(f"    requires b = {b_eff:.4f}  (one-loop beta coefficient). Not a SM-natural")
print(f"    integer; and transmutation gives a SCALE, not a dimensionless ratio")
print(f"    relative to M_Pl without inputting M_Pl. So this RELABELS, not derives.")
# D2 instanton: ratio ~ exp(-S), S = ? compare to alpha^21 = exp(21 ln alpha)
S_inst = -math.log(target)
print(f"\n D2 instanton  ratio ~ exp(-S):  S = -ln(aGe) = {S_inst:.3f}")
print(f"    Compare alpha^21 = exp(21 ln alpha) = exp({21*math.log(alpha):.3f}) ->")
print(f"    S = 21*|ln alpha| = {21*abs(math.log(alpha)):.3f}.  So the v59 form IS")
print(f"    an instanton-like e^{{-S}} with S = 21 ln(1/alpha) = (dim Spin(7)) ln(1/alpha).")
print(f"    i.e. 21 instanton-units, one per Spin(7) generator. This is the SAME")
print(f"    'product over 21' statement in exponent language -- self-consistent,")
print(f"    but still posits the count 21, not derives a single-instanton action.")
# 8pi^2/g^2 instanton action with g^2 = 5 sqrt(alpha)?
S_8pi2 = 8*math.pi**2 / (5*sa)
print(f"\n    (aside) a single SU(2)_L instanton S=8pi^2/g_W^2 = 8pi^2/(5 sqrt a) = {S_8pi2:.1f}")
print(f"    is ~{S_8pi2/S_inst:.1f}x too large -- a SINGLE gauge instanton does NOT")
print(f"    give 1e-45; you need exactly the 21-fold (per-generator) version.")

print("\n" + "=" * 76)
print("VERDICT on G8")
print("=" * 76)
print("- ROBUST: the EXPONENT is forced to 21 = dim Spin(7) (to 0.004%). The data")
print("  demands exactly 21; no neighbouring dim works. f_g ~ alpha^(21/2).")
print("- MECHANISM (partial): alpha^(21/2) IS the natural value of a determinant /")
print("  top-form / coherent product over the 21 Spin(7) generators (B2/B3 exact),")
print("  and equivalently an e^{-S} with S = 21 ln(1/alpha) (D2). The product vs")
print("  additive-gauge contrast (det/volume vs sum/index) is structurally clean.")
print("- NOT DERIVED: no Lagrangian maps the graviton coupling to that determinant/")
print("  top-form; the Haar volume term (B5) is the best candidate but unbuilt.")
print("  The PREFACTOR is not uniquely pinned (21/16 vs 17/13 vs 4/3).")
print("- BOTTOM LINE: G8 is a VALUE conjecture with a STRONG (theorem-adjacent)")
print("  exponent and a PLAUSIBLE-but-unbuilt product mechanism. The exponent is")
print("  the most law-like gravity number in v59; the prefactor and Lagrangian")
print("  remain open.")
