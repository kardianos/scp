#!/usr/bin/env python3
"""
v59/synthesis/gravity_magnitude_test.py

The gravity COUPLING MAGNITUDE -- the historically fatal problem (V6 was ~1e40 too
strong). In the OBE, gravity is  box Omega_grav = f_g * rho_grav, rho_grav = Tr(M^dag M)
(the second-moment / mass charge, established live in gravity_charge_test.py). The
point-mass interaction energy is  f_g^2 m1 m2 / (4 pi r), so

        G_N = f_g^2 / (4 pi)   and   alpha_G(e) = G_N m_e^2 / (hbar c) = (m_e/M_Pl)^2 .

So the magnitude question = "what fixes f_g (equivalently G_N) at the ~1e-45 level?"

v59 already has a structural conjecture: G_e = alpha_G(electron) = (21/16) alpha^21,
with 21 = dim Spin(7) and 21/16 = dim Spin(7)/dim Cl(3,1). alpha is the IR value
alpha(0) (gravity is long-range/IR). This script:
  1. pins the exponent and prefactor against CODATA;
  2. tests it as a UNIVERSAL G_N across particles;
  3. tests a 'combination' MECHANISM: gravity = product over the 21 Spin(7)
     generators each ~sqrt(alpha) (multiplicative), vs gauge forces = additive
     subsets (g_W^2 = 5 sqrt(alpha), 5 = 21-16 generators);
  4. an HONEST numerology check: is (21/16, 21) uniquely clean, or one of many fits?
"""

import math

# --- constants (CODATA-ish, MeV) ---
mPl   = 1.220890e22            # Planck mass
alpha = 1 / 137.035999084      # alpha(0), IR
m = {"electron": 0.51099895, "muon": 105.6583755, "tau": 1776.86,
     "proton": 938.272, "W": 80377.0, "Z": 91187.6, "top": 172500.0}

aG = {p: (mm / mPl)**2 for p, mm in m.items()}   # alpha_G(particle) = (m/M_Pl)^2

dimSpin7, dimCl31, dimG2 = 21, 16, 14

print("=" * 74)
print("GRAVITY MAGNITUDE: is G_e = (21/16) alpha^21  (21 = dim Spin(7)) ?")
print("=" * 74)

aGe = aG["electron"]
pred = (dimSpin7 / dimCl31) * alpha**dimSpin7
n_fit = math.log(aGe / (dimSpin7 / dimCl31)) / math.log(alpha)
print(f"\n alpha_G(e) (CODATA)        = {aGe:.5e}")
print(f" (21/16) alpha^21           = {pred:.5e}   ratio = {pred/aGe:.4f}  ({abs(pred/aGe-1)*100:.2f}%)")
print(f" exponent giving prefactor 21/16:  n = {n_fit:.4f}   (dim Spin(7) = 21)")
print(f" => exponent is integer 21 to {abs(n_fit-21)/21*100:.3f}%  -- not 'near 21', it IS 21.")
print(f"\n Hierarchy form:  m_e/M_Pl = sqrt(21/16) * alpha^(21/2)")
print(f"   m_e/M_Pl = {m['electron']/mPl:.4e} ;  sqrt(21/16) alpha^10.5 = "
      f"{math.sqrt(dimSpin7/dimCl31)*alpha**(dimSpin7/2):.4e}")

# alpha(M_Z) check: gravity is IR, so it must use alpha(0), NOT alpha(M_Z).
aMZ = 1/127.951
print(f"\n IR consistency: (21/16) alpha(M_Z)^21 = {(21/16)*aMZ**21:.3e} "
      f"(off by {(21/16)*aMZ**21/aGe:.1f}x) -- gravity correctly needs alpha(0), not alpha(M_Z).")

# ---------------------------------------------------------------------------
# 2. Universality: G_N is one constant; alpha_G(particle) = (m/m_e)^2 * G_e.
# ---------------------------------------------------------------------------
print("\n" + "-" * 74)
print("UNIVERSAL G_N test:  alpha_G(particle) =? (m/m_e)^2 * (21/16) alpha^21")
print("-" * 74)
for p, mm in m.items():
    predp = (mm / m["electron"])**2 * (dimSpin7/dimCl31) * alpha**dimSpin7
    print(f"  {p:8s}: alpha_G = {aG[p]:.3e}   pred = {predp:.3e}   ratio {predp/aG[p]:.4f}")
print("  => one G_N fits all (alpha_G ~ m^2). The alpha^21 is specifically the")
print("     ELECTRON-electron coupling; other particles follow by mass^2. Consistent.")

# ---------------------------------------------------------------------------
# 3. Combination MECHANISM hypothesis: multiplicative over 21 generators.
#    gauge: g_W^2 = 5 sqrt(alpha)   (5 = 21-16 generators, ADDITIVE subset)
#    gravity: g_grav = product over ALL 21 Spin(7) gens of sqrt(alpha)
#             => g_grav^2 = (sqrt(alpha))^(2*21) = alpha^21   (MULTIPLICATIVE)
# ---------------------------------------------------------------------------
print("\n" + "-" * 74)
print("COMBINATION mechanism: additive gauge vs multiplicative gravity")
print("-" * 74)
gW2 = 5 * math.sqrt(alpha)
print(f"  gauge   g_W^2 = 5*sqrt(alpha)    = {gW2:.4f}   (5 = 21-16 gens, ADDITIVE)")
print(f"  gravity g_grav^2 = (sqrt(alpha))^(2*21) = alpha^21 = {alpha**21:.3e}")
print(f"          (MULTIPLICATIVE over all 21 Spin(7) generators)")
print(f"  prefactor 21/16 = dim Spin(7)/dim Cl(3,1) (the embedding index base).")
print("  HYPOTHESIS: gravity is the unique coupling requiring ALL 21 generators")
print("  coherently (a product), so it is alpha^(21/2)-suppressed; gauge forces use")
print("  additive subsets (5), so stay O(sqrt(alpha)). This DISTINGUISHES why gravity")
print("  is uniquely weak. It is a hypothesis -- no Lagrangian derives the product.")

# ---------------------------------------------------------------------------
# 4. HONEST numerology check: scan integer exponents, see which gives a clean
#    O(1) low-denominator-rational prefactor. Is (21/16, 21) uniquely clean?
# ---------------------------------------------------------------------------
print("\n" + "-" * 74)
print("NUMEROLOGY CHECK: for each integer n, prefactor c=alpha_G(e)/alpha^n")
print("-" * 74)
from fractions import Fraction
print(f"  {'n':>3s} {'c=aGe/alpha^n':>14s}  nearest simple rational (den<=20)")
for nn in range(18, 25):
    c = aGe / alpha**nn
    fr = Fraction(c).limit_denominator(20)
    note = ""
    if 0.3 < c < 3:
        note = f"<- O(1); {fr} = {float(fr):.4f} (err {abs(float(fr)/c-1)*100:.1f}%)"
        if nn == 21:
            note += "  [21=dimSpin7, 21/16=dimSpin7/dimCl31]"
    print(f"  {nn:3d} {c:14.4e}  {note}")
print("  => only n=21 gives an O(1) prefactor near a structural rational (21/16);")
print("     n!=21 give c far from 1 or with ugly denominators. The EXPONENT is the")
print("     real signal (pinned to 21.000); the prefactor 21/16 is plausible but, as")
print("     a free O(1) rational, weaker evidence. NOTE 21 is over-determined in v59")
print("     (21 = dimSpin7 = 28-7 = 14+7 = 35-14), so the integer alone is not unique.")

print("\n" + "=" * 74)
print("VERDICT on the gravity magnitude")
print("=" * 74)
print("- LIVE numerically: G_e = (21/16) alpha^21 matches at 0.25%, and the exponent")
print("  is 21.000 = dim Spin(7) to 0.002% (using alpha(0), as IR gravity requires).")
print("  So f_g = sqrt(4 pi G_N) ~ alpha^(21/2): the V6 1e40 overshoot is EXACTLY the")
print("  missing alpha^21 ~ 1e-45 suppression. The hierarchy m_e/M_Pl = sqrt(21/16)")
print("  alpha^(21/2) is the magnitude statement.")
print("- MECHANISM: still a hypothesis. The 'multiplicative over 21 generators' picture")
print("  explains why gravity (all gens, product) is uniquely weak vs gauge (additive")
print("  subset, 5=21-16), but NO Lagrangian derives the product or the 21/16 prefactor.")
print("- HONEST: this is a VALUE conjecture (like alpha itself), not a derivation. The")
print("  exponent 21 is a strong signal; the prefactor and mechanism are not established.")
