#!/usr/bin/env python3
"""
rg_invariance_test.py  --  the highest-value falsifiable test in the cluster (B2).

CLAIM under test: Koide Q is invariant under a COMMON (flavor-universal) rescaling
m_i -> eta*m_i, because Q = sum m / (sum sqrt m)^2 -> eta*sum m / (eta*(sum sqrt m)^2)
= Q. So if QCD running were flavor-universal at one loop, Q would be RG-invariant and
the 5-7% observed drift would be due ONLY to flavor-NON-universal effects (different
anomalous dimensions are NOT the issue at LO QCD -- gamma_m is flavor-universal; the
drift comes from threshold matching, EW/Yukawa, and the mixed-scale inputs).

This script:
  (1) proves the common-rescaling invariance analytically (numeric check).
  (2) shows that one-loop QCD gamma_m IS flavor-universal => the LO-QCD-run Q is
      EXACTLY the input-scale Q. Hence the scale drift in quark_koide_rg.py is NOT
      from universal QCD running; it is from the mismatched INPUT scales of the
      datasets (each table defines its masses with different threshold/scheme
      choices) plus genuinely flavor-dependent (Yukawa/EW) running of the top.
  (3) isolates how much of the drift the top quark's large Yukawa (flavor-NON-
      universal) contributes to Q_u.

Conclusion preview: Q is protected against the universal part of running but NOT
against the top-Yukawa (non-universal) part, which is exactly why Q_u drifts ~5%.
This makes Q_u especially un-predictive; Q_d (no large Yukawa) drifts less from
running and more from the s-quark mass uncertainty.
"""

import numpy as np

def koide_Q(masses):
    m = np.asarray(masses, float)
    return float(np.sum(m) / np.sum(np.sqrt(m))**2)

print("=" * 78)
print("(1) Common-rescaling invariance of Koide Q  (analytic, numeric check)")
print("=" * 78)
base_down = [4.67e-3, 93.4e-3, 4.18]
base_up   = [2.16e-3, 1.27, 172.69]
for eta in [0.5, 1.0, 2.0, 10.0]:
    Qd = koide_Q([eta*m for m in base_down])
    Qu = koide_Q([eta*m for m in base_up])
    print(f"  eta={eta:5.1f}:  Q_d={Qd:.6f}  Q_u={Qu:.6f}   (identical => invariant)")
print("  => Q is EXACTLY invariant under a common (flavor-universal) rescaling. [proved]")

print()
print("=" * 78)
print("(2) One-loop QCD mass anomalous dimension is FLAVOR-UNIVERSAL")
print("=" * 78)
print("""
  m_i(mu) = m_i(mu0) * [alpha_s(mu)/alpha_s(mu0)]^(gamma0/2b0),
  gamma0 = 8 (=2*C_F*3? -> gamma_m0 = 6 C_F = 8 in common norm), SAME for all flavors.
  Therefore the LO-QCD running multiplies ALL same-type quark masses by the SAME
  factor eta_QCD(mu). By (1), Q is UNCHANGED by pure LO-QCD running.

  COROLLARY: the 2-7% scale drift seen in quark_koide_rg.py is NOT one-loop QCD
  running. It is:
     - mismatched INPUT scales / thresholds between the literature datasets, and
     - flavor-NON-universal running: the TOP Yukawa contributes a large, flavor-
       specific anomalous dimension to m_t (and EW corrections), which DOES change Q_u.
""")

# Demonstrate: apply a common QCD factor to one consistent input set and confirm Q fixed.
mZ_down = [2.75e-3, 0.0535, 2.89]
mZ_up   = [1.29e-3, 0.627, 171.7]
for eta in [1.0, 1.8]:  # ~ running 2GeV<->MZ common factor magnitude
    print(f"  common factor eta={eta}:  Q_d {koide_Q([eta*m for m in mZ_down]):.6f}, "
          f"Q_u {koide_Q([eta*m for m in mZ_up]):.6f}  (unchanged)")

print()
print("=" * 78)
print("(3) How much does the TOP (flavor-non-universal) drive Q_u drift?")
print("=" * 78)
print("""
  Vary ONLY m_t (the flavor-non-universal runner), hold m_u,m_c fixed at MZ values,
  and watch Q_u. This isolates the non-universal contribution.
""")
mu_uc = [1.29e-3, 0.627]
print(f"  {'m_t (GeV)':>10} {'Q_u':>10} {'gap to 23/27 %':>16}")
for mt in [150, 162, 171.7, 172.7, 384]:
    Qu = koide_Q(mu_uc + [mt])
    gap = (23/27 - Qu)/Qu*100
    print(f"  {mt:10.1f} {Qu:10.6f} {gap:16.3f}")
print("""
  => Q_u is dominated by m_t (it's by far the largest sqrt-mass). Because m_t runs
     flavor-non-universally (top Yukawa), Q_u is NOT protected -- its drift is
     intrinsic, not removable by a common rescaling. The 0.34% match at the PDG
     mixed scale is therefore a coincidence of one particular m_t convention
     (the pole-ish 172.69), NOT an RG-invariant statement.
""")

print("=" * 78)
print("(4) phi_d / phi_l hierarchy check (G6, avenue C2: radiative O(alpha) size)")
print("=" * 78)
alphaMZ = 1/127.951
alpha0 = 1/137.036
phi_l = 2/9
phi_d = 0.10859
phi_u = 0.07251
print(f"  phi_d/phi_l = {phi_d/phi_l:.4f}")
print(f"  14*alpha(MZ)/(2/9) = {14*alphaMZ/phi_l:.4f}  -> gap {abs(14*alphaMZ/phi_l - phi_d/phi_l)/(phi_d/phi_l)*100:.2f}%")
print(f"  phi_u/phi_l = {phi_u/phi_l:.4f}")
print(f"  10*alpha(0)/(2/9) = {10*alpha0/phi_l:.4f}  -> gap {abs(10*alpha0/phi_l - phi_u/phi_l)/(phi_u/phi_l)*100:.2f}%")
print("""
  The HIERARCHY phi_quark ~ O(alpha) * phi_lepton is the defensible content:
  the quark phases are an order-alpha effect (radiative?), the lepton phase is O(1)
  (tree, = 2/9). The specific integers (14, 10) remain conjectural; the SIZE (O(alpha))
  is the robust observation.
""")
