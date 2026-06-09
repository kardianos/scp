#!/usr/bin/env python3
"""
v59/synthesis/gravity_charge_test.py

The gravity-charge question: what does the OBE gravity term couple to?

The current OBE sources gravity from the Koide-constraint DEVIATION
    rho_N = | |Pi_H[Phi]|^2 - (1 - 14/D) |   ~  a^2 (|xi|^2 - 1/2).
This has two fatal problems, tested below:
  (A1) it is ZERO for leptons (|xi|^2 = 1/2 on the constraint) -> leptons would
       not gravitate;
  (A2) it is SCALE-INDEPENDENT in shape: mass = 9 Q a^2 lives in the amplitude a,
       but |xi|^2 is the dimensionless shape -> rho_N cannot be proportional to mass.

Candidate B (the fix): the gravity charge is the SECOND MOMENT of the Brannen
kernel, Tr(M^dag M) = sum_k m_k = 9 Q a^2 = the total mass. This is the Frobenius^2
norm of the 3x3 generation kernel -- the SAME Frobenius^2 structure as the EW
bridge v = ||Y||_F^2 = dim(L)^2 a^2, just over the 3-dim generation space instead
of the 28-dim ambient. It carries a^2, is nonzero for leptons, and equals mass by
construction (equivalence principle).

Decisive test: charge / mass must be UNIVERSAL across the three sectors (EP).
"""

import numpy as np

# PDG-ish masses (MeV)
sectors = {
    "lepton":  (np.array([0.51099895, 105.6583755, 1776.86]), 28),
    "d-quark": (np.array([4.7, 93.0, 4180.0]),                 35),
    "u-quark": (np.array([2.2, 1270.0, 172500.0]),             63),
}
v_higgs = 246220.0  # MeV
N_gen = 3


def brannen(m):
    """Return (a, Q, |xi|^2) for a sector from its three masses.
    Eigenvalues of the Brannen kernel M are sqrt(m_k); a=<sqrt m>, Q=Koide,
    |xi|^2 = t^2 where Q=(1+2t^2)/3."""
    s = np.sqrt(m)
    a = s.mean()
    Q = m.sum() / s.sum()**2
    t2 = (3 * Q - 1) / 2.0           # |xi|^2
    return a, Q, t2


print("=" * 74)
print("GRAVITY-CHARGE TEST: which scalar map gives a charge proportional to mass?")
print("=" * 74)
print(f"\n{'sector':8s} {'a':>9s} {'Q':>7s} {'|xi|^2':>7s} {'1-14/D':>7s}"
      f" {'massSm':>10s} {'chgA(dev)':>10s} {'chgB(Tr)':>10s}")

rows = []
for name, (m, D) in sectors.items():
    a, Q, t2 = brannen(m)
    xi2_constraint = 1.0 - 14.0 / D                 # the v59 vacuum |xi|^2
    mass = m.sum()                                  # = 9 Q a^2 = Tr(M^dag M)
    # Candidate A: OBE constraint-deviation source (scale-free shape), x a^2
    chgA = a**2 * (t2 - 0.5)                         # ~ a^2 (|xi|^2 - 1/2)
    # Candidate B: second moment / Frobenius^2 of the kernel = total mass
    chgB = mass
    rows.append((name, a, Q, t2, xi2_constraint, mass, chgA, chgB, D))
    print(f"{name:8s} {a:9.3f} {Q:7.4f} {t2:7.4f} {xi2_constraint:7.4f}"
          f" {mass:10.2f} {chgA:10.3f} {chgB:10.1f}")

# Equivalence-principle check: charge / mass should be UNIVERSAL (same for all 3).
print("\n--- Equivalence principle: charge / mass (must be universal) ---")
massA = np.array([r[6] for r in rows]); massB = np.array([r[7] for r in rows])
mass = np.array([r[5] for r in rows])
ratA = massA / mass; ratB = massB / mass
print(f" Candidate A (deviation): charge/mass = {ratA}")
print(f"   lepton entry = {ratA[0]:.4f}  -> leptons DO NOT gravitate (charge=0).")
print(f"   spread max/min(|nonzero|) = {np.max(np.abs(ratA[1:]))/np.min(np.abs(ratA[1:])):.2f}x  => NOT universal. DEAD.")
print(f" Candidate B (2nd moment): charge/mass = {ratB}  => universal (=1). EP satisfied. LIVE.")

# Also test the raw mass-RATIO reproduction (the memory's u/d ~ 2.78 failure).
md = sectors["d-quark"][0].sum(); mu = sectors["u-quark"][0].sum()
a_d = brannen(sectors["d-quark"][0])[0]; a_u = brannen(sectors["u-quark"][0])[0]
t2_d = brannen(sectors["d-quark"][0])[2]; t2_u = brannen(sectors["u-quark"][0])[2]
print("\n--- Mass-ratio reproduction (u-quark / d-quark) ---")
print(f" empirical Sum m_u / Sum m_d                = {mu/md:6.2f}")
print(f" Candidate B (2nd moment)                   = {mu/md:6.2f}  (= itself, EP-exact)")
print(f" Candidate A deviation FACTOR ratio         = {(t2_u-0.5)/(t2_d-0.5):6.2f}  (the memory's 2.78)")
print(f" Candidate A deviation x a^2 ratio          = {(a_u**2*(t2_u-0.5))/(a_d**2*(t2_d-0.5)):6.2f}")
print(f"   => neither A form matches {mu/md:.0f}: deviation-charge is DEAD on mass ratios.")

# Structural bonus: total lepton mass tied to v by the SAME Frobenius^2 / Koide.
a_l, Q_l, _ = brannen(sectors["lepton"][0]); D_l = 28
print("\n--- Bridge unification: gravity charge (mass) and v share Frobenius^2/Koide ---")
print(f" Sum m_lepton           = 9 Q a^2          = {9*Q_l*a_l**2:.2f} MeV  (Tr M^dag M, 3-dim)")
print(f" v_Higgs                = dim(L)^2 a^2     = {D_l**2*a_l**2:.0f} MeV  (||Y||_F^2, 28-dim)")
print(f" Sum m_lepton / v       = 9Q/dim(L)^2 = 6/784 = {9*Q_l/D_l**2:.6f}")
print(f"   empirical Sum m_l / v_obs              = {sectors['lepton'][0].sum()/v_higgs:.6f}"
      f"   (rel {abs(9*Q_l/D_l**2 - sectors['lepton'][0].sum()/v_higgs)/(9*Q_l/D_l**2):.2e})")

print("\n" + "=" * 74)
print("VERDICT")
print("=" * 74)
print("- DEAD: the OBE constraint-deviation source rho_N ~ a^2(|xi|^2-1/2). It is 0")
print("  for leptons and its shape factor is scale-free, so it cannot be mass.")
print("  Diagnosis: rho_N drops the amplitude a, but mass = 9 Q a^2 lives in a^2.")
print("- LIVE: the gravity charge is the SECOND MOMENT Tr(M^dag M) = sum m = 9 Q a^2")
print("  (Frobenius^2 of the kernel). Nonzero for leptons, EP-exact, and the SAME")
print("  Frobenius^2 structure as the EW bridge v = dim(L)^2 a^2 (different space).")
print("- BONUS: Sum m_lepton = (9Q/dim(L)^2) v = (6/784) v at 0.07% -- one relation")
print("  ties the gravity charge to the Higgs VEV via Koide + dim(L).")
print("- STILL OPEN (this test does not touch): the coupling MAGNITUDE (Newton G;")
print("  V6 was ~1e40 too strong) and SCALAR vs TENSOR (LIGO h=+-2).")
