#!/usr/bin/env python3
"""
quark_phases_scales.py  --  v59 gap attack, quark-flavour cluster (G6).

Test the conjectured structural forms for the quark Brannen phases phi_d, phi_u
and the scale ratios a_u^2/a_d^2, a_u^2/a_l^2 -- WITH overfitting controls.

G6 facts to scrutinize:
  * phi_d ~ +0.110 rad, phi_u ~ -2.02 rad (or fundamental-domain -0.0725 rad).
  * phi_q != Q_q/3 (the EXACT lepton relation phi_l = Q_l/3 = 2/9 does NOT extend).
  * prior numerology: phi_d ~ 14*alpha(MZ) (dim G2 * alpha), phi_u ~ -10*alpha(0)
    (dim Spin(5) * alpha) at ~0.7%.
  * alt suggestive: 3*phi_d ~ 1/3, 3*phi_u ~ -3/14.
  * a_u^2/a_d^2 ~ 35, a_u^2/a_l^2 ~ 72.

The phase is defined only mod the generation S_3: only cos(3 phi) is an observable.
So the HONEST test is on cos(3 phi), not phi itself. We test all candidate forms
against that, and we count how many "structural integer x alpha" forms land within
a comparable tolerance (the overfitting control).
"""

import numpy as np
import json

alpha0 = 1/137.035999084     # alpha(0), Thomson limit
alphaMZ = 1/127.951          # alpha(M_Z), MS-bar

# --------------------------------------------------------------------------
# Empirical Brannen phases (from FINDINGS_brannen_dynamics.md high-precision fit)
# These come from fitting sqrt(m_k) = a(1 + 2t cos(phi + 2pi k/3)) to the masses.
# Note: the phase depends on the mass SET (hence scale) used, just like Q.
# We use the same PDG-native set as the prior v59 fit for continuity, but flag it.
# --------------------------------------------------------------------------
phi_d_emp = 0.10859     # rad  (down sector, fundamental domain)
phi_u_emp = -0.07251    # rad  (up sector, fundamental domain; = -2.02 + 2pi/3 region)
phi_l_emp = 0.22222     # rad  (lepton, = 2/9 exactly to 1e-5)

print("=" * 78)
print("Quark Brannen phases: candidate structural forms")
print("=" * 78)
print(f"  phi_l (lepton)  = {phi_l_emp:.5f}  matches 2/9 = {2/9:.5f}  (EXACT, phi_l=Q_l/3)")
print(f"  phi_d (down)    = {phi_d_emp:.5f}")
print(f"  phi_u (up)      = {phi_u_emp:.5f}")
print()
print("  Does phi_q = Q_q/3 extend?  (the lepton relation)")
print(f"    Q_d/3 = (11/15)/3 = {(11/15)/3:.5f}   vs phi_d = {phi_d_emp:.5f}  "
      f"-> gap {abs((11/15)/3 - phi_d_emp)/phi_d_emp*100:.1f}%  (NO)")
print(f"    Q_u/3 = (23/27)/3 = {(23/27)/3:.5f}   vs phi_u = {phi_u_emp:.5f}  "
      f"-> WRONG SIGN/magnitude (NO)")

# --------------------------------------------------------------------------
# Candidate forms. We test BOTH the raw phi and cos(3 phi) (the true observable).
# --------------------------------------------------------------------------
candidates_d = {
    "1/9 (=Q_l/6)":           1/9,
    "14*alpha(MZ) [dimG2]":   14*alphaMZ,
    "14*alpha(0)":            14*alpha0,
    "15*alpha(0)":            15*alpha0,
    "Q_d/3 = 11/45":          (11/15)/3,
}
candidates_u = {
    "-1/14 (=-1/dimG2... )":  -1/14,
    "-10*alpha(0) [dimSp5]":  -10*alpha0,
    "-10*alpha(MZ)":          -10*alphaMZ,
    "-9*alpha(0)":            -9*alpha0,
    "-Q_u/3":                 -(23/27)/3,
}

def report(name, emp, cands):
    print()
    print(f"  --- {name}: phi_emp = {emp:.5f}, observable cos(3 phi) = {np.cos(3*emp):.5f} ---")
    print(f"    {'candidate':<28} {'value':>10} {'gap_phi%':>9} {'cos3phi':>9} {'gap_cos%':>9}")
    for label, val in cands.items():
        gphi = abs(val - emp)/abs(emp)*100
        c_emp = np.cos(3*emp); c_val = np.cos(3*val)
        gcos = abs(c_val - c_emp)/abs(c_emp)*100 if c_emp != 0 else float('nan')
        print(f"    {label:<28} {val:10.5f} {gphi:9.2f} {c_val:9.5f} {gcos:9.2f}")

report("DOWN", phi_d_emp, candidates_d)
report("UP", phi_u_emp, candidates_u)

# --------------------------------------------------------------------------
# OVERFITTING CONTROL for the "integer * alpha" pattern.
# How many small integers n give n*alpha within 1% of phi_d (or phi_u)?
# If many do, the "14*alpha" / "10*alpha" reading is not special.
# --------------------------------------------------------------------------
print()
print("=" * 78)
print("Overfitting control: how many 'n*alpha' forms fit each phase?")
print("=" * 78)
for nm, emp in [("phi_d", phi_d_emp), ("phi_u", abs(phi_u_emp))]:
    print(f"  {nm} = {emp:.5f}:")
    for tol in [1, 2, 5]:
        hits_a0 = [n for n in range(1, 40) if abs(n*alpha0 - emp)/emp*100 < tol]
        hits_aMZ = [n for n in range(1, 40) if abs(n*alphaMZ - emp)/emp*100 < tol]
        print(f"    within {tol}%:  n*alpha(0) -> n in {hits_a0};  "
              f"n*alpha(MZ) -> n in {hits_aMZ}")
    # density of the n*alpha grid near the phase
    spacing0 = alpha0/emp*100
    print(f"    [grid spacing of n*alpha(0) near phi is {spacing0:.2f}% of phi -> "
          f"a 1%-tolerance window catches ~{2*1.0/spacing0:.1f} integers on average]")

# --------------------------------------------------------------------------
# Scale ratios a_u^2/a_d^2 and a_u^2/a_l^2.
# a_X = (1/3) sum sqrt(m_k) for sector X (since sum sqrt(m_k) = 3 a).
# So a_X^2 = (sum sqrt m)^2 / 9.  Ratios are scheme/scale dependent too.
# --------------------------------------------------------------------------
print()
print("=" * 78)
print("Scale ratios a_X^2 = (sum sqrt m_k)^2 / 9  (claimed a_u^2/a_d^2~35, a_u^2/a_l^2~72)")
print("=" * 78)

# lepton masses (GeV) -- charged leptons, essentially scale-independent (QED tiny)
m_lep = [0.5109989e-3, 0.1056584, 1.77686]

mass_sets = {
    "PDG_native_mixed": {
        "up":   [2.16e-3, 1.27, 172.69],
        "down": [4.67e-3, 93.4e-3, 4.18],
    },
    "MSbar_MZ": {
        "up":   [1.29e-3, 0.627, 171.7],
        "down": [2.75e-3, 0.0535, 2.89],
    },
    "MSbar_2GeV": {
        "up":   [2.16e-3, 1.24, 384.0],
        "down": [4.67e-3, 93.4e-3, 7.4],
    },
}

def a_sq(masses):
    return (np.sqrt(np.array(masses, float)).sum())**2 / 9.0

a_l2 = a_sq(m_lep)
print(f"  lepton a_l^2 = {a_l2:.6e} GeV   (scale-stable; QED running negligible)")
print()
print(f"  {'set':<20} {'a_u^2/a_d^2':>12} {'vs 35':>8} {'a_u^2/a_l^2':>12} {'vs 72':>8}")
print("  " + "-" * 64)
ratios = {}
for label, ms in mass_sets.items():
    au2 = a_sq(ms["up"]); ad2 = a_sq(ms["down"])
    r_ud = au2/ad2; r_ul = au2/a_l2
    ratios[label] = {"a_u2_over_a_d2": r_ud, "a_u2_over_a_l2": r_ul}
    print(f"  {label:<20} {r_ud:12.2f} {r_ud/35:8.2f} {r_ul:12.2f} {r_ul/72:8.2f}")

print()
print("  Interpretation: the scale ratios are NOT scale-stable (they move with")
print("  the heavy-quark masses, esp. m_t and m_b which dominate sum sqrt m).")
print("  '35' and '72' match only at the PDG mixed-scale convention -- the same")
print("  convention that makes Q match. They are not single-scale predictions.")

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
out = {
    "phi_emp": {"d": phi_d_emp, "u": phi_u_emp, "l": phi_l_emp},
    "phi_q_eq_Q_over_3": False,
    "scale_ratios": ratios,
    "a_l2_GeV": a_l2,
    "note": "All phases/ratios inherit RG/scheme dependence; integer*alpha is overfit-prone.",
}
with open("quark_phases_scales.json", "w") as f:
    json.dump(out, f, indent=2)
print()
print("Saved quark_phases_scales.json")
