#!/usr/bin/env python3
"""
v59/furey_construction/05_quark_sector.py

Variant E: Quark sector via three-generation extension and triality.

The lepton-sector analysis (steps 1-11) gives three generations from the
Z_3 cyclic structure on Cl(3,1) (or equivalently, the Z_3 triality of
Spin(8) acting on the octonion factor).  In the Furey construction, the
QUARK sector arises from the COLOR SU(3) acting on the octonion imaginary
units.

Quark masses (PDG 2024, MS-bar at characteristic scales):
  u: 2.16  MeV       d: 4.67  MeV       s: 93.4  MeV
  c: 1.27  GeV       b: 4.18  GeV       t: 172.69 GeV

In contrast to the clean Koide structure for leptons, the quark mass
spectrum is messier (running, QCD effects).  But the same algebraic
machinery should apply.

Hypothesis: each quark flavor (u, d, s, c, b, t) corresponds to a triality
branch × color triplet structure, with mass eigenvalues derived from
Brannen-like expressions in the extended algebra.

This script:
  1. Applies the Brannen formula to quark mass ratios.
  2. Checks if any quark-triplet satisfies Koide-like identities.
  3. Identifies what specific algebraic content is needed for quark masses.
"""

import numpy as np
import math
from itertools import combinations, permutations

print("="*72)
print("Variant E: Quark sector analysis")
print("="*72)


# =========================================================================
# Part A: Quark masses and basic ratios
# =========================================================================
print()
print("-"*72)
print("Part A: Quark masses")
print("-"*72)

quark_masses = {
    'u': 2.16e-3,    # GeV
    'd': 4.67e-3,
    's': 93.4e-3,
    'c': 1.27,
    'b': 4.18,
    't': 172.69,
}
# Generation pairs:
gens_up = ['u', 'c', 't']
gens_down = ['d', 's', 'b']

print("\nQuark masses (PDG 2024, GeV):")
for q, m in quark_masses.items():
    print(f"  m_{q} = {m:.5g} GeV")


# =========================================================================
# Part B: Koide ratio for up-type and down-type quarks
# =========================================================================
print()
print("-"*72)
print("Part B: Koide ratios for quark generations")
print("-"*72)

def koide_Q(masses):
    nu = np.sqrt(np.array(masses))
    return float(np.sum(nu**2) / np.sum(nu)**2)

# Up-type: (u, c, t)
m_up = [quark_masses[q] for q in gens_up]
Q_up = koide_Q(m_up)
print(f"\nUp-type quarks (u, c, t):")
print(f"  masses: {m_up}")
print(f"  Q_up = {Q_up:.6f}")
print(f"  2/3 = {2/3:.6f}")
print(f"  Δ from 2/3: {Q_up - 2/3:+.6f}")

# Down-type: (d, s, b)
m_down = [quark_masses[q] for q in gens_down]
Q_down = koide_Q(m_down)
print(f"\nDown-type quarks (d, s, b):")
print(f"  masses: {m_down}")
print(f"  Q_down = {Q_down:.6f}")
print(f"  Δ from 2/3: {Q_down - 2/3:+.6f}")

# Combined Koide:
m_all = [quark_masses[q] for q in ['u', 'd', 's', 'c', 'b', 't']]
print(f"\nQuark mass ratios (relative to up):")
for q in ['u', 'd', 's', 'c', 'b', 't']:
    print(f"  m_{q}/m_u = {quark_masses[q]/quark_masses['u']:.4f}")


# =========================================================================
# Part C: Brannen parametrization for quarks
# =========================================================================
print()
print("-"*72)
print("Part C: Brannen parametrization for up-type and down-type quarks")
print("-"*72)

def brannen_fit(masses):
    sqm = np.sqrt(np.array(masses))
    a = np.sum(sqm) / 3
    if a <= 0:
        return None
    z = (sqm - a) / (a * np.sqrt(2))
    A = np.array([[1.0, 0], [-0.5, -np.sqrt(3)/2], [-0.5, np.sqrt(3)/2]])
    best = None
    for perm in permutations(range(3)):
        c = z[list(perm)]
        sol, _, _, _ = np.linalg.lstsq(A, c, rcond=None)
        cp, sp = sol
        norm = np.hypot(cp, sp)
        phi = np.arctan2(sp, cp)
        c_pred = A @ sol
        inv_perm = np.argsort(perm)
        sqm_pred = a * (1 + np.sqrt(2) * c_pred[inv_perm])
        resid = np.linalg.norm(sqm - sqm_pred) / np.linalg.norm(sqm)
        if best is None or resid < best[0]:
            best = (resid, a, phi, norm)
    return best

print(f"\nUp-type quarks Brannen fit:")
fit_up = brannen_fit(m_up)
print(f"  a (scale) = {fit_up[1]:.5g} sqrt(GeV)")
print(f"  phi       = {fit_up[2]:.6f} rad")
print(f"  b/a eff   = {fit_up[3]:.6f}  (1 = perfect Brannen)")
print(f"  residual  = {fit_up[0]:.3e}")

print(f"\nDown-type quarks Brannen fit:")
fit_down = brannen_fit(m_down)
print(f"  a (scale) = {fit_down[1]:.5g} sqrt(GeV)")
print(f"  phi       = {fit_down[2]:.6f} rad")
print(f"  b/a eff   = {fit_down[3]:.6f}")
print(f"  residual  = {fit_down[0]:.3e}")


# =========================================================================
# Part D: Cross-quark identities
# =========================================================================
print()
print("-"*72)
print("Part D: cross-quark / lepton identities")
print("-"*72)

# Various conjectured relations:
# Gatto-Sartori-Tonin (1968): tan(theta_C) ~ sqrt(m_d / m_s)
import math
theta_C_exp = 13.04 * math.pi / 180  # Cabibbo angle
sqrt_ratio = math.sqrt(quark_masses['d'] / quark_masses['s'])

print(f"\nGatto-Sartori-Tonin relation:")
print(f"  tan(theta_C) (empirical) = {math.tan(theta_C_exp):.6f}")
print(f"  sqrt(m_d/m_s)            = {sqrt_ratio:.6f}")
print(f"  ratio                     = {math.tan(theta_C_exp) / sqrt_ratio:.4f}")

# Koide-like relation across all 6 quarks?
sqm_all = [math.sqrt(m) for m in m_all]
Q_all = sum(s**2 for s in sqm_all) / sum(sqm_all)**2
print(f"\nKoide-like ratio over all 6 quarks:")
print(f"  Q = sum(sqrt(m))^2 / sum(sqrt(m))^2 ... wait this is wrong")
print(f"  Q = sum(m) / sum(sqrt(m))^2 = {sum(m_all) / sum(sqm_all)**2:.6f}")


# =========================================================================
# Part E: Structural interpretation
# =========================================================================
print()
print("="*72)
print("Summary")
print("="*72)
print(f"""
Quark sector findings:

  Up-type Koide Q_up = {Q_up:.4f}        (cf. 2/3 = {2/3:.4f})
  Down-type Koide Q_down = {Q_down:.4f}   (cf. 2/3 = {2/3:.4f})

  Q_up is NOT 2/3 (gap = {Q_up - 2/3:.3f}).
  Q_down is NOT 2/3 (gap = {Q_down - 2/3:.3f}).

  Up-type Brannen phi = {fit_up[2]:.4f} rad
  Down-type Brannen phi = {fit_down[2]:.4f} rad
  Charged-lepton Brannen phi = 0.2222 rad
  These are all distinct.

Implication: the quark sectors do NOT satisfy Koide Q = 2/3.  The Brannen
parametrization formally applies but with different phi values for up-type
and down-type.

The lepton Koide identity is a special feature of the lepton sector.  In
the Furey program, this is because leptons are SU(3)_c SINGLETS (no
color), while quarks are TRIPLETS.  The Koide identity reflects the
single-generation-singlet structure; for triplets, color mixing modifies
the eigenvalue spectrum.

Quark masses also include large QCD running corrections (MS-bar masses at
different scales).  A clean algebraic relation may require running to a
common scale or restricting to pole masses (for c, b, t).

CONCLUSION: the quark sector is NOT a clean structural test of the
G_2/Spin(7) framework as the lepton sector is.  Quark mass predictions
require additional content (color SU(3), running, mixing with the Higgs
sector).

Variant E: COMPLETE (negative for tight quark Koide identities).
""")
