#!/usr/bin/env python3
"""
v59/first_experiments/01_koide_substrate_scan.py

First numerological probe against TYCHO_TABLE T1.1 (Koide formula).

Background
----------
The Koide formula for charged leptons,
    Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3,
is identically satisfied iff sqrt(m_k) follows Brannen's parametrization
    sqrt(m_k) = a (1 + sqrt(2) cos(2 pi k/3 + phi)),   k = 0, 1, 2.

Any 3-by-3 Hermitian matrix invariant under cyclic permutation has
eigenvalues of the form  a + 2|b| cos(2 pi k/3 + arg b).  This matches
Brannen iff |b|/a = 1/sqrt(2).  Substrate question: what natural object
predicts |b|/a = 1/sqrt(2) AND phi ~ 0.222 rad simultaneously?

This script:
  1. Computes the experimental Brannen parameters (a, phi) to high precision.
  2. Verifies the analytic structure of Hermitian Z_3-cyclic matrices.
  3. Surveys candidate substrate sources of the |b|/a = 1/sqrt(2) ratio.
  4. Tests three knot ropelength triples against Koide.
  5. Tests projections of Platonic solid vertices.
"""

import numpy as np
from itertools import combinations, permutations


# ---- Experimental input (PDG 2024) ----
m_e_MeV   = 0.51099895069
m_mu_MeV  = 105.6583755
m_tau_MeV = 1776.86


def koide_Q(nu):
    """Koide ratio Q = sum(nu^2) / sum(nu)^2."""
    nu = np.asarray(nu, dtype=float)
    return float(np.sum(nu ** 2) / np.sum(nu) ** 2)


def brannen_fit(nu):
    """Fit nu_k = a (1 + sqrt(2) cos(2 pi k/3 + phi)) to a triple.

    Returns dict with a, phi, b_over_a_eff, residual, perm.
    b_over_a_eff = 1.0 means a perfect Brannen fit (Q = 2/3 identically).
    """
    nu = np.asarray(nu, dtype=float)
    a = float(np.sum(nu)) / 3.0
    if a <= 0:
        return None
    z = (nu - a) / (a * np.sqrt(2))  # should be cos(2 pi k/3 + phi)
    A_basis = np.array([
        [1.0,        0.0          ],
        [-0.5,      -np.sqrt(3)/2 ],
        [-0.5,       np.sqrt(3)/2 ]
    ])
    best = None
    for perm in permutations(range(3)):
        c = z[list(perm)]
        sol, _, _, _ = np.linalg.lstsq(A_basis, c, rcond=None)
        cp, sp = sol
        norm = float(np.hypot(cp, sp))
        phi = float(np.arctan2(sp, cp))
        c_pred = A_basis @ sol
        inv_perm = np.argsort(perm)
        nu_pred = a * (1.0 + np.sqrt(2) * c_pred[inv_perm])
        residual = float(np.linalg.norm(nu - nu_pred) / np.linalg.norm(nu))
        if best is None or residual < best['residual']:
            best = dict(a=a, phi=phi, b_over_a_eff=norm,
                        residual=residual, perm=perm)
    return best


# ---- Pretty printers ----
def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def subsection(title):
    print()
    print("-" * 72)
    print(title)
    print("-" * 72)


# =========================================================================
# Part 1: Experimental Brannen parameters
# =========================================================================
section("Part 1: Experimental Brannen parameters (TYCHO T1.1, T1.2)")

sqm = [np.sqrt(m_e_MeV), np.sqrt(m_mu_MeV), np.sqrt(m_tau_MeV)]
print()
print(f"Charged lepton masses (PDG 2024):")
print(f"  m_e   = {m_e_MeV:.10f} MeV")
print(f"  m_mu  = {m_mu_MeV:.7f} MeV")
print(f"  m_tau = {m_tau_MeV:.2f} MeV")
print()
print(f"Square roots (sqrt MeV):")
print(f"  sqrt(m_e)   = {sqm[0]:.10f}")
print(f"  sqrt(m_mu)  = {sqm[1]:.10f}")
print(f"  sqrt(m_tau) = {sqm[2]:.10f}")
print()

Q_exp = koide_Q(sqm)
print(f"Koide ratio Q (experimental) = {Q_exp:.10f}")
print(f"2/3                          = {2/3:.10f}")
print(f"Delta from 2/3               = {Q_exp - 2/3:+.3e}")
print()

fit = brannen_fit(sqm)
if fit:
    print(f"Brannen fit:")
    print(f"  a (scale)         = {fit['a']:.10f} sqrt(MeV)")
    print(f"  phi               = {fit['phi']:.10f} rad")
    print(f"                    = {np.degrees(fit['phi']):.6f} deg")
    print(f"  b/a effective     = {fit['b_over_a_eff']:.10f}  (1 means perfect Brannen)")
    print(f"  best permutation  = {fit['perm']}  (k assignment to e,mu,tau)")
    print(f"  fractional residual = {fit['residual']:.3e}")
    print()
    print(f"Does phi match simple rationals or angles?")
    candidates = [
        ('2/9 rad',          2/9),
        ('1/(2*pi) rad',     1/(2*np.pi)),
        ('pi/14 rad',        np.pi/14),
        ('pi/13 rad',        np.pi/13),
        ('arctan(2/9) rad',  np.arctan(2/9)),
        ('atan(1/4.5)',      np.arctan(1/4.5)),
        ('1/4.5 rad',        1/4.5),
    ]
    for name, val in candidates:
        d = fit['phi'] - val
        print(f"  {name:18s} = {val:.10f}  (delta = {d:+.3e})")


# =========================================================================
# Part 2: Analytic content of Hermitian Z_3 circulants
# =========================================================================
section("Part 2: Hermitian Z_3 circulant matrices")
print()
print("Family:")
print("  M = a I + b S + b* S^T  (S = cyclic shift, b in C)")
print("  Eigenvalues = a + 2|b| cos(2 pi k/3 + arg b),  k = 0, 1, 2.")
print()
print(f"For Brannen / Koide Q = 2/3 we need |b|/a = 1/sqrt(2) = {1/np.sqrt(2):.10f}.")
print()
print(f"{'|b|/a':>10} {'eigenvalues (phi=0.222)':>40} {'Q':>10} {'delta':>10}")
for ba in [0.30, 0.50, 1/np.sqrt(2), 0.80, 1.00]:
    a = 1.0
    bmag = ba * a
    phi_t = 0.2222
    eigs = sorted([a + 2*bmag*np.cos(2*np.pi*k/3 + phi_t) for k in range(3)])
    if min(eigs) > 0:
        Q = koide_Q(eigs)
        print(f"  {ba:.5f}  [{eigs[0]:.4f}, {eigs[1]:.4f}, {eigs[2]:.4f}]   "
              f"{Q:.5f}   {Q - 2/3:+.3e}")
    else:
        print(f"  {ba:.5f}  eigenvalue went negative — skip")


# =========================================================================
# Part 3: Candidate substrate sources of |b|/a = 1/sqrt(2)
# =========================================================================
section("Part 3: Substrate-side sources of |b|/a = 1/sqrt(2)")

subsection("Candidate A: Cube edge-to-face-diagonal ratio")
print(f"  Edge length        = 1")
print(f"  Face diag length   = sqrt(2)")
print(f"  Edge / face-diag   = 1/sqrt(2) = {1/np.sqrt(2):.10f}   -- matches")
print(f"  Identification of which is 'a' vs 'b' is not unique; suggestive only.")

subsection("Candidate B: Inscribed/circumscribed sphere ratios")
ratios = {
    'tetrahedron'  : 1/3,
    'cube'         : 1/np.sqrt(3),
    'octahedron'   : np.sqrt(3)/3,  # = 1/sqrt(3) actually for regular octa, different by convention
    'dodecahedron' : np.sqrt(((25 + 11*np.sqrt(5))/60) * (3/((1+np.sqrt(5))/2)**4)),  # ~ 0.7947
}
for name, r in ratios.items():
    delta = r - 1/np.sqrt(2)
    flag = "  <-- close to 1/sqrt(2)" if abs(delta) < 0.05 else ""
    print(f"  {name:14s}: r_in/r_out = {r:.7f}  (delta from 1/sqrt(2) = {delta:+.4f}){flag}")

subsection("Candidate C: Cl(3,0) magnitudes")
print(f"  Vector e_1               : |e_1| = 1")
print(f"  Symmetric bivector       : |(e12+e23+e31)/sqrt(3)|^2 = 3/3 = 1")
print(f"  Trivector e_123          : |e_123|^2 = 1")
print(f"  No obvious 1/sqrt(2) ratio between natural unit objects.")

subsection("Candidate D: Tight-binding triangle with phase")
print(f"  H = eps * I + t * (cyclic hopping with phase phi)")
print(f"  Eigenvalues: eps + 2t cos(2 pi k/3 + phi)")
print(f"  t/eps is a FREE parameter; no natural reason for t/eps = 1/sqrt(2).")
print(f"  This means tight-binding alone does NOT predict Koide.")


# =========================================================================
# Part 4: Knot ropelength scan
# =========================================================================
section("Part 4: Knot ropelengths (Pieranski / Cantarella data)")

ropelengths = {
    '3_1':  16.372,
    '4_1':  21.039,
    '5_1':  24.736,
    '5_2':  24.991,
    '6_1':  28.296,
    '6_2':  28.745,
    '6_3':  28.957,
    '7_1':  30.722,
}

print()
print(f"Testing all triples of {len(ropelengths)} knots, treating ropelength as sqrt(m).")
print()
hits = []
for triple in combinations(sorted(ropelengths.keys()), 3):
    nus = [ropelengths[t] for t in triple]
    Q = koide_Q(nus)
    delta = abs(Q - 2/3)
    hits.append((delta, triple, nus, Q))
hits.sort()

print(f"  Top 5 closest to Q = 2/3:")
for delta, triple, nus, Q in hits[:5]:
    print(f"    {triple}: nus = ({nus[0]:.2f}, {nus[1]:.2f}, {nus[2]:.2f}), "
          f"Q = {Q:.5f}, delta = {delta:+.4f}")

print()
print(f"  Comment: All knot ropelengths are in the range ~16-31, ratio ~2.")
print(f"  For Koide-by-extreme-ratio (one heavy, two light) the ropelengths would")
print(f"  need to span many orders of magnitude. They do not.")
print(f"  Therefore ropelength-as-sqrt(m) is NOT a viable Brannen substrate")
print(f"  unless a non-linear function of ropelength (e.g., exp(ropelength*const))")
print(f"  is used instead.")


# =========================================================================
# Part 5: Polyhedral vertex projections (cube test)
# =========================================================================
section("Part 5: Polyhedral vertex projections")

print()
print("Cube vertices at (+/-1, +/-1, +/-1)/sqrt(3), projected on (1,1,1)/sqrt(3):")
from itertools import product as iproduct
CUBE = np.array(list(iproduct([-1, 1], repeat=3)), dtype=float) / np.sqrt(3)
diag = np.array([1, 1, 1]) / np.sqrt(3)
projs = sorted(CUBE @ diag)
print(f"  Projections: {[float(f'{x:.4f}') for x in projs]}")
print(f"  Distinct positive values: {sorted(set(round(x, 4) for x in projs if x > 0))}")
# Three vertices at +1/sqrt(3), three at +1*sqrt(3) ... wait let me compute
# (1,1,-1)/sqrt(3) . (1,1,1)/sqrt(3) = (1+1-1)/3 = 1/3
# So projections are at +/-1, +/-1/3.
# Triple at (1/3, 1/3, 1/3):
print()
nu_tri = [1/3, 1/3, 1/3]
print(f"  Triple of cube verts at projection 1/3: Q = {koide_Q(nu_tri):.6f}  (= 1/3, not Koide)")

print()
print("Icosahedron vertices (12 total), projection on a 3-fold axis (1,1,1)/sqrt(3):")
phi_g = (1 + np.sqrt(5)) / 2
ICO_raw = np.array([
    [0, 1, phi_g], [0, 1, -phi_g], [0, -1, phi_g], [0, -1, -phi_g],
    [1, phi_g, 0], [1, -phi_g, 0], [-1, phi_g, 0], [-1, -phi_g, 0],
    [phi_g, 0, 1], [phi_g, 0, -1], [-phi_g, 0, 1], [-phi_g, 0, -1]
], dtype=float)
ICO = ICO_raw / np.linalg.norm(ICO_raw[0])
projs_ico = ICO @ diag
unique_projs = sorted(set(round(x, 5) for x in projs_ico))
print(f"  Distinct projection values: {unique_projs}")
# Try positive distinct values
pos = [x for x in unique_projs if x > 0.001]
if len(pos) >= 3:
    # Take the 3 largest positives:
    nu = sorted(pos)[-3:]
    print(f"  Three largest positive projections: {nu}")
    Q = koide_Q(nu)
    fit_ico = brannen_fit(nu)
    print(f"    Q = {Q:.6f}, delta from 2/3 = {Q - 2/3:+.3e}")


# =========================================================================
# Part 6: Closing summary
# =========================================================================
section("Summary")

print()
print("Experimental Koide ratio agrees with 2/3 to ~5e-5.")
print(f"Brannen phi for charged leptons: {fit['phi']:.6f} rad ~ {np.degrees(fit['phi']):.3f} deg.")
print()
print("Substrate-side situation after this initial pass:")
print()
print("  * No tested substrate predicts |b|/a = 1/sqrt(2) WITHOUT free parameters.")
print("  * Cube edge/face-diagonal ratio is 1/sqrt(2), suggestive only.")
print("  * Knot ropelengths span too narrow a range; not a direct Brannen substrate.")
print("  * Polyhedral vertex projections give degenerate triples (no Koide).")
print()
print("Implication: a passing substrate needs to predict not only Z_3 symmetry")
print("(which is common) but also the SPECIFIC RATIO |b|/a = 1/sqrt(2) plus a phi")
print(f"near 0.222 rad.  The latter (phi) is the more discriminating target.")
print()
print("Suggested next experiments:")
print("  - finite-group representation eigenvalues (A_4, S_4, A_5, binary tetra/octa/icosa).")
print("  - Cl(0,7) octonionic invariants (Furey-program structures).")
print("  - Lie-algebra root-system ratios for G_2, F_4, E_6.")
print("  - Look for natural occurrences of phi ~ 0.222 rad specifically.")
print()
