#!/usr/bin/env python3
"""
v59/synthesis/obe_geometric_6_7.py

Test OBE bridge variants 6 and 7 (the Yukawa-free, geometric readings).

Option 7 (Kaluza-Klein separation): a MASSLESS field on R^3 x M_int (compact
  internal manifold) gives, at long range, a single 4D zero mode -> 1/r force,
  while the massive KK tower gives short-range (Yukawa) corrections that build up
  the full higher-D 1/r^{d+1} at short range. The zero-mode coupling is DILUTED
  by 1/Vol_int. => the radial law (1/r^2) and the bridge count (Vol_int) live in
  SEPARATE factors of K. So 1/r^2 and the 28-dilution are compatible.

Option 6 (rectangular generation x ambient overlap): leptons = 3 modes (rows)
  overlapping D_ambient internal directions (columns). Masses = singular values
  of a 3 x D overlap matrix W; no Yukawa coupling. The bridge becomes
      sqrt(v) = (D / N_gen) * sum(sqrt(m)),
  which holds at 0.03% ONLY for the lepton ambient D=28, and fails for the quark
  ambients D=35, 63 -- the lepton-specificity of 03_higgs_bridge_result.md,
  in SVD language.
"""

import numpy as np

# ===========================================================================
# OPTION 7 : Kaluza-Klein separation  (radial law vs internal-volume count)
# ===========================================================================
print("=" * 70)
print("OPTION 7: KK separation -- 1/r zero mode + Yukawa tower, 1/Vol dilution")
print("=" * 70)


def kk_green(r, R, d, nmax=40):
    """Massless propagator on R^3 x T^d (flat internal torus, circumference R),
    same internal point:  G(r) = (1/Vol) sum_n exp(-|n|/R * r) / (4 pi r),
    Vol = R^d, n in Z^d. Returns G(r)."""
    Vol = R**d
    ns = np.arange(-nmax, nmax + 1)
    grids = np.meshgrid(*([ns] * d), indexing="ij")
    absn = np.sqrt(sum(g.astype(float)**2 for g in grids)).ravel()
    m = absn / R                                  # KK masses |n|/R
    # sum over modes; broadcast over r
    G = np.zeros_like(r)
    for mi in m:
        G += np.exp(-mi * r) / (4 * np.pi * r)
    return G / Vol


def slope(r, f, rlo, rhi):
    msk = (r >= rlo) & (r <= rhi) & (f > 0)
    return np.polyfit(np.log(r[msk]), np.log(f[msk]), 1)[0]


r = np.linspace(0.02, 40.0, 8000)
R = 1.0                                            # internal size
for d in (1, 2):
    G = kk_green(r, R, d)
    s_long = slope(r, G, 14.0, 30.0)               # r >> R: expect -1 (4D, zero mode)
    s_short = slope(r, G, 0.03, 0.2)               # r << R: expect -(d+1) (full dim)
    print(f"\n internal dim d={d}, size R={R}:")
    print(f"   long-range  (r>>R) slope = {s_long:+.3f}   (expect -1: 4D 1/r zero mode)")
    print(f"   short-range (r<<R) slope = {s_short:+.3f}   (expect -{d+1}: full (3+{d})D)")

# 1/Vol dilution of the zero-mode coupling: vary R (=> Vol=R^d), measure G*4pi*r at
# r >> R (r=35, so even R=4 has e^{-r/R}=e^{-8.75}~1e-4 negligible).
print("\n zero-mode coupling ~ 1/Vol_int  (the 'spread over the ambient' dilution):")
for R in (1.0, 2.0, 4.0):
    G = kk_green(r, R, d=1)
    coup = (G * 4 * np.pi * r)[np.argmin(np.abs(r - 35.0))]   # -> 1/Vol at r >> R
    print(f"   R={R} (Vol={R}):  long-range coupling 4pi r G = {coup:.4f}   (1/Vol = {1/R:.4f})")

print("\n => radial law (1/r) is the 4D zero mode; the internal volume only sets the")
print("    COUPLING (1/Vol per leg). Bridge count (Vol~28) and 1/r^2 are independent.")

# ===========================================================================
# OPTION 6 : rectangular generation x ambient overlap  (Yukawa-free)
# ===========================================================================
print("\n" + "=" * 70)
print("OPTION 6: rectangular 3 x D overlap -- masses = singular values")
print("=" * 70)

# PDG-ish masses (MeV)
m_lep = np.array([0.51099895, 105.6583755, 1776.86])           # e, mu, tau
m_dq = np.array([4.7, 93.0, 4180.0])                            # d, s, b  (MS-bar)
m_uq = np.array([2.2, 1270.0, 172500.0])                        # u, c, t
v_higgs = 246220.0                                              # MeV
sqrt_v = np.sqrt(v_higgs)
N_gen = 3

def koide(m):
    return m.sum() / (np.sqrt(m).sum())**2

# Brannen shape for leptons
s = np.sqrt(m_lep)
a_l = s.mean()
Q = koide(m_lep)
# phase: sqrt(m)/a - 1 = sqrt(2) cos(phi + 2pi k/3); take the largest (tau, k=0)
phi = np.arccos(((s.max() / a_l) - 1.0) / np.sqrt(2.0))
print(f"\n leptons: Q = {Q:.5f} (2/3={2/3:.5f}),  phi = {phi:.4f} (2/9={2/9:.4f}),"
      f"  |phi-2/9| = {abs(phi-2/9):.4f}")

print("\n bridge in overlap form:  sqrt(v) =? (D/N_gen) * sum(sqrt(m))")
for name, m, D in [("lepton", m_lep, 28), ("d-quark", m_dq, 35), ("u-quark", m_uq, 63)]:
    pred = (D / N_gen) * np.sqrt(m).sum()
    ratio = sqrt_v / (np.sqrt(m).sum() / N_gen)        # = sqrt(v)/a  -> should be D
    print(f"   {name:8s} D={D:2d}:  sqrt(v)/a = {ratio:6.2f}  (law wants {D})"
          f"   {'MATCH' if abs(ratio-D)/D < 0.01 else 'fails (%.1fx)' % (max(ratio,D)/min(ratio,D))}")

# SVD realization: build a 3 x 28 overlap W with singular values = sqrt(m_lep),
# democratic over the 28 ambient columns (3 orthonormal Z3-phase column patterns).
D = 28
rng = np.random.default_rng(0)
# 3 orthonormal column patterns over 28 directions
X = rng.standard_normal((D, 3))
Qc, _ = np.linalg.qr(X)                 # Qc: 28x3 orthonormal columns
W = (np.sqrt(m_lep)[:, None] * Qc.T)    # 3x28, rows = s_k * pattern_k
sv = np.linalg.svd(W, compute_uv=False)
print("\n SVD of a democratic 3x28 overlap W (columns = ambient L-directions):")
print(f"   singular values  = {np.sort(sv)[::-1]}")
print(f"   target sqrt(m)   = {np.sort(np.sqrt(m_lep))[::-1]}")
print(f"   max rel. err     = {np.max(np.abs(np.sort(sv)-np.sort(np.sqrt(m_lep)))/np.sqrt(m_lep)):.2e}")
print(f"   nuclear norm ||W||_* = {sv.sum():.4f} = sum sqrt(m) = {np.sqrt(m_lep).sum():.4f}")
print(f"   => sqrt(v) = (28/3)||W||_* = {(28/3)*sv.sum():.2f}  vs  sqrt(v_obs) = {sqrt_v:.2f}"
      f"   ({100*abs((28/3)*sv.sum()-sqrt_v)/sqrt_v:.3f}%)")

print("\n HONEST NOTE: the SVD recovers whatever singular values are inserted, so")
print(" Option 6 RELOCATES the Brannen shape (t, phi) into the generation rows --")
print(" it does not derive them. What it DOES make structural: the ambient leg")
print(" count D and the lepton-specific bridge sqrt(v)=(D/3)*sum(sqrt m), D=28 only.")
