#!/usr/bin/env python3
"""
integration_v58/bridge_vhiggs_cl7.py

THE v_Higgs BRIDGE COMPUTATION (requested 2026-05-24).

Question: in the *actual* Furey Cl(7)_even 8x8 representation, does the lepton
color-singlet mass scale a_l relate to the 28-dim L-grade (= Lambda^2 + Lambda^6
= so(8)) norm by the STRUCTURAL factor 28 = dim(L), as the conjecture
    v_Higgs = dim(L)^2 * a_l^2   <=>   sqrt(v) = 28 * a_l
requires?  Or does the rep's natural norm give sqrt(28) (quadrature) instead of
28 (linear) -- in which case the extra sqrt(28) is the Koide sqrt-mass linearity,
not something Cl(7) supplies on its own.

We build Cl(7) as 8x8 matrices (gamma_k^2 = -I, anticommuting), form the 28
L-grade blades (21 bivectors Lambda^2 + 7 hexads Lambda^6), confirm grade<->square
sign, check Frobenius orthogonality, and compute the two candidate norms:
  (A) vacuum-as-VECTOR-in-L : democratic sum over 28 directions, scale a_l
  (B) vacuum-as-BILINEAR-on-L : 28x28 = 784 components, each ~ a_l  (Yukawa = matrix)
and the lepton-specificity check (do D=35 d-quarks / D=63 u-quarks follow the
same rule -- the integration plan says they do NOT).
"""

import numpy as np

np.set_printoptions(linewidth=200, suppress=True)

# ---------------------------------------------------------------------------
# 1. Build Cl(7): seven anticommuting 8x8 real-structured matrices, gamma^2 = -I.
#    Standard tensor-of-Pauli construction on (C^2)^3 = C^8.
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X  = np.array([[0, 1], [1, 0]], dtype=complex)
Y  = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z  = np.array([[1, 0], [0, -1]], dtype=complex)

def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)

# Hermitian anticommuting set K_1..K_7 (K_a^2 = +I, {K_a,K_b}=2 delta).
K = [
    kron3(X, I2, I2),
    kron3(Y, I2, I2),
    kron3(Z, X, I2),
    kron3(Z, Y, I2),
    kron3(Z, Z, X),
    kron3(Z, Z, Y),
    kron3(Z, Z, Z),
]
# gamma_a = i K_a  =>  gamma_a^2 = -I, still anticommuting.
gamma = [1j * Ka for Ka in K]

def is_anticomm(a, b):
    return np.allclose(a @ b + b @ a, 0)

# sanity: gamma^2 = -I, anticommuting
assert all(np.allclose(g @ g, -np.eye(8)) for g in gamma), "gamma^2 != -I"
for a in range(7):
    for b in range(a + 1, 7):
        assert is_anticomm(gamma[a], gamma[b]), f"gamma_{a},gamma_{b} commute"
print("Cl(7) rep OK: 7 gammas, gamma^2 = -I, mutually anticommuting (8x8).")

def blade(idxs):
    M = np.eye(8, dtype=complex)
    for k in idxs:
        M = M @ gamma[k]
    return M

from itertools import combinations

# ---------------------------------------------------------------------------
# 2. Grades.  L = Lambda^2 (C(7,2)=21) + Lambda^6 (C(7,6)=7) = 28.  F = Lambda^4 (35).
# ---------------------------------------------------------------------------
biv  = [blade(c) for c in combinations(range(7), 2)]   # 21
hexa = [blade(c) for c in combinations(range(7), 6)]    # 7
four = [blade(c) for c in combinations(range(7), 4)]    # 35
Lgrade = biv + hexa
print(f"dim Lambda^2={len(biv)}, Lambda^6={len(hexa)}, L=dim(L)={len(Lgrade)}; F=Lambda^4={len(four)}")
assert len(Lgrade) == 28 and len(four) == 35

# grade -> square sign:  blades in L square to -I (complex structure), F to +I.
def sq_sign(M):
    s = M @ M
    if np.allclose(s, np.eye(8)):   return +1
    if np.allclose(s, -np.eye(8)):  return -1
    return 0
nL_neg = sum(1 for M in Lgrade if sq_sign(M) == -1)
nF_pos = sum(1 for M in four   if sq_sign(M) == +1)
print(f"L blades squaring to -I (complex structures): {nL_neg}/28")
print(f"F blades squaring to +I (real structures):    {nF_pos}/35")

# ---------------------------------------------------------------------------
# 3. Frobenius inner product on 8x8: <A,B> = Re Tr(A^dag B).  Are the 28
#    L-blades ORTHOGONAL?  (If yes -> vacuum-as-vector L2 norm = sqrt(28)*a_l.)
# ---------------------------------------------------------------------------
def fro_ip(A, B):
    return np.real(np.trace(A.conj().T @ B))

G = np.array([[fro_ip(A, B) for B in Lgrade] for A in Lgrade])
diag = np.diag(G)
offmax = np.max(np.abs(G - np.diag(diag)))
print(f"\nL-blade Gram matrix: diag entries all = {diag[0]:.1f} (=Tr I_8=8)?  "
      f"{np.allclose(diag, diag[0])}; max |off-diagonal| = {offmax:.2e}")
print("=> the 28 L-blades are Frobenius-ORTHOGONAL and equal-norm." if offmax < 1e-9
      else "=> L-blades NOT orthogonal.")

# ---------------------------------------------------------------------------
# 4. The two candidate vacua and their norms, in units of the per-direction
#    scale a_l.  Set a_l = 1; read off the factor multiplying a_l.
# ---------------------------------------------------------------------------
a_l = 1.0
# normalize each L-blade to unit Frobenius norm:
e = [M / np.sqrt(fro_ip(M, M)) for M in Lgrade]

# (A) vacuum as a VECTOR in L: democratic, each direction amplitude a_l.
#     V = a_l * sum_i e_i.  L2 (energy/Frobenius) norm:
V = a_l * sum(e)
normA_L2 = np.sqrt(fro_ip(V, V))
normA_L1 = a_l * len(e)            # L1 sum of amplitudes = coherent/linear sum
print(f"\n(A) vacuum-as-vector-in-L, scale a_l=1:")
print(f"    L2 (Frobenius) norm  sqrt(v)/a_l = {normA_L2:.4f}   (= sqrt(28) = {np.sqrt(28):.4f})")
print(f"    L1 (linear)   norm   sqrt(v)/a_l = {normA_L1:.4f}   (= 28)")

# (B) vacuum/mass as a BILINEAR (matrix) on the 28-dim L: 28x28 = 784 components,
#     each of democratic magnitude a_l.  Frobenius^2 = 784 a_l^2 => sqrt(v)/a_l = 28.
D = len(Lgrade)
Yuk = a_l * np.ones((D, D))         # democratic 28x28 (magnitudes a_l)
v_from_bilinear = np.sum(np.abs(Yuk) ** 2)   # Frobenius^2 = sum of |components|^2
print(f"\n(B) vacuum-as-bilinear-on-L (Yukawa is a {D}x{D} matrix), each comp a_l=1:")
print(f"    Frobenius^2 = #components * a_l^2 = {v_from_bilinear:.0f} = dim(L)^2 = {D*D}")
print(f"    => sqrt(v)/a_l = {np.sqrt(v_from_bilinear):.4f}  (= 28, the conjecture)")
print(f"    NOTE: this uses the natural L2 norm but on the {D}x{D}={D*D} ENTRIES,")
print(f"          not on the {D} directions. The '28^2' is dim(L)^2 = #matrix entries.")

# ---------------------------------------------------------------------------
# 5. Empirical anchor + lepton-specificity (does sqrt(v)=D*a hold for quarks?).
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("Empirical anchor and lepton-specificity")
print("=" * 70)
# charged leptons (GeV)
me, mmu, mtau = 0.51099895e-3, 105.6583755e-3, 1776.86e-3
S_lep = np.sqrt(me) + np.sqrt(mmu) + np.sqrt(mtau)
a_lep = S_lep / 3
vH = 246.21965
print(f"a_lepton = (sum sqrt m)/3 = {a_lep:.5f} GeV^1/2;  sqrt(vH) = {np.sqrt(vH):.4f}")
print(f"  sqrt(vH)/a_lepton = {np.sqrt(vH)/a_lep:.4f}   vs dim(L)=28   "
      f"(err {abs(np.sqrt(vH)/a_lep-28)/28*100:.3f}%)")

# If the rule were universal sqrt(v)=D*a, quark scales would be sqrt(vH)/D:
for name, D_sec in [("d-quark (D=35)", 35), ("u-quark (D=63)", 63)]:
    a_pred = np.sqrt(vH) / D_sec
    print(f"  predicted a_{name} = sqrt(vH)/{D_sec} = {a_pred:.5f} GeV^1/2 "
          f"(=> m-scale {a_pred**2*1e3:.3f} MeV)")
# quark Brannen scales (rough, running masses ~2 GeV): a_q = (sum sqrt m_q)/3
# light-quark MSbar @2GeV: m_d~4.7 MeV, m_s~93 MeV, m_b~4180 MeV (d-type);
md, ms, mb = 4.7e-3, 93.4e-3, 4180e-3
mu_, mc, mt = 2.16e-3, 1270e-3, 172500e-3
a_d = (np.sqrt(md)+np.sqrt(ms)+np.sqrt(mb))/3
a_u = (np.sqrt(mu_)+np.sqrt(mc)+np.sqrt(mt))/3
print(f"  actual a_d (running) ~ {a_d:.5f};  sqrt(vH)/a_d = {np.sqrt(vH)/a_d:.2f}  (vs 35)")
print(f"  actual a_u (running) ~ {a_u:.5f};  sqrt(vH)/a_u = {np.sqrt(vH)/a_u:.2f}  (vs 63)")
print("=> sqrt(v)=D*a is LEPTON-SPECIFIC; quarks do NOT follow it (as integration plan noted).")
