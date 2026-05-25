#!/usr/bin/env python3
"""
integration_v58/g2_koide_derivation.py

THE KEY BUILD (2026-05-24): derive the Koide split t² = (D - dimG₂)/D from the
G₂-content of the maximally-symmetric vacuum, in the ACTUAL algebra.

Chain to test:
  - maximal symmetry is the PRIMITIVE selection principle (v58 energy is flat on the
    vacuum manifold, `05_v58_vacuum_alignment.md`);
  - the lepton sector is L = so(8) (D=28), the mass-bearing grade (proven L=skew=so(8));
  - G₂ = Aut(𝕆) is the universal "inert core" (octonion automorphisms);
  - decompose L under G₂; the mass-splitting (Brannen phase ξ) lives in the NON-G₂
    complement; the maximally-mixed (uniform) vacuum gives it weight
        t² = dim(L \\ G₂)/dim(L) = (D - dimG₂)/D = (28-14)/28 = 1/2  ⇒ Koide Q = 2/3.

This script builds the octonions, constructs G₂ as the derivation algebra (verify dim=14),
decomposes so(7)/so(8) under G₂, and checks where the Brannen ξ-bivectors sit.
Self-contained (numpy).
"""
import numpy as np
from itertools import combinations
np.set_printoptions(linewidth=170, suppress=True)

# ===========================================================================
# 1. Octonions 𝕆 = ℝ ⊕ Im𝕆 (e0=1 real; e1..e7 imaginary).  Fano triples (cyclic).
# ===========================================================================
# Standard cyclic convention: line (i, i+1, i+3) mod 7 (on 1..7), e_i e_{i+1}=e_{i+3}.
TRIPLES = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
mult = np.zeros((8,8,8))  # mult[a,b,c] = coeff of e_c in e_a*e_b ; indices 0..7, e0=1
for a in range(8): mult[0,a,a]=1; mult[a,0,a]=1           # 1 is identity
for a in range(1,8): mult[a,a,0]=-1                        # e_a^2 = -1
for (i,j,k) in TRIPLES:
    for (x,y,z) in [(i,j,k),(j,k,i),(k,i,j)]:              # cyclic: positive
        mult[x,y,z]=1; mult[y,x,z]=-1                      # anti-cyclic: negative

def omul(u,v):  # octonion product of two 8-vectors
    return np.einsum('a,b,abc->c', u, v, mult)

# --- verify it's a valid (alternative, normed) octonion algebra ---
E=[np.eye(8)[a] for a in range(8)]
ok_norm=True; ok_alt=True
for a in range(1,8):
    for b in range(1,8):
        p=omul(E[a],E[b])
        if a!=b and abs(np.dot(p,p)-1)>1e-9: ok_norm=False    # |e_a e_b|=1
    # alternativity: e_a(e_a e_b) = (e_a e_a) e_b = -e_b
    for b in range(8):
        if np.linalg.norm(omul(E[a],omul(E[a],E[b])) + E[b])>1e-9: ok_alt=False
print(f"Octonion algebra: norm-basis OK={ok_norm}, alternative OK={ok_alt}")

# ===========================================================================
# 2. so(7) acting on Im𝕆 = ℝ⁷ (span e1..e7): the 21 antisymmetric generators.
# ===========================================================================
# Basis of so(7): A_{pq} (p<q in 1..7), (A_{pq})_{rs}=δ_pr δ_qs - δ_ps δ_qr  on the 7-dim Im.
idx7=list(range(1,8))                      # octonion imaginary labels 1..7
def Apq(p,q):                              # 7x7 antisymmetric (acting on Im, basis e1..e7)
    M=np.zeros((7,7)); ip=idx7.index(p); iq=idx7.index(q)
    M[ip,iq]=1; M[iq,ip]=-1; return M
so7_basis=[Apq(p,q) for p,q in combinations(idx7,2)]
print(f"dim so(7) = {len(so7_basis)} (expect 21)")

# extend a 7x7 Im-rotation to act as a derivation candidate on full 𝕆 (fix e0):
def lift(D7):                              # 7x7 on Im -> 8x8 on 𝕆 (e0 fixed)
    D8=np.zeros((8,8)); D8[1:,1:]=D7; return D8

# ===========================================================================
# 3. G₂ = derivations of 𝕆: D(uv) = (Du)v + u(Dv).  Solve for D in so(7).
# ===========================================================================
# Linear conditions on the 21 coefficients c_k (D = Σ c_k so7_basis[k]):
#   for all basis pairs (a,b): D(e_a e_b) - (D e_a) e_b - e_a (D e_b) = 0.
rows=[]
for k,B in enumerate(so7_basis):
    D8=lift(B)
    # derivation defect for this basis generator, over all (a,b) in Im:
    pass
# Build the linear map: defect_k(a,b,:) = D8(e_a e_b) - (D8 e_a)e_b - e_a (D8 e_b)
A_big=[]  # each row: contribution of generator k to one scalar equation component
for a in range(1,8):
    for b in range(1,8):
        # vector equation (8 comps) linear in c_k:
        for comp in range(8):
            row=[]
            for B in so7_basis:
                D8=lift(B)
                Dab = omul(D8@E[a], E[b]) + omul(E[a], D8@E[b])  # (Du)v+u(Dv)
                Dprod = D8 @ omul(E[a],E[b])                      # D(uv)
                row.append((Dprod - Dab)[comp])
            A_big.append(row)
A_big=np.array(A_big)                       # (Nrows, 21)
# Null space = derivation algebra G₂
u,s,vt=np.linalg.svd(A_big)
tol=1e-9
null_dim=int(np.sum(s<tol)) + (vt.shape[0]-len(s))
G2_coeffs=vt[len(s)-null_dim:] if null_dim>0 else vt[np.sum(s>=tol):]
# robust null-space:
ns_mask = np.concatenate([s<tol, np.ones(vt.shape[0]-len(s),bool)])
G2_coeffs = vt[ns_mask]
print(f"dim(derivation algebra G₂) = {G2_coeffs.shape[0]}  (expect 14)")

# ===========================================================================
# 4. Decompose so(7)=21 under G₂: G₂(14) ⊕ complement(7). Then so(8)=14⊕7⊕7.
# ===========================================================================
dim_so7=21; dim_G2=G2_coeffs.shape[0]; dim_comp_so7=dim_so7-dim_G2
print(f"so(7) = G₂({dim_G2}) ⊕ complement({dim_comp_so7})  [complement = the G₂ '7' (vector)]")
dim_L=28  # = dim so(8) = Λ²⊕Λ⁶ = 21+7
# so(8) under G₂ = so(7)⊕ℝ⁷ = [G₂(14)⊕7]⊕7 = 14 ⊕ 7 ⊕ 7
dim_nonG2_in_L = dim_L - dim_G2
print(f"L = so(8) = {dim_L} = G₂({dim_G2}) ⊕ nonG₂({dim_nonG2_in_L})  [nonG₂ = 7⊕7, the two vectors]")

# ===========================================================================
# 5. Where does the Brannen phase ξ live?  ξ-bivectors = {e01,e02,e12} (the ℍ-slice
#    imaginary units, LeptonComplexStructure).  In octonion-index terms the lepton
#    triple is an ASSOCIATIVE triple (a quaternion ⊂ 𝕆), e.g. {1,2,4} (a Fano line).
#    Test: is the su(2) of an associative triple inside G₂, or in the complement?
# ===========================================================================
def coeffvec(B7):  # express a 7x7 so(7) element in the so7_basis coords (21-vector)
    return np.array([np.sum(B7*so7_basis[k])/2.0 for k in range(len(so7_basis))])
# projector onto G₂ within so(7) (orthonormalize G2_coeffs first)
Gq,_=np.linalg.qr(G2_coeffs.T)             # 21 x 14 orthonormal columns spanning G₂
def frac_in_G2(B7):
    v=coeffvec(B7); v=v/np.linalg.norm(v)
    proj=Gq @ (Gq.T @ v)
    return np.dot(proj,proj)               # fraction of |v|^2 lying in G₂
# the Fano-line (associative) bivectors for line (1,2,4):
for (p,q) in [(1,2),(1,4),(2,4)]:
    print(f"  bivector A_{{{p}{q}}} (associative-triple): fraction in G₂ = {frac_in_G2(Apq(p,q)):.4f}")
# a NON-associative pair (not in any common simple structure), e.g. (1,3) [3 not with 1,2 in a line w/ 4]
for (p,q) in [(1,6),(2,7),(3,1)]:
    print(f"  bivector A_{{{p}{q}}} (generic):            fraction in G₂ = {frac_in_G2(Apq(p,q)):.4f}")

# ===========================================================================
# 6. The maximal-mixing result: uniform weight over L ⇒ t² = nonG₂ fraction.
# ===========================================================================
print()
print("=== Koide split from maximal mixing over L, G₂ as inert core ===")
for name,D in [("lepton (L)",28),("d-quark (F=Λ⁴)",35),("u-quark (L⊕F)",63)]:
    t2=(D-14)/D; Q=(1+2*t2)/3
    print(f"  {name:16s} D={D:2d}: t²=(D-dimG₂)/D=(D-14)/{D}={t2:.4f}  ⇒  Q=(1+2t²)/3 = {Q}")
print(f"  lepton t²=1/2 ⇒ Q=2/3 EXACTLY (D_lepton=28=2·dimG₂={2*dim_G2}).")
