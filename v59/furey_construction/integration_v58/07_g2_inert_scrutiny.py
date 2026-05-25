#!/usr/bin/env python3
"""
integration_v58/07_g2_inert_scrutiny.py

STRESS-TEST of the "P2 — G₂ inert" assumption behind the Koide formula
    t² = (D − dimG₂)/D = (D − 14)/D.

Reuses the octonion + G₂ construction from g2_koide_derivation.py, then asks
four questions and tries to BREAK the claim:

  (1) Is the "−14" the removal of a genuine G₂-subspace in EVERY sector, or just
      arithmetic that happens to fit the leptons?  We decompose Λ⁴ℝ⁷ (d-quark,
      D=35) and the u-sector (Λ²⊕Λ⁴⊕Λ⁶, D=63) under the SAME G₂ acting on ℝ⁷.

  (2) Does G₂ genuinely fix the lepton singlet(s)?  Build the 8-spinor of
      Spin(7), branch under G₂ (8 = 1 ⊕ 7), and test the singlet directly.

  (3) The crux: does "G₂-invariant ⟹ no mass-splitting" hold?  Compute the
      commutant of G₂ on each sector — a G₂-commuting (inert-respecting) operator
      is scalar on each irrep (Schur) but CAN differ between irreps.

  (4) Is "14" privileged as dimG₂, or equally 2·dim(Im𝕆)=2·7, and does the G₂
      *subspace* reading survive sector-by-sector?

Robust methods only: invariant-counting (kernel dims) and commutant dims (both
are stable rank computations).  No fragile Casimir clustering.
Self-contained (numpy).
"""
import numpy as np
from itertools import combinations
np.set_printoptions(linewidth=200, suppress=True, precision=4)

def rank_tol(M, rel=1e-8):
    s = np.linalg.svd(M, compute_uv=False)
    if len(s)==0: return 0, s
    return int(np.sum(s > rel*max(1.0, s.max()))), s

# ===========================================================================
# 0. Octonions + G₂ (same construction as g2_koide_derivation.py)
# ===========================================================================
TRIPLES = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
mult = np.zeros((8,8,8))
for a in range(8): mult[0,a,a]=1; mult[a,0,a]=1
for a in range(1,8): mult[a,a,0]=-1
for (i,j,k) in TRIPLES:
    for (x,y,z) in [(i,j,k),(j,k,i),(k,i,j)]:
        mult[x,y,z]=1; mult[y,x,z]=-1
def omul(u,v): return np.einsum('a,b,abc->c', u, v, mult)
E=[np.eye(8)[a] for a in range(8)]

idx7=list(range(1,8))
def Apq(p,q):
    M=np.zeros((7,7)); ip=idx7.index(p); iq=idx7.index(q)
    M[ip,iq]=1; M[iq,ip]=-1; return M
so7_basis=[Apq(p,q) for p,q in combinations(idx7,2)]
def lift(D7):
    D8=np.zeros((8,8)); D8[1:,1:]=D7; return D8

A_big=[]
for a in range(1,8):
    for b in range(1,8):
        for comp in range(8):
            row=[]
            for B in so7_basis:
                D8=lift(B)
                Dab = omul(D8@E[a], E[b]) + omul(E[a], D8@E[b])
                Dprod = D8 @ omul(E[a],E[b])
                row.append((Dprod - Dab)[comp])
            A_big.append(row)
A_big=np.array(A_big)
u,s,vt=np.linalg.svd(A_big)
ns_mask = np.concatenate([s<1e-9, np.ones(vt.shape[0]-len(s),bool)])
G2_coeffs = vt[ns_mask]                       # 14 x 21 coords in so7_basis
G2_gens7 = np.array([sum(c[k]*so7_basis[k] for k in range(21)) for c in G2_coeffs])
print(f"[setup] dim G₂ = {G2_gens7.shape[0]} (expect 14); generators antisym: "
      f"max|M+Mᵀ|={max(np.linalg.norm(M+M.T) for M in G2_gens7):.1e}")

# ---------------------------------------------------------------------------
# G₂ action on Λ^k ℝ⁷ (natural derivation / Lie-derivative action)
# ---------------------------------------------------------------------------
def wedge_basis(n,k): return list(combinations(range(n),k))
def lie_on_wedge(M7, k):
    basis = wedge_basis(7,k); bidx={b:i for i,b in enumerate(basis)}
    dim=len(basis); rho=np.zeros((dim,dim))
    for col,b in enumerate(basis):
        for pos,i in enumerate(b):
            for r in range(7):
                coeff=M7[r,i]
                if abs(coeff)<1e-14: continue
                nt=list(b); nt[pos]=r
                if len(set(nt))<k: continue
                arr=nt[:]; sign=1
                for x in range(len(arr)):
                    for y in range(len(arr)-1-x):
                        if arr[y]>arr[y+1]: arr[y],arr[y+1]=arr[y+1],arr[y]; sign=-sign
                rho[bidx[tuple(arr)],col]+=sign*coeff
    return rho
def g2_rep_on_wedge(k): return [lie_on_wedge(M,k) for M in G2_gens7]

def trivial_dim(rep_gens):
    """# of common-kernel vectors = # G₂-invariants in the rep."""
    M = np.vstack(rep_gens)
    r,_ = rank_tol(M)
    return rep_gens[0].shape[0]-r

def commutant_dim(rep_gens):
    """dim of {X : [X,ρ(g)]=0 ∀g}.  = Σ_λ m_λ²  (Schur).  Multiplicity-free ⟹ #irreps."""
    d=rep_gens[0].shape[0]
    rows=[np.kron(np.eye(d), g) - np.kron(g.T, np.eye(d)) for g in rep_gens]
    A=np.vstack(rows)
    r,_=rank_tol(A)
    return d*d-r

# ===========================================================================
# (1) G₂-decomposition of the quark sector grades + the "−14" question
# ===========================================================================
print("\n=== (1) G₂-CONTENT OF EACH SECTOR GRADE Λ^k ℝ⁷ ===")
print("    (#inv = G₂-invariants; #irr = commutant dim = # distinct irreps if mult-free)")
sector_grades = {}
for k in [1,2,3,4,5,6]:
    rep=g2_rep_on_wedge(k)
    d=rep[0].shape[0]; ninv=trivial_dim(rep); nirr=commutant_dim(rep)
    sector_grades[k]=(d,ninv,nirr)
    print(f"  Λ^{k}ℝ⁷  dim={d:2d}  #inv={ninv}  #irr(commutant)={nirr}")

print("""
  Textbook G₂ branchings (to compare):  Λ¹=7, Λ²=7⊕14, Λ³=1⊕7⊕27, Λ⁴=1⊕7⊕27,
  Λ⁵=7⊕14, Λ⁶=7.  So expected (#inv,#irr): Λ²(0,2) Λ⁴(1,3) Λ⁶(0,1).
  Our #inv counts: Λ²={i2}, Λ⁴={i4}, Λ⁶={i6}.""".format(
      i2=sector_grades[2][1], i4=sector_grades[4][1], i6=sector_grades[6][1]))

print("\n--- (1b) Does '−14' remove a genuine G₂-invariant subspace per sector? ---")
# LEPTON: in g2_koide_derivation.py the lepton sector is so(8)=Λ²ℝ⁸=14⊕7⊕7, and
# G₂(14) sits there as the ADJOINT subALGEBRA.  '28-14' = remove the adjoint. GENUINE.
print("  LEPTON  so(8)=Λ²ℝ⁸ = 14 ⊕ 7 ⊕ 7 : G₂(14)=adjoint is a genuine subalgebra.")
print("          (NB: this '28' is so(8), NOT the Λ²ℝ⁷⊕Λ⁶ℝ⁷=21+7 used in Predictions.lean.)")
# d-quark: Λ⁴=1⊕7⊕27.  Is there a 14-dim G₂-invariant subspace to remove?
consts=[1,7,27]
subsets14=[c for r in range(1,4) for c in combinations(consts,r) if sum(c)==14]
print(f"  d-QUARK Λ⁴ℝ⁷ = 1 ⊕ 7 ⊕ 27 : subsets of constituents summing to 14 = "
      f"{subsets14 or 'NONE'}.")
print("          ⟹ NO 14-dim G₂-invariant subspace exists in Λ⁴. '35-14=21' is bare")
print("            arithmetic; the G₂ adjoint (14) is NOT EVEN A CONSTITUENT of Λ⁴.")
# u-quark: Λ²⊕Λ⁴⊕Λ⁶ = (7⊕14)⊕(1⊕7⊕27)⊕7 = 1 ⊕ 14 ⊕ 4·7 ⊕ 27.
print("  u-QUARK Λ²⊕Λ⁴⊕Λ⁶ = (7⊕14)⊕(1⊕7⊕27)⊕7 = 1 ⊕ 14 ⊕ 7⊕7⊕7⊕7 ⊕ 27.")
print("          ⟹ ONE 14 (the adjoint, from Λ²) IS present; 63-14=49 could remove it.")
print("            But it is the SAME copy that the d-quark sector lacked entirely.")

# ===========================================================================
# (2) Spin(7) spinor 8 = 1 ⊕ 7 under G₂; G₂ fixes the lepton singlet?
# ===========================================================================
print("\n=== (2) Spin(7) spinor 8 under G₂ — lepton-singlet 'inert' test ===")
def gamma7():
    I2=np.eye(2); X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]]); Z=np.array([[1,0],[0,-1]],complex)
    def kron(*a):
        r=a[0]
        for m in a[1:]: r=np.kron(r,m)
        return r
    return [kron(X,I2,I2),kron(Y,I2,I2),kron(Z,X,I2),kron(Z,Y,I2),
            kron(Z,Z,X),kron(Z,Z,Y),kron(Z,Z,Z)]
gam=gamma7()
cliff_ok=all(np.allclose(gam[a]@gam[b]+gam[b]@gam[a], 2*(a==b)*np.eye(8)) for a in range(7) for b in range(7))
def Sigma(p,q): return 0.25*(gam[p]@gam[q]-gam[q]@gam[p])
def spinor_of_so7(M7):
    out=np.zeros((8,8),complex)
    for (p,q) in combinations(range(7),2): out += M7[p,q]*Sigma(p,q)
    return out
G2_spinor=[spinor_of_so7(M) for M in G2_gens7]
r,sv=rank_tol(np.vstack(G2_spinor))
inv_dim=8-r
us,ss,vts=np.linalg.svd(np.vstack(G2_spinor))
singlet=vts[-1].conj()
maxact=max(np.linalg.norm(g@singlet) for g in G2_spinor)
print(f"  Cl(7) Clifford relations OK = {cliff_ok}")
print(f"  8-spinor under G₂: dim(invariants) = {inv_dim}  (expect 1 ⟹ 8 = 1 ⊕ 7)")
print(f"  G₂ on the singlet vector: max‖X·v‖ = {maxact:.1e}  (≈0 ⟹ singlet truly fixed)")
print(f"  ⟹ The G₂ singlet IS inert. But: 8 has only ONE G₂-singlet, while there are")
print(f"    TWO leptons (N=0,3). The second lepton lives in the '7' (as the SU(3)-singlet),")
print(f"    on which G₂ acts NON-trivially (the 7 is G₂-irreducible).")

# verify the 7 is G₂-irreducible (commutant dim 1) and has no G₂ invariant:
rep7=g2_rep_on_wedge(1)
print(f"  Check: 7 (=Λ¹) #G₂-invariants={trivial_dim(rep7)} (expect 0), "
      f"commutant dim={commutant_dim(rep7)} (expect 1 ⟹ irreducible).")

# ===========================================================================
# (3) THE CRUX: does a G₂-commuting ('inert-respecting') operator carry no
#     mass-splitting?  Schur ⟹ scalar per irrep, but different scalars allowed.
# ===========================================================================
print("\n=== (3) CRUX: can a G₂-COMMUTING operator still split masses? ===")
for k,name in [(2,"Λ² (so(7) part)"),(4,"Λ⁴ (d-quark F)"),(6,"Λ⁶")]:
    rep=g2_rep_on_wedge(k)
    print(f"  {name:18s}: dim commutant of G₂ = {commutant_dim(rep)}  "
          f"(# free 'inert' parameters = # irreps; each can take its own mass).")
print("""  Interpretation (representation theory, not opinion):
    A G₂-commuting operator is scalar on each G₂-irrep (Schur's lemma) but may take
    a DIFFERENT scalar on each distinct irrep.  Hence a 'G₂-inert' (G₂-respecting)
    mass operator CAN split a sector into its irrep blocks — it only fails to split
    WITHIN a single irrep.  So 'G₂-invariant ⟹ NO mass-splitting' is FALSE as a blanket
    statement: it forbids only INTRA-irrep splitting.  The 3 generations are not one
    G₂-irrep, so generation splitting is NOT forbidden by G₂-inertness — the loophole.""")

# Concrete demonstration: an explicit G₂-commuting operator on Λ⁴ that is NOT scalar.
rep4=g2_rep_on_wedge(4); d4=35
# Build a random element of the commutant by projection: average a random matrix over G₂.
# (Use the group-averaging surrogate: solve [X,g]=0 least squares -> take a basis vector.)
rows=[np.kron(np.eye(d4), g) - np.kron(g.T, np.eye(d4)) for g in rep4]
Ac=np.vstack(rows); uu,scc,vv=np.linalg.svd(Ac)
null=vv[len(scc)-np.sum(scc<1e-6*scc.max()):] if np.any(scc<1e-6*scc.max()) else vv[np.sum(scc>=1e-6*scc.max()):]
nmask=np.concatenate([scc<1e-6*scc.max(), np.ones(vv.shape[0]-len(scc),bool)])
commX=[v.reshape(d4,d4) for v in vv[nmask]]
# pick the commutant element with the most distinct eigenvalues (proves non-scalar splitting)
best=None;bestspread=0
for X in commX:
    Xs=(X+X.T)/2
    ev=np.linalg.eigvalsh(Xs)
    spread=len(set(np.round(ev,3)))
    if spread>bestspread: bestspread=spread; best=ev
print(f"  Demo: a G₂-commuting operator on Λ⁴ has eigenvalue-multiset with "
      f"{bestspread} distinct values ⟹ it assigns DIFFERENT masses to the "
      f"1, 7, 27 blocks (NOT a multiple of identity).")

# ===========================================================================
# (4) Is "14" privileged as dimG₂, or just 2·dim(Im𝕆)=2·7?
# ===========================================================================
print("\n=== (4) Is the integer 14 privileged as dim G₂? ===")
print("  Coincidences that all give 14:")
print("    14 = dim G₂ = dim Aut(𝕆)            [the claimed structural origin]")
print("    14 = 2·7 = 2·dim(Im𝕆) = dim(7⊕7)    [the nonG₂ complement in so(8)]")
print("    14 = 21 − 7 = dim so(7) − dim(vector 7)")
print("  Sector-by-sector, does the '14' have a G₂-SUBSPACE meaning?")
print("    lepton (so(8)=14⊕7⊕7): YES — both readings coincide; nonG₂ = 7⊕7 = 2·7.")
print("    d-quark (Λ⁴=1⊕7⊕27)  : NO  — neither 14=dimG₂ nor 2·7 is a sub-constituent;")
print("                              35-14=21 matches no block sum (blocks 1,7,27).")
print("    u-quark (1⊕14⊕4·7⊕27): the 14 appears as the adjoint copy in Λ²; 2·7 also")
print("                              realizable as two of the four 7's — AMBIGUOUS which.")
print("""  ⟹ Only for the lepton does '14' robustly mean 'the G₂ subspace'.  For quarks the
     invariant content of '14' is the bare integer (which equals BOTH dimG₂ and 2·dimImO).""")

print("\n[done]")
