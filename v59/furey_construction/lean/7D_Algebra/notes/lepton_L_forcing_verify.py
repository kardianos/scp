#!/usr/bin/env python3
"""Round 5: the lepton 2x2 block from L vs F, full picture, and the honest verdict.

We separate two notions of 'lepton mass':
  (D) DIAGONAL (Majorana-type) entries on indices 0 and 7 individually.
  (O) OFF-DIAGONAL (Dirac-type) {0,7} coupling joining N=0 <-> N=3.

For the Koide/Brannen kernel the relevant object is the 2x2 (in fact the full
flavor) mass matrix on the lepton SECTOR, i.e. the 2x2 block on {0,7}. Let's get
the FULL 2x2 blocks each grade supplies and their ranks, and decide.
"""
import itertools
import numpy as np
octMultTable = [
  [(1,0),(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7)],
  [(1,1),(-1,0),(1,3),(-1,2),(1,5),(-1,4),(-1,7),(1,6)],
  [(1,2),(-1,3),(-1,0),(1,1),(1,6),(1,7),(-1,4),(-1,5)],
  [(1,3),(1,2),(-1,1),(-1,0),(1,7),(-1,6),(1,5),(-1,4)],
  [(1,4),(-1,5),(-1,6),(-1,7),(-1,0),(1,1),(1,2),(1,3)],
  [(1,5),(1,4),(-1,7),(1,6),(-1,1),(-1,0),(-1,3),(1,2)],
  [(1,6),(1,7),(1,4),(-1,5),(-1,2),(1,3),(-1,0),(-1,1)],
  [(1,7),(-1,6),(1,5),(1,4),(-1,3),(-1,2),(1,1),(-1,0)],
]
def gamma(k):
    M=[[0]*8 for _ in range(8)]
    for p in range(8):
        for m in range(8):
            s,t=octMultTable[k+1][m]
            if t==p: M[p][m]=s
    return M
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(8)) for j in range(8)] for i in range(8)]
def diag(M): return [M[i][i] for i in range(8)]
ID8=[[1 if i==j else 0 for j in range(8)] for i in range(8)]
def prod(ks):
    M=ID8
    for k in ks: M=mm(M,gamma(k))
    return M
def blk(M,inds): return [[M[a][b] for b in inds] for a in inds]

L_biv=[mm(gamma(i),gamma(j)) for i in range(7) for j in range(i+1,7)]
hexads=[[k for k in range(7) if k!=d] for d in range(7)]
L_six=[prod(h) for h in hexads]
L_all=L_biv+L_six
F_all=[prod(list(q)) for q in itertools.combinations(range(7),4)]
lep=[0,7]

def block_span_rank(ops,inds):
    rows=[]
    for M in ops:
        B=blk(M,inds)
        rows.append([B[i][j] for i in range(len(inds)) for j in range(len(inds))])
    R=np.array(rows)
    return np.linalg.matrix_rank(R)

print("Full lepton 2x2 {0,7} block span rank:")
print("  from L-grade:", block_span_rank(L_all,lep), "(out of 4)")
print("  from F-grade:", block_span_rank(F_all,lep), "(out of 4)")

# enumerate distinct 2x2 lepton blocks from L
def blocks(ops,inds):
    s=set()
    for M in ops:
        B=blk(M,inds)
        s.add(tuple(tuple(r) for r in B))
    return s
print("\nDistinct lepton 2x2 blocks from L:", sorted(blocks(L_all,lep)))
print("Distinct lepton 2x2 blocks from F:", sorted(blocks(F_all,lep)))

# So which structural component can the L-grade lepton block be? It is purely
# antisymmetric (epsilon-like): [[0,1],[-1,0]]. That is the SO(2)~U(1) rotation
# joining the two singlets. A Brannen/Koide MASS is a SYMMETRIC (Hermitian after
# i-factor) 2x2. Let's classify symmetric vs antisymmetric content per grade.
def classify(ops,inds,name):
    sym=set(); antisym=set()
    for M in ops:
        B=np.array(blk(M,inds))
        if (B==B.T).all() and B.any(): sym.add(tuple(B.flatten()))
        if (B==-B.T).all() and B.any(): antisym.add(tuple(B.flatten()))
    print(f"{name}: symmetric lepton blocks {sorted(sym)}; antisymmetric {sorted(antisym)}")
classify(L_all,lep,"L-grade")
classify(F_all,lep,"F-grade")

# Interpretation check: in Cl(7)_even, the lepton singlets are N=0 (e_R) and
# N=3 (nu_R) -> they are CHARGE CONJUGATES / opposite-parity Weyl components.
# A real mass connecting them is the Dirac mass. L gives the antisymmetric
# (=i*symmetric after complexification, i.e. a genuine Hermitian Dirac mass);
# F gives diagonal (Majorana) + symmetric off-diag. Both nonzero.

# FINAL honest tally
print("\n================ VERDICT INPUTS ================")
print("L lepton block: rank", block_span_rank(L_all,lep), "; content = antisymmetric only (epsilon).")
print("F lepton block: rank", block_span_rank(F_all,lep), "; content = diagonal + symmetric off-diag.")
print("Both grades CAN furnish a lepton mass term. Neither is excluded by")
print("availability alone. => lepton=L is NOT forced by the 8x8 mass-availability")
print("computation. The G2-singlet-in-F argument does NOT force lepton=F either;")
print("F's only color-neutral diagonal is the GLOBAL chirality operator, shared")
print("with quarks, not a lepton-isolating channel.")
