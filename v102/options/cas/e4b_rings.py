#!/usr/bin/env python3
"""E4(b), as run: exact spin-S Heisenberg rings (J = 1), S = 1, 3/2, 2,
n = 4, 6, 8; gap = lowest S^z = 1 level minus the S^z = 0 ground state.
A gapless chain has gap * n roughly constant; a gapped one has gap -> Delta_S.
(Replaces the open-chain DMRG of e4_ceiling.py part b, which was too slow
for the question it answers.)"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G
t0 = time.time()
print("# E4(b) exact spin-S Heisenberg rings, J = 1")
for S in (1.0, 1.5, 2.0):
    d = int(round(2*S + 1)); Sz, Sp = G.spin_ops(S)
    for n in (4, 6, 8):
        at = lambda op, i: sp.kron(sp.kron(sp.identity(d**i), sp.csr_matrix(op)), sp.identity(d**(n - i - 1)), format='csr')
        Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
        H = sp.csr_matrix((d**n, d**n))
        for i in range(n):
            j = (i + 1) % n
            H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
        tot = sum(Z).diagonal()
        E = {}
        for m in (0, 1):
            keep = np.where(np.abs(tot - m) < 1e-9)[0]
            Hm = H[keep][:, keep]
            E[m] = (np.linalg.eigvalsh(Hm.toarray())[0] if Hm.shape[0] < 3000
                    else eigsh(Hm, k=1, which='SA', tol=1e-10)[0][0])
        gap = E[1] - E[0]
        print(f"  S={S} n={n}: E0={E[0]:+.6f}  gap={gap:.6f}  gap*n={gap*n:.4f}  [{time.time()-t0:.0f}s]", flush=True)
print("  literature (open/periodic chains, n -> infinity): Delta_1 = 0.4105 J; S = 3/2 gapless; "
      "Delta_2 = 0.0876-0.0892 J (xi ~ 49); Delta_3 ~ 0.010 J (xi ~ 640)")
print(f"done [{time.time()-t0:.0f}s]")
