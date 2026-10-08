#!/usr/bin/env python3
"""L1 post-registration convergence check: n = 8 without truncation, and
sum c <= 4, 3, 2 (S^z = 0), g/Delta = 0.1.  Output: out/l1_exact8.out."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import numpy as np
from scipy.sparse.linalg import eigsh
from slaw import DELTA, basis_ring, hamiltonian, total_S, heis_ring, fit_J
g = 0.1*DELTA
print("# L1 post-registration convergence check: n=8 without truncation (S^z=0), g/Delta=0.1")
H, Z, P = heis_ring(8); szd = sum(Z).diagonal(); keep = np.where(np.abs(szd) < 1e-9)[0]
hw = np.linalg.eigvalsh(H[keep][:, keep].toarray())[:13]
for cmax in (None, 4, 3, 2):
    st, ix = basis_ring(8, cmax=cmax, sz=0); Hm = hamiltonian(8, g, st, ix)
    w, v = eigsh(Hm, k=20, which='SA', tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
    J, miss = fit_J(w[1:13]-w[0], hw[1:13]-hw[0])
    print(f"  cmax={str(cmax):>4} dim={Hm.shape[0]:6d} E0={w[0]:+.7f} J_fit={J:.6f} max miss={miss.max():.3e} S0={total_S(v[:,0],st,ix):.4f}")
