#!/usr/bin/env python3
"""Post-registration diagnostic: bracket the transport threshold on the
n = 6 ring (S^z = 0, sum c <= 2), g/Delta = 0.1."""
import sys
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from p1_micro import micro, S2_sz0
for t in (0.12, 0.14, 0.16, 0.18):
    H, space, basis = micro(6, 0.1, t, cmax=2)
    w, v = eigsh(H, k=4, which='SA', tol=1e-12)
    o = np.argsort(w); v = v[:, o]
    print(f"  n=6 g/D=0.1 t/D={t:.2f}: S_ground={S2_sz0(v[:, 0], basis):.3f}  <space>={float(v[:, 0]**2@space):.3f}", flush=True)
