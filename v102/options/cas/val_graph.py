#!/usr/bin/env python3
"""Validate graph.py against ../../lock/cas/slaw.py (the validated S-law ring)."""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0]); sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/../../lock/cas')
import graph as G
import slaw
t0 = time.time()
g = 0.1*G.DELTA
print(f"EF match: {np.abs(G.EF[:3]-slaw.EF).max():.1e}; DELTA {G.DELTA:.6f} vs {slaw.DELTA:.6f}")
for (n, t, cmax) in ((6, 0.0, None), (6, 0.05*G.DELTA, 2), (4, 0.9, None)):
    gr = G.Graph(n, G.ring_links(n), g, t=[t]*n, rule='T2', ceil='fixed', caps=[2]*n, cmax=cmax)
    H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
    H = H - n*G.EF[2]*__import__('scipy.sparse', fromlist=['x']).identity(H.shape[0])
    st, ix = slaw.basis_ring(n, cmax=cmax, sz=0, local=(t == 0))
    Hs = slaw.hamiltonian(n, g, st, ix, t=t)
    k = 8
    w1 = np.sort(eigsh(H, k=k, which='SA', tol=1e-13)[0]) if H.shape[0] > 2000 else np.linalg.eigvalsh(H.toarray())[:k]
    w2 = np.sort(eigsh(Hs, k=k, which='SA', tol=1e-13)[0]) if Hs.shape[0] > 2000 else np.linalg.eigvalsh(Hs.toarray())[:k]
    herm = abs(H - H.T).max()
    print(f"n={n} t={t:.4f} cmax={cmax}: dim graph {H.shape[0]} (reachable) vs slaw {Hs.shape[0]};  max|dE| lowest {k} = {np.abs(w1-w2).max():.1e};  hermitian {herm:.1e}  [{time.time()-t0:.0f}s]")
# twist at theta = 0 equals the wedge; theta = 2 pi is a gauge (A -> -A)
for th in (0.0, 2*np.pi):
    M = G.twist_matrix(th)
    print(f"twist({th:.3f}) = {M.tolist()}")
