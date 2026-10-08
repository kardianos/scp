#!/usr/bin/env python3
"""E1 — transport inside the ring (../PREREG.md §3 E1).
n = 6 ring, caps 2 (C-lock), general basis, sum c <= 3, S^z = 0, g/Delta = 0.1.
Rules: T0 (t = 0), T1 (T2 rule at t = 0.05 Delta, 0.10 Delta), T2, T3, T4 (t = 0.9)."""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
n = 6; g = 0.1*G.DELTA
hx = G.heis_spectrum(n)[:13]; hx = hx[1:] - hx[0]
print(f"# E1 transport inside the ring: n = {n}, caps 2, sum c <= 3, S^z = 0, g/Delta = 0.1, Delta = {G.DELTA:.5f}")
print(f"  K1 bars: S_ground = 0, <sum c>/n <= 0.1, L1 miss < 5%")
cases = [("T0", 'T2', 0.0), ("T1 0.05D", 'T2', 0.05*G.DELTA), ("T1 0.10D", 'T2', 0.10*G.DELTA),
         ("T2", 'T2', 0.9), ("T3", 'T3', 0.9), ("T4", 'T4', 0.9)]
for (name, rule, t) in cases:
    gr = G.Graph(n, G.ring_links(n), g, t=[t]*n, rule=rule, ceil='fixed', caps=[2]*n, cmax=3)
    H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
    H = H - n*G.EF[2]*sp.identity(H.shape[0])
    k = 14
    if H.shape[0] > 3000:
        w, v = eigsh(H, k=k, which='SA', tol=1e-12)
    else:
        w, v = np.linalg.eigh(H.toarray()); w, v = w[:k], v[:, :k]
    o = np.argsort(w); w = w[o]; v = v[:, o]
    S = G.total_S(v[:, 0], gr.states)
    N, C = G.obs_arrays(gr.states, n, n)
    p0 = np.abs(v[:, 0])**2
    sc = float(p0@C.sum(1)); nvar = float(p0@((N - 2)**2).sum(1))
    J, miss = G.fit_J(w[1:13] - w[0], hx)
    ok = abs(S) < 1e-3 and sc/n <= 0.1 and miss.max() < 0.05
    print(f"  {name:9s} rule={rule} t={t:.4f} (t/Delta={t/G.DELTA:.3f}): dim={H.shape[0]:7d}  E0={w[0]:+.6f}  S_ground={S:.3f}  "
          f"<sum c>={sc:.4f} (/n {sc/n:.4f})  sum<(N_i-2)^2>={nvar:.4f}  J_fit={J:.6f}  gap={w[1]-w[0]:.6f} (={(w[1]-w[0])/J:.3f} J)  "
          f"L1 max miss={miss.max()*100:.2f}%  K1 {'PASS' if ok else 'FAIL'}  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
