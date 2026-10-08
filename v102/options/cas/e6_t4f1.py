#!/usr/bin/env python3
"""E6(a) combination, small: F1 along the fiber x T4 (t = 0.9) on the n = 4 ring
(caps 2, no truncation of space, S^z = 0), against F1 x T0 (t = 0).
Space flow moves a converted quantum off its link, so R6's local law
N_i + sum c = 2 can fail at a cell and T4's predicate can open."""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G
t0 = time.time(); n = 4; g = 0.1*G.DELTA
hx = G.heis_spectrum(n)[:13]; hx = hx[1:] - hx[0]
print(f"# E6(a') F1 along the fiber x transport rule, n = {n}, caps 2, S^z = 0, g/Delta = 0.1")
for ts in (0.0, 0.03, 0.06, 0.12):
    for (lab, rule, tf) in (("T0 (t=0)", 'T2', 0.0), ("T4 (t=0.9)", 'T4', 0.9)):
        shop = [(l, (l + 1) % n, ts) for l in range(n)] if ts > 0 else None
        gr = G.Graph(n, G.ring_links(n), g, t=[tf]*n, rule=rule, ceil='fixed', caps=[2]*n, shop=shop)
        H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
        H = H - n*G.EF[2]*sp.identity(H.shape[0])
        w, v = (np.linalg.eigh(H.toarray()) if H.shape[0] < 4000 else eigsh(H, k=14, which='SA', tol=1e-12))
        o = np.argsort(w); w = w[o][:14]; v = v[:, o][:, :14]
        S = G.total_S(v[:, 0], gr.states)
        N, C = G.obs_arrays(gr.states, n, n)
        p0 = np.abs(v[:, 0])**2; sc = float(p0@C.sum(1))
        # weight of states where some cell is NOT at N + sum c = 2 (T4's predicate open)
        open_ = np.array([any(cells[i][0] + cells[i][1] + cs[(i - 1) % n] + cs[i] != 2 for i in range(n)) for cells, cs in gr.states])
        J, miss = G.fit_J(w[1:13] - w[0], hx)
        print(f"  t_s={ts:.2f} {lab:10s}: dim={H.shape[0]:6d}  S={S:.3f}  <sum c>/n={sc/n:.4f}  ground weight with T4 open={p0[open_].sum():.4f}  "
              f"J_fit={J:.6f}  L1 miss={miss.max()*100:.2f}%  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
