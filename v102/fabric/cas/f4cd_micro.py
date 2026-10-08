#!/usr/bin/env python3
"""F4(c, d) — the microscopic S-law ring with one hole (real link quanta
allowed; the check on reading V), and its translation under transport t.
n = 6, caps 2 (C-lock), T2 rule (the hole's only way to move: a full
neighbour's quantum hops in), sum c <= 3, S^z = 1/2.  The hole is seeded
at every position; at t = 0 the six positions are disconnected (exact
degeneracy).  Reported per t: the lowest 12 levels with total spin, the
band width W of the lowest six (the hole's positional states), and the
ground spin."""
import sys, time, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../options/cas')
import graph as G

t0 = time.time()
n = 6; gr = 0.1; g = gr*G.DELTA; J = 0.009231
def total_S_half(v, states):
    acc = {}
    for k, (cells, cs) in enumerate(states):
        if abs(v[k]) < 1e-13: continue
        for i, (a, b) in enumerate(cells):
            if b == 0: continue
            nc = list(cells); nc[i] = (a + 1, b - 1)
            key = (tuple(nc), cs)
            acc[key] = acc.get(key, 0.0) + v[k]*np.sqrt(b*(a + 1))
    s2 = sum(abs(x)**2 for x in acc.values())
    return -0.5 + np.sqrt(1.0 + s2)             # S(S+1) = <S-S+> + Sz(Sz+1), Sz = 1/2

seeds = []
for h in range(n):
    for combo in itertools.product([(2, 0), (1, 1), (0, 2)], repeat=n - 1):
        cells = list(combo); cells.insert(h, (1, 0))
        seeds.append((tuple(cells), (0,)*n))
E2 = (n - 1)*G.EF[2] + G.EF[1]
print(f"# F4(c,d) microscopic one-hole ring n = {n}, g/Delta = {gr}, J = {J}; T2 transport t; sum c <= 3; S^z = 1/2")
for (lab, t) in (("0", 0.0), ("0.1J", 0.1*J), ("0.3J", 0.3*J), ("J", J), ("3J", 3*J), ("10J", 10*J), ("0.1Delta", 0.1*G.DELTA), ("0.17Delta", 0.17*G.DELTA)):
    gm = G.Graph(n, G.ring_links(n), g, t=[t]*n, rule='T2', ceil='fixed', caps=[2]*n, cmax=3)
    H = gm.build(seeds, keep=lambda s: G.sz_of(s) == 1)
    H = H - E2*sp.identity(H.shape[0])
    k = 14
    w, v = eigsh(H, k=k, which='SA', tol=1e-11); o = np.argsort(w); w = w[o]; v = v[:, o]
    S = [total_S_half(v[:, q], gm.states) for q in range(k)]
    N, C = G.obs_arrays(gm.states, n, n)
    sc = float(np.abs(v[:, 0])**2@C.sum(1))
    W = w[5] - w[0]
    print(f"  t={lab:9s} ({t:.5f}, t/J={t/J:.2f}): dim={H.shape[0]:7d}  ground S={S[0]:.2f}  <sum c>/n={sc/n:.4f}  "
          f"W(lowest 6)={W/J:.5f} J  levels (E-E0)/J: " + " ".join(f"{(e-w[0])/J:.3f}(S={s:.1f})" for e, s in zip(w[:12], S[:12])) +
          f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
