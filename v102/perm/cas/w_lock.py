#!/usr/bin/env python3
"""W1 — is the widened (3+4) fiber a lock?  g/Delta = 0.1.
(a) The two-cell W cluster: its effective spectrum by (S, I).
(b) Effective rings n = 4, 6 (h2 + r3), sectors (2S^z, 2I^z): ground state, its
    (S, I), the lowest excitations by (S, I); microscopic W ring n = 4 as the check.
(c) The same for the 3+2 fiber (M = 2) through the same code, as a control."""
import sys, time, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import graphm as GM, weff
t0 = time.time()
g = 0.1*GM.DELTA
def low(H, k):
    if H.shape[0] <= 3000:
        w = np.linalg.eigvalsh(H.toarray()); return w[:k]
    return np.sort(eigsh(H, k=k, which='SA', tol=1e-11)[0])
for (name, M, P) in (("3+2 (control)", 2, GM.P_2), ("W (3+4)", 4, GM.P_W)):
    h2, r3, info = weff.effective(M, P, g)
    print(f"# {name}: cluster separations/overlaps: " + ", ".join(f"{k}:{v[0]:.3f}/{v[1]:.3f}" for k, v in info.items()) + f"  [{time.time()-t0:.0f}s]", flush=True)
    w2 = np.linalg.eigvalsh(h2[(2, 2)])
    print(f"  two-cell effective levels: " + " ".join(f"{x:+.6f}" for x in np.unique(np.round(w2, 8))))
    for n in (4, 6):
        caps = [2]*n
        T = weff.chain_terms(caps, h2, r3, ring=True)
        targets = [(0,), (2,), (4,)] if M == 2 else [(0, 0), (2, 0), (0, 2), (2, 2), (4, 0), (0, 4)]
        res = {}
        for tg in targets:
            H, *_ = weff.sector_hamiltonian(M, caps, T, tg)
            res[tg] = low(H, 6)
        E0 = res[targets[0]][0]
        print(f"  n={n} effective ring: E0 = {E0:+.6f}; lowest per sector (E-E0): " +
              "  ".join(f"{tg}: " + ",".join(f"{x-E0:.5f}" for x in res[tg][:3]) for tg in targets) + f"  [{time.time()-t0:.0f}s]", flush=True)
    # microscopic n = 4 ring check (sector 0)
    n = 4
    gm = GM.GraphM(M, n, [(i, (i + 1) % n) for i in range(n)], g, P, [2]*n)
    seeds = [(tuple(c), (0,)*n) for c in itertools.product(GM.cell_states(M, 2), repeat=n)]
    tg0 = targets[0]
    Hm = gm.build(seeds, keep=lambda s: weff.labels_of(M, tuple(map(sum, zip(*s[0])))) == tg0 if False else
                  (lambda cells: tuple(sum(x) for x in zip(*[weff.labels_of(M, c) for c in cells])) == tg0)(s[0]))
    Hm = Hm - n*GM.EF[2]*sp.identity(Hm.shape[0])
    wm = low(Hm, 6)
    print(f"  n=4 microscopic check, sector {tg0}: dim {Hm.shape[0]}, lowest: " + " ".join(f"{x:+.6f}" for x in wm[:4]) + f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
