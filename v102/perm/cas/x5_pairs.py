#!/usr/bin/env python3
"""X5 (reported, not gated) — two holes: do they bind?  Microscopic S-law
(real link quanta allowed, so two adjacent holes may convert: a (1,1) pair
sits at Delta_11 = -0.024, below its space quantum).
Binding B = E(two holes) + E(no hole) - 2 E(one hole), same ring, same
truncation; B < 0 means the pair is lower than two separate holes.
(a) 3+2, n = 8, t = 0, holes at separation d = 1..4 (pinned), S^z = 0;
(b) 3+2, n = 6, transport t = J (T2), holes free, S^z = 0;
(c) W, n = 4, t = 0, d = 1, 2, sector (0,0);
(d) W, n = 4, t = J, holes free, sector (0,0)."""
import sys, time, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../options/cas')
import graph as G, graphm as GM

t0 = time.time()
g = 0.1*G.DELTA; J = 0.009231
def e0(H):
    return (np.linalg.eigvalsh(H.toarray())[0] if H.shape[0] < 3000 else eigsh(H, k=1, which='SA', tol=1e-11)[0][0])

def ring32(n, holes, t, cmax, free=False):
    seeds = []
    for combo in itertools.product([(2, 0), (1, 1), (0, 2)], repeat=n - len(holes)):
        it = iter(combo)
        if free:
            for hs in itertools.combinations(range(n), len(holes)):
                it = iter(combo)
                cells = [(1, 0) if i in hs else next(it) for i in range(n)]
                seeds.append((tuple(cells), (0,)*n))
        else:
            cells = [(1, 0) if i in holes else next(it) for i in range(n)]
            seeds.append((tuple(cells), (0,)*n))
    gr = G.Graph(n, G.ring_links(n), g, t=[t]*n, rule='T2', ceil='fixed', caps=[2]*n, cmax=cmax)
    tgt = 0 if len(holes) % 2 == 0 else 1
    H = gr.build(seeds, keep=lambda s: G.sz_of(s) == tgt)
    return e0(H), H.shape[0]

print("# X5 two holes, binding B = E2 + E0 - 2 E1 (J = J_h2 = 0.009231)")
n = 8
E0, d0 = ring32(n, [], 0.0, 3); E1, d1 = ring32(n, [0], 0.0, 3)
print(f"  (a) 3+2 n={n} t=0: E0={E0:+.6f} E1={E1:+.6f}  [{time.time()-t0:.0f}s]", flush=True)
for d in (1, 2, 3, 4):
    E2, d2 = ring32(n, [0, d], 0.0, 3)
    print(f"      separation {d}: E2={E2:+.6f}  B = {(E2+E0-2*E1)/J:+.4f} J  (dim {d2})  [{time.time()-t0:.0f}s]", flush=True)
n = 6
E0, _ = ring32(n, [], J, 3); E1, _ = ring32(n, [0], J, 3, free=True); E2, d2 = ring32(n, [0, 1], J, 3, free=True)
print(f"  (b) 3+2 n={n} t=J, holes free: B = {(E2+E0-2*E1)/J:+.4f} J  (two-hole dim {d2})  [{time.time()-t0:.0f}s]", flush=True)

def ringW(n, holes, t, free=False):
    full = GM.cell_states(4, 2); one = GM.cell_states(4, 1)
    seeds = []
    hsets = list(itertools.combinations(range(n), len(holes))) if free else [tuple(holes)]
    for hs in hsets:
        for combo in itertools.product(*[(one if i in hs else full) for i in range(n)]):
            seeds.append((tuple(combo), (0,)*n))
    gm = GM.GraphM(4, n, [(i, (i + 1) % n) for i in range(n)], g, GM.P_W, [2]*n, t=[t]*n)
    tg = (0, 0) if len(holes) % 2 == 0 else (1, 1)
    H = gm.build(seeds, keep=lambda s: GM.charges_W(s[0]) == tg)
    return e0(H), H.shape[0]

n = 4
E0, _ = ringW(n, [], 0.0); E1, _ = ringW(n, [0], 0.0)
print(f"  (c) W n={n} t=0: E0={E0:+.6f} E1={E1:+.6f}  [{time.time()-t0:.0f}s]", flush=True)
for d in (1, 2):
    E2, d2 = ringW(n, [0, d], 0.0)
    print(f"      separation {d}: B = {(E2+E0-2*E1)/J:+.4f} J  (dim {d2})  [{time.time()-t0:.0f}s]", flush=True)
E0, _ = ringW(n, [], J); E1, _ = ringW(n, [0], J, free=True); E2, d2 = ringW(n, [0, 1], J, free=True)
print(f"  (d) W n={n} t=J, holes free: B = {(E2+E0-2*E1)/J:+.4f} J  (dim {d2})  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
