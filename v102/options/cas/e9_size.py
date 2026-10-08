#!/usr/bin/env python3
"""E9 — size (../PREREG.md §3 E9).

usage: e9_size.py a    Z1 coalescence: cube (two 4-rings + 4 rungs, exact
                       S-law, sum c <= 3) vs two 4-rings; hexagonal prism vs
                       two 6-rings at leading order (spin-1 Heisenberg, J_h2)
       e9_size.py b    Z2 twist: E0(n, theta), n = 4, 6 exact, 8 (sum c <= 3)
       e9_size.py c    Z3 energy per quantum: rings, open chains with cap-1
                       ends, free light
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../lock/cas')
import graph as G

t0 = time.time()
part = sys.argv[1]
g = 0.1*G.DELTA

def low(H, k=6):
    if H.shape[0] > 3000:
        w, v = eigsh(H, k=k, which='SA', tol=1e-12)
    else:
        w, v = np.linalg.eigh(H.toarray()); w, v = w[:k], v[:, :k]
    o = np.argsort(w); return w[o], v[:, o]

def heis_graph(n, links, sz=0):
    Sz, Sp = G.spin_ops(1)
    at = lambda op, i: sp.kron(sp.kron(sp.identity(3**i), sp.csr_matrix(op)), sp.identity(3**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((3**n, 3**n))
    for (i, j) in links:
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    keep = np.where(np.abs(sum(Z).diagonal() - sz) < 1e-9)[0]
    return H[keep][:, keep]

if part == 'a':
    print("# E9(a) Z1 coalescence (X1 between full cells; every cell keeps 2 quanta, R6 local law over its links)")
    ring4 = G.ring_links(4)
    gr = G.Graph(4, ring4, g, ceil='fixed', caps=[2]*4)
    H = gr.build(G.ring_seeds(4, 4, sz=0), keep=lambda s: G.sz_of(s) == 0); H = H - 4*G.EF[2]*sp.identity(H.shape[0])
    w4, _ = low(H)
    cube = ring4 + G.ring_links(4, 4) + [(i, i + 4) for i in range(4)]
    gc = G.Graph(8, cube, g, ceil='fixed', caps=[2]*8, cmax=3)
    seeds = G.ring_seeds(8, len(cube), sz=0)
    Hc = gc.build(seeds, keep=lambda s: G.sz_of(s) == 0)
    Hc = Hc - 8*G.EF[2]*sp.identity(Hc.shape[0])
    wc, vc = low(Hc, 8)
    Sc = G.total_S(vc[:, 0], gc.states)
    print(f"  4-ring (exact): E0 = {w4[0]:+.6f}  per cell {w4[0]/4:+.6f}  gap {w4[1]-w4[0]:.6f}")
    print(f"  cube (exact S-law, sum c <= 3, dim {Hc.shape[0]}): E0 = {wc[0]:+.6f}  per cell {wc[0]/8:+.6f}  S_ground {Sc:.3f}  "
          f"first gap (S^z=0 sector) {wc[1]-wc[0]:.6f}")
    print(f"  coalescence energy E(cube) - 2 E(4-ring) = {wc[0] - 2*w4[0]:+.6f}  ({'joining is favoured' if wc[0] < 2*w4[0] else 'apart is favoured'})  [{time.time()-t0:.0f}s]", flush=True)
    # leading order: spin-1 Heisenberg with J per converting link
    for (name, n, links) in (("4-ring", 4, ring4), ("cube", 8, cube), ("6-ring", 6, G.ring_links(6)),
                             ("hexagonal prism", 12, G.ring_links(6) + G.ring_links(6, 6) + [(i, i + 6) for i in range(6)])):
        Hh = heis_graph(n, links)
        wh = np.sort(eigsh(Hh, k=4, which='SA', tol=1e-10)[0]) if Hh.shape[0] > 2000 else np.linalg.eigvalsh(Hh.toarray())[:4]
        print(f"  leading order (spin-1 Heisenberg, units of J): {name:16s} E0 = {wh[0]:+.5f}  per cell {wh[0]/n:+.5f}  gap {wh[1]-wh[0]:.5f}  [{time.time()-t0:.0f}s]", flush=True)

if part == 'b':
    print("# E9(b) Z2 frame twist on one link: A_l -> b_i^T eps U(theta) b_j; all S^z (the twist breaks S^z)")
    for (n, cmax) in ((4, None), (6, None), (8, 3)):
        row = []
        for th in (0.0, np.pi/4, np.pi/2, np.pi):
            gr = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n, twist={0: G.twist_matrix(th)}, cmax=cmax)
            H = gr.build(G.ring_seeds(n, n)); H = H - n*G.EF[2]*sp.identity(H.shape[0])
            w, _ = low(H, 2)
            row.append(w[0])
        print(f"  n={n}: E0(theta = 0, pi/4, pi/2, pi) = " + " ".join(f"{x:+.6f}" for x in row) +
              f";  per cell: " + " ".join(f"{x/n:+.6f}" for x in row) + f";  cost of theta = pi: {row[3]-row[0]:+.6f}  [{time.time()-t0:.0f}s]", flush=True)

if part == 'c':
    print("# E9(c) Z3 energy per quantum (absolute: sum E_f(N_i) + conversion energy, over the number of quanta)")
    for n in range(4, 11):
        cmax = None if n <= 8 else 3
        gr = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n, cmax=cmax)
        H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
        w, _ = low(H, 2)
        Q = 2*n
        print(f"  ring  n={n:2d}: quanta {Q:2d}  E/quantum = {w[0]/Q:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    for n in range(3, 11):
        cmax = None if n <= 8 else 3
        caps = [1] + [2]*(n - 2) + [1]
        links = [(i, i + 1) for i in range(n - 1)]
        gr = G.Graph(n, links, g, ceil='fixed', caps=caps, cmax=cmax)
        import itertools
        seeds = []
        for combo in itertools.product(*[[(a, c - a) for a in range(c + 1)] for c in caps]):
            if sum(a - b for a, b in combo) in (0, 1, -1):
                seeds.append((tuple(combo), (0,)*len(links)))
        H = gr.build(seeds, keep=lambda s: G.sz_of(s) == 0)
        w, _ = low(H, 2)
        Q = sum(caps)
        print(f"  open  n={n:2d}: quanta {Q:2d}  E/quantum = {w[0]/Q:.6f}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"  a lone full cell, no conversion: E_f(2)/2 = {G.EF[2]/2:.6f}")
    print(f"  free light at its band bottom, t = 0.9: 1-D E_f(1) - 2t = {G.EF[1]-1.8:+.6f};  cubic lattice E_f(1) - 6t = {G.EF[1]-5.4:+.6f}")
print(f"done [{time.time()-t0:.0f}s]")
