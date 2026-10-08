#!/usr/bin/env python3
"""E6 — space flow on lock links (../PREREG.md §3 E6).

(a) F1 along the fiber: each ring link's oscillator hops to the adjacent
    links' oscillators (t_s), n = 6 ring, caps 2, field t = 0 (T0), general
    basis, sum c <= 3, S^z = 0.  Combination: the same with T4 at t = 0.9
    (space flow breaks R6's local law, so T4's predicate can open).
(c) F2 rate law: first-jump rate gamma <sum c> on ring ground states
    (n = 4, 6 exact; 8 with sum c <= 3); lifetime in 1/J; and the energy
    that one jump leaves unaccounted if the departing quantum carries Omega.
usage: e6_space.py a | c
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
part = sys.argv[1]
g = 0.1*G.DELTA

def lowest(H, k=14):
    if H.shape[0] > 3000:
        w, v = eigsh(H, k=k, which='SA', tol=1e-12)
    else:
        w, v = np.linalg.eigh(H.toarray()); w, v = w[:k], v[:, :k]
    o = np.argsort(w); return w[o], v[:, o]

if part == 'a':
    n = 6
    hx = G.heis_spectrum(n)[:13]; hx = hx[1:] - hx[0]
    print(f"# E6(a) F1 along the fiber: n = {n}, caps 2, sum c <= 3, S^z = 0, g/Delta = 0.1; on-ring space band half-width 2 t_s vs Delta = {G.DELTA:.4f}")
    for (rule, tf) in (('T2', 0.0), ('T4', 0.9)):
        for ts in (0.0, 0.03, 0.06, 0.12, 0.24):
            shop = [(l, (l + 1) % n, ts) for l in range(n)] if ts > 0 else None
            gr = G.Graph(n, G.ring_links(n), g, t=[tf]*n, rule=rule, ceil='fixed', caps=[2]*n, shop=shop, cmax=3)
            H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
            H = H - n*G.EF[2]*sp.identity(H.shape[0])
            w, v = lowest(H)
            S = G.total_S(v[:, 0], gr.states)
            N, C = G.obs_arrays(gr.states, n, n)
            p0 = np.abs(v[:, 0])**2
            sc = float(p0@C.sum(1))
            J, miss = G.fit_J(w[1:13] - w[0], hx)
            ok = abs(S) < 1e-3 and sc/n <= 0.1 and miss.max() < 0.05
            lab = 'T0 (t=0)' if tf == 0 else 'T4 (t=0.9)'
            print(f"  {lab:10s} t_s={ts:.2f} (2t_s/Delta={2*ts/G.DELTA:.3f}): dim={H.shape[0]:7d}  E0={w[0]:+.6f}  S={S:.3f}  "
                  f"<sum c>/n={sc/n:.4f}  J_fit={J:.6f}  gap/J={(w[1]-w[0])/J:.3f}  L1 miss={miss.max()*100:.2f}%  "
                  f"K1 {'PASS' if ok else 'FAIL'}  [{time.time()-t0:.0f}s]", flush=True)

if part == 'c':
    print(f"# E6(c) F2 rate-law space flow: a jump c_l at rate gamma per link; ring ground states, g/Delta = 0.1")
    print(f"  readings: gamma = s_k = 0.06; gamma = s_k*s_disp = 0.018")
    for (n, cmax) in ((4, None), (6, None), (8, 3)):
        gr = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n, cmax=cmax)
        H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
        H = H - n*G.EF[2]*sp.identity(H.shape[0])
        w, v = lowest(H, 14)
        hx = G.heis_spectrum(n)[:13]; hx = hx[1:] - hx[0]
        J, miss = G.fit_J(w[1:13] - w[0], hx)
        N, C = G.obs_arrays(gr.states, n, n)
        p0 = np.abs(v[:, 0])**2; csum = float(p0@C.sum(1))
        # one jump c_0 on the ground state: the depleted state's energy (relative, field part only)
        st = gr.states; idx = gr.idx
        out = {}
        for k, (cells, cs) in enumerate(st):
            if cs[0] == 0 or abs(v[k, 0]) < 1e-14: continue
            cs2 = list(cs); cs2[0] -= 1
            key = (cells, tuple(cs2))
            out[key] = out.get(key, 0.0) + np.sqrt(cs[0])*v[k, 0]
        nrm = sum(x*x for x in out.values())
        Eafter = sum(x*x*(gr.diag(key) - n*G.EF[2]) for key, x in out.items())/nrm
        dE = Eafter - w[0]
        for gam in (0.06, 0.018):
            tau = 1/(gam*csum)
            print(f"  n={n}: <sum c>={csum:.4f}  J={J:.6f}  gamma={gam}: first-jump rate {gam*csum:.3e}, lifetime {tau:.1f} = {tau*J:.3f}/J")
        print(f"        one jump on link 0: system energy changes by {dE:+.4f} (diagonal part); the quantum leaving carries Omega = {G.OMEGA}; "
              f"unaccounted = {dE + G.OMEGA:+.4f} per jump (= {(dE + G.OMEGA)/G.DELTA:.2f} Delta)")
print(f"done [{time.time()-t0:.0f}s]")
