#!/usr/bin/env python3
"""E8 — the admissible window of g (../PREREG.md §3 E8).

usage: e8_g.py static GR      cluster coefficients + n = 6 exact ring + hold
       e8_g.py margin GR      t_c(g) on the n = 6 ring (sum c <= 2, general basis)
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../lock/cas')
import graph as G
from clusters import effective

t0 = time.time()
part, gr_ = sys.argv[1], float(sys.argv[2])
g = gr_*G.DELTA
n = 6

if part == 'static':
    h2, r3, co, info = effective(g)
    sepmin = min(v[0] for v in info.values()); evmin = min(v[1] for v in info.values())
    r3max = max(np.abs(x).max() for x in r3.values())
    print(f"# E8 static g/Delta={gr_} (g={g:.5f})")
    print(f"  clusters: J={co['J']:.6f} = {co['J']/(2*g*g/G.DELTA):.4f} x 2g^2/Delta;  K/J={co['K']/co['J']:+.4f};  "
          f"max|r3|/J={r3max/co['J']:.4f};  J_end/J={co['J_end']/co['J']:.3f};  min separation to converted levels={sepmin:.4f};  "
          f"min overlap eigenvalue={evmin:.4f}  (Haldane phase of the bilinear-biquadratic chain: -1 < K/J < 1)")
    gr = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n)
    H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0).toarray() - n*G.EF[2]*np.eye(len(gr.states))
    w, v = np.linalg.eigh(H)
    hx = G.heis_spectrum(n)[:13]; hx = hx[1:] - hx[0]
    J, miss = G.fit_J(w[1:13] - w[0], hx)
    S = G.total_S(v[:, 0], gr.states)
    N, C = G.obs_arrays(gr.states, n, n)
    sc = float(np.abs(v[:, 0])**2@C.sum(1))
    print(f"  n=6 exact ring: S={S:.3f}  <sum c>/n={sc/n:.4f}  J_fit={J:.6f}  gap/J={(w[1]-w[0])/J:.4f}  L1 miss={miss.max()*100:.2f}%  "
          f"K1(static) {'PASS' if (abs(S) < 1e-3 and sc/n <= 0.1 and miss.max() < 0.05) else 'FAIL'}  [{time.time()-t0:.0f}s]", flush=True)
    # hold: sudden start, all S^z, T = 100/J
    gra = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n)
    Ha = gra.build(G.ring_seeds(n, n)).toarray() - n*G.EF[2]*np.eye(len(gra.states))
    wa, Va = np.linalg.eigh(Ha)
    Na, Ca = G.obs_arrays(gra.states, n, n)
    T = 100/J; ts = np.linspace(0, T, 401)
    rng = np.random.default_rng(1001)
    starts = {"neel": [(2, 0) if i % 2 == 0 else (0, 2) for i in range(n)]}
    for name, cells in starts.items():
        v0 = np.zeros(len(gra.states)); v0[gra.idx[(tuple(cells), (0,)*n)]] = 1
        a0 = Va.T@v0
        traj = (Va@(a0[:, None]*np.exp(-1j*np.outer(wa, ts)))).T
        P = np.abs(traj)**2
        cl = P@Ca; ff = (P@Na).sum(1)/(2*n)
        print(f"  hold ({name}, T=100/J={T:.0f}): max<c_l>={cl.max():.4f}  field fraction min={ff.min():.4f}  2nd-half mean={ff[200:].mean():.4f}  "
              f"K1(hold >= 0.85) {'PASS' if ff[200:].mean() >= 0.85 else 'FAIL'}  [{time.time()-t0:.0f}s]", flush=True)

if part == 'margin':
    st = None
    print(f"# E8 margin g/Delta={gr_}: n=6, sum c <= 2, general basis, S^z = 0")
    gr = G.Graph(n, G.ring_links(n), g, t=[0.0]*n, rule='T2', ceil='fixed', caps=[2]*n, cmax=2)
    prev_ok = None
    for tr in (0.05, 0.10, 0.15, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50):
        gr.t = [tr*G.DELTA]*n
        H = gr.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0)
        H = H - n*G.EF[2]*sp.identity(H.shape[0])
        w, v = eigsh(H, k=4, which='SA', tol=1e-11); o = np.argsort(w); v = v[:, o]
        S = G.total_S(v[:, 0], gr.states)
        ok = abs(S) < 1e-3
        print(f"  t/Delta={tr:.2f}: S_ground={S:.3f} {'lock' if ok else 'NOT the lock'}  [{time.time()-t0:.0f}s]", flush=True)
        if not ok: break
print(f"done [{time.time()-t0:.0f}s]")
