#!/usr/bin/env python3
"""E5 — cooling a hot lock by light scattering (../PREREG.md §3 E5, K-b).

n = 4 ring (cap 2, t = 0 inside: T0 = T4), one boundary link b = (cell 0, o);
o is the end site of a semi-infinite photon lead (hop t_out = 0.9).  One
photon.  Impurity space P = P1 (photon at o, any ring state, polarization s)
+ P2 (photon converted with a ring-cell-0 quantum into the boundary
oscillator c_b at Omega).  The lead is integrated out exactly as the
self-energy Sigma(w) = (w - i sqrt(4t^2 - w^2))/2 of a semi-infinite chain,
one per channel (ring eigenstate f, polarization s), with w the photon's
kinetic energy in that channel.  F1: c_b also hops into a semi-infinite
space lead (hop t_s = s_k), the same self-energy at space kinetic energy.

S_ba = delta_ba - i sqrt(G_a G_b) Gfull_ba,  G = -2 Im Sigma (open channels);
unitarity (plus absorption under F1) is printed.

Per photon: P_{f<-i}(eps) = 1/2 sum_{s,s'} |S_{(f,s'),(i,s)}|^2 at E = e_i + eps.
Sources: cold (eps uniform in the lowest delta of the band), resonant
(eps_res +- 0.05), broadband (whole band).  Master equation over photons
from the hot Neel start (populations; photons arrive one at a time).

usage: e5_cool.py gr F       (gr = g/Delta, F = F0 | F1 | X0)
"""
import sys, time, itertools
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
gr_ = float(sys.argv[1]); F = sys.argv[2]
X1 = F != 'X0'
n = 4; g = gr_*G.DELTA; TOUT = 0.9; TS = 0.06
E2ref = n*G.EF[2]

# ---------------------------------------------------------------- ring (P1 ring factor) and depleted ring (P2)
ring = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n)
ring.build(G.ring_seeds(n, n))
rst = ring.states; rid = {s: k for k, s in enumerate(rst)}
import scipy.sparse as sp
Hr = ring.build(G.ring_seeds(n, n)).toarray() - E2ref*np.eye(len(rst))
assert ring.states == rst
szr = np.array([G.sz_of(s) for s in rst])
# eigenbasis, block-diagonal in S^z
er = np.zeros(len(rst)); Ur = np.zeros((len(rst), len(rst)))
order = []
for m in sorted(set(szr)):
    ix = np.where(szr == m)[0]
    w, v = np.linalg.eigh(Hr[np.ix_(ix, ix)])
    for k in range(len(ix)):
        order.append((w[k], m, ix, v[:, k]))
order.sort(key=lambda x: x[0])
for c, (w, m, ix, v) in enumerate(order):
    er[c] = w; Ur[ix, c] = v
szf = np.array([o[1] for o in order])
NR = len(rst)

# depleted ring: cell 0 holds one quantum less (N_0 + c_3 + c_0 = 1)
dseeds = []
for combo in itertools.product([(2, 0), (1, 1), (0, 2)], repeat=n - 1):
    for c0 in ((1, 0), (0, 1)):
        dseeds.append(((c0,) + combo, (0,)*n))
dep = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[1] + [2]*(n - 1))
Hd = dep.build(dseeds).toarray()
dst = dep.states; did = {s: k for k, s in enumerate(dst)}
# P2 energies relative to (ring full + photon at o): cells + Omega (c_b = 1) - E2ref - E_f(1)
H22 = Hd + (G.OMEGA - E2ref - G.EF[1])*np.eye(len(dst))
lam, Ud = np.linalg.eigh(H22)
ND = len(dst)

# ---------------------------------------------------------------- V: P1 (ring state r, pol s) -> P2 (depleted state d)
V = np.zeros((ND, NR*2))
if X1:
    for k, (cells, cs) in enumerate(rst):
        for s in (0, 1):                       # photon at o: (1,0) or (0,1)
            for al in (0, 1):
                be = s
                m = G.WEDGE[al, be]
                if m == 0 or cells[0][al] == 0: continue
                nc = list(cells); x = list(nc[0]); x[al] -= 1; nc[0] = tuple(x)
                d = did.get((tuple(nc), cs))
                if d is None: raise RuntimeError("depleted state missing")
                V[d, 2*k + s] += g*m*np.sqrt(cells[0][al]*1)
# to eigenbases: P1 columns (f, s) and P2 rows (eigen m)
U1 = np.kron(Ur, np.eye(2))                     # (r,s) -> (f,s)
Vt = Ud.T@V@U1                                  # (m) x (f,s)
e1 = np.repeat(er, 2)                           # channel energies (ring part)

def sigma(w, t):
    w = np.asarray(w, dtype=float)
    inside = np.abs(w) < 2*t
    out = np.empty(w.shape, complex)
    out[inside] = (w[inside] - 1j*np.sqrt(4*t*t - w[inside]**2))/2
    wo = w[~inside]
    out[~inside] = (wo - np.sign(wo)*np.sqrt(wo**2 - 4*t*t))/2
    return out

# conserved total 2S^z: channel (f, s) carries szf[f] + (+1 if s == 0 else -1); a P2 state its own
Mch = np.repeat(szf, 2) + np.tile([1, -1], NR)
szd_basis = np.array([G.sz_of(s) for s in dst])
lam = np.zeros(ND); Ud = np.zeros((ND, ND)); szd = np.zeros(ND, int)
c = 0
for m in sorted(set(szd_basis)):
    ix = np.where(szd_basis == m)[0]
    w, v = np.linalg.eigh(H22[np.ix_(ix, ix)])
    for k in range(len(ix)):
        lam[c] = w[k]; Ud[ix, c] = v[:, k]; szd[c] = m; c += 1
Vt = Ud.T@V@U1
blocks = {}
for M in sorted(set(Mch)):
    chs = np.where(Mch == M)[0]; ps = np.where(szd == M)[0]
    blocks[M] = (chs, ps, Vt[np.ix_(ps, chs)])
leak = max(abs(Vt[np.ix_(np.where(szd == M1)[0], np.where(Mch == M2)[0])]).max(initial=0)
           for M1 in set(szd) for M2 in set(Mch) if M1 != M2)
print(f"  S^z blocks: {len(blocks)}; largest channel block {max(len(b[0]) for b in blocks.values())}; cross-sector V {leak:.1e}")

def smatrix(E, a):
    """S column for incoming channel a (if open): (open channel ids, |S|^2 column)."""
    chs, ps, Vb = blocks[Mch[a]]
    w = E - e1[chs]
    ka = np.where(chs == a)[0][0]
    if not abs(w[ka]) < 2*TOUT: return None, None
    S1 = sigma(w, TOUT)
    S2 = sigma(E - lam[ps], TS) if F == 'F1' else np.zeros(len(ps), complex)
    D2 = 1.0/(E - lam[ps] - S2 + 1e-13j)
    A = np.diag(E - e1[chs] - S1) - (Vb.T*D2)@Vb
    Gam = -2*S1.imag
    r = np.zeros(len(chs), complex); r[ka] = 1
    Gc = np.linalg.solve(A, r)
    op = np.where(np.abs(w) < 2*TOUT)[0]
    col = -1j*np.sqrt(Gam[op])*Gc[op]*np.sqrt(Gam[ka])
    col[np.where(op == ka)[0][0]] += 1
    return chs[op], np.abs(col)**2

eres = G.OMEGA - G.EF[2]
print(f"# E5 cooling by light scattering: n={n} ring, g/Delta={gr_}, {'X0 (no boundary conversion)' if not X1 else 'X1 + ' + F}; "
      f"t_out={TOUT}, t_s={TS}, eps_res={eres:.5f}; ring dim {NR} (all S^z), depleted dim {ND}  [{time.time()-t0:.0f}s]")
ground = 0
print(f"  ring: E0={er[0]:+.6f}  first gaps: " + " ".join(f"{x:.5f}" for x in (er[1:8] - er[0])) +
      f"  spin manifold top (level 80) = {er[80]-er[0]:.5f}, level 81 = {er[81]-er[0]:.5f}")

def trans_matrix(eps_list):
    """Average over eps of P_{f<-i}; returns M (NR x NR) and absorption per i; unitarity worst."""
    M = np.zeros((NR, NR)); absb = np.zeros(NR); worst = 0.0
    for eps in eps_list:
        E_by_i = er + eps
        for i in range(NR):
            E = E_by_i[i]
            for a in (2*i, 2*i + 1):
                op, col = smatrix(E, a)
                if op is None: continue
                tot = col.sum(); worst = max(worst, abs(tot - 1) if F != 'F1' else max(0, tot - 1))
                absb[i] += (1 - tot)/2/len(eps_list)
                np.add.at(M[:, i], op//2, col/2/len(eps_list))
    return M, absb, worst

def neel_pop():
    v = np.zeros(NR)
    v[rid[(((2, 0), (0, 2), (2, 0), (0, 2)), (0,)*n)]] = 1
    return (Ur.T@v)**2

sources = []
for delta in (0.005, 0.01, 0.02, 0.05):
    sources.append((f"cold delta={delta}", np.linspace(-2*TOUT + delta/20, -2*TOUT + delta, 10)))
sources.append(("resonant eps_res+-0.05", np.linspace(eres - 0.05, eres + 0.05, 11)))
sources.append(("broadband", np.linspace(-2*TOUT + 0.01, 2*TOUT - 0.01, 37)))
p0 = neel_pop()
print(f"  hot start (Neel): ground pop {p0[0]:.4f}, spin-manifold weight {p0[:81].sum():.4f}, mean energy above E0 {(p0@er)-er[0]:.5f}")
for (name, eps) in sources:
    M, absb, worst = trans_matrix(eps)
    # unitarity is exact in the physics; numerically it slips near the infinitely
    # narrow in-band resonances of F0/X0.  Where there is no absorption channel,
    # renormalize each column to 1; where there is (F1), cap column + absorption at 1.
    cs = M.sum(0) + absb
    corr = np.abs(cs - 1).max() if F != 'F1' else max(0.0, (cs - 1).max())
    if F != 'F1':
        absb = np.zeros(NR); M = M/M.sum(0)[None, :]
    else:
        over = cs > 1; M[:, over] /= cs[over][None, :]; absb[over] /= cs[over]
    # photon-number dependence by repeated squaring of the augmented map
    # [[M, 0], [absb, 1]] (last state = destroyed); stable for a substochastic M.
    Aug = np.zeros((NR + 1, NR + 1)); Aug[:NR, :NR] = M; Aug[NR, :NR] = absb; Aug[NR, NR] = 1.0
    q0 = np.concatenate([p0, [0.0]])
    Pk = Aug.copy(); hit = None; hist = []
    for k in range(0, 44):                      # N = 2^k photons
        q = Pk@q0; Nph = 2**k
        if hit is None and q[0] >= 0.9: hit = (Nph, q[NR])
        if k in (0, 3, 7, 10, 13, 17, 20, 23, 27, 30, 33, 40, 43):
            hist.append((Nph, q[0], q[NR], (q[:NR]@er)/max(q[:NR].sum(), 1e-300) - er[0]))
        Pk = Pk@Pk
    lw, Vr = np.linalg.eig(M)
    j = np.argmax(np.abs(lw)); st_ = np.abs(np.real(Vr[:, j])); st_ /= st_.sum()
    print(f"  source {name:24s}: unitarity worst {worst:.1e} (column correction {corr:.1e}); per-photon destruction from ground {absb[0]:.3e}, "
          f"max over levels {absb.max():.3e}; per-photon ground->other {1-M[0,0]-absb[0]:.3e}; "
          f"[{time.time()-t0:.0f}s]")
    print("     " + "  ".join(f"N={k:.0e}: p0={a:.4f} dest={b:.2e} E-E0={c:.5f}" for k, a, b, c in hist))
    print(f"     reach p0 >= 0.9: {'N_ph <= %.3g (first power of 2), cumulative destruction %.3e' % hit if hit else 'never within 2^43 = 8.8e12 photons'};  "
          f"stationary of the map: p0 = {st_[0]:.4f}, mean E-E0 = {st_@er - er[0]:.5f}", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
