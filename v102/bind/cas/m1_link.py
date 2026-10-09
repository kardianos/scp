#!/usr/bin/env python3
"""M1 (../PREREG.md) — the (1,1) link in the exact S-law against the effective
pricings of lo.py (P11 = zero, core, score, upper, lower).

Pinned holes (t = 0) on rings, the microscopic S-law (graphm, real link quanta
allowed, total space quanta <= cmax) and the leading-order model (lo) with each
pricing.  Binding B = E(two holes) + E(no hole) - 2 E(one hole), in J0.
For the two-hole sector at separation d = 1, every exact eigenstate is
characterized by its energy (as B), the space quanta on the shared link, the
total space quanta, and its overlap with each lo ground state placed in the
space-free configurations.
usage: m1_link.py M n cmax"""
import sys, time, itertools, functools
print = functools.partial(print, flush=True)
import numpy as np
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo

t0 = time.time()
M, n, cmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA
links = [(i, (i + 1) % n) for i in range(n)]
def charges(cells):
    return GM.charges_W(cells) if M == 4 else (sum(c[0] - c[1] for c in cells),)
def exact(holes, tgt, dense=False):
    full = GM.cell_states(M, 2); one = GM.cell_states(M, 1)
    seeds = [(tuple(c), (0,)*n) for c in itertools.product(*[(one if i in holes else full) for i in range(n)])]
    gm = GM.GraphM(M, n, links, g, P, [2]*n, cmax=cmax)
    H = gm.build(seeds, keep=lambda s: charges(s[0]) == tgt)
    if dense or H.shape[0] < 9000:
        w, V = np.linalg.eigh(H.toarray())
    else:
        w, V = eigsh(H, k=60, which='SA', tol=1e-11); o = np.argsort(w); w, V = w[o], V[:, o]
    return w, V, gm
def lo_ground(holes, tgt, p11):
    lo.P11 = p11
    locs = [lo.Local(M, 'hole' if i in holes else 'full') for i in range(n)]
    T = lo.chain_terms(locs, P, g, ring=True)
    H, basis, dims, strides = lo.sector(locs, T, tgt, lo.labfun_factory(M, False))
    w, V = np.linalg.eigh(H.toarray()) if H.shape[0] < 6000 else eigsh(H, k=4, which='SA', tol=1e-12)
    o = np.argsort(w); w, V = w[o], V[:, o]
    lo.P11 = 'zero'
    occ = [[locs[i].occ[(c//strides[i]) % dims[i]] for i in range(n)] for c in basis]
    return w, V, occ
def efsum(holes):
    return sum(GM.EF[1] if i in holes else GM.EF[2] for i in range(n))
q1 = (1, 1) if M == 4 else (1,)
z = (0, 0) if M == 4 else (0,)
print(f"# M1: the (1,1) link, exact S-law (cmax={cmax}) against lo pricings; {'W' if M == 4 else '3+2'} ring n={n}, pinned holes; B in J0")
w0, _, _ = exact([], z); w1, _, _ = exact([0], q1)
E0x, E1x = w0[0], w1[0]
print(f"  exact: E0 = {E0x:.9f}, E1 = {E1x:.9f}  [{time.time()-t0:.0f}s]")
pr = ['zero', 'core', 'score', 'upper', 'lower']
lo0 = {p: lo_ground([], z, p)[0][0] for p in pr}
lo1 = {p: lo_ground([0], q1, p)[0][0] for p in pr}
for d in (1, 2, 3):
    if d > n//2: break
    holes = [0, d]
    for tgt in ([z] + ([(2, 2)] if M == 4 else [(2,)])):
        wx, Vx, gm = exact(holes, tgt)
        Bx = (wx + E0x - 2*E1x)/J0       # B of each exact level (the E_f sums cancel)
        line = f"  d={d} sector {tgt}: exact ground B = {Bx[0]:+.4f}"
        if d == 1:
            cl = np.array([s[1][0] for s in gm.states]); ct = np.array([sum(s[1]) for s in gm.states])
            wl = (np.abs(Vx)**2*cl[:, None]).sum(0); wt = (np.abs(Vx)**2*ct[:, None]).sum(0)
            line += f" (shared-link space {wl[0]:.3f}, total space {wt[0]:.3f})"
        print(line + f"  [{time.time()-t0:.0f}s]")
        for p in pr:
            wl_, Vl, occ = lo_ground(holes, tgt, p)
            Bl = (wl_[0] + lo0[p] - 2*lo1[p])/J0
            s = f"      lo {p:5s}: B = {Bl:+.4f}"
            if d == 1:
                vec = np.zeros(len(gm.states))
                for k, oc in enumerate(occ):
                    key = (tuple(oc), (0,)*n)
                    if key in gm.idx: vec[gm.idx[key]] = Vl[k, 0]
                ov = np.abs(Vx.T@vec)**2
                top = np.argsort(-ov)[:4]
                s += "   overlap with exact levels: " + ", ".join(f"[B={Bx[k]:+.3f}: {ov[k]:.3f}, link {wl[k]:.2f}]" for k in top)
                s += f"  (sum over the levels kept {ov.sum():.3f})"
            print(s)
        if d == 1:
            print("      lowest exact levels (B, shared-link space): " + " ".join(f"({Bx[k]:+.2f},{wl[k]:.2f})" for k in range(min(16, len(wx)))))
print(f"done [{time.time()-t0:.0f}s]")
