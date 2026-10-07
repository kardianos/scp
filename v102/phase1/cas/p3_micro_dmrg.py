#!/usr/bin/env python3
"""P3/P5 — xi from the microscopic chain (DMRG).
Usage: p3_micro_dmrg.py g t L chi [linkcap]     (Delta = 1)
Sites alternate cell, link, cell, ...: L cells (end cells cap 1, bulk
cap 2, R6), L-1 links (scalar space mode, c <= linkcap, default 2).
A chemical potential mu = 0.1 per quantum pins full occupation without
making space free (Delta - 2 mu > 0); <Q> is checked against the sector.  xi from |<S_i.S_i+r>| sqrt(r) ~ exp(-r/xi), r in [6, 20],
pairs centred on the middle cell.  'eff' as g runs the effective
spin-1 chain (spin-1/2 ends) with the same L and window.
"""
import sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dmrg import DMRG, build_mpo, heisenberg_chain

t0 = time.time()
gs, ts, L, chi = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
linkcap = int(sys.argv[5]) if len(sys.argv) > 5 else 2
mu = 0.1        # a space quantum still costs Delta - 2 mu = 0.8 Delta: the pin cannot fill links

def cell_space(cap):
    st = [(a, N - a) for N in range(cap + 1) for a in range(N + 1)]
    pos = {s: k for k, s in enumerate(st)}; d = len(st)
    b1 = np.zeros((d, d)); b2 = np.zeros((d, d))
    for (a, b), k in pos.items():
        if a > 0: b1[pos[(a - 1, b)], k] = np.sqrt(a)
        if b > 0: b2[pos[(a, b - 1)], k] = np.sqrt(b)
    Nop = np.diag([a + b for a, b in st]).astype(float)
    Sz = np.diag([0.5*(a - b) for a, b in st])
    Sp = b1.T@b2                     # b1^+ b2, stays inside the cap
    return d, b1, b2, Nop, Sz, Sp

def link_space(cap):
    d = cap + 1
    c = np.diag(np.sqrt(np.arange(1, d)), 1)
    return d, c

if gs == 'eff':
    W, const, dims, ops = heisenberg_chain([0.5] + [1.0]*(L - 2) + [0.5])
    cells = list(range(L)); Szs = [o[0] for o in ops]; Sps = [o[1] for o in ops]
    sites_of_cell = cells
else:
    g = float(gs)
    dims, onsite, Szs, Sps, b1s, b2s, Ns = [], [], [], [], [], [], []
    sites_of_cell = []
    for c in range(L):
        cap = 1 if c in (0, L - 1) else 2
        d, b1, b2, Nop, Sz, Sp = cell_space(cap)
        sites_of_cell.append(len(dims))
        dims.append(d); onsite.append(-mu*Nop); Szs.append(Sz); Sps.append(Sp)
        b1s.append(b1); b2s.append(b2); Ns.append(Nop)
        if c < L - 1:
            dl, cop = link_space(linkcap)
            dims.append(dl); onsite.append((1.0 - 2*mu)*(cop.T@cop))
            Szs.append(None); Sps.append(None); b1s.append(cop); b2s.append(None); Ns.append(2*cop.T@cop)
    terms = []
    for c in range(L - 1):
        i, l, j = sites_of_cell[c], sites_of_cell[c] + 1, sites_of_cell[c + 1]
        cop = b1s[l]; Il = np.eye(dims[l])
        B1i, B2i, B1j, B2j = b1s[i], b2s[i], b1s[j], b2s[j]
        terms += [(i, [B1i, cop.T, B2j], g), (i, [B2i, cop.T, B1j], -g),
                  (i, [B1i.T, cop, B2j.T], g), (i, [B2i.T, cop, B1j.T], -g)]
        if ts != 0:
            for (Bi, Bj) in ((B1i, B1j), (B2i, B2j)):
                terms += [(i, [Bi.T, Il, Bj], ts), (i, [Bi, Il, Bj.T], ts)]
    W, const = build_mpo(dims, onsite, terms)

d = DMRG(W, dims, chi=chi, seed=7)
E = d.run(sweeps=10, etol=1e-10) + const
out = [f"{sys.argv[1:]}", f"E={E:.10f}"]
if gs != 'eff':
    Q = sum(d.expect1(Ns[s], s) for s in range(len(dims)))
    space = sum(d.expect1(b1s[s].T@b1s[s], s) for s in range(len(dims)) if Szs[s] is None)
    Qmax = 2*(L - 2) + 2
    out += [f"<Q>={Q:.6f} (sector {Qmax}: {'OK' if abs(Q - Qmax) < 1e-4 else 'WRONG SECTOR'})", f"<space>/link={space/(L-1):.5f}"]
mid = L//2
rs = np.arange(2, 25)
C = []
for r in rs:
    ci = mid - r//2; cj = ci + r
    si, sj = sites_of_cell[ci], sites_of_cell[cj]
    full_z = [None]*len(dims); full_m = [None]*len(dims)
    full_z[sj] = Szs[sj]; full_m[sj] = Sps[sj].T
    zz = d.corr_row(Szs[si], si, full_z, [sj])[sj]
    pm = d.corr_row(Sps[si], si, full_m, [sj])[sj]
    C.append(zz + pm)
C = np.array(C)
win = (rs >= 6) & (rs <= 20)
xi = -1/np.polyfit(rs[win], np.log(np.abs(C[win])*np.sqrt(rs[win])), 1)[0]
out += [f"xi={xi:.4f}", f"C(2)={C[0]:+.5f} C(10)={C[8]:+.5e} C(20)={C[18]:+.5e}", f"wall={time.time()-t0:.0f}s"]
print("  ".join(out), flush=True)
