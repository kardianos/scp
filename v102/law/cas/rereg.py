#!/usr/bin/env python3
"""Part R — the two defective Phase 1 bars, re-registered on new cases
(../PREREG.md Part R).

R-P1'  n = 6 microscopic ring (sum c <= 2, S^z = 0), g/Delta = 0.1,
       t/Delta in {0.05, 0.15}; the first 12 excitation energies WITH
       multiplicity vs the spin-1 Heisenberg ring; one fitted J_eff.
R-P6c' new slit geometry: separation +-16, half-width 2, source
       momentum 1.0, centred 4 rows off axis.  The eraser equals the
       coherent single-slit sum.
"""
import sys
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
import scipy.sparse as sp
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/../../phase1/cas')
from p1_micro import micro, S2_sz0, spin1, sector

print("# R-P1'  n=6 ring, levels with multiplicity (S^z = 0 sector)")
H, szd = spin1(6)
wh = np.linalg.eigvalsh(sector(H, szd, 0).toarray())
hx = wh[1:13] - wh[0]
ok_all = True
for t in (0.05, 0.15):
    Hm, space, basis = micro(6, 0.1, t, cmax=2)
    w, v = eigsh(Hm, k=20, which='SA', tol=1e-12)
    o = np.argsort(w); w = w[o]; v = v[:, o]
    ex = w[1:13] - w[0]
    J = float(ex@hx/(hx@hx))
    miss = np.abs(ex/(J*hx) - 1)
    S0 = S2_sz0(v[:, 0], basis)
    ok = miss.max() < 0.05 and abs(S0) < 1e-6
    ok_all &= ok
    print(f"  t/D={t:.2f}: J_eff={J:.6f}  max miss over 12 levels={miss.max():.3e}  S_ground={S0:.4f}  {'PASS' if ok else 'MISS'}", flush=True)
print(f"R-P1' {'PASS' if ok_all else 'MISS'}")

print("\n# R-P6c'  new slit geometry")
NX, NY = 140, 161; xw, xs = 60, 125; yc = NY//2; sep, hw = 16, 2
def open_cells(blocked):
    ok = np.ones((NX, NY), bool); ok[xw, :] = False
    for (c, on) in ((yc - sep, 'S1'), (yc + sep, 'S2')):
        if on not in blocked: ok[xw, c - hw:c + hw + 1] = True
    return ok
def run(blocked):
    ok = open_cells(blocked)
    idx = -np.ones((NX, NY), int); cells = np.argwhere(ok); idx[ok] = np.arange(len(cells))
    r, c = [], []
    for (x, y) in cells:
        for (dx, dy) in ((1, 0), (0, 1)):
            X, Y = x + dx, y + dy
            if X < NX and Y < NY and ok[X, Y]:
                r += [idx[x, y], idx[X, Y]]; c += [idx[X, Y], idx[x, y]]
    Hs = -sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(cells), len(cells)))
    X, Y = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    psi0 = np.exp(-((X - 32)**2 + (Y - yc - 4)**2)/(2*7.0**2) + 1j*1.0*X)
    v = psi0[ok]/np.linalg.norm(psi0[open_cells(())])
    traj = expm_multiply(-1j*Hs, v, start=0, stop=80, num=161, endpoint=True)
    return traj[:, idx[xs, :]]
F1, F2, Fb = run(('S2',)), run(('S1',)), run(())
I_er = (np.abs((F1/np.sqrt(2) + F2/np.sqrt(2))/np.sqrt(2))**2).sum(0)
I_coh = (np.abs(F1 + F2)**2).sum(0)/4
d = np.max(np.abs(I_er - I_coh))/np.max(I_coh)
print(f"  eraser vs coherent single-slit sum: max rel diff = {d:.2e}  {'PASS' if d < 1e-10 else 'MISS'}  (bar 1e-10)")
I_b = (np.abs(Fb)**2).sum(0)
print(f"  reported, not gated: both-open vs coherent sum (scattering between slits) = {np.max(np.abs(I_b - (np.abs(F1+F2)**2).sum(0)))/np.max(I_b):.2e}")
