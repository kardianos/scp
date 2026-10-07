#!/usr/bin/env python3
"""P6 — the slit from the unlocked, linearized field (a lattice calculator).

Under R1 the linear transport is the same hop for both plane components,
so a single quantum's two-component wavefunction evolves component by
component (R2: the one-quantum sector of the free boson field).  A lone
quantum cannot convert: the conversion operator needs two quanta, one
per link end (gates/WEDGE.md W1).

Square lattice, continuous-time hop, hard wall with two slits; a
Gaussian packet aimed at the wall; screen intensity integrated in time.
Interference fraction V = max|I - I_1 - I_2| / max(I_1 + I_2), where
I_1, I_2 are the single-slit runs.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

NX, NY = 140, 161
xw, xs = 60, 125                       # wall column, screen column
yc = NY//2; sep, w = 12, 3             # slit centres yc +- sep, half-width w

def open_cells(blocked):
    ok = np.ones((NX, NY), bool)
    ok[xw, :] = False
    for (c, on) in ((yc - sep, 'S1'), (yc + sep, 'S2')):
        if on not in blocked:
            ok[xw, c - w:c + w + 1] = True
    return ok

def hamiltonian(ok):
    idx = -np.ones((NX, NY), int); cells = np.argwhere(ok)
    idx[ok] = np.arange(len(cells))
    r, c = [], []
    for (x, y) in cells:
        for (dx, dy) in ((1, 0), (0, 1)):
            X, Y = x + dx, y + dy
            if X < NX and Y < NY and ok[X, Y]:
                r += [idx[x, y], idx[X, Y]]; c += [idx[X, Y], idx[x, y]]
    H = -sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(cells), len(cells)))
    return H, idx

def run(blocked, amp=1.0):
    ok = open_cells(blocked); H, idx = hamiltonian(ok)
    X, Y = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    psi0 = amp*np.exp(-((X - 30)**2 + (Y - yc)**2)/(2*6.0**2) + 1j*1.2*X)
    v = psi0[ok]; v = v/np.linalg.norm(psi0[open_cells(())])
    ts = np.linspace(0, 70, 141)
    traj = expm_multiply(-1j*H, v, start=0, stop=70, num=141, endpoint=True)
    scr = idx[xs, :]
    return (np.abs(traj[:, scr])**2).sum(0), traj, idx

I_both, T_both, idx_b = run(())
I_1, T_1, idx_1 = run(('S2',))
I_2, T_2, idx_2 = run(('S1',))
V = lambda I: np.max(np.abs(I - I_1 - I_2))/np.max(I_1 + I_2)
print(f"# P6 slit calculator  (lattice {NX}x{NY}, slits at y={yc-sep},{yc+sep}, half-width {w})")
print(f"  scalar: interference fraction V_scalar = {V(I_both):.4f}")

# (a) two-component field, arbitrary polarization, no polarizers: same hop per component
rng = np.random.default_rng(6); al, ph = rng.uniform(0, np.pi), rng.uniform(0, 2*np.pi)
p = np.array([np.cos(al), np.sin(al)*np.exp(1j*ph)])
I_T = np.abs(p[0])**2*I_both + np.abs(p[1])**2*I_both     # component-wise identical evolution
# do it honestly: evolve each component separately
I_T2 = sum(np.abs(pa)**2*run(())[0] for pa in p)
print(f"P6a two-component vs scalar: max|I_T - I_S|/max I_S = {np.max(np.abs(I_T2 - I_both))/np.max(I_both):.2e}  (bar 1e-10)")

# (b) orthogonal polarizers: slit 1 passes component 1 only, slit 2 component 2 only; input (1,1)/sqrt2
psi1 = T_1/np.sqrt(2); psi2 = T_2/np.sqrt(2)     # component 1 sees S2 blocked; component 2 sees S1 blocked
scr1, scr2 = idx_1[xs, :], idx_2[xs, :]
I_pol = (np.abs(psi1[:, scr1])**2 + np.abs(psi2[:, scr2])**2).sum(0)
V_pol = np.max(np.abs(I_pol - (I_1 + I_2)/2))/np.max((I_1 + I_2)/2)
print(f"P6b polarizers: interference fraction V_pol = {V_pol:.2e}  (bar < 0.01)")

# (c) eraser: project on (1,1)/sqrt2 at the screen
I_er = (np.abs((psi1[:, scr1] + psi2[:, scr2])/np.sqrt(2))**2).sum(0)
V_er = np.max(np.abs(I_er - (I_1 + I_2)/4))/np.max((I_1 + I_2)/4)
print(f"P6c eraser: V_eraser = {V_er:.6f}   V_scalar = {V(I_both):.6f}   |diff| = {abs(V_er - V(I_both)):.2e}  (bar 1e-10 as registered)")
I_coh = (np.abs(T_1[:, scr1] + T_2[:, scr2])**2).sum(0)/4
print(f"    identity check: eraser = coherent sum of the single-slit fields, max rel diff = {np.max(np.abs(I_er - I_coh))/np.max(I_coh):.2e}")
I_S_vs_sum = np.max(np.abs(I_both - (np.abs(T_1[:, scr1] + T_2[:, scr2])**2).sum(0)))/np.max(I_both)
print(f"    both-open field vs sum of single-slit fields (multiple scattering between slits): {I_S_vs_sum:.2e}")
