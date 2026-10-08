#!/usr/bin/env python3
"""L2 — the lock holds its field (real time), and the R-c control decays
(../PREREG.md §3 L2).  Also the L6 numbers (internal energy and canonical
temperature after the sudden start, n = 6).

(a) S-law rings n = 6 (exact, all S^z) and n = 8 (sum c <= 3, all S^z);
    initial states: Neel product and three seeded random Jones products;
    evolve to T = 100/J (J fitted on the same ring, L1 procedure).
(b) Control R-c: n = 4, sum c <= 1, each link's mode replaced by M = 400
    modes uniform over +-W = +-1.0 around the pair resonance, couplings
    g/sqrt(M); Neel start; golden rule Gamma = sum_l (pi g^2/W) <A_l^+A_l>.
    Same-ring S1 check: n = 4, sum c <= 1, one oscillator per link.
"""
import sys, time, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import (DELTA, EF, OMEGA, basis_ring, hamiltonian, diag_obs, product_state,
                  neel_cells, random_jones_cells, AdagA_two_full, heis_ring, fit_J)

t0 = time.time()
def stamp(): return f"[{time.time()-t0:5.0f}s]"
gr = 0.1; g = gr*DELTA
lam = np.linalg.eigvalsh(AdagA_two_full()).max()
BOUND_C = 5*gr*gr*lam
BOUND_F = 1 - 2.5*gr*gr*lam
print(f"# L2 hold.  g/Delta = {gr}, g = {g:.5f}, lambda_max(A^+A) = {lam:.4f}")
print(f"  bars: <c_l^+c_l> <= 5 (g/D)^2 lam = {BOUND_C:.4f} at all links and times; "
      f"second-half mean field fraction >= {BOUND_F:.4f}; relative energy drift < 1e-8")

def ring_J(n, cmax):
    st, ix = basis_ring(n, cmax=cmax, sz=0); H = hamiltonian(n, g, st, ix)
    w = np.sort(eigsh(H, k=14, which='SA', tol=1e-12)[0]) if H.shape[0] > 3000 else np.linalg.eigvalsh(H.toarray())[:14]
    Hh, Z, P = heis_ring(n); keep = np.where(np.abs(sum(Z).diagonal()) < 1e-9)[0]
    hw = np.linalg.eigvalsh(Hh[keep][:, keep].toarray())[:13]
    J, miss = fit_J(w[1:13] - w[0], hw[1:13] - hw[0])
    return J, w[0], w[1] - w[0]

def initial_states(n):
    out = [("neel", neel_cells(n))]
    for seed in (1, 2, 3):
        out.append((f"jones{seed}", random_jones_cells(n, np.random.default_rng(1000 + seed))))
    return out

ok_a = True
L6 = {}
for (n, cmax) in ((6, None), (8, 3)):
    J, E0, gap = ring_J(n, cmax)
    T = 100/J; NT = 401; ts = np.linspace(0, T, NT)
    st, ix = basis_ring(n, cmax=cmax, sz=None)
    H = hamiltonian(n, g, st, ix)
    C, N, Z = diag_obs(n, st)
    Etot0 = n*EF[2]
    print(f"\n  ring n={n} cmax={cmax}: dim(all S^z)={H.shape[0]}  J_fit={J:.6f}  E0={E0:+.6f}  gap={gap:.6f}  "
          f"T=100/J={T:.0f}  {stamp()}", flush=True)
    dense = H.shape[0] <= 8000
    if dense:
        w, V = np.linalg.eigh(H.toarray())
    for name, cells in initial_states(n):
        v0 = product_state(cells, st, ix, n); v0 /= np.linalg.norm(v0)
        if dense:
            a0 = V.T@v0
            traj = (V@(a0[:, None]*np.exp(-1j*np.outer(w, ts)))).T
            if n == 6:
                p = np.abs(a0)**2
                L6[name] = (p, w)
        else:
            traj = expm_multiply(-1j*H, v0, start=0, stop=T, num=NT, endpoint=True)
        P = np.abs(traj)**2
        cl = P@C                       # (NT, n)
        ff = (P@N).sum(1)/(2*n)
        En = np.real(np.einsum('ti,ti->t', traj.conj(), (H@traj.T).T))
        nrm = P.sum(1)
        drift_abs = np.max(np.abs(En - En[0])); drift_rel = drift_abs/Etot0
        cmax_t = cl.max(); ff2 = ff[NT//2:].mean()
        ok = drift_rel < 1e-8 and cmax_t <= BOUND_C and ff2 >= BOUND_F
        ok_a &= ok
        print(f"    {name:7s}: <H'>(0)={En[0]:+.2e}  drift abs={drift_abs:.1e} (={drift_abs/J:.1e} J) rel to total "
              f"{Etot0:.3f}: {drift_rel:.1e}  norm err={np.max(np.abs(nrm-1)):.1e}  max<c_l>={cmax_t:.4f}  "
              f"mean<c_l>={cl.mean():.4f}  field fraction: min={ff.min():.4f} 2nd-half mean={ff2:.4f}  "
              f"{'ok' if ok else 'MISS'}  {stamp()}", flush=True)
print(f"L2(a) {'PASS' if ok_a else 'MISS'}")

# ---------------------------------------------------------------- (b) control R-c
print(f"\n# L2(b) control R-c: pitchless space (S1 refused)  {stamp()}", flush=True)
n = 4; M = 400; Wd = 1.0
rel = 2*(EF[2] - EF[1])                       # pair release energy; resonance
eps = np.linspace(-Wd, Wd, M)                 # detunings of the M continuum modes
gm = g/np.sqrt(M)
fullc = [(a, 2 - a) for a in range(3)]; halfc = [(a, 1 - a) for a in range(2)]
states = []; idx = {}
for cells in itertools.product(fullc, repeat=n):
    idx[(None, None, cells)] = len(states); states.append((None, None, cells))
for l in range(n):
    i, j = l, (l + 1) % n
    for m in range(M):
        for cells in itertools.product(*[(halfc if k in (i, j) else fullc) for k in range(n)]):
            idx[(l, m, cells)] = len(states); states.append((l, m, cells))
rows, cols, vals = [], [], []
for k, (l, m, cells) in enumerate(states):
    if l is not None:
        rows.append(k); cols.append(k); vals.append(eps[m]); continue
    for ll in range(n):
        i, j = ll, (ll + 1) % n
        for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):
            if cells[i][ai] == 0 or cells[j][aj] == 0: continue
            amp = sgn*np.sqrt(cells[i][ai]*cells[j][aj])
            nc = list(cells)
            x = list(nc[i]); x[ai] -= 1; nc[i] = tuple(x)
            y = list(nc[j]); y[aj] -= 1; nc[j] = tuple(y)
            for m in range(M):
                k2 = idx[(ll, m, tuple(nc))]
                rows += [k2, k]; cols += [k, k2]; vals += [gm*amp, gm*amp]
Hc = sp.csr_matrix((vals, (rows, cols)), shape=(len(states), len(states)))
v0 = np.zeros(len(states), complex)
v0[idx[(None, None, tuple((2, 0) if q % 2 == 0 else (0, 2) for q in range(n)))]] = 1
AA = 4.0*n                                    # <A_l^+A_l> on a Neel bond (2,0)|(0,2) is 4; n links
Gam = np.pi*g*g*AA/Wd
Tc = 5/Gam
ts = np.linspace(0, Tc, 201)
traj = expm_multiply(-1j*Hc, v0, start=0, stop=Tc, num=201, endpoint=True)
isconv = np.array([s[0] is not None for s in states])
Pc = (np.abs(traj)**2)[:, isconv].sum(1)
print(f"  dim={len(states)}  Gamma_golden={Gam:.5f}  1/Gamma={1/Gam:.2f}  5/Gamma={Tc:.2f}  "
      f"recurrence 2pi M/(2W)={2*np.pi*M/(2*Wd):.0f}")
for f in (0.2, 0.5, 1.0):
    k = int(round(f*(len(ts) - 1)))
    print(f"    t = {f*5:.1f}/Gamma: P(converted) = {Pc[k]:.4f}   golden rule 1-exp(-Gamma t) = {1-np.exp(-Gam*ts[k]):.4f}")
ok_b1 = Pc[-1] > 0.9
print(f"  R-c: P(converted) at 5/Gamma = {Pc[-1]:.4f}  (bar > 0.9)  {'ok' if ok_b1 else 'MISS'}")
# same-ring S1 check
st, ix = basis_ring(n, cmax=1, sz=None); H1 = hamiltonian(n, g, st, ix); C1, N1, _ = diag_obs(n, st)
v1 = product_state(neel_cells(n), st, ix, n)
J4 = ring_J(4, None)[0]
T1 = max(Tc, 100/J4)
traj1 = expm_multiply(-1j*H1, v1, start=0, stop=T1, num=401, endpoint=True)
cl1 = (np.abs(traj1)**2)@C1
ok_b2 = cl1.max() <= BOUND_C
print(f"  S1 same ring (n=4, sum c<=1), t in [0, {T1:.0f}]: max <c_l> = {cl1.max():.4f}  (bar {BOUND_C:.4f})  "
      f"{'ok' if ok_b2 else 'MISS'}")
print(f"L2(b) {'PASS' if (ok_b1 and ok_b2) else 'MISS'}")
print(f"L2 {'PASS' if (ok_a and ok_b1 and ok_b2) else 'MISS'}")

# ---------------------------------------------------------------- L6 numbers (n = 6)
print(f"\n# L6 (report) internal energy after the sudden start, n = 6  {stamp()}")
p, w = next(iter(L6.values()))[0], next(iter(L6.values()))[1]
st, ix = basis_ring(6, sz=None); C6, _, _ = diag_obs(6, st)
J6, E06, gap6 = ring_J(6, None)
# spin manifold = the 3^6 = 729 lowest eigenstates (all S^z); check the separation
sep = w[729] - w[728]
print(f"  spin manifold: lowest 729 levels span [{w[0]:+.5f}, {w[728]:+.5f}]; next level {w[729]:+.5f} "
      f"(separation {sep:.4f} = {sep/DELTA:.2f} Delta)")
wl = w[:729]
def mean_E(beta):
    x = np.exp(-beta*(wl - wl[0])); return (wl*x).sum()/x.sum()
from scipy.optimize import brentq
for name, (p, w_) in L6.items():
    wlow = p[:729].sum(); El = (p[:729]*w_[:729]).sum()/wlow
    Ehigh = (p[729:]*w_[729:]).sum()
    if El < wl.mean():
        beta = brentq(lambda b: mean_E(b) - El, 1e-6, 1e6/J6)
        Tn = 1/beta
    else:
        Tn = np.inf
    print(f"  {name:7s}: weight in spin manifold = {wlow:.4f}; its mean energy above E0 = {El-w[0]:.5f} "
          f"= {(El-w[0])/gap6:.2f} gaps;  canonical T = {Tn:.5f} = {Tn/gap6:.2f} gaps = {Tn/J6:.2f} J;  "
          f"energy in space-carrying levels = {Ehigh:+.5f}")
print(f"done {stamp()}")
