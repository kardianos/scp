#!/usr/bin/env python3
"""L3 — the lock's quantum walks by conversion (../PREREG.md §3 L3, as
amended: n = 10, sum c <= 3).

Kick |phi> = S^+_0 |0> / norm on the S-law ring; spin profile
m_r(t) = <phi(t)| S^z_r |phi(t)> (sums to 1).  Compare with the same kick
on the spin-1 Heisenberg ring n = 10 at tau = J t, J fitted on this ring
(L1 procedure).  Gates: sup |m^S - m^H| < 0.05 over tau in [0, 4] at
g/Delta = 0.05 and 0.1; field never travels (printed); g = 0: nothing
moves (printed)."""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import DELTA, basis_ring, hamiltonian, diag_obs, splus, heis_ring, fit_J

t0 = time.time()
def stamp(): return f"[{time.time()-t0:5.0f}s]"
n, cmax = 10, 3
TAU = np.linspace(0, 4, 81)

# ---------------------------------------------------------------- Heisenberg reference
H, Z, P = heis_ring(n)
szd = sum(Z).diagonal()
k0 = np.where(np.abs(szd) < 1e-9)[0]; k1 = np.where(np.abs(szd - 1) < 1e-9)[0]
H0h = H[k0][:, k0]; H1h = H[k1][:, k1]
wh, vh = eigsh(H0h, k=14, which='SA', tol=1e-12); o = np.argsort(wh); wh = wh[o]; vh = vh[:, o]
Sp0 = P[0][k1][:, k0]
ph = Sp0@vh[:, 0]; ph /= np.linalg.norm(ph)
trh = expm_multiply(-1j*H1h, ph, start=0, stop=TAU[-1], num=len(TAU), endpoint=True)
Zd = np.stack([Z[r].diagonal()[k1] for r in range(n)], 1)
mH = (np.abs(trh)**2)@Zd
print(f"# L3 walk by conversion, n = {n}, sum c <= {cmax}")
print(f"  Heisenberg ring: dim S^z=0 {len(k0)}, S^z=1 {len(k1)}; E0 = {wh[0]:.6f}; sum_r m_r(0) = {mH[0].sum():.6f}  {stamp()}", flush=True)

# ---------------------------------------------------------------- S-law sectors
st0, ix0 = basis_ring(n, cmax=cmax, sz=0)
st1, ix1 = basis_ring(n, cmax=cmax, sz=2)
D0 = hamiltonian(n, 0.0, st0, ix0); X0 = hamiltonian(n, 1.0, st0, ix0) - D0
D1 = hamiltonian(n, 0.0, st1, ix1); X1 = hamiltonian(n, 1.0, st1, ix1) - D1
Sp = splus(0, st0, ix1)
C1, N1, Z1 = diag_obs(n, st1)
print(f"  S-law: dim S^z=0 {len(st0)}, S^z=1 {len(st1)}  {stamp()}", flush=True)

ok_all = True
res = {}
for gr in (0.05, 0.1):
    g = gr*DELTA
    H0 = (D0 + g*X0).tocsr(); H1 = (D1 + g*X1).tocsr()
    w, v = eigsh(H0, k=14, which='SA', tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
    J, miss = fit_J(w[1:13] - w[0], wh[1:13] - wh[0])
    phi = Sp@v[:, 0]; phi /= np.linalg.norm(phi)
    tr = expm_multiply(-1j*H1, phi, start=0, stop=TAU[-1]/J, num=len(TAU), endpoint=True)
    Pr = np.abs(tr)**2
    m = Pr@Z1
    dev = np.abs(m - mH).max()
    Nd = np.abs(Pr@N1 - 2).max(); cmx = (Pr@C1).max()
    ok = dev < 0.05
    ok_all &= ok
    res[gr] = (m, J)
    print(f"  g/Delta={gr:.2f}: J_fit={J:.6f} (L1-style miss {miss.max():.2e});  sum_r m_r = {m.sum(1).min():.6f}..{m.sum(1).max():.6f};"
          f"  sup|m^S - m^H| over tau in [0,4] = {dev:.4f}  {'ok' if ok else 'MISS'};  field: max|<N_i>-2| = {Nd:.4f} "
          f"(<= 2 max<c> = {2*cmx:.4f})  {stamp()}", flush=True)
    if gr == 0.1:
        tr0 = expm_multiply(-1j*D1, phi, start=0, stop=TAU[-1]/J, num=9, endpoint=True)
        m0 = (np.abs(tr0)**2)@Z1
        print(f"    g = 0 (same kicked state, same times): max |m_r(t) - m_r(0)| = {np.abs(m0 - m0[0]).max():.1e}  (nothing moves)")
m5, J5 = res[0.05]; m1, J1 = res[0.1]
print(f"  collapse between g/Delta = 0.05 and 0.1 in units of 1/J: sup|m(0.05) - m(0.1)| = {np.abs(m5 - m1).max():.4f}; "
      f"J ratio {J1/J5:.4f} (2g^2/Delta ratio 4)")
print("  profile at tau = 0, 1, 2, 3, 4 (r = 0..5; S-law g/D=0.1 | Heisenberg):")
dist = np.minimum(np.arange(n), n - np.arange(n))
for k in (0, 20, 40, 60, 80):
    fS = [m1[k][dist == r].sum() for r in range(n//2 + 1)]; fH = [mH[k][dist == r].sum() for r in range(n//2 + 1)]
    print(f"   tau={TAU[k]:.0f}: " + " ".join(f"{a:+.3f}" for a in fS) + "  |  " + " ".join(f"{b:+.3f}" for b in fH))
w2S = [(m1[k]*dist**2).sum() for k in range(len(TAU))]; w2H = [(mH[k]*dist**2).sum() for k in range(len(TAU))]
print("  spread <d^2>(tau): " + " ".join(f"{TAU[k]:.1f}:{w2S[k]:.3f}/{w2H[k]:.3f}" for k in range(0, 81, 10)))
print(f"L3 {'PASS' if ok_all else 'MISS'}  {stamp()}")
