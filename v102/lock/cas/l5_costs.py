#!/usr/bin/env python3
"""L5 — reports, not gates (../PREREG.md §4 L5).

(a) Space transport: one bond (two full cells, field singlet) whose link
    oscillator is the first site of a space lead (site pitch Omega, hop
    t_s).  Single-conversion sector, which reduces exactly to
    |field singlet> --g sqrt(6)--> |pair singlet, space at c> + lead.
    A z-coordinated tree lead maps onto a chain with first hop sqrt(z) t_s
    and later hops sqrt(z-1) t_s (half-bandwidth 2 sqrt(z-1) t_s).
    Reported: the bound-state weight on the field level, the late-time
    field survival, and the band-edge condition.
(b) The t > 0 margin under the S-law: n = 6 ring, general basis, sum c <= 2
    (registered) and 3 (check), S^z = 0; total spin of the ground state.
"""
import sys, time
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import DELTA, basis_ring, hamiltonian, total_S

t0 = time.time()
def stamp(): return f"[{time.time()-t0:5.0f}s]"
gr = 0.1; g = gr*DELTA; geff = g*np.sqrt(6.0)

def bond_lead(ts, z, M=6000):
    """Tridiagonal chain: site 0 = field singlet (energy 0), site 1 = space at
    the link (Delta), sites 2.. = lead (Delta).  Returns (bound weight on
    site 0, band bottom, late-time survival)."""
    d = np.full(M + 2, DELTA); d[0] = 0.0
    e = np.empty(M + 1); e[0] = geff
    if z == 1:                       # semi-infinite 1-D lead attached at its end
        e[1:] = ts
        half = 2*ts
    else:
        e[1] = np.sqrt(z)*ts; e[2:] = np.sqrt(z - 1)*ts
        half = 2*np.sqrt(z - 1)*ts
    w, V = eigh_tridiagonal(d, e)
    bottom = DELTA - half
    a0 = V[0]**2
    below = w < bottom - 1e-9
    Zb = a0[below].sum()
    # late-time survival, averaged over t in [T/2, T] with T well before the lead-end echo
    vmax = 2*half if half > 0 else 1.0
    T = 0.4*M/max(vmax, 1e-9)
    ts_ = np.linspace(T/2, T, 400)
    amp = (a0[None, :]*np.exp(-1j*np.outer(ts_, w))).sum(1)
    return Zb, bottom, np.mean(np.abs(amp)**2), T

print(f"# L5(a) space transport vs the S1 oscillator; one bond, g/Delta = {gr}, g sqrt(lambda_max) = {geff:.5f}, Delta = {DELTA:.5f}")
print("  lead            t_s     t_s/Delta  band bottom - pair   bound weight Z_b   late survival   (T)")
for ts_ratio in (0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6):
    ts = ts_ratio*DELTA
    Zb, bottom, surv, T = bond_lead(ts, 1)
    print(f"  1-D             {ts:.4f}  {ts_ratio:5.2f}      {bottom:+.4f}             {Zb:.4f}            {surv:.4f}        ({T:.0f})")
s_k = 0.06; zl = 2*(7.28 - 1)
print(f"\n  reading t_s = s_k = {s_k} (s_k is a diffusive conductance, not an amplitude; this is a reading)")
for (name, z) in (("1-D fiber lead", 1), (f"tree, z_l = {zl:.1f}", zl)):
    Zb, bottom, surv, T = bond_lead(s_k, z)
    print(f"  {name:22s}: band bottom - pair = {bottom:+.4f}   Z_b = {Zb:.4f}   late survival = {surv:.4f}")
print(f"  lattice edge estimate z_l t_s = {zl*s_k:.4f} vs Delta = {DELTA:.4f}: band bottom - pair = {DELTA - zl*s_k:+.4f} "
      f"({'holds' if DELTA - zl*s_k > 0 else 'band reaches the pair: decays'})")
print(f"  condition: the lock keeps a bound field state iff the space band's lower edge stays above the pair, "
      f"i.e. half-bandwidth < Delta = {DELTA:.4f}")

# ---------------------------------------------------------------- (b)
print(f"\n# L5(b) t > 0 margin under the S-law, n = 6 ring, general basis, S^z = 0, g/Delta = {gr}  {stamp()}", flush=True)
for cmax in (2, 3):
    st, ix = basis_ring(6, cmax=cmax, sz=0, local=False)
    print(f"  sum c <= {cmax}: dim = {len(st)}  {stamp()}", flush=True)
    prev = None
    for tr in (0.0, 0.05, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25):
        H = hamiltonian(6, g, st, ix, t=tr*DELTA)
        w, v = eigsh(H, k=6, which='SA', tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
        S = total_S(v[:, 0], st, ix)
        print(f"    t/Delta = {tr:.2f} (t = {tr*DELTA:.4f}): E0 = {w[0]:+.6f}  S_ground = {S:.3f}  "
              f"{'lock (S=0)' if abs(S) < 1e-3 else 'NOT the lock'}  {stamp()}", flush=True)
print(f"done {stamp()}")
