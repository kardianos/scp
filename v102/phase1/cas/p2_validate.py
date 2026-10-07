#!/usr/bin/env python3
"""P2 — validate the DMRG instrument on the effective chain.
Open spin-1 Heisenberg chain with spin-1/2 end cells (R6: a degree-1
end cell holds one quantum).  J = 1.  Literature: e_inf = -1.401484,
gap 0.41050, xi 6.03 (White & Huse 1993)."""
import sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dmrg import DMRG, heisenberg_chain, spin_ops

def chain(Lb): return [0.5] + [1.0]*Lb + [0.5]

def ground(Lb, chi, Sz=None, seed=1):
    W, c, dims, ops = heisenberg_chain(chain(Lb), Sz_target=Sz, lam=1.0)
    d = DMRG(W, dims, chi=chi, seed=seed)
    E = d.run(sweeps=12, etol=1e-11) + c
    return E, d, ops

t0 = time.time()
print("# P2 instrument validation (effective chain, spin-1/2 end cells)")
res = {}
for Lb, chi in ((40, 64), (60, 64), (80, 64)):
    E0, d0, ops = ground(Lb, chi, Sz=0.0)
    E1, _, _ = ground(Lb, chi, Sz=1.0, seed=2)
    E2, _, _ = ground(Lb, chi, Sz=2.0, seed=3)
    res[(Lb, chi)] = (E0, E1, E2, d0, ops)
    print(f"  L_bulk={Lb:3d} chi={chi:3d}  E0={E0:.9f}  gap1={E1-E0:.6f}  E(Sz=2)-E0={E2-E0:.6f}"
          f"  ratio={(E2-E0)/(E1-E0):.4f}   [{time.time()-t0:.0f}s]", flush=True)
e_inf = (res[(80, 64)][0] - res[(60, 64)][0])/20
print(f"P2a energy per bulk bond (E0(80)-E0(60))/20 = {e_inf:.6f}   literature -1.401484   miss {abs(e_inf+1.401484):.1e}  (bar 1e-4)")
Ls = np.array([40, 60, 80.]); gs = np.array([res[(L, 64)][1] - res[(L, 64)][0] for L in (40, 60, 80)])
A = np.stack([np.ones(3), 1/(Ls + 1)**2], 1)
coef = np.linalg.lstsq(A, gs, rcond=None)[0]
print(f"P2b gap(L) = D_inf + a/(L+1)^2 fit: D_inf = {coef[0]:.5f}   literature 0.41050   miss {abs(coef[0]-0.4105)/0.4105:.2%}  (bar 2%)")
print("    chi convergence: see p3_eff (L=48, chi=48) for the same window")
# xi from bulk correlations, L=80 chi=64
E0, d0, ops = res[(80, 64)][0], res[(80, 64)][3], res[(80, 64)][4]
N = len(d0.dims); mid = N//2
Zs = [o[0] for o in ops]; Ps = [o[1] for o in ops]
rs = np.arange(2, 31)
C = []
for r in rs:
    i = mid - r//2; j = i + r
    zz = d0.corr_row(Zs[i], i, Zs, [j])[j]
    pm = d0.corr_row(Ps[i], i, [P.T for P in Ps], [j])[j]
    C.append(zz + pm)        # S.S = zz + (S+S- + S-S+)/2 = zz + pm for a real symmetric state
C = np.array(C)
win = (rs >= 6) & (rs <= 24)
y = np.log(np.abs(C[win])*np.sqrt(rs[win]))
slope = np.polyfit(rs[win], y, 1)[0]
xi = -1/slope
print(f"P2c xi from |<S_i.S_i+r>| sqrt(r) ~ e^(-r/xi), r in [6,24]: xi = {xi:.3f}   literature 6.03   miss {abs(xi-6.03)/6.03:.2%}  (bar 5%)")
print("    r   <S.S>")
for r, c in zip(rs[::4], C[::4]): print(f"   {r:3d}  {c:+.6e}")
ratios = [(res[k][2]-res[k][0])/(res[k][1]-res[k][0]) for k in ((40,64),(60,64),(80,64))]
print(f"P4b(open chain) (E(Sz=2)-E0)/(2 gap1) = " + ", ".join(f"{r/2:.4f}" for r in ratios) + "   (bar >= 0.98)")
print(f"done [{time.time()-t0:.0f}s]")
