#!/usr/bin/env python3
"""L4(a) — collect the effective-chain DMRG runs (l4_mass.py) and
extrapolate the gap as D + a/(L+1)^2 (P2's form).  Gate: Delta_inf/J_h2 > 0.2.
Reported: Delta_inf/J against 0.4105, and the Heisenberg reference with
the same L and procedure."""
import glob
import numpy as np
here = __file__.rsplit('/', 1)[0]
E = {}
for f in glob.glob(here + '/out/l4_*_*_*.out'):
    for line in open(f):
        p = line.split()
        if len(p) >= 6 and p[0] in ('slaw', 'heis'):
            E[(p[0], int(p[1]), int(p[2]))] = (float(p[4]), float(p[5]))
print("# L4(a) effective-chain gap, chi = 64; chain [cap1] + [cap2]*L + [cap1]")
res = {}
for mode in ('heis', 'slaw'):
    Ls = sorted({k[1] for k in E if k[0] == mode and (mode, k[1], 0) in E and (mode, k[1], 1) in E})
    gaps = np.array([E[(mode, L, 1)][0] - E[(mode, L, 0)][0] for L in Ls])
    J = E[(mode, Ls[0], 0)][1] if Ls else float('nan')
    for L, gp in zip(Ls, gaps):
        print(f"  {mode}  L_bulk={L:3d}  E0={E[(mode, L, 0)][0]:+.8f}  gap/J={gp:.6f}")
    if len(Ls) >= 2:
        A = np.stack([np.ones(len(Ls)), 1/(np.array(Ls) + 1.0)**2], 1)
        coef, *_ = np.linalg.lstsq(A, gaps, rcond=None)
        resid = gaps - A@coef
        res[mode] = (Ls, gaps, coef)
        print(f"  {mode}: Delta_inf/J = {coef[0]:.5f}  (fit residual max {np.abs(resid).max():.1e}; J = {J:.6f})")
if 'heis' in res:
    print(f"  Heisenberg reference vs literature 0.4105: miss {abs(res['heis'][2][0] - 0.4105)/0.4105:.2%} "
          f"(P2 at L = 40, 60, 80 gave 0.41317)")
if 'slaw' in res and 'heis' in res:
    Lh = dict(zip(res['heis'][0], res['heis'][1]))
    print("  ratio S-law/Heisenberg at equal L: " + "  ".join(f"L={L}:{g/Lh[L]:.4f}" for L, g in zip(res['slaw'][0], res['slaw'][1]) if L in Lh))
    D = res['slaw'][2][0]
    print(f"  S-law Delta_inf/J_h2 = {D:.5f}; vs 0.4105: {D/0.4105 - 1:+.2%}; vs the same-procedure Heisenberg "
          f"{res['heis'][2][0]:.5f}: {D/res['heis'][2][0] - 1:+.2%}")
    print(f"L4(a) {'PASS' if D > 0.2 else 'MISS'}  (bar Delta_inf/J_h2 > 0.2)")
