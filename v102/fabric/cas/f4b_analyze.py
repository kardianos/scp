#!/usr/bin/env python3
"""F4(b) analysis: the hole's spin cloud on the long chains.
For each L_bulk: the S^z = 1/2 profile summed over both sides at distance d
from the hole; its decay length from a log-linear fit of |a_d| over
d = 2 .. Lh - 8 (clear of the end cells); the moment within d <= 3, 6, 12;
the residual moment within 4 cells of each end; and E(3/2) - E(1/2) against
the no-hole chain's gap at the same length (../../lock/ L4: 0.4951, 0.4713
at L_bulk = 24, 32 and 0.4506 at 48; Delta_inf = 0.4357)."""
import sys, glob
import numpy as np
here = __file__.rsplit('/', 1)[0]
gap_nohole = {24: 0.4951, 32: 0.4713, 48: 0.4506}
print("# F4(b) the hole's cloud on effective open chains (reading V), g/Delta = 0.1, DMRG chi = 64")
res = {}
for f in sorted(glob.glob(here + '/out/f4b_*_1.out')):
    Lb = int(f.split('_')[-2])
    lines = open(f).read().splitlines()
    if len(lines) < 2: continue
    E1 = float(lines[0].split('E/J=')[1].split()[0])
    p = np.array([float(x) for x in lines[1].split()[1:]])
    L = len(p); h = L//2; Lh = h - 1
    a = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
    ds = np.arange(2, Lh - 7)
    sl, ic = np.polyfit(ds, np.log(np.abs(a[ds])), 1)
    lam = -1/sl
    cum = {r: a[:r + 1].sum() for r in (3, 6, 12) if r < len(a)}
    endm = p[:4].sum() + p[-4:].sum()
    f3 = f.replace('_1.out', '_3.out')
    E3 = None
    try:
        E3 = float(open(f3).read().splitlines()[0].split('E/J=')[1].split()[0])
    except Exception:
        pass
    res[Lb] = lam
    print(f"  L_bulk={Lb:3d}: decay length {lam:.3f} grains (fit d = 2..{Lh-8});  moment within d<=3: {cum.get(3, float('nan')):+.4f}, "
          f"d<=6: {cum.get(6, float('nan')):+.4f}, d<=12: {cum.get(12, float('nan')):+.4f};  residual at the ends {endm:+.4f};  "
          + (f"E(3/2) - E(1/2) = {E3-E1:.4f} J" if E3 is not None else "E(3/2) pending")
          + (f" (no-hole gap at this length {gap_nohole[Lb]:.4f} J)" if Lb in gap_nohole else ""))
if len(res) >= 2:
    ks = sorted(res)
    a, b = res[ks[-2]], res[ks[-1]]
    print(f"  N2: decay length {a:.3f} (L_bulk={ks[-2]}) -> {b:.3f} (L_bulk={ks[-1]}): change {abs(b-a)/a*100:.1f}%  "
          f"({'PASS' if abs(b-a)/a < 0.10 else 'MISS'}: bar < 10%)")
