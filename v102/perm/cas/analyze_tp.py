#!/usr/bin/env python3
"""Analysis of the TeNPy chains (tp_chain.py), energies in J0 = 2 g^2/Delta.
W: the gap from the (1,1) member of the first multiplet (and (1,0) at
L = 24 as a check), extrapolated D + a/(L+1)^2; the (1/2,1/2) hole's
cloud and its decay length over d = 2..Lh-3 and over a common d = 2..6;
E(3/2,3/2) - E(1/2,1/2) against the gap (N5).
3+2: the no-hole gap; the mobile hole's density, hole-centred correlation
C(d) and its decay length; E(3/2) - E(1/2) against the gap (N5)."""
import re
import numpy as np
here = __file__.rsplit('/', 1)[0]

def read(f):
    try:
        L = open(f).read().splitlines()
    except FileNotFoundError:
        return None
    L = [x for x in L if x.strip()]
    if not L or 'E/J0=' not in L[0]: return None
    out = {'E': float(re.search(r'E/J0=([-0-9.]+)', L[0]).group(1)),
           'trunc': float(re.search(r'trunc=([-0-9.e+]+)', L[0]).group(1))}
    for line in L[1:]:
        k, *v = line.split()
        if k in ('Sz', 'Iz', 'nh', 'C'): out[k] = np.array([float(x) for x in v])
    return out

def decay(a, lo, hi):
    ds = np.arange(lo, hi + 1); ds = ds[(ds < len(a))]
    ds = ds[np.abs(a[ds]) > 1e-12]
    if len(ds) < 3: return float('nan')
    return -1/np.polyfit(ds, np.log(np.abs(a[ds])), 1)[0]

def f(M, mode, L, tg, t, chi):
    return read(f"{here}/out/tp_M{M}_{mode}_{L}_{tg}_t{t}_chi{chi}.out")

print("# TeNPy chains, leading order, energies in J0")
gap = {}
for L in (16, 24, 32):
    a, b = f(4, 'nohole', L, '0,0', 0, 128), f(4, 'nohole', L, '2,2', 0, 128)
    if a and b:
        gap[L] = b['E'] - a['E']
        print(f"  W L_bulk={L}: gap {gap[L]:.5f} J0 (trunc {max(a['trunc'], b['trunc']):.1e})")
c = f(4, 'nohole', 24, '2,0', 0, 128); a = f(4, 'nohole', 24, '0,0', 0, 128)
if c and a: print(f"  W L_bulk=24: the (1,0) member gives {c['E']-a['E']:.5f} J0 (same multiplet?)")
if len(gap) >= 2:
    Ls = np.array(sorted(gap), float); y = np.array([gap[int(x)] for x in Ls])
    (D, aa), *_ = np.linalg.lstsq(np.stack([np.ones_like(Ls), 1/(Ls + 1)**2], 1), y, rcond=None)
    print(f"  W gap extrapolated: {D:.5f} J0  K: {'PASS (> 0.1 J)' if D > 0.1 else 'FAIL'}")
lam = {}
for L in (16, 24, 32):
    r = f(4, 'hole', L, '1,1', 0, 128)
    if not r: continue
    p, q = r['Sz'], r['Iz']; h = len(p)//2
    prof = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
    lam[L] = (decay(prof, 2, h - 3), decay(prof, 2, 6))
    line = (f"  W hole L_bulk={L}: hole {p[h]:+.4f}, neighbours {prof[1]:+.4f}; decay length own range {lam[L][0]:.3f} (d=2..{h-3}), "
            f"common d=2..6 {lam[L][1]:.3f}; max|Sz-Iz| {np.abs(p-q).max():.1e}; trunc {r['trunc']:.1e}")
    r3 = f(4, 'hole', L, '3,3', 0, 128)
    if r3 and L in gap:
        e = r3['E'] - r['E']; line += f";  E(3/2,3/2)-E(1/2,1/2) = {e:.5f} vs gap {gap[L]:.5f} -> {'BOUND' if e < gap[L] - 1e-4 else 'above'}"
    print(line)
ks = sorted(lam)
if len(ks) >= 2:
    a, b = ks[-2], ks[-1]
    for k, name in ((0, 'own range'), (1, 'common d=2..6')):
        x, y = lam[a][k], lam[b][k]
        print(f"  W N2 ({name}): {x:.3f} -> {y:.3f}  change {abs(y-x)/x*100:.1f}%  ({'PASS' if abs(y-x)/x < 0.10 else 'MISS'})")
print()
g2 = {}
for L in (24, 40):
    a, b = f(2, 'nohole', L, '0', 0, 64), f(2, 'nohole', L, '2', 0, 64)
    if a and b: g2[L] = b['E'] - a['E']; print(f"  3+2 L_bulk={L}: no-hole gap {g2[L]:.5f} J0")
lamT = {}
for t in ('0.1', '1.0', '3.0'):
    for L in (24, 40):
        r = f(2, 'mob', L, '1', t, 64)
        if not r: continue
        C = r['C']; l = decay(C, 2, (L + 1)//2 - 3)
        # conditional (hole-centred) profile: divide by the hole weight that has a partner at distance d
        nh = r['nh']; n = len(nh)
        w = np.array([sum(nh[i]*((i + d < n) + (i - d >= 0 and d > 0)) for i in range(n)) for d in range(n)])
        c = C/np.maximum(w, 1e-12)
        lc = decay(c, 2, (L + 1)//2 - 3)
        lamT[(t, L)] = l; lamT[('c', t, L)] = lc
        print(f"  conditional profile t={t} L_bulk={L}: " + " ".join(f"{x:+.3f}" for x in c[:8]) + f"; decay length {lc:.3f}")
        line = f"  3+2 mobile hole t={t} J0 L_bulk={L}: C(0..4) " + " ".join(f"{x:+.3f}" for x in C[:5]) + f"; decay length {l:.3f}; trunc {r['trunc']:.1e}"
        r3 = f(2, 'mob', L, '3', t, 64)
        if r3 and L in g2:
            e = r3['E'] - r['E']; line += f";  E(3/2)-E(1/2) = {e:.5f} vs gap {g2[L]:.5f} -> {'BOUND' if e < g2[L] - 1e-4 else 'above'}"
        print(line)
    if (t, 24) in lamT and (t, 40) in lamT:
        x, y = lamT[(t, 24)], lamT[(t, 40)]
        print(f"  N2 under T (t={t}), raw C(d) as registered: {x:.3f} -> {y:.3f}  change {abs(y-x)/x*100:.1f}%  ({'PASS' if abs(y-x)/x < 0.10 else 'MISS'})")
        x, y = lamT[('c', t, 24)], lamT[('c', t, 40)]
        print(f"  N2 under T (t={t}), conditional profile: {x:.3f} -> {y:.3f}  change {abs(y-x)/x*100:.1f}%  ({'PASS' if abs(y-x)/x < 0.10 else 'MISS'})")

def r90(cum_profile):
    """Smallest d with |cumulative moment within d| >= 0.9 x 1/2 (range-free size)."""
    c = np.cumsum(cum_profile)
    for d in range(len(c)):
        if abs(c[d]) >= 0.45 and all(abs(c[k]) >= 0.45 for k in range(d, len(c))): return d
    return float('nan')
print()
print("# range-free size: the radius holding 90% of the hole's spin-1/2 (the cumulative moment stays >= 0.45 from there on)")
for L in (16, 24, 32):
    r = f(4, 'hole', L, '1,1', 0, 128)
    if r:
        p = r['Sz']; h = len(p)//2
        prof = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
        print(f"  W static hole L_bulk={L}: r90 = {r90(prof)}")
for t in ('0.1', '1.0', '3.0'):
    for L in (24, 40):
        r = f(2, 'mob', L, '1', t, 64)
        if r: print(f"  3+2 mobile hole t={t} L_bulk={L}: r90 = {r90(r['C'])}")
