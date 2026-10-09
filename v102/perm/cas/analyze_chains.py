#!/usr/bin/env python3
"""Analysis of the DMRG chains (x12: W; x3c: 3+2 mobile hole), energies in J0.
K(W): gap(L) = E(2S^z=2) - E(0) without a hole; fit D + a/(L+1)^2.
Hole (W static; 3+2 static and mobile): the decay length of the staggered
spin cloud (static: S^z profile about the middle; mobile: the hole-centred
correlation C(d)); E(higher multiplet) - E(ground) against the no-hole gap
at the same length (N5)."""
import glob, re
import numpy as np
here = __file__.rsplit('/', 1)[0]

def read(f):
    try:
        L = open(f).read().splitlines()
    except FileNotFoundError:
        return None
    if not L: return None
    m = re.search(r'E/J0=([-0-9.]+)', L[0])
    out = {'E': float(m.group(1))}
    for line in L[1:]:
        k, *v = line.split()
        out[k] = np.array([float(x) for x in v])
    return out

def decay(a, dmin, dmax):
    ds = np.arange(dmin, dmax + 1)
    ds = ds[np.abs(a[ds]) > 1e-12]
    if len(ds) < 3: return float('nan')
    sl = np.polyfit(ds, np.log(np.abs(a[ds])), 1)[0]
    return -1/sl

print("# chains (leading order, energies in J0 = 2 g^2/Delta)")
# ---- W
gaps = {}
for L in (8, 12, 16, 20):
    r0 = read(f"{here}/out/x12_W_nohole_{L}_0,0_chi64.out"); r2 = read(f"{here}/out/x12_W_nohole_{L}_2,2_chi64.out")
    if r0 and r2:
        gaps[L] = r2['E'] - r0['E']
        print(f"  W no hole L_bulk={L}: gap (the (1,1) member of the first multiplet) = {gaps[L]:.5f} J0;  target sums {r2['Sz'].sum():+.3f}/{r2['Iz'].sum():+.3f}")
if len(gaps) >= 2:
    Ls = np.array(sorted(gaps)); y = np.array([gaps[L] for L in Ls])
    A = np.stack([np.ones_like(Ls, dtype=float), 1/(Ls + 1.0)**2], 1)
    (D, a), *_ = np.linalg.lstsq(A, y, rcond=None)
    print(f"  W gap extrapolated D + a/(L+1)^2: D = {D:.5f} J0  ({'gapped' if D > 0.1 else 'NOT clearly gapped'})")
for L in (8, 12, 16, 20):
    rh = read(f"{here}/out/x12_W_hole_{L}_1,1_chi64.out")
    if not rh: continue
    p = rh['Sz']; n = len(p); h = n//2
    a = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
    lam = decay(a, 2, h - 3)
    q = rh['Iz']; ai = np.array([q[h]] + [q[h - d] + q[h + d] for d in range(1, h + 1)])
    line = (f"  W hole L_bulk={L}: hole S^z {p[h]:+.4f}, neighbours {a[1]:+.4f}; S^z decay length {lam:.3f} (fit d=2..{h-3}); "
            f"I^z decay length {decay(ai, 2, h - 3):.3f}; max |S^z - I^z| profile {np.abs(p - q).max():.1e}")
    for tg in ('3,3',):
        rx = read(f"{here}/out/x12_W_hole_{L}_{tg}_chi64.out")
        if rx: line += f";  E({tg}) - E(1,1) = {rx['E']-rh['E']:.5f}"
    if L in gaps: line += f";  no-hole gap {gaps[L]:.5f}"
    print(line)
wl = {}
for L in (8, 12, 16, 20):
    rh = read(f"{here}/out/x12_W_hole_{L}_1,1_chi64.out")
    if rh:
        p = rh['Sz']; h = len(p)//2
        wl[L] = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
ks = sorted(wl)
if len(ks) >= 2:
    a, b = ks[-2], ks[-1]
    own = (decay(wl[a], 2, len(wl[a]) - 4), decay(wl[b], 2, len(wl[b]) - 4))
    com = (decay(wl[a], 2, 4), decay(wl[b], 2, 4))
    print(f"  W N2: two longest L_bulk={a},{b}: own fit ranges {own[0]:.3f} -> {own[1]:.3f} ({abs(own[1]-own[0])/own[0]*100:.1f}%); "
          f"common range d=2..4: {com[0]:.3f} -> {com[1]:.3f} ({abs(com[1]-com[0])/com[0]*100:.1f}%)")
# ---- 3+2 mobile hole
print()
g2 = {}
for L in (24, 40):
    r0 = read(f"{here}/out/x3c_M2_nohole_{L}_0_t0.out"); r1 = read(f"{here}/out/x3c_M2_nohole_{L}_2_t0.out")
    if r0 and r1: g2[L] = r1['E'] - r0['E']; print(f"  3+2 no hole L_bulk={L}: gap {g2[L]:.5f} J0")
    rs = read(f"{here}/out/x3c_M2_hole_{L}_1_t0.out"); rs3 = read(f"{here}/out/x3c_M2_hole_{L}_3_t0.out")
    if rs:
        p = rs['Sz']; n = len(p); h = n//2
        a = np.array([p[h]] + [p[h - d] + p[h + d] for d in range(1, h + 1)])
        print(f"  3+2 static hole L_bulk={L}: decay length {decay(a, 2, h - 5):.3f}" +
              (f";  E(3/2)-E(1/2) = {rs3['E']-rs['E']:.5f} (no-hole gap {g2.get(L, float('nan')):.5f})" if rs3 else ""))
    for t in ('0.1', '1.0', '3.0'):
        rm = read(f"{here}/out/x3c_M2_mob_{L}_1_t{t}.out"); rm3 = read(f"{here}/out/x3c_M2_mob_{L}_3_t{t}.out")
        if not rm: continue
        C = rm['C']; lam = decay(C, 2, len(C)//2 - 2)
        line = f"  3+2 mobile hole L_bulk={L} t={t} J0: holes {rm['nh'].sum():.4f}, C(0..4) " + " ".join(f"{x:+.3f}" for x in C[:5]) + f"; decay length of C {lam:.3f}"
        if rm3:
            e = rm3['E'] - rm['E']; gg = g2.get(L, float('nan'))
            line += f";  E(3/2)-E(1/2) = {e:.5f} vs no-hole gap {gg:.5f} -> {'BOUND' if e < gg - 1e-4 else 'above'}"
        print(line)
