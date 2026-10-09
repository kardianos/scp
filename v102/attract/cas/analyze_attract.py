#!/usr/bin/env python3
"""Analysis for ../PREREG.md: B1 (the moving W hole), A1 (ladder from
transport), A2 (ladder from the all-order terms), A3 (two holes), A4 (W's
(3/2,1/2) channel).  Energies in J0 (TeNPy, leading order) or in the bulk J
(A2, all orders)."""
import re, glob
import numpy as np
here = __file__.rsplit('/', 1)[0]
PERM = here + '/../../perm/cas/out'

def read(f):
    try:
        L = [x for x in open(f).read().splitlines() if x.strip()]
    except FileNotFoundError:
        return None
    if not L or 'E/J' not in L[0]: return None
    out = {'E': float(re.search(r'E/J0?=([-0-9.]+)', L[0]).group(1))}
    m = re.search(r'trunc=([-0-9.e+]+)', L[0])
    if m: out['trunc'] = float(m.group(1))
    for line in L[1:]:
        k, *v = line.split()
        out[k] = np.array([float(x) for x in v])
    return out

def tp(M, mode, L, tg, t, chi):
    r = read(f"{here}/out/tp_M{M}_{mode}_{L}_{tg}_t{t}_chi{chi}.out")
    return r if r else read(f"{PERM}/tp_M{M}_{mode}_{L}_{tg}_t{t}_chi{chi}.out")

def decay(a, lo, hi):
    ds = np.arange(lo, min(hi, len(a) - 1) + 1); ds = ds[np.abs(a[ds]) > 1e-12]
    if len(ds) < 3: return float('nan')
    return -1/np.polyfit(ds, np.log(np.abs(a[ds])), 1)[0]

def extrap(Ls, ms):
    Ls = np.array(Ls, float); ms = np.array(ms)
    if len(Ls) < 2: return float('nan')
    A = np.stack([np.ones_like(Ls), 1/(Ls + 1)**2], 1)
    return np.linalg.lstsq(A, ms, rcond=None)[0][0]

def r90(prof):
    c = np.cumsum(prof)
    for d in range(len(c)):
        if all(abs(c[k]) >= 0.45 for k in range(d, len(c))): return d
    return float('nan')

gapW = {}
for L in (16, 24, 32):
    a, b = tp(4, 'nohole', L, '0,0', 0, 128), tp(4, 'nohole', L, '2,2', 0, 128)
    if a and b: gapW[L] = b['E'] - a['E']
gap2 = {}
for L in (24, 40):
    a, b = tp(2, 'nohole', L, '0', 0, 64), tp(2, 'nohole', L, '2', 0, 64)
    if a and b: gap2[L] = b['E'] - a['E']
print(f"# gaps (J0): W {gapW};  3+2 {gap2}")

print("\n# B1 / A1 (W): the moving W hole")
for t in ('0.1', '1.0', '3.0'):
    lam = {}; marg = {}
    for L in (16, 24):
        r = tp(4, 'mob', L, '1,1', t, 128)
        if not r: continue
        C = r['C']; nh = r['nh']; n = len(nh)
        w = np.array([sum(nh[i]*((i + d < n) + (i - d >= 0 and d > 0)) for i in range(n)) for d in range(n)])
        c = C/np.maximum(w, 1e-12)
        lam[L] = decay(C, 2, (L + 1)//2 - 3)
        line = f"  t={t} L_bulk={L}: trunc {r.get('trunc', float('nan')):.1e}; C(0..4) " + " ".join(f"{x:+.3f}" for x in C[:5]) + \
               f"; decay (registered) {lam[L]:.3f}; conditional decay {decay(c, 2, (L + 1)//2 - 3):.3f}; r90 {r90(C)}"
        for tg in ('3,3', '3,1'):
            r3 = tp(4, 'mob', L, tg, t, 128)
            if r3 and L in gapW:
                m = r3['E'] - r['E'] - gapW[L]; marg.setdefault(tg, {})[L] = m
                line += f";  E({tg})-E(1,1)-gap = {m:+.4f}"
        print(line)
    if len(lam) == 2:
        x, y = lam[16], lam[24]
        print(f"   N2 (moving W hole, t={t}): {x:.3f} -> {y:.3f}  change {abs(y-x)/x*100:.1f}%  ({'PASS' if abs(y-x)/x < 0.10 else 'MISS'})")
    for tg, d in marg.items():
        if len(d) == 2:
            ex = extrap(list(d), list(d.values()))
            print(f"   ladder ({tg}) at t={t}: margins {d}; extrapolated {ex:+.4f}  -> {'BOUND' if ex < -0.01 and d[max(d)] < 0 else 'not bound'}")

print("\n# A1 (3+2): ladder from transport")
for t in ('0.1', '1.0', '3.0', '5.0', '7.0'):
    d = {}; gs = {}
    for L in (24, 40):
        r1, r3 = tp(2, 'mob', L, '1', t, 64), tp(2, 'mob', L, '3', t, 64)
        if r1 and r3 and L in gap2: d[L] = r3['E'] - r1['E'] - gap2[L]
    if d:
        ex = extrap(list(d), list(d.values())) if len(d) == 2 else float('nan')
        print(f"  t={t}: margins {{{', '.join(f'{k}: {v:+.4f}' for k, v in d.items())}}}; extrapolated {ex:+.4f} -> "
              f"{'BOUND' if ex < -0.01 and d[max(d)] < 0 else 'not bound'}")

print("\n# A2 (3+2, all orders): ladder from the fourth-order terms (energies in the bulk J)")
for g in ('0.1', '0.2', '0.3'):
    d = {}
    for L in (24, 40):
        f = lambda mode, tw: read(f"{here}/out/a2_g{g}_{mode}_{L}_{tw}.out")
        h1, h3, n0, n2 = f('hole', 1), f('hole', 3), f('nohole', 0), f('nohole', 2)
        if g == '0.1' and L == 40:
            h1 = h1 or read(here + '/../../fabric/cas/out/f4b_40_1.out')
            h3 = h3 or read(here + '/../../fabric/cas/out/f4b_40_3.out')
        if h1 and h3 and n0 and n2:
            d[L] = (h3['E'] - h1['E']) - (n2['E'] - n0['E'])
    if d:
        ex = extrap(list(d), list(d.values())) if len(d) == 2 else float('nan')
        print(f"  g/Delta={g}: margins {{{', '.join(f'{k}: {v:+.4f}' for k, v in d.items())}}}; extrapolated {ex:+.4f} -> "
              f"{'BOUND' if ex < -0.01 and d[max(d)] < 0 else 'not bound'}")

print("\n# A3 (3+2): two holes")
far = {tg: tp(2, 'pair10', 24, tg, 0, 64) for tg in ('0', '2')}
for d in (2, 3, 4, 5, 6):
    line = f"  static pair d={d}:"
    for tg in ('0', '2'):
        r = tp(2, f'pair{d}', 24, tg, 0, 64)
        if r and far[tg]: line += f"  S^z={int(tg)//2}: E(d)-E(10) = {r['E']-far[tg]['E']:+.4f} J0"
    print(line)
for t in ('0.1', '1.0'):
    for L in (24, 40):
        e0 = tp(2, 'nohole', L + 2, '0', 0, 64); e1 = tp(2, 'mob', L + 1, '1', t, 64)
        for tg in ('0', '2'):
            e2 = tp(2, 'mob2', L, tg, t, 64)
            if e0 and e1 and e2:
                B = e2['E'] + e0['E'] - 2*e1['E']
                Pr = e2.get('Pr'); mean = float((np.arange(len(Pr))*Pr).sum()/max(Pr.sum(), 1e-12)) if Pr is not None else float('nan')
                pk = int(np.argmax(Pr)) if Pr is not None else -1
                print(f"  mobile pair t={t} L_bulk={L} 2S^z={tg}: B = {B:+.4f} J0;  pair distance: mean {mean:.2f}, peak {pk} "
                      f"(P(peak) {Pr[pk]:.3f})")

print("\n# A3 (W): two holes, static pairs, sectors (2S^z, 2I^z)")
farW = {tg: tp(4, 'pair10', 24, tg, 0, 128) for tg in ('0,0', '2,0', '2,2')}
for d in (2, 3, 4, 6):
    line = f"  d={d}:"
    for tg in ('0,0', '2,0', '2,2'):
        r = tp(4, f'pair{d}', 24, tg, 0, 128)
        if r and farW[tg]: line += f"  {tg}: E(d)-E(10) = {r['E']-farW[tg]['E']:+.4f}"
        elif r: line += f"  {tg}: E = {r['E']:.4f}"
    print(line)
print("  (a level in (0,0) absent from (2,0) and (2,2) is (S,I) = (0,0); in (2,0) absent from (2,2) is (1,0); (2,2) is (1,1))")
for L in (16, 24):
    e0 = tp(4, 'nohole', L + 2, '0,0', 0, 128); e1 = tp(4, 'mob', L + 1, '1,1', '0.1', 128)
    for tg in ('0,0', '2,0', '2,2'):
        e2 = tp(4, 'mob2', L, tg, '0.1', 128)
        if e0 and e1 and e2:
            Pr = e2.get('Pr'); mean = float((np.arange(len(Pr))*Pr).sum()/max(Pr.sum(), 1e-12)) if Pr is not None else float('nan')
            print(f"  W mobile pair t=0.1 L_bulk={L} sector {tg}: B = {e2['E']+e0['E']-2*e1['E']:+.4f} J0; mean pair distance {mean:.2f}")

print("\n# A4 (W): (3/2,1/2) against (3/2,3/2), static hole")
for L in (24, 32):
    a, b, c = tp(4, 'hole', L, '1,1', 0, 128), tp(4, 'hole', L, '3,1', 0, 128), tp(4, 'hole', L, '3,3', 0, 128)
    if a and b and c and L in gapW:
        print(f"  L_bulk={L}: E(3,1)-E(1,1) = {b['E']-a['E']:.5f}; E(3,3)-E(1,1) = {c['E']-a['E']:.5f}; gap {gapW[L]:.5f}")
