#!/usr/bin/env python3
"""bind M2/M3b/M3c analysis, by the bars of ../PREREG.md §4."""
import glob, re, sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
OUT = here + '/out'
BAR = {4: dict(bound=-0.10, free=-0.03, W1=0.388, gap=1.32, lens=(18, 26), Nmax=6),
       2: dict(bound=-0.05, free=-0.02, W1=0.112, gap=0.44, lens=(26, 42), Nmax=6)}
PR = ['zero', 'core', 'score', 'upper', 'lower']
def load(M, nb, nh, tg, p):
    f = f"{OUT}/m2_M{M}_n{nb}_h{nh}_{tg}_{p}.out"
    if not os.path.exists(f) or os.path.getsize(f) == 0: return None
    L = open(f).read().split('\n')
    m = re.search(r'E/J0=(\S+).*trunc=(\S+)', L[0])
    if not m: return None
    r = dict(E=float(m.group(1)), trunc=float(m.group(2)))
    for ln in L[1:]:
        if ln.startswith('Pr '): r['Pr'] = np.array([float(x) for x in ln.split()[1:]])
        if ln.startswith('nh '): r['nh'] = np.array([float(x) for x in ln.split()[1:]])
    return r
def sec(M, nh): return ('0,0' if nh % 2 == 0 else '1,1') if M == 4 else ('0' if nh % 2 == 0 else '1')
def meand(r, nh):
    if nh < 2 or 'Pr' not in r: return float('nan')
    Pr = r['Pr']; d = np.arange(len(Pr)); return float((d*Pr).sum()/Pr.sum())
def adj(r, nh):
    if nh < 2 or 'Pr' not in r: return float('nan')
    return float(r['Pr'][1]/(nh*(nh - 1)/2))
for M in (2, 4):
    b = BAR[M]; name = 'W' if M == 4 else '3+2'
    print(f"\n=== {name}: bars bound b <= {b['bound']}, free b >= {b['free']} (longer length); W1 = {b['W1']}, gap = {b['gap']} J0 ===")
    E0 = {nb: (load(M, nb, 0, sec(M, 0), 'zero') or {}).get('E') for nb in b['lens']}
    res = {}
    for p in PR:
        for nb in b['lens']:
            if E0[nb] is None: continue
            rows = {}
            r1 = load(M, nb, 1, sec(M, 1), p)
            if r1 is None: continue
            e1 = r1['E'] - E0[nb]
            for nh in range(1, b['Nmax'] + 1):
                r = load(M, nb, nh, sec(M, nh), p)
                if r is None: continue
                B = r['E'] - E0[nb] - nh*e1
                rows[nh] = dict(B=B, md=meand(r, nh), adj=adj(r, nh), tr=r['trunc'])
            for nh in rows:
                rows[nh]['b'] = rows[nh]['B'] - rows[nh - 1]['B'] if nh - 1 in rows else (0.0 if nh == 1 else float('nan'))
            res[(p, nb)] = rows
            if rows:
                print(f"  {p:5s} n_bulk={nb}: " + "  ".join(f"N={nh}: B={v['B']:+.3f} b={v['b']:+.3f} <d>={v['md']:.2f} adj={v['adj']:.2f} tr={v['tr']:.0e}" for nh, v in sorted(rows.items()) if nh >= 2))
    # verdicts
    print(f"  --- verdicts ({name}) ---")
    for p in PR:
        la, lb = b['lens']
        ra, rb = res.get((p, la), {}), res.get((p, lb), {})
        if not ra: continue
        two = bool(rb)
        st = []
        for nh in range(2, b['Nmax'] + 1):
            if nh not in ra: break
            ba = ra[nh]['b']; bb = rb[nh]['b'] if two and nh in rb else None
            if bb is None:
                v = 'bound*' if ba <= b['bound'] else ('free*' if ba >= b['free'] else 'marg*')
                comp = None
            else:
                if ba <= b['bound'] and bb <= b['bound']: v = 'bound'
                elif bb >= b['free']: v = 'free'
                else: v = 'marginal'
                ma, mb = ra[nh]['md'], rb[nh]['md']
                comp = abs(mb - ma)/ma
            st.append((nh, v, comp))
        nstar = 1
        for nh, v, comp in st:
            if v == 'bound' and comp is not None and comp < 0.15: nstar = nh
            else: break
        forced = None
        if nstar >= 2:
            nxt = [x for x in st if x[0] == nstar + 1]
            if nxt:
                nh, v, comp = nxt[0]
                forced = (v == 'free') or (comp is not None and comp >= 0.15)
        p3 = [(nh, abs(ra[nh]['b']) < b['W1']) for nh, v, c in st if v == 'bound']
        p4 = [(nh, abs(ra[nh]['b']) > b['gap']) for nh, v, c in st if v == 'bound']
        print(f"  {p:5s}: " + ", ".join(f"N={nh} {v}" + (f" (d-change {c:.0%})" if c is not None else "") for nh, v, c in st)
              + f"  => N* = {nstar}" + ("" if nstar < 2 else f", size forced (i): {forced}")
              + f"; P3 open at {[nh for nh, o in p3 if o]}, P4 open at {[nh for nh, o in p4 if o]}")
    # multiplets
    la = b['lens'][0]
    for p in ('zero', 'upper'):
        out = []
        for nh in range(2, b['Nmax'] + 1):
            r = load(M, la, nh, sec(M, nh), p)
            if r is None: continue
            if M == 4:
                alts = ('2,0', '0,2') if nh % 2 == 0 else ('3,1', '1,3')
            else:
                alts = ('2',) if nh % 2 == 0 else ('3',)
            ds = []
            for a in alts:
                ra_ = load(M, la, nh, a, p)
                ds.append(None if ra_ is None else ra_['E'] - r['E'])
            out.append(f"N={nh}: " + " ".join(f"E({a})-E(min)={'?' if x is None else f'{x:+.4f}'}" for a, x in zip(alts, ds)))
        if out: print(f"  multiplets {p} n_bulk={la}: " + "; ".join(out))
