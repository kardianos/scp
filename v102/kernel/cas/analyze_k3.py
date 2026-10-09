#!/usr/bin/env python3
"""K3 analysis (../K3_PREREG.md): kernel bars per run, the replay with projection (R12), and the
registered readouts. Part A: A_M{2,4}_{d,b}_s{1..5}.out (light). Part B: B1mid_M*.out, B1end_M*.out,
B2_M*.out (seeded, no light).
usage: analyze_k3.py dir"""
import sys, re, glob, subprocess, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1]

def run(f):
    txt = open(f).read()
    r = re.search(r'# RESULT lockq events=(\d+) Nh=(\d+) window_viol=(\d+)', txt)
    if not r: return None
    cv = re.search(r'check_violations=(\d+)', txt); dr = re.search(r'# RESULT drift_rel (\S+)', txt)
    sw = re.search(r'# RESULT lockq_proj max_sweeps=(\d+)', txt)
    n = int(re.search(r'K2a: M=\d+ n=(\d+)', txt).group(1))
    rep = subprocess.run([sys.executable, '-I', here + '/k2d_replay.py', f, 'score', 'proj'], capture_output=True, text=True).stdout
    ab = re.search(r'\|<H>/J0 diff\| (\S+), max\|nh diff\| (\S+), max\|sz diff\| (\S+)', rep)
    fin = re.search(r'final: N_h=(\d+), nh = ([0-9. ]+); adjacency (\S+); mean B (\S+) J0(?:; bound weight \(B_k < -0.05\) (\S+))?', rep)
    aft = re.search(r'after the last event \(t=(\S+)\): nh = ([0-9. ]+); profile change to tL, half sum \|nh\(tL\) - nh\|: (\S+)', rep)
    pair = re.search(r'pair: <d> of the bound part (\S+) -> (\S+) \(weight (\S+) -> (\S+)\); <d> of the unbound part (\S+) -> (\S+)', rep)
    d = dict(n=n, ev=int(r.group(1)), Nh=int(r.group(2)), wv=int(r.group(3)), cv=int(cv.group(1)) if cv else -1,
             drift=float(dr.group(1)) if dr else float('nan'), sweeps=int(sw.group(1)) if sw else -1,
             ab=max(float(x) for x in ab.groups()) if ab else float('nan'), rep=rep)
    d['nh'] = np.array([float(x) for x in fin.group(2).split()]) if fin else np.zeros(n)
    d['adj'] = float(fin.group(3)) if fin else 0.0
    d['B'] = float(fin.group(4)) if fin else 0.0
    d['bw'] = float(fin.group(5)) if fin and fin.group(5) else float('nan')
    d['nh_e'] = np.array([float(x) for x in aft.group(2).split()]) if aft else d['nh']
    d['Dw'] = float(aft.group(3)) if aft else 0.0
    d['pair'] = tuple(float(x) for x in pair.groups()) if pair else None
    # events, with the N = 1 neighbours each had (holes so far, or an end cell)
    holes, evs = set(), []
    for t, site, nu, de in re.findall(r'# LQ event t=(\S+) site=(\d+) .*?nu=(\d+) .*?delta_exact=(\S+)', txt):
        s = int(site); k = sum(1 for q in (s - 1, s + 1) if q in holes or q in (0, n - 1))
        evs.append((float(t), s, int(nu), float(de), k)); holes.add(s)
    d['evs'] = evs
    hs = [i for i in range(n) if d['nh'][i] > 0.5]
    blocks, cur = [], []
    for i in hs:
        if cur and i == cur[-1] + 1: cur.append(i)
        else:
            if cur: blocks.append(cur)
            cur = [i]
    if cur: blocks.append(cur)
    d['sites'], d['blocks'] = hs, [len(b) for b in blocks if len(b) >= 2]
    d['ok'] = abs(d['drift']) <= 1e-12 and d['cv'] == 0 and d['wv'] == 0 and d['ab'] <= 1e-8
    return d

def cls(d):
    if not np.isnan(d['bw']) and d['bw'] >= 0.10: return 'nucleus'
    if d['adj'] >= 0.5 and d['B'] > 0: return 'block'
    if d['Nh'] >= 1: return 'singles'
    return 'empty'

def prof(nh):
    w = nh.sum()
    if w <= 0: return float('nan'), float('nan')
    x = np.arange(len(nh)); m = (x*nh).sum()/w
    return m, np.sqrt(((x - m)**2*nh).sum()/w)

bad = 0
print("=== Part A: light (mid placement) ===")
print("run          ev Nh wv cv  drift     A/B     sweeps  adj    B/J0     bound   D      sites            blocks  class")
summ = {}
for f in sorted(glob.glob(D + '/A_M*_*_s*.out')):
    tag = os.path.basename(f)[:-4]
    M, light, sd = re.match(r'A_M(\d)_(\w)_s(\d)', tag).groups()
    d = run(f)
    if d is None: print(f"{tag:12s} MISSING RESULT"); bad += 1; continue
    bad += not d['ok']; c = cls(d)
    summ.setdefault((M, light), []).append((d, c))
    print(f"{tag:12s} {d['ev']:2d} {d['Nh']:2d} {d['wv']:2d} {d['cv']:2d}  {d['drift']:+.1e} {d['ab']:.1e}  {d['sweeps']:4d}  {d['adj']:5.2f} {d['B']:+8.3f}  {d['bw']:6.3f}  {d['Dw']:.3f}  {str(d['sites']):16s} {str(d['blocks']):7s} {c}" + ("" if d['ok'] else "   <-- BAR"))
print("\n  events (t, site, mode, delta, N=1 neighbours):")
for (M, light), rows in sorted(summ.items()):
    for k, (d, c) in enumerate(rows):
        print(f"  A_M{M}_{light}_#{k + 1}: " + " ".join(f"({t:.2f},{s},{nu},{de:+.3f},{nb})" for t, s, nu, de, nb in d['evs']))
print("\n  Q1/Q4 per configuration:")
for (M, light), rows in sorted(summ.items()):
    cnt = {k: sum(1 for d, c in rows if c == k) for k in ('nucleus', 'block', 'singles', 'empty')}
    nhs = [d['Nh'] for d, c in rows]; bls = [tuple(d['blocks']) for d, c in rows]
    forced = len(set(nhs)) == 1 and len(set(bls)) == 1
    print(f"  {'W' if M == '4' else '3+2'} {'bright' if light == 'b' else 'dim':6s}: classes {cnt}; N_h {nhs}; blocks {[list(b) for b in bls]}; "
          f"bound weights {['%.3f' % d['bw'] for d, c in rows]}; mean B {['%+.2f' % d['B'] for d, c in rows]}; forces its build: {forced}")
print("\n=== Part B: seeded (no light) ===")
for f in sorted(glob.glob(D + '/B*_M*.out')):
    tag = os.path.basename(f)[:-4]
    d = run(f)
    if d is None: print(f"{tag:10s} MISSING RESULT"); bad += 1; continue
    bad += not d['ok']
    m0, s0 = prof(d['nh_e']); m1, s1 = prof(d['nh'])
    line = (f"{tag:10s} drift {d['drift']:+.1e} A/B {d['ab']:.1e} cv {d['cv']} wv {d['wv']} sweeps {d['sweeps']}: N_h {d['Nh']}; "
            f"<x> {m0:.3f} -> {m1:.3f}, spread {s0:.3f} -> {s1:.3f}; D {d['Dw']:.4f}")
    if tag.startswith('B1'):
        line += "  => " + ('walks' if d['Dw'] >= 0.5 else 'pinned' if d['Dw'] < 0.1 else 'slow')
    if tag.startswith('B2'):
        line += f"; bound weight {d['bw']:.4f}; mean B {d['B']:+.3f} J0"
        if d['pair']:
            b0, b1, w0, w1, u0, u1 = d['pair']
            holds = w1 >= 0.10 and abs(b1 - b0) < 0.15*b0
            line += f"; bound part <d> {b0:.3f} -> {b1:.3f} (weight {w0:.4f} -> {w1:.4f}), unbound part <d> {u0:.3f} -> {u1:.3f}  => pair {'holds' if holds else 'does not hold'}"
    print(line + ("" if d['ok'] else "   <-- BAR"))
    print("    final nh: " + " ".join(f"{x:.3f}" for x in d['nh']))
print(f"\nbars (drift <= 1e-12, check 0, window 0, replay A/B with projection <= 1e-8): {'ALL PASS' if bad == 0 else f'{bad} FAIL'}")
