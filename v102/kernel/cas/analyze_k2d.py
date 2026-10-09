#!/usr/bin/env python3
"""K2d D3 analysis: kernel bars per run, the replay A/B, and the readouts.
usage: analyze_k2d.py dir_of_d3_outputs  (files M{2,4}_{b,d}_p{0..4}_s{1,2,3}.out)"""
import sys, re, glob, subprocess, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
NAMES = ['zero', 'core', 'score', 'upper', 'lower']
rows = []
for f in sorted(glob.glob(sys.argv[1] + '/M*_p*_s*.out')):
    tag = os.path.basename(f)[:-4]
    M, light, p, sd = re.match(r'M(\d)_(\w)_p(\d)_s(\d)', tag).groups()
    txt = open(f).read()
    r = re.search(r'# RESULT lockq events=(\d+) Nh=(\d+) window_viol=(\d+) max\|delta_exact\|=(\S+)', txt)
    cv = re.search(r'check_violations=(\d+)', txt)
    dr = re.search(r'# RESULT drift_rel (\S+)', txt)
    if not r: rows.append((tag, None)); continue
    rep = subprocess.run([sys.executable, '-I', here + '/k2d_replay.py', f, NAMES[int(p)]], capture_output=True, text=True).stdout
    ab = re.search(r'\|<H>/J0 diff\| (\S+), max\|nh diff\| (\S+), max\|sz diff\| (\S+)', rep)
    nb = re.search(r'events \(site, nu, hole density on its neighbours\): (.*); next to an end: (\d+)', rep)
    fin = re.search(r'final: N_h=(\d+), nh = ([0-9. ]+); adjacency (\S+); mean B (\S+) J0(?:; bound weight \(B_k < -0.05\) (\S+))?', rep)
    nbr = [float(x.split(',')[2].rstrip(')')) for x in nb.group(1).split()] if nb and nb.group(1).strip() else []
    rows.append((tag, dict(M=int(M), light=light, p=int(p), sd=int(sd), ev=int(r.group(1)), Nh=int(r.group(2)), wv=int(r.group(3)),
                           dmax=float(r.group(4)), cv=int(cv.group(1)) if cv else -1, drift=float(dr.group(1)) if dr else float('nan'),
                           ab=max(float(x) for x in ab.groups()) if ab else float('nan'),
                           adjform=sum(1 for x in nbr if x > 0.5), endform=int(nb.group(2)) if nb else 0,
                           nh=fin.group(2).strip() if fin else '', adj=float(fin.group(3)) if fin else 0.0,
                           B=float(fin.group(4)) if fin else 0.0, bw=float(fin.group(5)) if fin and fin.group(5) else float('nan'))))
bad = 0
print("run            ev Nh  wv cv  drift     A/B     adjform endform  adj    B/J0     bound  final nh")
for tag, d in rows:
    if d is None: print(f"{tag:14s} MISSING RESULT"); bad += 1; continue
    ok = abs(d['drift']) <= 1e-12 and d['cv'] == 0 and d['wv'] == 0 and d['ab'] <= 1e-8
    bad += not ok
    print(f"{tag:14s} {d['ev']:2d} {d['Nh']:2d}  {d['wv']}  {d['cv']}  {d['drift']:+.1e} {d['ab']:.1e}  {d['adjform']:3d}    {d['endform']:3d}   {d['adj']:5.2f} {d['B']:+8.3f}  {d['bw']:6.3f}  {d['nh']}" + ("" if ok else "   <-- BAR"))
print(f"\nbars (drift <= 1e-12, check 0, window 0, replay A/B <= 1e-8): {'ALL PASS' if bad == 0 else f'{bad} FAIL'}")
print("\nsummary per fabric x light x pricing (means over seeds): events, final N_h, events next to a hole, final adjacency, bound weight")
for M in (2, 4):
    for light in ('d', 'b'):
        for p in range(5):
            ds = [d for t, d in rows if d and d['M'] == M and d['light'] == light and d['p'] == p]
            if not ds: continue
            bw = [d['bw'] for d in ds if not np.isnan(d['bw'])]
            print(f"  {'W' if M == 4 else '3+2'} {'bright' if light == 'b' else 'dim':6s} {NAMES[p]:5s}: events {np.mean([d['ev'] for d in ds]):.1f}, N_h {[d['Nh'] for d in ds]}, next to a hole {sum(d['adjform'] for d in ds)}, adjacency {np.mean([d['adj'] for d in ds]):.2f}, bound weight {('%.3f' % np.mean(bw)) if bw else '-'} {['%.2f' % x for x in bw]}")
