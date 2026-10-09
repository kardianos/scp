#!/usr/bin/env python3
"""K2e E3 analysis: singlet core (lockq_p11=2, lockq_proj=1) in light. Kernel bars per run,
the replay with projection, and D3's readouts plus the final (1,1) singlet weight.
usage: analyze_k2e.py dir_of_e3_outputs  (files M{2,4}_{b,d}_s{1,2,3}.out)"""
import sys, re, glob, subprocess, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(sys.argv[1] + '/M*_*_s*.out')):
    tag = os.path.basename(f)[:-4]
    M, light, sd = re.match(r'M(\d)_(\w)_s(\d)', tag).groups()
    txt = open(f).read()
    r = re.search(r'# RESULT lockq events=(\d+) Nh=(\d+) window_viol=(\d+) max\|delta_exact\|=(\S+)', txt)
    cv = re.search(r'check_violations=(\d+)', txt)
    dr = re.search(r'# RESULT drift_rel (\S+)', txt)
    sw = re.search(r'# RESULT lockq_proj max_sweeps=(\d+)', txt)
    if not r: rows.append((tag, None)); continue
    rep = subprocess.run([sys.executable, '-I', here + '/k2d_replay.py', f, 'score', 'proj'], capture_output=True, text=True).stdout
    ab = re.search(r'\|<H>/J0 diff\| (\S+), max\|nh diff\| (\S+), max\|sz diff\| (\S+)', rep)
    nb = re.search(r'events \(site, nu, hole density on its neighbours\): (.*); next to an end: (\d+)', rep)
    fin = re.search(r'final: N_h=(\d+), nh = ([0-9. ]+); adjacency (\S+); mean B (\S+) J0(?:; bound weight \(B_k < -0.05\) (\S+))?', rep)
    sg = re.search(r'\(1,1\) singlet weight (\S+)', rep)
    nbr = [float(x.split(',')[2].rstrip(')')) for x in nb.group(1).split()] if nb and nb.group(1).strip() else []
    rows.append((tag, dict(M=int(M), light=light, sd=int(sd), ev=int(r.group(1)), Nh=int(r.group(2)), wv=int(r.group(3)),
                           cv=int(cv.group(1)) if cv else -1, drift=float(dr.group(1)) if dr else float('nan'),
                           sweeps=int(sw.group(1)) if sw else -1,
                           ab=max(float(x) for x in ab.groups()) if ab else float('nan'),
                           adjform=sum(1 for x in nbr if x > 0.5), endform=int(nb.group(2)) if nb else 0,
                           nbr=nbr, nh=fin.group(2).strip() if fin else '', adj=float(fin.group(3)) if fin else 0.0,
                           B=float(fin.group(4)) if fin else 0.0, bw=float(fin.group(5)) if fin and fin.group(5) else float('nan'),
                           sg=sg.group(1) if sg else '-', ev_txt=nb.group(1).strip() if nb else '')))
bad = 0
print("run       ev Nh  wv cv  drift     A/B     sweeps adjform endform  adj    B/J0     bound  singlet   final nh")
for tag, d in rows:
    if d is None: print(f"{tag:9s} MISSING RESULT"); bad += 1; continue
    ok = abs(d['drift']) <= 1e-12 and d['cv'] == 0 and d['wv'] == 0 and d['ab'] <= 1e-8
    bad += not ok
    print(f"{tag:9s} {d['ev']:2d} {d['Nh']:2d}  {d['wv']}  {d['cv']}  {d['drift']:+.1e} {d['ab']:.1e}  {d['sweeps']:3d}   {d['adjform']:3d}    {d['endform']:3d}   {d['adj']:5.2f} {d['B']:+8.3f}  {d['bw']:6.3f}  {d['sg']:>8s}  {d['nh']}" + ("" if ok else "   <-- BAR"))
print(f"\nbars (drift <= 1e-12, check 0, window 0, replay A/B with projection <= 1e-8): {'ALL PASS' if bad == 0 else f'{bad} FAIL'}")
print("\nevents per run (site, nu, hole density on its neighbours just before):")
for tag, d in rows:
    if d: print(f"  {tag}: {d['ev_txt']}")
