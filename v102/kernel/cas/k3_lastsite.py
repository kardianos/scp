#!/usr/bin/env python3
"""K3 check (reported, not gated). usage: k3_lastsite.py kernel.out"""
# For a K3 run: the final lock state (replayed), and every possible next jump (vacant bulk site, mode):
# projected jump weight and the window detuning delta = w1e - Omega - (Ef1 - Ef2) - (<H>_next - <H>) (rate units).
# Check: the same formula on the state just before the last fired event reproduces its printed delta_exact.
import os, sys, runpy, re, numpy as np
from scipy.sparse.linalg import expm_multiply
f = sys.argv[1]
sys.argv = ['k2d_replay.py', f, 'score', 'proj']
G = runpy.run_path(os.path.dirname(os.path.abspath(__file__)) + '/k2d_replay.py', run_name='k2d')
txt = open(f).read()
m = re.search(r'Ef1=(\S+) Ef2=(\S+) Omega=(\S+)', txt); Ef1, Ef2, Om = map(float, m.groups())
w1e = float(re.findall(r'w1e_o=(\S+)', txt)[-1])
locs, P, M, n, strides, dims, sect, project = G['locs'], G['P'], G['M'], G['n'], G['strides'], G['dims'], G['sect'], G['project']
def cand(psi, basis, tgt, site, nu):
    L = locs[site]; dig = (basis//strides[site]) % dims[site]; out = {}
    for kk, c in enumerate(basis):
        if abs(psi[kk]) < 1e-300 or L.N[dig[kk]] != 2: continue
        occ = L.occ[dig[kk]]
        for mu in range(M):
            p = P[mu, nu]
            if p == 0 or occ[mu] == 0: continue
            o2 = list(occ); o2[mu] -= 1; kd = L.idx[tuple(o2)]
            c2 = c + (kd - dig[kk])*strides[site]
            out[c2] = out.get(c2, 0.0) + p*np.sqrt(occ[mu])*psi[kk]
    if not out: return 0.0, 0.0, float('nan')
    mu = [mu for mu in range(M) if P[mu, nu] != 0][0]
    q = np.array([1 if mu == 0 else -1]) if M == 2 else np.array(G['GM'].MODES_W[mu])
    nt = tuple(int(x) for x in np.array(tgt[:-1]) - q) + (int(tgt[-1]) + 1,)
    H, b, _, _ = sect(nt); pos = {c: i for i, c in enumerate(b)}
    x = np.zeros(len(b), complex)
    for c, a in out.items(): x[pos[c]] += a
    w0 = np.vdot(x, x).real; x = project(x, nt); w = np.vdot(x, x).real
    if w < 1e-28: return w0, w, float('nan')
    En = np.vdot(x, H@x).real/w
    return w0, w, En
def scan(psi, basis, tgt, label):
    H = sect(tgt)[0]; Ex = np.vdot(psi, H@psi).real/np.vdot(psi, psi).real
    print(label)
    for s in range(1, n - 1):
        for nu in range(2 if M == 2 else 4):
            w0, w, En = cand(psi, basis, tgt, s, nu)
            if w0 < 1e-6: continue
            d = w1e - Om - (Ef1 - Ef2) - (En - Ex) if w > 1e-28 else float('nan')
            print(f"  site {s} mode {nu}: weight before projection {w0:.4f}, after {w:.3e}, delta {d:+.4f} ({'in' if abs(d) < 0.12 else 'OUT of'} the window)")
# check on the last fired event
ev = G['ev']; te, site, nu = ev[-1]
print(f"{f.rsplit('/',1)[1]}: w1e={w1e} Ef1={Ef1} Ef2={Ef2} Omega={Om}")
psi, tgt, basis, tnow = G['v'][:, 0].astype(complex), G['z0'], G['b0'], 0.0
for (t1, s1, n1) in ev[:-1]:
    if t1 > tnow: psi = expm_multiply(-1j*(t1 - tnow)*sect(tgt)[0], psi); tnow = t1
    psi, tgt, basis = G['jump'](psi, basis, tgt, s1, n1)
if te > tnow: psi = expm_multiply(-1j*(te - tnow)*sect(tgt)[0], psi)
H = sect(tgt)[0]; Ex = np.vdot(psi, H@psi).real
w0, w, En = cand(psi, basis, tgt, site, nu)
printed = float(re.findall(r'delta_exact=(\S+)', txt)[-1])
print(f"  check, last event (site {site} mode {nu}): computed delta {w1e - Om - (Ef1 - Ef2) - (En - Ex):+.6f}, printed {printed:+.6f}")
scan(G['psi'], G['basis'], G['tgt'], "  final state, every next jump:")
