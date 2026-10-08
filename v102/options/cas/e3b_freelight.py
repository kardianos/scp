#!/usr/bin/env python3
"""E3b — what each (transport rule x ceiling) does to free light, far from any
lock (added during the campaign; reported under K3, not a registered gate).

An open chain of 24 cells, no lock anywhere.  Two photons, polarization 1
from the left and polarization 2 from the right, Gaussian packets at
kinetic energy eps = 0 (k = pi/2), moving toward each other, t = 0.9.  No
conversion (S3 needs a pair across a link; the photons are on different
cells most of the time; conversion is off here to isolate transport).
Free bosons pass through each other.  Reported: P(pol-1 photon ends right
of the pol-2 photon), P(the two photons are adjacent at the end), and the
largest probability that both photons stay frozen in place.
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
n = 24; t = 0.9
links = [(i, i + 1) for i in range(n - 1)]
xa, xb, sig, k0 = 6, 17, 1.5, np.pi/2
xs = np.arange(n)
pa = np.exp(-(xs - xa)**2/(4*sig**2))*np.exp(-1j*k0*xs); pa /= np.linalg.norm(pa)   # moves right (v = -2t sin k, k = -pi/2)
pb = np.exp(-(xs - xb)**2/(4*sig**2))*np.exp(1j*k0*xs); pb /= np.linalg.norm(pb)    # moves left
print(f"# E3b free light: two photons on an open chain of {n} cells, t = {t}, packets at {xa} (pol 1, rightward) and {xb} (pol 2, leftward)")
for rule in ('T2', 'T4'):
    for ceil in ('lock', 'geo', 'held'):
        kw = dict(ceil='fixed', caps=None) if ceil == 'lock' else dict(ceil=ceil)
        gr = G.Graph(n, links, 0.0, conv=[False]*len(links), t=[t]*len(links), rule=rule, **kw)
        seeds = []
        for x in range(n):
            for y in range(n):
                if x == y: continue
                cells = [(0, 0)]*n; cells[x] = (1, 0); cells[y] = (0, 1)
                seeds.append((tuple(cells), (0,)*len(links)))
        H = gr.build(seeds)
        D = H.shape[0]
        v0 = np.zeros(D, complex)
        for x in range(n):
            for y in range(n):
                if x == y: continue
                cells = [(0, 0)]*n; cells[x] = (1, 0); cells[y] = (0, 1)
                k = gr.idx.get((tuple(cells), (0,)*len(links)))
                if k is not None: v0[k] += pa[x]*pb[y]
        lost = 1 - np.linalg.norm(v0)**2
        v0 /= np.linalg.norm(v0)
        T = (xb - xa)/(2*t)*1.0 + (xb - xa)/(4*t)   # meet at (xb-xa)/(4t), read after they have separated by about the start gap
        vt = expm_multiply(-1j*H, v0, start=0, stop=T, num=3, endpoint=True)[-1]
        P = np.abs(vt)**2
        pos1 = np.zeros(D); pos2 = np.zeros(D)
        for k, (cells, cs) in enumerate(gr.states):
            for i, (a, b) in enumerate(cells):
                if a: pos1[k] += i*a/max(a, 1)
                if b: pos2[k] += i*b/max(b, 1)
        passed = P[pos1 > pos2].sum(); adj = P[np.abs(pos1 - pos2) == 1].sum()
        frozen0 = abs(np.vdot(v0, vt))**2
        print(f"  rule={rule} ceil={ceil:4s}: dim={D:5d}  initial weight on forbidden states {lost:.2e}  "
              f"P(passed)={passed:.4f}  P(adjacent at end)={adj:.4f}  overlap with the start={frozen0:.4f}  [{time.time()-t0:.0f}s]", flush=True)
# adjacent start: do two photons on neighbouring cells move at all?
print("  two photons started on neighbouring cells 11, 12 (pol 1, pol 2), no packet:")
for rule in ('T2', 'T4'):
    for ceil in ('lock', 'geo', 'held'):
        kw = dict(ceil='fixed', caps=None) if ceil == 'lock' else dict(ceil=ceil)
        gr = G.Graph(n, links, 0.0, conv=[False]*len(links), t=[t]*len(links), rule=rule, **kw)
        cells = [(0, 0)]*n; cells[11] = (1, 0); cells[12] = (0, 1)
        s0 = (tuple(cells), (0,)*len(links))
        H = gr.build([s0])
        print(f"    rule={rule} ceil={ceil:4s}: states reachable from the adjacent pair = {H.shape[0]}"
              f"{'  (frozen: the pair never moves)' if H.shape[0] == 1 else ''}")
print(f"done [{time.time()-t0:.0f}s]")
