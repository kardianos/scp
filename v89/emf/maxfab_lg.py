#!/usr/bin/env python3
"""maxfab-lg — CAN EM5 AVOID DELAUNAY? The link-graph complex, weighted.

v1 (unweighted link-graph triangles) gave a flat resonance; the DEC on
Delaunay+barycentric gave a cone (dec_ladder_*: bar met). This probe
asks whether the production kernel could build EM5 faces from its OWN
link rule (triangles of linked triples) with cheap diagonal weights:

    *1_e = (sum_{f∋e} A_f/3) / len(e)         (lumped dual area)
    *2_f = (dbar / kappa) / A_f               (lumped dual length,
           kappa = (faces/edge)/4 — multiplicity correction vs the
           cubic-lattice reference where 4 plaquettes/edge gives c=C)

PRE-REGISTERED READING (asymmetric by design): a CONE (ladder fit,
same bar as dec) is DECISIVE — production EM5 can use link-graph
faces, no Delaunay needed. A flat/band result is NOT decisive against
the complex (could blame the heuristic weights); the fallback is the
proven in-kernel Delaunay route.

usage: maxfab_lg.py foam.tsv [kx ...]
"""
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree

FOAM = sys.argv[1] if len(sys.argv) > 1 else "foam.tsv"
KXS = [float(a) for a in sys.argv[2:]] or \
      [0.261799, 0.523599, 0.785398, 1.047198, 1.308997, 1.570796,
       1.832596, 2.094395]

d = np.loadtxt(FOAM, skiprows=1)
P, R = d[:, 1:4], d[:, 4]
NC = len(P)

# links: the production rule d < 1.15*(ri+rj)
pairs = cKDTree(P).query_pairs(r=2.1, output_type="ndarray")
dd = np.linalg.norm(P[pairs[:, 0]] - P[pairs[:, 1]], axis=1)
E = pairs[dd < 1.15 * (R[pairs[:, 0]] + R[pairs[:, 1]])]
E = E[np.lexsort((E[:, 1], E[:, 0]))]
NE = len(E)
edge_id = {(int(a), int(b)): i for i, (a, b) in enumerate(E)}
adj = [set() for _ in range(NC)]
for a, b in E:
    adj[a].add(int(b)); adj[b].add(int(a))

# faces: all linked triples i<j<k
tris = []
for a, b in E:
    for k in sorted(adj[a] & adj[b]):
        if k > b:
            tris.append((int(a), int(b), int(k)))
faces = np.array(tris)
NF = len(faces)
print(f"# maxfab_lg: NC={NC} NE={NE} NF={NF} (faces/edge={3*NF/NE:.2f})",
      flush=True)

fa, fb, fc = faces[:, 0], faces[:, 1], faces[:, 2]
e_ab = np.array([edge_id[(x, y)] for x, y in zip(fa, fb)])
e_bc = np.array([edge_id[(x, y)] for x, y in zip(fb, fc)])
e_ac = np.array([edge_id[(x, y)] for x, y in zip(fa, fc)])
rows = np.repeat(np.arange(NF), 3)
cols = np.stack([e_ab, e_bc, e_ac], axis=1).ravel()
vals = np.tile([1.0, 1.0, -1.0], NF)
d1 = coo_matrix((vals, (rows, cols)), shape=(NF, NE)).tocsr()
rows = np.repeat(np.arange(NE), 2)
d0 = coo_matrix((np.tile([-1.0, 1.0], NE), (rows, E.ravel())),
                shape=(NE, NC)).tocsr()

emid = 0.5 * (P[E[:, 0]] + P[E[:, 1]])
evec = P[E[:, 1]] - P[E[:, 0]]
elen = np.linalg.norm(evec, axis=1)
eu = evec / elen[:, None]
farea = 0.5 * np.linalg.norm(np.cross(P[fb] - P[fa], P[fc] - P[fa]), axis=1)
dbar = float(np.mean(elen))
kappa = (3.0 * NF / NE) / 4.0

# lumped weights
adual = np.zeros(NE)
for esel in (e_ab, e_bc, e_ac):
    np.add.at(adual, esel, farea / 3.0)
star1 = adual / elen
star2 = (dbar / kappa) / farea
inv1 = 1.0 / star1

lo, hi = P.min(axis=0), P.max(axis=0)
cen = 0.5 * (lo + hi)
rint = 0.5 * float(np.min(hi - lo)) - 2.0
rr = np.linalg.norm(emid - cen, axis=1)
env = np.clip((rint - rr) / 2.0, 0.0, 1.0)
env = 0.5 - 0.5 * np.cos(np.pi * env)

v = np.random.default_rng(1).standard_normal(NE)
for _ in range(60):
    v = inv1 * (d1.T @ (star2 * (d1 @ v)))
    v /= np.linalg.norm(v)
lam = v @ (inv1 * (d1.T @ (star2 * (d1 @ v))))
dt = min(0.02, 1.5 / np.sqrt(lam))
print(f"# lambda_max={lam:.2f} dt={dt:.4f}", flush=True)

div_op = (d0.T @ diags(star1) @ d0).tocsr()
dg = div_op.diagonal().copy()
dg[dg <= 0] = 1.0
Mpre = diags(1.0 / dg)
results = []
for kx in KXS:
    T = 180.0 if kx <= 0.55 else 120.0
    Ef = env * np.cos(kx * emid[:, 0]) * eu[:, 1]
    B = np.zeros(NF)
    rhs = d0.T @ (star1 * Ef)
    phi, info = cg(div_op, rhs, tol=1e-10, maxiter=4000, M=Mpre)
    Ef -= d0 @ phi
    pw = star1 * env * eu[:, 1] * np.exp(-1j * kx * emid[:, 0])
    ns = int(T / dt)
    rec = np.empty(ns, dtype=complex)
    H0 = 0.5 * (Ef @ (star1 * Ef) + B @ (star2 * B))
    for s in range(ns):
        rec[s] = pw @ Ef
        B -= dt * (d1 @ Ef)
        Ef += dt * (inv1 * (d1.T @ (star2 * B)))
    H1 = 0.5 * (Ef @ (star1 * Ef) + B @ (star2 * B))
    w = np.hanning(ns)
    ts = np.arange(ns) * dt
    grid = np.arange(0.02, 3.5, 0.005)
    amp2 = np.abs(np.exp(1j * np.outer(grid, ts)) @ (rec * w)) ** 2
    wpk = grid[int(np.argmax(amp2))]
    print(f"# RESULT lg kx={kx:.2f} omega_peak={wpk:.3f} "
          f"Hdrift={(H1 - H0) / H0:.2e}", flush=True)
    results.append((kx, wpk))

ks = np.array([k for k, _ in results if k <= 1.05])
ws = np.array([w_ for k, w_ in results if k <= 1.05])
if len(ks) >= 3:
    m, bfit = np.polyfit(ks, ws, 1)
    pred = m * ks + bfit
    sst = np.sum((ws - ws.mean()) ** 2)
    r2 = 1 - np.sum((ws - pred) ** 2) / sst if sst > 0 else 1.0
    ok = abs(bfit) <= 0.1 * abs(m) and r2 >= 0.99
    print(f"# CONE fit (k<=1): omega = {m:.3f}*k {bfit:+.3f}, R2={r2:.4f} "
          f"-> {'CONE — link-graph complex viable' if ok else 'not a cone here (weights not exonerated; Delaunay route stands)'}")
