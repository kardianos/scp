#!/usr/bin/env python3
"""maxfab v2/v3 — DEC Maxwell on a Delaunay complex over the foam cells.

v1 (maxfab.c, all-triangles complex, unweighted) FALSIFIED: flat
resonance ~1.48, no k-dependence — over-complete (~7 faces/edge),
metric-less. This is the fix, per EMF.md 2026-07-31:

  v2: Delaunay tetrahedralization -> true 2-complex (interior faces
      shared by exactly 2 tets) + BARYCENTRIC hodge weights
        *1_e = A_dual(e)/len(e)     (dual polygon area / primal length)
        *2_f = l_dual(f)/A_f        (dual edge length / primal area)
      Faraday (metric-free):  dB_f/dt = -(d1 E)_f
      Ampere:                 dE_e/dt = +(*1^-1 d1^T *2 B)_e
      On a cubic lattice these weights give omega = C k exactly.

  v3: TRANSVERSE PROJECTION of the initial packet — remove the
      gradient (curl-kernel) component by solving the weighted Poisson
      (d0^T *1 d0) phi = d0^T *1 E, E <- E - d0 phi — plus a Hann
      window on the probe DFT. v2's low-k reading was masked by the
      frozen gradient kernel's spectral leakage; after projection the
      whole spectrum is propagating modes.

PRE-REGISTERED CONE BAR (before scoring): omega_peak(k) over the small-k
half (k <= 1.0) fits omega = m*k + b with |b| <= 0.1*m and R^2 >= 0.99
=> the operator class supports a linear cone on THIS foam geometry
(what the production scalar-hop band could not do: EM1 measured
v_g = 0.610*sin(1.410 k)). Verdict feeds the EM5 law design.

usage: maxfab_dec.py foam.tsv [kx ...]
"""
import sys

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg

FOAM = sys.argv[1] if len(sys.argv) > 1 else "foam.tsv"
KXS = [float(a) for a in sys.argv[2:]] or [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.4, 1.7, 2.0]

pts = []
with open(FOAM) as f:
    f.readline()
    for ln in f:
        p = ln.split()
        if len(p) >= 5:
            pts.append((float(p[1]), float(p[2]), float(p[3])))
P = np.array(pts)
NC = len(P)

tri = Delaunay(P)
tets = tri.simplices
# trim hull slivers: cut = 1.6 x median Delaunay edge (geometry-scaled;
# reproduces the hand cut 2.2 on the foam and 2.6 on the annealed set.
# NN-median failed here: dmin contact pairs bias it low on dense foam),
# plus real volume
a, b, c, d = (P[tets[:, i]] for i in range(4))
edge_all = np.concatenate([np.linalg.norm(x - y, axis=1) for x, y in
                           [(a, b), (a, c), (a, d), (b, c), (b, d), (c, d)]])
cut = 1.6 * float(np.median(edge_all))
emax = np.max([np.linalg.norm(x - y, axis=1) for x, y in
               [(a, b), (a, c), (a, d), (b, c), (b, d), (c, d)]], axis=0)
vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6.0
keep = (emax < cut) & (vol > 1e-4)
print(f"# trim: cut={cut:.3f} kept={keep.sum()}/{len(tets)}",
      flush=True)
tets = tets[keep]
NT = len(tets)

edge_id, faces_id = {}, {}
def eid(i, j):
    k = (i, j) if i < j else (j, i)
    if k not in edge_id:
        edge_id[k] = len(edge_id)
    return edge_id[k]
def fid(i, j, k):
    key = tuple(sorted((i, j, k)))
    if key not in faces_id:
        faces_id[key] = len(faces_id)
    return faces_id[key]

tet_faces = np.empty((NT, 4), dtype=np.int64)
for t, (i, j, k, l) in enumerate(tets):
    tet_faces[t] = (fid(j, k, l), fid(i, k, l), fid(i, j, l), fid(i, j, k))
faces = np.array(sorted(faces_id, key=faces_id.get))
NF = len(faces)
for i, j, k in faces:
    eid(i, j); eid(j, k); eid(i, k)
edges = np.array(sorted(edge_id, key=edge_id.get))
NE = len(edges)
print(f"# maxfab_dec: NC={NC} NT={NT} NF={NF} NE={NE} "
      f"(faces/edge={3*NF/NE:.2f})", flush=True)

# d1 (NF x NE): face (a<b<c) boundary = +(a,b) +(b,c) -(a,c)
fa, fb, fc = faces[:, 0], faces[:, 1], faces[:, 2]
e_ab = np.array([edge_id[(x, y)] for x, y in zip(fa, fb)])
e_bc = np.array([edge_id[(x, y)] for x, y in zip(fb, fc)])
e_ac = np.array([edge_id[(x, y)] for x, y in zip(fa, fc)])
rows = np.repeat(np.arange(NF), 3)
cols = np.stack([e_ab, e_bc, e_ac], axis=1).ravel()
vals = np.tile([1.0, 1.0, -1.0], NF)
d1 = coo_matrix((vals, (rows, cols)), shape=(NF, NE)).tocsr()

# d0 (NE x NC): edge (a<b) = phi_b - phi_a
rows = np.repeat(np.arange(NE), 2)
cols = edges.ravel()
vals = np.tile([-1.0, 1.0], NE)
d0 = coo_matrix((vals, (rows, cols)), shape=(NE, NC)).tocsr()

# geometry
emid = 0.5 * (P[edges[:, 0]] + P[edges[:, 1]])
evec = P[edges[:, 1]] - P[edges[:, 0]]
elen = np.linalg.norm(evec, axis=1)
eu = evec / elen[:, None]
fbar = (P[fa] + P[fb] + P[fc]) / 3.0
farea = 0.5 * np.linalg.norm(np.cross(P[fb] - P[fa], P[fc] - P[fa]), axis=1)
tbar = P[tets].mean(axis=1)

# *2_f = (sum over incident tets |bary(T)-bary(f)|) / A_f
ldual = np.zeros(NF)
for col in range(4):
    fi = tet_faces[:, col]
    ldual += np.zeros(NF)  # keep shape
    np.add.at(ldual, fi, np.linalg.norm(tbar - fbar[fi], axis=1))
star2 = ldual / farea

# *1_e = A_dual(e)/len(e); A_dual = sum over (tet, face) flags around e of
# area(mid(e), bary(f), bary(T))
adual = np.zeros(NE)
tf_e = ((e_ab, e_bc, e_ac))
for col in range(4):
    fi = tet_faces[:, col]
    for esel in tf_e:
        ei = esel[fi]
        v1 = fbar[fi] - emid[ei]
        v2 = tbar - emid[ei]
        np.add.at(adual, ei, 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1))
star1 = adual / elen
inv1 = 1.0 / star1

# mode frame: box-wide envelope (smooth radial taper), k-projected probe.
# A localized packet at small kx spans dk ~ 1/sig >> kx and the local-probe
# DFT latches onto the wrong spectral line (measured: non-monotonic
# omega(k) at kx<=0.5). Seeding the whole interior with cos(kx x) and
# recording the complex projection onto that exact spatial mode isolates
# omega(kx) — the standard dispersion instrument.
lo, hi = P.min(axis=0), P.max(axis=0)
cen = 0.5 * (lo + hi)
rint = 0.5 * float(np.min(hi - lo)) - 2.0
rr = np.linalg.norm(emid - cen, axis=1)
env = np.clip((rint - rr) / 2.0, 0.0, 1.0)
env = 0.5 - 0.5 * np.cos(np.pi * env)   # cos taper, 2-unit skin
amp = 1.0

# CFL from power iteration on L = inv1 d1^T star2 d1
v = np.random.default_rng(1).standard_normal(NE)
for _ in range(60):
    v = inv1 * (d1.T @ (star2 * (d1 @ v)))
    v /= np.linalg.norm(v)
lam = v @ (inv1 * (d1.T @ (star2 * (d1 @ v))))
dt = min(0.02, 1.5 / np.sqrt(lam))
print(f"# lambda_max={lam:.2f} dt={dt:.4f}", flush=True)

div_op = (d0.T @ diags(star1) @ d0).tocsr()
dg = div_op.diagonal().copy()
dg[dg <= 0] = 1.0          # isolated nodes: zero rhs, keep phi=0
Mpre = diags(1.0 / dg)
results = []
for kx in KXS:
    T = 180.0 if kx <= 0.55 else 120.0
    E = amp * env * np.cos(kx * emid[:, 0]) * eu[:, 1]
    B = np.zeros(NF)

    # v3: transverse projection (remove weighted-gradient component)
    rhs = d0.T @ (star1 * E)
    r0 = np.linalg.norm(rhs)
    phi, info = cg(div_op, rhs, tol=1e-10, maxiter=4000, M=Mpre)
    E -= d0 @ phi
    r1 = np.linalg.norm(d0.T @ (star1 * E))
    print(f"# kx={kx:.2f} div: {r0:.3e} -> {r1:.3e} (cg info={info})",
          flush=True)

    # complex projection onto the seeded spatial mode (weighted)
    pw = star1 * env * eu[:, 1] * np.exp(-1j * kx * emid[:, 0])
    ns = int(T / dt)
    rec = np.empty(ns, dtype=complex)
    H0 = 0.5 * (E @ (star1 * E) + B @ (star2 * B))
    for s in range(ns):
        rec[s] = pw @ E
        B -= dt * (d1 @ E)
        E += dt * (inv1 * (d1.T @ (star2 * B)))
    H1 = 0.5 * (E @ (star1 * E) + B @ (star2 * B))

    w = np.hanning(ns)
    wr = rec * w
    ts = np.arange(ns) * dt
    grid = np.arange(0.02, 3.5, 0.005)
    amp2 = np.abs(np.exp(1j * np.outer(grid, ts)) @ wr) ** 2
    wpk = grid[int(np.argmax(amp2))]
    print(f"# RESULT dec kx={kx:.2f} omega_peak={wpk:.3f} "
          f"Hdrift={(H1 - H0) / H0:.2e}", flush=True)
    results.append((kx, wpk))

ks = np.array([k for k, _ in results if k <= 1.05])
ws = np.array([w_ for k, w_ in results if k <= 1.05])
if len(ks) >= 3:
    m, bfit = np.polyfit(ks, ws, 1)
    pred = m * ks + bfit
    ssr = np.sum((ws - pred) ** 2)
    sst = np.sum((ws - ws.mean()) ** 2)
    r2 = 1 - ssr / sst if sst > 0 else 1.0
    ok = abs(bfit) <= 0.1 * abs(m) and r2 >= 0.99
    print(f"# CONE fit (k<=1): omega = {m:.3f}*k {bfit:+.3f}, R2={r2:.4f} "
          f"-> {'CONE (bar met)' if ok else 'NOT a cone (bar failed)'}")
    print("# bar: |b| <= 0.1*m and R2 >= 0.99 (pre-registered)")
