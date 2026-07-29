#!/usr/bin/env python3
"""
v89 MORPHOLOGICAL SEARCH — structural backtracer over the real foam.

Enumerates and searches the topology space of consonant networks on
foam_s20260727.tsv (the simulation-standard foam), scoring every candidate
with a lightweight independent scorer that mirrors the kernel's NETGATE
report (cellfab.c init=net) plus the static-lock theory:

  gate(psi)      = ((1+cos psi)/2)^8                      [gate_of, cellfab.c:489]
  psi_f          = th_i - w_i d / C - th_j                [cellfab.c:1811]
  pitch          w = w2/(1 + q_detune * x), x = Em/cap    [cellfab.c:1715]
  pair ladder    (w_i + w_j) d / C = 2 pi m               [B1 design note]
  link rule      d < 1.15 (r_i + r_j)                     [cellfab.c:701]
  static lock    uniform w per connected component  =>  struts: d = pi/w (m=1)
  cable          forward-tuned, back gate = gate(-2 w d)  [comp12 anchor 0.100]
  cycle closure  w * L_cycle = 2 pi m  for every independent cycle

Laws: laws_V2g.cfg — w2=2.9 q_detune=1.2 cap=2.5 C=1 p_gate=8
      gamma_res_m=0.10 mob_floor=0.004 s_pull=0.15 r0=0.85 rjit=0.06.

Leak model (ranking device; the sim is the judge):
  intended strut  : (1-gf) + (1-gb)          want both gates open
  intended cable  : (1-g_fwd) + W_BACK*g_back  back-gate openness = desert-
                    protected frustration (comp12: gb=0.100, longest-lived)
  parasite        : gf + gb                   off-rung openness radiates
  vacuum port     : res_11(w) * <G>_rand * W_VAC per unlisted foam neighbour
                    (res_11 = gm^2/(gm^2 + (w2-w)^2): the skirt trickle —
                    reproduces the measured lifetime-rises-with-load law)
Phases th are solved by weighted least squares over INTENDED edges only
(locks form on the strong channels; parasites are scored, not optimized
away) from a BFS-tree init that fixes the winding basin.

Outputs: ledger.tsv (every candidate one line), shortlist.json (top 8 + 2
negative controls, concrete cell ids / edges / labels / omega / cycle
integers), nets/*.net seedfiles (frozen V/E format), report.json.

Usage:  python3 search.py [--quick] [--seed N] [--tries N] [--restarts N]
"""

import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
FOAM = os.path.join(HERE, "..", "foam", "foam_s20260727.tsv")

# ---------------------------------------------------------------- laws (V2g)
W2 = 2.9          # base pitch
QDET = 1.2        # q_detune
CAP = 2.5
CC = 1.0          # C
PGATE = 8
GAMMA_M = 0.10    # gamma_res_m (dense sector acceptance)
GAMMA_ROUGH = 0.5
ROUGH_K = 0.35
MOB_FLOOR = 0.004
S_PULL = 0.15
COMB_LIMIT = 6
X_SKIRT = 0.0617  # measured vacuum-skirt death boundary
X_SEED_MIN = 0.25  # campaign floor for seeded load
X_MAX = 0.90
LINK_MARGIN = 1.15

W_OMEGA_MIN = W2 / (1.0 + QDET * X_MAX)       # 1.3942  (x=0.90)
W_OMEGA_MAX = W2 / (1.0 + QDET * X_SEED_MIN)  # 2.2308  (x=0.25)

# scoring weights (documented in MORPHO.md; ranking only)
W_BACK = 0.5
W_PAR = 1.0
W_VAC = 0.11          # sympathetic-floor conductance ratio sqrt(mob_floor*cap/Em)
G_RAND = 0.19638      # <gate>_psi~U = C(16,8)/2^16 (time-avg gate of a beating pair)
PEN_OPEN = 0.25       # per vertex with intended degree < 2 (chain end: MASS 1.7x)
PEN_BRIDGE = 0.12     # per bridge edge (edge carrying no closed cycle)


def gate(psi):
    return (0.5 * (1.0 + np.cos(psi))) ** PGATE


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def omega_of_x(x):
    return W2 / (1.0 + QDET * x)


def x_of_omega(w):
    return (W2 / w - 1.0) / QDET


def res_11(w):
    """1:1 comb resonance of a structure voice against the x~0 room (w2)."""
    det = W2 - w
    return GAMMA_M * GAMMA_M / (GAMMA_M * GAMMA_M + det * det)


# ---------------------------------------------------------------- foam
class Foam:
    def __init__(self, path):
        d = np.loadtxt(path, skiprows=1)
        assert np.all(d[:, 0].astype(int) == np.arange(len(d))), "ids not sequential"
        self.xyz = d[:, 1:4].copy()
        self.r = d[:, 4].copy()
        self.n = len(d)
        self.tree = cKDTree(self.xyz)
        pairs = self.tree.query_pairs(LINK_MARGIN * 2 * self.r.max(),
                                      output_type="ndarray")
        dd = np.linalg.norm(self.xyz[pairs[:, 0]] - self.xyz[pairs[:, 1]], axis=1)
        cut = LINK_MARGIN * (self.r[pairs[:, 0]] + self.r[pairs[:, 1]])
        m = dd < cut
        self.pairs = pairs[m]
        self.pdist = dd[m]
        self.adj = [[] for _ in range(self.n)]
        for (a, b), dist in zip(self.pairs, self.pdist):
            self.adj[a].append((b, dist))
            self.adj[b].append((a, dist))
        self.deg = np.array([len(a) for a in self.adj])

    def linked(self, i, j):
        d = np.linalg.norm(self.xyz[i] - self.xyz[j])
        return d < LINK_MARGIN * (self.r[i] + self.r[j]), d

    def links_within(self, cells):
        """All foam links internal to the cell set: (i, j, d) with i,j cell ids."""
        cs = set(int(c) for c in cells)
        out = []
        for i in cs:
            for j, d in self.adj[i]:
                if j > i and j in cs:
                    out.append((i, j, d))
        return out


# ---------------------------------------------------------------- scorer
class Score:
    __slots__ = ("ok", "omega", "x", "leak_edge", "leak_par", "leak_vac",
                 "per_voice", "labels", "psi", "gates", "npar", "gpar_max",
                 "n", "ne", "theta", "cycles", "nopen", "nbridge", "penalty",
                 "total", "broken", "par_pairs", "wind")


def solve_theta(cells, edges, dvec, omega):
    """Phases by weighted LSQ over intended edges from a BFS-tree init.
    Returns theta dict (basin = natural phase drops, winding geometric)."""
    idx = {c: k for k, c in enumerate(cells)}
    n = len(cells)
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(edges):
        adj[idx[a]].append((idx[b], e, +1))
        adj[idx[b]].append((idx[a], e, -1))
    theta = np.zeros(n)
    seen = np.zeros(n, bool)
    # BFS tree init: th_j = th_i - w*d  (forward drop along the tree)
    for root in range(n):
        if seen[root]:
            continue
        seen[root] = True
        stack = [root]
        while stack:
            u = stack.pop()
            for v, e, s in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    theta[v] = theta[u] - s * omega * dvec[e]
                    stack.append(v)
    # unwrapped targets from the init basin
    drop = np.zeros(len(edges))
    for e, (a, b) in enumerate(edges):
        i, j = idx[a], idx[b]
        psi0 = wrap(theta[i] - theta[j] - omega * dvec[e])
        drop[e] = (theta[i] - theta[j]) - psi0  # = w*d within the basin, exactly
    # LSQ: minimize sum (th_i - th_j - drop_e)^2   (gauge: mean zero per comp)
    A = np.zeros((len(edges) + 1, n))
    rhs = np.zeros(len(edges) + 1)
    for e, (a, b) in enumerate(edges):
        A[e, idx[a]] = 1.0
        A[e, idx[b]] = -1.0
        rhs[e] = drop[e]
    A[len(edges), :] = 1e-3  # gauge pin
    sol, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    return {c: sol[idx[c]] for c in cells}


def fundamental_cycles(cells, edges):
    """Fundamental cycle basis (edge index lists with signs) via BFS tree."""
    idx = {c: k for k, c in enumerate(cells)}
    n = len(cells)
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(edges):
        adj[idx[a]].append((idx[b], e))
        adj[idx[b]].append((idx[a], e))
    parent = [-1] * n
    pedge = [-1] * n
    depth = [0] * n
    seen = [False] * n
    tree = set()
    order = []
    for root in range(n):
        if seen[root]:
            continue
        seen[root] = True
        stack = [root]
        while stack:
            u = stack.pop()
            order.append(u)
            for v, e in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    parent[v] = u
                    pedge[v] = e
                    depth[v] = depth[u] + 1
                    tree.add(e)
                    stack.append(v)
    cycles = []
    for e, (a, b) in enumerate(edges):
        if e in tree:
            continue
        u, v = idx[a], idx[b]
        pu, pv = [], []
        uu, vv = u, v
        while depth[uu] > depth[vv]:
            pu.append(pedge[uu]); uu = parent[uu]
        while depth[vv] > depth[uu]:
            pv.append(pedge[vv]); vv = parent[vv]
        while uu != vv:
            pu.append(pedge[uu]); uu = parent[uu]
            pv.append(pedge[vv]); vv = parent[vv]
        cycles.append([e] + pu + pv)
    return cycles


def bridges_and_open(cells, edges):
    """Count intended-graph bridge edges and vertices with degree < 2."""
    deg = {c: 0 for c in cells}
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    nopen = sum(1 for c in cells if deg[c] < 2)
    cyc = fundamental_cycles(cells, edges)
    on_cycle = set()
    for cy in cyc:
        on_cycle.update(cy)
    nbridge = len(edges) - len(on_cycle)
    return nopen, nbridge


def is_bipartite(cells, edges):
    color = {}
    adj = {c: [] for c in cells}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    for root in cells:
        if root in color:
            continue
        color[root] = 0
        stack = [root]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    stack.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def score_state(foam, cells, edges, omega, theta=None, want_detail=False):
    """Score an embedded labeled state at pitch omega. edges: list of (i, j)
    foam-cell pairs (intended). Everything else linked in the set = parasite."""
    sc = Score()
    cells = [int(c) for c in cells]
    x = x_of_omega(omega)
    sc.omega, sc.x, sc.n, sc.ne = omega, x, len(cells), len(edges)
    links = foam.links_within(cells)
    ldict = {}
    for i, j, d in links:
        ldict[(min(i, j), max(i, j))] = d
    eset = set((min(a, b), max(a, b)) for a, b in edges)
    sc.broken = [e for e in eset if e not in ldict]
    if sc.broken:
        sc.ok = False
        sc.total = sc.per_voice = 99.0
        sc.leak_edge = sc.leak_par = sc.leak_vac = 99.0
        sc.npar = 0
        sc.gpar_max = 0.0
        sc.labels = []
        sc.cycles = []
        sc.nopen = sc.nbridge = 0
        sc.penalty = 0.0
        sc.wind = 0.0
        sc.par_pairs = []
        return sc
    sc.ok = True
    dvec = np.array([ldict[e] for e in sorted(eset)])
    elist = sorted(eset)
    if theta is None:
        theta = solve_theta(cells, elist, dvec, omega)
    sc.theta = theta
    # intended edges: auto-label strut / cable-fwd / cable-bwd
    leak_e = 0.0
    labels = []
    psis = []
    gs = []
    for (a, b), d in zip(elist, dvec):
        pf = wrap(theta[a] - omega * d - theta[b])
        pb = wrap(theta[b] - omega * d - theta[a])
        gf, gb = float(gate(pf)), float(gate(pb))
        c_strut = (1 - gf) + (1 - gb)
        c_cf = (1 - gf) + W_BACK * gb
        c_cb = (1 - gb) + W_BACK * gf
        cost = min(c_strut, c_cf, c_cb)
        lab = ("strut", "cable_fwd", "cable_bwd")[
            int(np.argmin([c_strut, c_cf, c_cb]))]
        leak_e += cost
        labels.append(lab)
        psis.append((float(pf), float(pb)))
        gs.append((gf, gb))
    # parasites
    leak_p = 0.0
    npar = 0
    gpmax = 0.0
    par_pairs = []
    for i, j, d in links:
        key = (min(i, j), max(i, j))
        if key in eset:
            continue
        pf = wrap(theta[i] - omega * d - theta[j])
        pb = wrap(theta[j] - omega * d - theta[i])
        gf, gb = float(gate(pf)), float(gate(pb))
        leak_p += W_PAR * (gf + gb)
        gpmax = max(gpmax, gf, gb)
        npar += 1
        par_pairs.append((i, j, round(d, 4), round(gf, 4), round(gb, 4)))
    # vacuum ports: unlisted foam neighbours of every structure cell
    cset = set(cells)
    ideg = {c: 0 for c in cells}
    for a, b in elist:
        ideg[a] += 1
        ideg[b] += 1
    nports = 0
    for c in cells:
        inpar = sum(1 for j, _ in foam.adj[c] if j in cset)
        nports += max(0, len(foam.adj[c]) - inpar)
    leak_v = nports * res_11(omega) * G_RAND * W_VAC
    # closure bookkeeping
    cyc = fundamental_cycles(cells, elist)
    cinfo = []
    wind = 0.0
    for cy in cyc:
        L = sum(dvec[e] for e in cy)
        m = omega * L / (2 * np.pi)
        cinfo.append((len(cy), round(float(L), 3), round(float(m), 4)))
        wind += abs(m - len(cy) / 2.0)
    sc.wind = round(float(wind), 3)
    sc.nopen, sc.nbridge = bridges_and_open(cells, elist)
    sc.penalty = PEN_OPEN * sc.nopen + PEN_BRIDGE * sc.nbridge
    sc.leak_edge = leak_e
    sc.leak_par = leak_p
    sc.leak_vac = leak_v
    sc.npar = npar
    sc.gpar_max = gpmax
    sc.labels = list(zip(elist, labels, psis, gs))
    sc.cycles = cinfo
    sc.par_pairs = par_pairs
    sc.per_voice = (leak_e + leak_p + leak_v) / max(1, len(cells))
    sc.total = sc.per_voice + sc.penalty / max(1, len(cells))
    return sc


def best_omega(foam, cells, edges, wgrid=None):
    """Scan pitch over the seedable window, return best score."""
    if wgrid is None:
        wgrid = np.linspace(W_OMEGA_MIN, W_OMEGA_MAX, 41)
    best = None
    for w in wgrid:
        s = score_state(foam, cells, edges, w)
        if best is None or s.total < best.total:
            best = s
    if best is not None and best.ok:
        # local refine
        for w in np.linspace(best.omega - 0.02, best.omega + 0.02, 21):
            if W_OMEGA_MIN <= w <= W_OMEGA_MAX:
                s = score_state(foam, cells, edges, w)
                if s.total < best.total:
                    best = s
    return best


# ---------------------------------------------------------------- embedder
def rand_rot(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def embed_once(foam, coords, edges, center, R, refine=3):
    """Snap ideal coords (rotated R, translated to center) to nearest free
    foam cells, then refine picks to minimize edge-length error vs the
    ideal edge lengths (mirrors the kernel's cube-uniformity passes)."""
    tgt = coords @ R.T + center
    nv = len(coords)
    picks = [-1] * nv
    used = set()
    for k in range(nv):
        dd, ii = foam.tree.query(tgt[k], k=8)
        got = False
        for cand in np.atleast_1d(ii):
            if int(cand) not in used:
                picks[k] = int(cand)
                used.add(int(cand))
                got = True
                break
        if not got:
            return None, None
    d_ideal = {e: np.linalg.norm(coords[e[0]] - coords[e[1]]) for e in edges}
    inc = [[] for _ in range(nv)]
    for e in edges:
        inc[e[0]].append(e)
        inc[e[1]].append(e)
    for _ in range(refine):
        for k in range(nv):
            cands = foam.tree.query_ball_point(tgt[k], 1.2)
            best_sc, best_c = None, picks[k]
            for cand in cands:
                cand = int(cand)
                if cand in used and cand != picks[k]:
                    continue
                s = 0.0
                for e in inc[k]:
                    o = e[1] if e[0] == k else e[0]
                    dd = np.linalg.norm(foam.xyz[cand] - foam.xyz[picks[o]])
                    s += (dd - d_ideal[e]) ** 2
                if best_sc is None or s < best_sc:
                    best_sc, best_c = s, cand
            if best_c != picks[k]:
                used.discard(picks[k])
                used.add(best_c)
                picks[k] = best_c
    cell_edges = [(picks[a], picks[b]) for a, b in edges]
    return picks, cell_edges


def embed_candidate(foam, cand, rng, tries=60, box_margin=4.0):
    """Try many centers/rotations; return (best_score, stats)."""
    coords, edges = cand["coords"], cand["edges"]
    lo, hi = box_margin, 24.0 - box_margin
    n_ok = 0
    spread = []
    best = None
    best_embed = None
    for _ in range(tries):
        center = rng.uniform(lo, hi, 3)
        R = rand_rot(rng)
        picks, cell_edges = embed_once(foam, coords, edges, center, R)
        if picks is None:
            continue
        # realized: every intended edge must be a foam link
        okall = True
        ds = []
        for a, b in cell_edges:
            lk, d = foam.linked(a, b)
            ds.append(d)
            if not lk:
                okall = False
        if not okall:
            continue
        n_ok += 1
        ds = np.array(ds)
        spread.append(ds.std())
        s = best_omega(foam, picks, cell_edges)
        if s is not None and s.ok and (best is None or s.total < best.total):
            best = s
            best_embed = (picks, cell_edges)
    stats = {
        "tries": tries,
        "embed_ok": n_ok,
        "realizability": n_ok / tries,
        "sigma_d_med": float(np.median(spread)) if spread else None,
    }
    return best, best_embed, stats


# ---------------------------------------------------------------- catalog
def ring_coords(n, d, z=0.0):
    R = d / (2 * math.sin(math.pi / n))
    return np.array([[R * math.cos(2 * math.pi * k / n),
                      R * math.sin(2 * math.pi * k / n), z] for k in range(n)])


def g_ring(n, d):
    c = ring_coords(n, d)
    e = [(k, (k + 1) % n) for k in range(n)]
    return {"coords": c, "edges": e}


def g_chain(n, d):
    c = np.array([[k * d, 0.15 * ((-1) ** k), 0] for k in range(n)])
    e = [(k, k + 1) for k in range(n - 1)]
    return {"coords": c, "edges": e}


def g_blob(n, spacing=1.2):
    rng = np.random.default_rng(7)
    c = rng.normal(scale=spacing, size=(n, 3))
    return {"coords": c, "edges": []}


def g_cube(a):
    c = np.array([[(k & 1), (k >> 1) & 1, (k >> 2) & 1] for k in range(8)],
                 float) * a
    e = [(i, i ^ (1 << b)) for i in range(8) for b in range(3) if i < (i ^ (1 << b))]
    return {"coords": c, "edges": e}


def g_tower(layers, a):
    """Stack of `layers` cubes: 2x2x(layers+1) grid points."""
    pts, idx = [], {}
    for z in range(layers + 1):
        for y in range(2):
            for x in range(2):
                idx[(x, y, z)] = len(pts)
                pts.append([x * a, y * a, z * a])
    e = []
    for (x, y, z), i in idx.items():
        for dx, dy, dz in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            j = idx.get((x + dx, y + dy, z + dz))
            if j is not None:
                e.append((i, j))
    return {"coords": np.array(pts, float), "edges": e}


def g_prism(ngon, d, h=None):
    h = h or d
    top = ring_coords(ngon, d, z=h)
    bot = ring_coords(ngon, d, z=0.0)
    c = np.vstack([bot, top])
    e = [(k, (k + 1) % ngon) for k in range(ngon)]
    e += [(ngon + k, ngon + (k + 1) % ngon) for k in range(ngon)]
    e += [(k, ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_tube(ngon, layers, d_ring, d_ax):
    """Open prism stack: `layers` rings of ngon cells, axial struts."""
    rings = [ring_coords(ngon, d_ring, z=i * d_ax) for i in range(layers)]
    c = np.vstack(rings)
    e = []
    for i in range(layers):
        e += [(i * ngon + k, i * ngon + (k + 1) % ngon) for k in range(ngon)]
    for i in range(layers - 1):
        e += [(i * ngon + k, (i + 1) * ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_trunc_oct(d):
    """Truncated octahedron: permutations of (0, +-1, +-2), edge sqrt(2)."""
    s = d / math.sqrt(2)
    pts = set()
    from itertools import permutations
    for p in permutations([0, 1, 2]):
        for s1 in (1, -1):
            for s2 in (1, -1):
                q = [p[0], p[1] * s1, p[2] * s2]
                pts.add(tuple(np.array(q, float)))
    pts = sorted(pts)
    c = np.array(pts) * s
    e = []
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-6:
                e.append((i, j))
    return {"coords": c, "edges": e}


def g_rhombic_dodeca(d):
    s = d / math.sqrt(3)
    pts = [np.array(v, float) for v in
           [(sx, sy, sz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]]
    pts += [np.array(v, float) for v in
            [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)]]
    c = np.array(pts) * s
    e = []
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-6:
                e.append((i, j))
    return {"coords": c, "edges": e}


def g_octahedron(d):
    s = d / math.sqrt(2)
    c = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                  [0, 0, 1], [0, 0, -1]], float) * s
    e = [(i, j) for i in range(6) for j in range(i + 1, 6)
         if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-6]
    return {"coords": c, "edges": e}


def g_icosahedron(d):
    phi = (1 + math.sqrt(5)) / 2
    v = []
    for a, b in [(1, phi), (-1, phi), (1, -phi), (-1, -phi)]:
        v += [[0, a, b], [a, b, 0], [b, 0, a]]
    c = np.array(v, float)
    edge0 = 2.0
    c *= d / edge0
    e = [(i, j) for i in range(12) for j in range(i + 1, 12)
         if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-4]
    return {"coords": c, "edges": e}


def g_antiprism(ngon, d):
    """Uniform antiprism, all edges d."""
    # top ring rotated by pi/n; solve height so slant edge = ring edge = d
    R = d / (2 * math.sin(math.pi / ngon))
    # slant^2 = (2 R sin(pi/2n))^2 + h^2 = d^2
    ch = 2 * R * math.sin(math.pi / (2 * ngon))
    h2 = d * d - ch * ch
    h = math.sqrt(max(h2, 1e-9))
    bot = ring_coords(ngon, d, 0.0)
    top = np.array([[R * math.cos(2 * math.pi * (k + 0.5) / ngon),
                     R * math.sin(2 * math.pi * (k + 0.5) / ngon), h]
                    for k in range(ngon)])
    c = np.vstack([bot, top])
    e = [(k, (k + 1) % ngon) for k in range(ngon)]
    e += [(ngon + k, ngon + (k + 1) % ngon) for k in range(ngon)]
    e += [(k, ngon + k) for k in range(ngon)]
    e += [(( k + 1) % ngon, ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_mobius(k, d):
    """Mobius quad strip = Mobius ladder M_k: rim C_2k + k rungs.
    Ribbon embed: rim step and rung both ~ d. k even => frustrated control."""
    n = 2 * k
    # ribbon: center circle radius R, half-width w, half-twist
    w = d / 2.0
    R = d * k / (2 * math.pi) * 1.05
    pts = []
    for j in range(k):
        t = 2 * math.pi * j / k
        er = np.array([math.cos(t), math.sin(t), 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        nvec = math.cos(t / 2) * er + math.sin(t / 2) * ez
        ccen = R * er
        pts.append(ccen + w * nvec)
        pts.append(ccen - w * nvec)
    c = np.array(pts)
    e = []
    for j in range(k):
        jn = (j + 1) % k
        a0, a1 = 2 * j, 2 * j + 1
        if jn == 0:
            b0, b1 = 1, 0  # half twist: rails swap at the seam
        else:
            b0, b1 = 2 * jn, 2 * jn + 1
        e.append((a0, b0))
        e.append((a1, b1))
        e.append((a0, a1))  # rung
    return {"coords": c, "edges": list(set(tuple(sorted(x)) for x in e))}


def g_theta(arms, d):
    """Theta graph: two poles joined by len(arms) paths (edge counts).
    Arms as circular arcs in planes fanned around the pole axis."""
    pole_sep_edges = max(arms)
    # poles on z axis; arc through angle pi in its own plane
    pts = [np.array([0, 0, 0.0]), None]
    coords = [np.array([0, 0, 0.0])]
    edges = []
    pole2 = None
    for ai, na in enumerate(arms):
        # arc of na edges from pole1 to pole2, chord length na*d*2/pi approx
        # place arc in plane rotated around z by ai*2pi/len(arms)
        phi = ai * 2 * math.pi / len(arms)
        # semicircle radius so that arc chord steps = d: step angle = pi/na
        Rарc = d / (2 * math.sin(math.pi / (2 * na)))
        prev = 0
        for s in range(1, na + 1):
            ang = math.pi * s / na
            local = np.array([Rарc * math.sin(ang),
                              0.0,
                              Rарc * (1 - math.cos(ang))])
            rot = np.array([[math.cos(phi), -math.sin(phi), 0],
                            [math.sin(phi), math.cos(phi), 0],
                            [0, 0, 1]])
            p = rot @ np.array([local[0], 0.02 * ai, 0]) + np.array([0, 0, local[2]])
            if s == na:
                if pole2 is None:
                    pole2 = len(coords)
                    coords.append(p)
                edges.append((prev, pole2))
            else:
                coords.append(p)
                edges.append((prev, len(coords) - 1))
                prev = len(coords) - 1
    return {"coords": np.array(coords), "edges": edges}


def g_cube_in_cube(a):
    inner = g_cube(a)
    b = a + 2 * a / math.sqrt(3)
    outer = g_cube(b)
    ci = inner["coords"] - a / 2
    co = outer["coords"] - b / 2
    c = np.vstack([ci, co])
    e = list(inner["edges"]) + [(i + 8, j + 8) for i, j in outer["edges"]]
    e += [(k, k + 8) for k in range(8)]  # radial
    return {"coords": c, "edges": e}


def g_hopf_rings(n, d, clearance_scale=1.0):
    """Two Hopf-linked rings (each n cells, spacing d): ring A in xy plane,
    ring B in xz plane centered on A's rim."""
    R = d / (2 * math.sin(math.pi / n))
    A = ring_coords(n, d)
    B = ring_coords(n, d)
    rot = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], float)  # xy -> xz
    B = B @ rot.T + np.array([R * clearance_scale, 0, 0])
    c = np.vstack([A, B])
    e = [(k, (k + 1) % n) for k in range(n)]
    e += [(n + k, n + (k + 1) % n) for k in range(n)]
    return {"coords": c, "edges": e}


def catalog(args):
    """(name, class, generator_dict, note) list — every entry gets a ledger row."""
    W_STRUT = 2.0944          # omega at d*=1.5 (x = 0.3205, comp12-class)
    cat = []
    # --- anchors / controls (calibration against measured lifetimes)
    cat.append(("blob12", "control", g_blob(12), "no intended edges; the measured -0.232%/tu class"))
    cat.append(("chain12", "control", g_chain(12, 1.5), "open chain; measured 1.7x worse than ring"))
    cat.append(("ring12_m6", "anchor", g_ring(12, 1.5), "unwound strut ring (m=6): died 2221"))
    cat.append(("ring12_m5", "anchor", g_ring(12, 1.25), "comp12 wound cable ring (m=5): alive t=5000"))
    cat.append(("ring6_m2", "anchor", g_ring(6, 1.10), "comp6 one-way (m=2): died 3836; 2nd-nbr 1.905 parasites?"))
    cat.append(("cube_a1.25", "anchor", g_cube(1.25), "H1 cube: died ~1200-1800; face diagonals 1.77 = parasites"))
    # --- bipartite strut shells (2D skins)
    cat.append(("cube_a1.5", "shell", g_cube(1.5), "P2: diagonals 2.12 past ceiling; edges at dbar"))
    cat.append(("tower2_a1.5", "shell", g_tower(2, 1.5), "2-cube stack = capped square tube, n=12"))
    cat.append(("tower3_a1.5", "shell", g_tower(3, 1.5), "3-cube stack, n=16"))
    cat.append(("hexprism_d1.5", "shell", g_prism(6, 1.5), "n=12 3-regular bipartite"))
    cat.append(("octprism_d1.5", "shell", g_prism(8, 1.5), "n=16 3-regular bipartite"))
    cat.append(("truncoct_d1.5", "shell", g_trunc_oct(1.5), "Kelvin cell n=24 bipartite 3-regular"))
    cat.append(("rhombdodec_d1.5", "shell", g_rhombic_dodeca(1.5), "PREDICT frustrated: face short-diagonal 1.73 parasite x12"))
    cat.append(("cubeincube_a1.5", "shell", g_cube_in_cube(1.5), "PREDICT broken: outer edges 3.23 beyond link ceiling"))
    # --- tubes (open homology carriers)
    cat.append(("tube6x3_d1.5", "tube", g_tube(6, 3, 1.5, 1.5), "hex tube 3 rings n=18"))
    cat.append(("tube8x3_d1.5", "tube", g_tube(8, 3, 1.5, 1.5), "oct tube 3 rings n=24"))
    cat.append(("tube6x4_d1.5", "tube", g_tube(6, 4, 1.5, 1.5), "hex tube 4 rings n=24"))
    # wound tube: ring edges at phi=5pi/6 (w=1), axial struts at pi; one omega
    w_ax = math.pi / 1.6
    d_ring_w = (math.pi - 2 * math.pi / 12) / w_ax
    cat.append(("woundtube12x2_w1", "tube", g_tube(12, 2, d_ring_w, 1.6),
                f"B1c co-rotating tube: d_ring={d_ring_w:.3f} d_ax=1.6, both rings m=5"))
    # --- theta / book (multi-cycle cable systems)
    cat.append(("theta_666_d1.5", "theta", g_theta([6, 6, 6], 1.5), "3 arms x6 edges; junction deg 3"))
    cat.append(("theta_668_d1.5", "theta", g_theta([6, 6, 8], 1.5), "unequal arms: relative cycle integers differ"))
    # --- negative controls (predicted frustrated)
    cat.append(("octahedron_d1.5", "control", g_octahedron(1.5), "NEG: odd cycles; also equator 4-cycles fail"))
    cat.append(("icosahedron_d1.5", "control", g_icosahedron(1.5), "NEG: triangulated shell"))
    cat.append(("antiprism4_d1.5", "control", g_antiprism(4, 1.5), "NEG: triangles; slant phi unreachable"))
    cat.append(("mobius6_d1.5", "control", g_mobius(6, 1.5), "NEG: Mobius ladder k=6 even => odd cycles"))
    # --- exotic-load ring (max load, foam-floor links)
    cat.append(("ring9_m2_x0.9", "ring", g_ring(9, 1.003), "phi=80deg one-way at x=0.90: max skirt margin, back gate ~3e-15"))
    cat.append(("ring8_m3", "ring", g_ring(8, 1.407), "phi=3pi/4 wound-1: back gate G(90)=0.0039, d=1.407"))
    cat.append(("ring10_m4", "ring", g_ring(10, 1.35), "phi=4pi/5: back gate G(72)=0.026"))
    cat.append(("ring16_m8", "ring", g_ring(16, 1.5), "large strut ring n=16"))
    # --- Hopf link molecule (topological bond, no channels between rings)
    cat.append(("hopf12x12_m5", "exotic", g_hopf_rings(12, 1.25),
                "two comp12 rings Hopf-linked; bond is topological only"))
    return cat


# ---------------------------------------------------------------- free search
def free_search(foam, rng, omega_class, n_restarts=6, steps=1200,
                nmin=8, nmax=32, seed_states=None, log=None):
    """Simulated annealing over foam subgraphs at fixed pitch class.
    State: (cells frozenset, intended edge frozenset). Labels are derived."""
    results = []
    lo, hi = 4.0, 20.0

    def rand_seed_state():
        while True:
            c0 = rng.integers(0, foam.n)
            p = foam.xyz[c0]
            if lo < p[0] < hi and lo < p[1] < hi and lo < p[2] < hi:
                break
        cells = [int(c0)]
        cset = {int(c0)}
        d_star = math.pi / omega_class
        while len(cells) < nmin + 4:
            # grow along links closest to the strut length
            frontier = []
            for c in cells:
                for j, d in foam.adj[c]:
                    if j not in cset:
                        frontier.append((abs(d - d_star), c, j))
            if not frontier:
                break
            frontier.sort()
            _, c, j = frontier[rng.integers(0, min(4, len(frontier)))]
            cells.append(int(j))
            cset.add(int(j))
        edges = set()
        for i, j, d in foam.links_within(cells):
            if abs(d - d_star) < 0.25:
                edges.add((min(i, j), max(i, j)))
        return cset, edges

    def state_score(cset, edges):
        n = len(cset)
        if n < 4:
            return None
        s = score_state(foam, sorted(cset), sorted(edges), omega_class)
        if not s.ok:
            return None
        extra = 0.0
        if n < nmin:
            extra += 0.3 * (nmin - n)
        if n > nmax:
            extra += 0.5 * (n - nmax)
        s.total += extra
        return s

    for restart in range(n_restarts):
        if seed_states and restart < len(seed_states):
            cset, edges = seed_states[restart]
            cset = set(int(c) for c in cset)
            edges = set((min(a, b), max(a, b)) for a, b in edges)
        else:
            cset, edges = rand_seed_state()
        cur = state_score(cset, edges)
        if cur is None:
            continue
        best = cur
        best_state = (frozenset(cset), frozenset(edges))
        T0, T1 = 0.08, 0.004
        for step in range(steps):
            T = T0 * (T1 / T0) ** (step / max(1, steps - 1))
            move = rng.random()
            nc, ne = set(cset), set(edges)
            if move < 0.30 and len(nc) < nmax + 2:
                # add a cell linked to the set, wire its best links
                cands = []
                d_star = math.pi / omega_class
                cl = list(nc)
                for _ in range(6):
                    c = cl[rng.integers(0, len(cl))]
                    for j, d in foam.adj[c]:
                        if j not in nc:
                            cands.append((abs(d - d_star), c, j, d))
                if not cands:
                    continue
                cands.sort()
                _, c, j, d = cands[rng.integers(0, min(5, len(cands)))]
                nc.add(int(j))
                for jj, dd in foam.adj[int(j)]:
                    if jj in nc and abs(dd - d_star) < 0.22:
                        ne.add((min(int(j), jj), max(int(j), jj)))
                if not any(int(j) in e for e in ne):
                    ne.add((min(int(j), c), max(int(j), c)))
            elif move < 0.55 and len(nc) > nmin - 2:
                # remove a low-degree cell
                degs = {}
                for a, b in ne:
                    degs[a] = degs.get(a, 0) + 1
                    degs[b] = degs.get(b, 0) + 1
                cl = sorted(nc, key=lambda c: degs.get(c, 0))
                victim = cl[rng.integers(0, min(4, len(cl)))]
                nc.discard(victim)
                ne = set(e for e in ne if victim not in e)
                # connectivity check
                if nc:
                    seen = set()
                    stack = [next(iter(nc))]
                    seen.add(stack[0])
                    while stack:
                        u = stack.pop()
                        for v, _ in foam.adj[u]:
                            if v in nc and v not in seen:
                                seen.add(v)
                                stack.append(v)
                    if seen != nc:
                        continue
            else:
                # toggle an internal link intended <-> parasite
                links = foam.links_within(nc)
                if not links:
                    continue
                i, j, d = links[rng.integers(0, len(links))]
                key = (min(i, j), max(i, j))
                if key in ne:
                    ne.discard(key)
                else:
                    ne.add(key)
            nxt = state_score(nc, ne)
            if nxt is None:
                continue
            dE = nxt.total - cur.total
            if dE < 0 or rng.random() < math.exp(-dE / max(T, 1e-9)):
                cset, edges, cur = nc, ne, nxt
                if cur.total < best.total:
                    best = cur
                    best_state = (frozenset(cset), frozenset(edges))
        results.append((best, best_state))
        if log is not None:
            log(f"  restart {restart}: best/voice={best.per_voice:.4f} "
                f"total={best.total:.4f} n={best.n} ne={best.ne} "
                f"npar={best.npar} open={best.nopen} bridges={best.nbridge}")
    return results


# ---------------------------------------------------------------- exotics
def exotic_fifth_bridge():
    """Interval-rung bridge (3:2): two unison clusters joined by fifth links.
    Paper numbers only; returns dict for the report."""
    out = {}
    # reachable fifth inside the seedable pitch window
    w_lo_min, w_lo_max = W_OMEGA_MIN, W_OMEGA_MAX / 1.5
    out["w_low_window"] = [w_lo_min, w_lo_max]
    w_lo = 1.45
    w_hi = 1.5 * w_lo
    out["example"] = {"w_low": w_lo, "w_high": w_hi,
                      "x_low": x_of_omega(w_lo), "x_high": x_of_omega(w_hi)}
    # bridge rung: (q*w_hi + p*w_lo) d = 2 pi m  with p=3 (low side), q=2
    Om = 2 * w_hi  # = 3*w_lo coincident partial
    for m in (1, 2):
        out[f"bridge_d_m{m}"] = math.pi * m / Om
    # acceptance: gw = gamma_m / (p q); tolerance in loads
    gw = GAMMA_M / 6.0
    ddet_dxhi = 2 * QDET * w_hi ** 2 / W2
    ddet_dxlo = 3 * QDET * w_lo ** 2 / W2
    out["acceptance_det"] = gw
    out["dx_tolerance_rms"] = gw / math.hypot(ddet_dxhi, ddet_dxlo)
    out["res_on_tooth"] = 1.0 / 6.0
    return out


def exotic_hopf_geometry(n=12, d=1.25):
    """Min distance between two Hopf-linked ring circles vs the link ceiling."""
    R = d / (2 * math.sin(math.pi / n))
    t = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    A = np.stack([R * np.cos(t), R * np.sin(t), 0 * t], 1)
    B = np.stack([R * np.cos(t) + R, 0 * t, R * np.sin(t)], 1)
    dmin = np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2))
    return {"R": R, "min_dist": float(dmin), "ceiling": 1.96,
            "clean": bool(dmin > 2.07)}


# ---------------------------------------------------------------- output
def edge_records(sc):
    recs = []
    for (a, b), lab, (pf, pb), (gf, gb) in sc.labels:
        recs.append({"i": int(a), "j": int(b), "label": lab,
                     "psi_f": round(pf, 4), "psi_b": round(pb, 4),
                     "gf": round(gf, 4), "gb": round(gb, 4)})
    return recs


def write_net(path, foam, sc, cells):
    """Frozen .net format: V x y z xk th2 / E a b."""
    idx = {c: k for k, c in enumerate(cells)}
    with open(path, "w") as f:
        f.write(f"# morpho export: omega={sc.omega:.5f} x={sc.x:.5f} "
                f"leak/voice={sc.per_voice:.4f} npar={sc.npar}\n")
        for c in cells:
            th = sc.theta[c] % (2 * math.pi)
            p = foam.xyz[c]
            f.write(f"V {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} {sc.x:.6f} {th:.6f}\n")
        for (a, b), lab, _, _ in sc.labels:
            f.write(f"E {idx[a]} {idx[b]}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--tries", type=int, default=None)
    ap.add_argument("--restarts", type=int, default=None)
    ap.add_argument("--steps", type=int, default=None)
    args = ap.parse_args()
    tries = args.tries or (12 if args.quick else 60)
    restarts = args.restarts or (2 if args.quick else 6)
    steps = args.steps or (250 if args.quick else 1200)

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    print(f"# loading foam {FOAM}")
    foam = Foam(FOAM)
    print(f"# foam: {foam.n} cells, {len(foam.pairs)} links, "
          f"dbar={foam.pdist.mean():.4f}, deg mean={foam.deg.mean():.1f}")

    ledger = []      # (name, class, verdict, key numbers)
    results = {}     # name -> dict for the report

    # ---------------- sweep 1: catalog
    print("\n# === SWEEP 1: catalog ===")
    for name, klass, cand, note in catalog(args):
        bip = is_bipartite(range(len(cand["coords"])), cand["edges"])
        best, emb, stats = embed_candidate(foam, cand, rng, tries=tries)
        row = {"name": name, "class": klass, "note": note,
               "n": len(cand["coords"]), "ne": len(cand["edges"]),
               "bipartite": bip, "stats": stats}
        if best is None or not best.ok:
            verdict = "UNREALIZABLE" if (stats["embed_ok"] == 0) else "NO-SCORE"
            row["verdict"] = verdict
            print(f"  {name:22s} [{klass}] n={row['n']:3d} -> {verdict} "
                  f"(embeds {stats['embed_ok']}/{stats['tries']})")
        else:
            labs = [l for _, l, _, _ in best.labels]
            row.update({
                "verdict": "scored", "omega": round(best.omega, 4),
                "x": round(best.x, 4),
                "leak_per_voice": round(best.per_voice, 4),
                "leak_edge": round(best.leak_edge, 4),
                "leak_par": round(best.leak_par, 4),
                "leak_vac": round(best.leak_vac, 4),
                "npar": best.npar, "gpar_max": round(best.gpar_max, 4),
                "nopen": best.nopen, "nbridge": best.nbridge,
                "wind": best.wind,
                "n_strut": labs.count("strut"),
                "n_cable": labs.count("cable_fwd") + labs.count("cable_bwd"),
                "cycles": best.cycles,
                "cells": [int(c) for c in emb[0]],
                "edges": edge_records(best),
                "penalty": round(best.penalty, 3),
                "total": round(best.total, 4),
            })
            row["score_obj"] = best
            print(f"  {name:22s} [{klass}] n={row['n']:3d} leak/voice="
                  f"{best.per_voice:.4f} (e={best.leak_edge:.3f} p={best.leak_par:.3f} "
                  f"v={best.leak_vac:.3f}) npar={best.npar} x={best.x:.3f} "
                  f"w={best.omega:.4f} strut/cable={row['n_strut']}/{row['n_cable']} "
                  f"realiz={stats['realizability']:.2f}")
        ledger.append(row)
        results[name] = row

    # ---------------- sweep 2: free search
    print("\n# === SWEEP 2: free search (backtracer) ===")
    free_finds = []
    omega_classes = [2.0944, 1.9635, 2.2308, 1.8000, 1.3963]
    seed_from = []
    for name in ("tube6x3_d1.5", "ring12_m5", "truncoct_d1.5"):
        r = results.get(name)
        if r and r.get("verdict") == "scored":
            cells = r["cells"]
            edges = [(e["i"], e["j"]) for e in r["edges"]]
            seed_from.append((set(cells), set((min(a, b), max(a, b))
                                              for a, b in edges)))
    for wc in omega_classes:
        print(f"  omega class {wc:.4f} (d*={math.pi/wc:.3f}, x={x_of_omega(wc):.3f})")
        res = free_search(foam, rng, wc, n_restarts=restarts, steps=steps,
                          seed_states=seed_from if wc == 2.0944 else None,
                          log=lambda s: print(s))
        for best, (cset, eset) in res:
            free_finds.append((best, sorted(cset), sorted(eset), wc))

    # polish top free finds with full omega scan
    free_finds.sort(key=lambda t: t[0].total)
    polished = []
    seen_sets = set()
    for best, cells, edges, wc in free_finds[:12]:
        key = frozenset(cells)
        if key in seen_sets:
            continue
        seen_sets.add(key)
        s = best_omega(foam, cells, edges)
        if s and s.ok:
            polished.append((s, cells, edges))
    polished.sort(key=lambda t: t[0].total)
    for k, (s, cells, edges) in enumerate(polished[:8]):
        name = f"free_{k}"
        labs = [l for _, l, _, _ in s.labels]
        row = {"name": name, "class": "free-search",
               "note": "backtracer find", "n": s.n, "ne": s.ne,
               "bipartite": is_bipartite(cells, edges),
               "verdict": "scored", "omega": round(s.omega, 4),
               "x": round(s.x, 4), "leak_per_voice": round(s.per_voice, 4),
               "leak_edge": round(s.leak_edge, 4),
               "leak_par": round(s.leak_par, 4),
               "leak_vac": round(s.leak_vac, 4),
               "npar": s.npar, "gpar_max": round(s.gpar_max, 4),
               "nopen": s.nopen, "nbridge": s.nbridge, "wind": s.wind,
               "n_strut": labs.count("strut"),
               "n_cable": labs.count("cable_fwd") + labs.count("cable_bwd"),
               "cycles": s.cycles, "cells": [int(c) for c in cells],
               "edges": edge_records(s), "penalty": round(s.penalty, 3),
               "total": round(s.total, 4),
               "stats": {"realizability": 1.0, "sigma_d_med": None,
                         "tries": 1, "embed_ok": 1},
               "score_obj": s}
        ledger.append(row)
        results[name] = row
        print(f"  {name}: n={s.n} ne={s.ne} leak/voice={s.per_voice:.4f} "
              f"npar={s.npar} open={s.nopen} bridges={s.nbridge} x={s.x:.3f}")

    # ---------------- sweep 3: exotics (paper)
    print("\n# === SWEEP 3: exotics (paper checks) ===")
    ex = {"fifth_bridge": exotic_fifth_bridge(),
          "hopf_geometry": exotic_hopf_geometry()}
    print(f"  fifth bridge: {ex['fifth_bridge']}")
    print(f"  hopf geometry: {ex['hopf_geometry']}")

    # ---------------- shortlist
    print("\n# === SHORTLIST ===")
    scored = [r for r in ledger if r.get("verdict") == "scored"
              and r["class"] not in ("control",)
              and r.get("x", 0) >= X_SEED_MIN - 1e-6
              and r.get("nopen", 9) == 0 and r.get("nbridge", 9) == 0]
    scored.sort(key=lambda r: r["leak_per_voice"])
    shortlist = []
    have = set()
    # top by leak, but guarantee >=2 wound (wind>0.3) candidates for diversity
    wound = [r for r in scored if r.get("wind", 0) > 0.3]
    for r in scored:
        if len(shortlist) >= 8:
            break
        if r["name"] in have:
            continue
        shortlist.append(r)
        have.add(r["name"])
    nw = sum(1 for r in shortlist if r.get("wind", 0) > 0.3)
    wi = 0
    while nw < 2 and wi < len(wound):
        r = wound[wi]
        wi += 1
        if r["name"] in have:
            continue
        shortlist[-1] = r
        have.add(r["name"])
        nw += 1
    negatives = [results.get("mobius6_d1.5"), results.get("octahedron_d1.5")]

    out = {"laws": {"w2": W2, "q_detune": QDET, "cap": CAP, "p_gate": PGATE,
                    "gamma_res_m": GAMMA_M, "x_skirt": X_SKIRT,
                    "seed_window_x": [X_SEED_MIN, X_MAX],
                    "omega_window": [W_OMEGA_MIN, W_OMEGA_MAX]},
           "score_weights": {"W_BACK": W_BACK, "W_PAR": W_PAR, "W_VAC": W_VAC,
                             "G_RAND": G_RAND, "PEN_OPEN": PEN_OPEN,
                             "PEN_BRIDGE": PEN_BRIDGE},
           "shortlist": [], "negative_controls": [], "exotics": ex}

    def strip(r):
        r2 = {k: v for k, v in r.items() if k not in ("score_obj",)}
        return r2

    for rank, r in enumerate(shortlist):
        rr = strip(r)
        rr["rank"] = rank + 1
        out["shortlist"].append(rr)
        print(f"  #{rank+1} {r['name']:22s} leak/voice={r['leak_per_voice']:.4f} "
              f"x={r['x']:.3f} npar={r['npar']} wind={r.get('wind',0)}")
        sc = r.get("score_obj")
        if sc is not None:
            write_net(os.path.join(HERE, "nets", f"{r['name']}.net"),
                      foam, sc, r["cells"])
    for r in negatives:
        if r:
            out["negative_controls"].append(strip(r))
            sc = r.get("score_obj")
            if sc is not None and r.get("verdict") == "scored":
                write_net(os.path.join(HERE, "nets", f"{r['name']}.net"),
                          foam, sc, r["cells"])

    with open(os.path.join(HERE, "shortlist.json"), "w") as f:
        json.dump(out, f, indent=1)

    # full ledger tsv
    with open(os.path.join(HERE, "ledger.tsv"), "w") as f:
        f.write("name\tclass\tn\tverdict\tbipartite\tomega\tx\tleak_per_voice\t"
                "leak_edge\tleak_par\tleak_vac\tnpar\tgpar_max\tnopen\tnbridge\t"
                "wind\tn_strut\tn_cable\trealizability\tnote\n")
        for r in ledger:
            f.write("\t".join(str(r.get(k, "")) for k in
                    ("name", "class", "n", "verdict", "bipartite", "omega", "x",
                     "leak_per_voice", "leak_edge", "leak_par", "leak_vac",
                     "npar", "gpar_max", "nopen", "nbridge", "wind",
                     "n_strut", "n_cable")) +
                    f"\t{r['stats']['realizability'] if r.get('stats') else ''}"
                    f"\t{r['note']}\n")

    with open(os.path.join(HERE, "report.json"), "w") as f:
        json.dump({"ledger": [strip(r) for r in ledger], "exotics": ex},
                  f, indent=1)
    print(f"\n# done in {time.time()-t0:.1f}s -> shortlist.json, ledger.tsv, "
          f"report.json, nets/")


if __name__ == "__main__":
    main()
