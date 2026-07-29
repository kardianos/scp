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
  static lock    uniform w per locked component  =>  struts: d = pi/w (m=1)
  cable          forward-tuned, back gate = gate(-2 w d)  [comp12 anchor 0.100]
  cycle closure  w * L_cycle = 2 pi m  for every independent cycle

Laws: laws_V2g.cfg — w2=2.9 q_detune=1.2 cap=2.5 C=1 p_gate=8
      gamma_res_m=0.10 mob_floor=0.004 s_pull=0.15 r0=0.85 rjit=0.06.

Leak model (a RANKING device; the sim is the judge):
  intended strut  : (1-gf) + (1-gb)            want both gates open
  intended cable  : (1-g_fwd) + W_BACK*g_back  back-gate openness = desert-
                    protected frustration (comp12: gb=0.100, longest-lived)
  parasite (same locked component): gf + gb at the solved phases — the lock
                    web pins psi, openness radiates off-rung
  parasite (cross-component / unlocked cell): 2*<G>_rand — relative gauge
                    drifts/slips, expectation of the gate over psi
  vacuum port     : res_11(w) * <G>_rand * W_VAC per unlisted foam neighbour
                    (res_11 = gm^2/(gm^2+(w2-w)^2): the skirt trickle — this
                    term reproduces the measured lifetime-rises-with-load law)
Phases are solved by least squares over INTENDED edges only (locks form on
the strong channels; parasites are scored, not optimised away) from a
BFS-tree init that fixes the winding basin — winding is geometric here:
it enters through the edge lengths and the pitch, exactly as in ring_m.

VALIDITY DOMAIN: anchored by measured lifetimes at per-voice x in
[0.128, 0.44] (skirt death ladder, MASS.md). Scores at x > 0.62 are
extrapolation and are flagged x_unmeasured; the pitch scan is clamped to
x <= 0.62 unless a candidate is explicitly heavy_ok.

Outputs: ledger.tsv (every candidate, one line), shortlist.json (top 8 +
2 negative controls: concrete cell ids, positions, edges with labels,
omega, per-cycle integers), nets/*.net (frozen V/E seedfiles), report.json.

Usage:  python3 search.py [--quick] [--seed N] [--tries N] [--restarts N]
"""

import argparse
import json
import math
import os
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
GAMMA_M = 0.10    # gamma_res_m (dense-sector acceptance)
MOB_FLOOR = 0.004
S_PULL = 0.15
COMB_LIMIT = 6
X_SKIRT = 0.0617  # measured vacuum-skirt death boundary
X_SEED_MIN = 0.25  # campaign floor for seeded load
X_MEASURED_MAX = 0.62  # heaviest measured survivor class ~0.44 + margin
X_MAX = 0.90
LINK_MARGIN = 1.15
R0 = 0.85

W_AT = lambda x: W2 / (1.0 + QDET * x)
W_OMEGA_MIN = W_AT(X_MAX)            # 1.3942
W_OMEGA_SAFE = W_AT(X_MEASURED_MAX)  # 1.6629  (x = 0.62 validity edge)
W_OMEGA_MAX = W_AT(X_SEED_MIN)       # 2.2308

# scoring weights (documented in MORPHO.md; ranking only)
W_BACK = 0.5
W_PAR = 1.0
W_VAC = 0.11          # sympathetic-floor conductance ratio ~ sqrt(mob_floor*cap/Em)
G_RAND = 0.19638      # <gate>_psi~U = C(16,8)/2^16 (drifting-pair expectation)
PEN_OPEN = 0.40       # per vertex with intended degree < 2 (chain end; MASS 1.7x)
PEN_BRIDGE = 0.20     # per bridge edge (an edge carrying no closed cycle)
PEN_FLUX = 1.00       # x normalized gate-flux imbalance (Kirchhoff: steady
                      # circulation needs per-vertex inflow = outflow; a
                      # source/sink pattern congests and unlocks in the sim)


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
    def __init__(self, xyz, r):
        self.xyz = np.asarray(xyz, float)
        self.r = np.asarray(r, float)
        self.n = len(self.r)
        self.tree = cKDTree(self.xyz)
        pairs = self.tree.query_pairs(LINK_MARGIN * 2 * self.r.max(),
                                      output_type="ndarray")
        if len(pairs) == 0:
            pairs = np.zeros((0, 2), int)
        dd = (np.linalg.norm(self.xyz[pairs[:, 0]] - self.xyz[pairs[:, 1]], axis=1)
              if len(pairs) else np.zeros(0))
        cut = (LINK_MARGIN * (self.r[pairs[:, 0]] + self.r[pairs[:, 1]])
               if len(pairs) else np.zeros(0))
        m = dd < cut
        self.pairs = pairs[m]
        self.pdist = dd[m]
        self.adj = [[] for _ in range(self.n)]
        for (a, b), dist in zip(self.pairs, self.pdist):
            self.adj[a].append((int(b), float(dist)))
            self.adj[b].append((int(a), float(dist)))
        self.deg = np.array([len(a) for a in self.adj])

    @classmethod
    def load(cls, path):
        d = np.loadtxt(path, skiprows=1)
        assert np.all(d[:, 0].astype(int) == np.arange(len(d))), "ids not sequential"
        return cls(d[:, 1:4], d[:, 4])

    def linked(self, i, j):
        d = float(np.linalg.norm(self.xyz[i] - self.xyz[j]))
        return d < LINK_MARGIN * (self.r[i] + self.r[j]), d

    def links_within(self, cells):
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
                 "per_voice", "labels", "npar", "gpar_max",
                 "n", "ne", "theta", "cycles", "nopen", "nbridge", "penalty",
                 "total", "broken", "par_pairs", "wind", "mass", "flux_imbal")


def components(cells, edges):
    comp = {c: -1 for c in cells}
    adj = {c: [] for c in cells}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    k = 0
    for root in cells:
        if comp[root] >= 0:
            continue
        comp[root] = k
        stack = [root]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if comp[v] < 0:
                    comp[v] = k
                    stack.append(v)
        k += 1
    return comp, k


def solve_theta(cells, edges, dvec, omega, drops=None, sweeps=250, eta=0.18):
    """Phases by the kernel's own lock relaxation: each edge pulls its
    receiver toward forward alignment at a rate weighted by the (open)
    gate — the rate-level fixed point of kappa_lock (cellfab.c:2034).
    Winding basins are respected exactly as in the sim: a dead gate pulls
    nothing (the desert). Init: BFS tree with per-edge drops (natural
    omega*d, or the candidate's target drops when supplied)."""
    idx = {c: k for k, c in enumerate(cells)}
    n = len(cells)
    if not edges:
        return {c: 0.0 for c in cells}
    if drops is None:
        drops = omega * dvec
    # drops may be SIGNED: value applies along storage direction (a -> b);
    # a negative value means the design transport runs b -> a.
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(edges):
        adj[idx[a]].append((idx[b], e, +1))
        adj[idx[b]].append((idx[a], e, -1))
    theta = np.zeros(n)
    seen = np.zeros(n, bool)
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
                    theta[v] = theta[u] - s * drops[e]
                    stack.append(v)
    ia = np.array([idx[a] for a, b in edges])
    ib = np.array([idx[b] for a, b in edges])
    wd = omega * dvec
    deg = np.zeros(n)
    np.add.at(deg, ib, 1.0)
    np.add.at(deg, ia, 1.0)
    deg = np.maximum(deg, 1.0)
    for it in range(sweeps):
        pf = wrap(theta[ia] - wd - theta[ib])   # error at receiver b
        pb = wrap(theta[ib] - wd - theta[ia])   # error at receiver a
        gf = gate(pf)
        gb = gate(pb)
        upd = np.zeros(n)
        np.add.at(upd, ib, eta * gf * pf)
        np.add.at(upd, ia, eta * gb * pb)
        upd /= deg
        theta += upd
        if it >= 50 and it % 10 == 0 and np.max(np.abs(upd)) < 1e-4:
            break
    return {c: theta[idx[c]] for c in cells}


def fundamental_cycles(cells, edges):
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
    for root in range(n):
        if seen[root]:
            continue
        seen[root] = True
        stack = [root]
        while stack:
            u = stack.pop()
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
        while depth[u] > depth[v]:
            pu.append(pedge[u]); u = parent[u]
        while depth[v] > depth[u]:
            pv.append(pedge[v]); v = parent[v]
        while u != v:
            pu.append(pedge[u]); u = parent[u]
            pv.append(pedge[v]); v = parent[v]
        cycles.append([e] + pu + pv)
    return cycles


def bridges_and_open(cells, edges):
    deg = {c: 0 for c in cells}
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    nopen = sum(1 for c in cells if deg[c] < 2)
    on_cycle = set()
    for cy in fundamental_cycles(cells, edges):
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


def score_state(foam, cells, edges, omega, count_ports=True, drops_map=None):
    """Score an embedded state at pitch omega. edges: intended (i, j) pairs.
    Every other foam link inside the set is a parasite. drops_map (optional):
    edge-key -> target retardation, fixing the winding basin of the seed."""
    sc = Score()
    cells = [int(c) for c in cells]
    sc.omega, sc.x, sc.n, sc.ne = omega, x_of_omega(omega), len(cells), len(edges)
    sc.mass = round(sc.n * sc.x * CAP, 3)
    links = foam.links_within(cells)
    ldict = {(min(i, j), max(i, j)): d for i, j, d in links}
    eset = set((min(a, b), max(a, b)) for a, b in edges)
    sc.broken = [e for e in eset if e not in ldict]
    if sc.broken:
        sc.ok = False
        sc.total = sc.per_voice = 99.0
        sc.leak_edge = sc.leak_par = sc.leak_vac = 99.0
        sc.npar = 0; sc.gpar_max = 0.0; sc.labels = []; sc.cycles = []
        sc.nopen = sc.nbridge = 0; sc.penalty = 0.0; sc.wind = 0.0
        sc.par_pairs = []; sc.theta = {}; sc.flux_imbal = 0.0
        return sc
    sc.ok = True
    elist = sorted(eset)
    dvec = np.array([ldict[e] for e in elist])
    comp, ncomp = components(cells, elist)
    drops = (np.array([drops_map.get(e, omega * ldict[e]) for e in elist])
             if drops_map else None)
    theta = solve_theta(cells, elist, dvec, omega, drops=drops)
    sc.theta = theta
    leak_e = 0.0
    labels = []
    for (a, b), d in zip(elist, dvec):
        pf = float(wrap(theta[a] - omega * d - theta[b]))
        pb = float(wrap(theta[b] - omega * d - theta[a]))
        gf, gb = float(gate(pf)), float(gate(pb))
        c_strut = (1 - gf) + (1 - gb)
        c_cf = (1 - gf) + W_BACK * gb
        c_cb = (1 - gb) + W_BACK * gf
        k = int(np.argmin([c_strut, c_cf, c_cb]))
        leak_e += (c_strut, c_cf, c_cb)[k]
        labels.append(((a, b), ("strut", "cable_fwd", "cable_bwd")[k],
                       (pf, pb), (gf, gb)))
    leak_p = 0.0
    npar = 0
    gpmax = 0.0
    par_pairs = []
    for i, j, d in links:
        key = (min(i, j), max(i, j))
        if key in eset:
            continue
        if comp[i] == comp[j]:
            pf = float(wrap(theta[i] - omega * d - theta[j]))
            pb = float(wrap(theta[j] - omega * d - theta[i]))
            gf, gb = float(gate(pf)), float(gate(pb))
        else:
            gf = gb = G_RAND  # relative gauge drifts/slips: expectation
        leak_p += W_PAR * (gf + gb)
        gpmax = max(gpmax, gf, gb)
        npar += 1
        par_pairs.append((i, j, round(d, 4), round(gf, 4), round(gb, 4)))
    nports = 0
    if count_ports:
        cset = set(cells)
        for c in cells:
            inside = sum(1 for j, _ in foam.adj[c] if j in cset)
            nports += max(0, len(foam.adj[c]) - inside)
    leak_v = nports * res_11(omega) * G_RAND * W_VAC
    cinfo = []
    wind = 0.0
    for cy in fundamental_cycles(cells, elist):
        L = float(sum(dvec[e] for e in cy))
        # ordered traversal of the cycle's edges -> phase winding number
        # m = (1/2pi) * sum of retardations (th_u - th_v mod 2pi): an exact
        # integer for any single-valued phase field — the homology label
        einc = {}
        for e in cy:
            a, b = elist[e]
            einc.setdefault(a, []).append((e, b))
            einc.setdefault(b, []).append((e, a))
        start = elist[cy[0]][0]
        walk = [start]
        used_e = set()
        cur = start
        while True:
            nxt = None
            for e, o in einc[cur]:
                if e not in used_e:
                    used_e.add(e)
                    nxt = o
                    break
            if nxt is None:
                break
            walk.append(nxt)
            cur = nxt
            if cur == start:
                break
        m = 0.0
        if len(walk) > 1 and walk[-1] == walk[0]:
            votes_fwd = 0
            dmap = {}
            for e in cy:
                a, b = elist[e]
                dmap[(a, b)] = dvec[e]
                dmap[(b, a)] = dvec[e]
            for k in range(len(walk) - 1):
                u, v = walk[k], walk[k + 1]
                m += (theta[u] - theta[v]) % (2 * np.pi)
                if abs(wrap(theta[u] - omega * dmap[(u, v)] - theta[v])) < np.pi / 2:
                    votes_fwd += 1
            m /= 2 * np.pi
            if votes_fwd < (len(walk) - 1) / 2.0:
                m = (len(walk) - 1) - m  # report along the transport direction
        cinfo.append((len(cy), round(L, 3), int(round(m))))
        wind += abs(m - len(cy) / 2.0)
    sc.wind = round(float(wind), 3)
    # Kirchhoff balance of open-gate flux: steady circulation requires
    # per-vertex inflow = outflow (a folded/pumped pattern congests)
    fin = {c: 0.0 for c in cells}
    fout = {c: 0.0 for c in cells}
    for (a, b), lab, (pf, pb), (gf, gb) in labels:
        fout[a] += gf; fin[b] += gf
        fout[b] += gb; fin[a] += gb
    tot_flux = sum(fin.values()) + sum(fout.values())
    sc.flux_imbal = (sum(abs(fin[c] - fout[c]) for c in cells)
                     / max(tot_flux, 1e-9))
    sc.nopen, sc.nbridge = bridges_and_open(cells, elist)
    sc.penalty = (PEN_OPEN * sc.nopen + PEN_BRIDGE * sc.nbridge
                  + PEN_FLUX * sc.flux_imbal * max(1, len(cells)))
    sc.leak_edge = leak_e
    sc.leak_par = leak_p
    sc.leak_vac = leak_v
    sc.npar = npar
    sc.gpar_max = gpmax
    sc.labels = labels
    sc.cycles = cinfo
    sc.par_pairs = par_pairs
    sc.per_voice = (leak_e + leak_p + leak_v) / max(1, len(cells))
    sc.total = sc.per_voice + sc.penalty / max(1, len(cells))
    return sc


def best_omega(foam, cells, edges, w_design=None, heavy_ok=False, span=0.07,
               drops_map=None):
    """Pitch selection: refine around the design pitch (ring_m-style retune,
    +-span), clamped to the seedable window and the measured-validity window
    unless heavy_ok. drops_map pins the winding basin during the scan."""
    wlo = W_OMEGA_MIN if heavy_ok else W_OMEGA_SAFE
    whi = W_OMEGA_MAX
    if w_design is None:
        grid = np.linspace(wlo, whi, 41)
    else:
        grid = w_design * np.linspace(1 - span, 1 + span, 29)
        grid = grid[(grid >= wlo - 1e-9) & (grid <= whi + 1e-9)]
        if len(grid) == 0:
            grid = np.array([min(max(w_design, wlo), whi)])
    best = None
    for w in grid:
        s = score_state(foam, cells, edges, w, drops_map=drops_map)
        if best is None or s.total < best.total:
            best = s
    if best is not None and best.ok:
        for w in np.linspace(best.omega - 0.01, best.omega + 0.01, 11):
            if wlo - 1e-9 <= w <= whi + 1e-9:
                s = score_state(foam, cells, edges, w, drops_map=drops_map)
                if s.total < best.total:
                    best = s
    return best


def design_basin(cand, cell_edges, ds):
    """For candidates carrying m_target: if every intended component is a
    simple cycle, return (omega_eff, drops_map) pinning each cycle's winding
    to m_target (ring_m semantics: omega = 2 pi m / L_actual)."""
    mt = cand.get("m_target")
    if not mt:
        return cand.get("omega"), None
    deg = {}
    for a, b in cell_edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    if any(v != 2 for v in deg.values()):
        return cand.get("omega"), None
    cells = sorted(deg)
    elist = [(min(a, b), max(a, b)) for a, b in cell_edges]
    comp, ncomp = components(cells, elist)
    Ls = {}
    for (a, b), d in zip(cell_edges, ds):
        Ls[comp[a]] = Ls.get(comp[a], 0.0) + d
    drops_map = {}
    for (a, b), d in zip(cell_edges, ds):
        L = Ls[comp[a]]
        drop = 2 * math.pi * mt * d / L
        # signed: design transport direction is the listed a -> b order
        drops_map[(min(a, b), max(a, b))] = drop if a < b else -drop
    w_eff = float(np.mean([2 * math.pi * mt / L for L in Ls.values()]))
    return w_eff, drops_map


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
    tgt = coords @ R.T + center
    nv = len(coords)
    picks = [-1] * nv
    used = set()
    for k in range(nv):
        _, ii = foam.tree.query(tgt[k], k=8)
        got = False
        for cand in np.atleast_1d(ii):
            if int(cand) not in used:
                picks[k] = int(cand)
                used.add(int(cand))
                got = True
                break
        if not got:
            return None, None
    d_ideal = {tuple(e): float(np.linalg.norm(coords[e[0]] - coords[e[1]]))
               for e in edges}
    inc = [[] for _ in range(nv)]
    for e in edges:
        inc[e[0]].append(tuple(e))
        inc[e[1]].append(tuple(e))
    for _ in range(refine):
        for k in range(nv):
            cands = foam.tree.query_ball_point(tgt[k], 1.35)
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
    return picks, [(picks[a], picks[b]) for a, b in edges]


def kabsch(P, Q):
    """Rigid transform (R, t) minimizing ||R P + t - Q||."""
    cp, cq = P.mean(0), Q.mean(0)
    H = (P - cp).T @ (Q - cq)
    U, S, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, dsign])
    R = Vt.T @ D @ U.T
    return R, cq - R @ cp


def embed_grow(foam, coords, edges, center, R0m, rng):
    """Sequential BFS growth: place vertices one at a time on foam cells,
    re-fitting the ideal frame to the placed cells (Kabsch) so local foam
    jitter does not accumulate. Handles closures better than rigid snap."""
    nv = len(coords)
    adj = [[] for _ in range(nv)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    root = int(rng.integers(0, nv))
    order = []
    seen = [False] * nv
    seen[root] = True
    queue = [root]
    while queue:
        u = queue.pop(0)
        order.append(u)
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True
                queue.append(v)
    if len(order) < nv:  # disconnected candidate: append the rest
        order += [k for k in range(nv) if not seen[k]]
    picks = {}
    used = set()
    Rm, tm = R0m, center - R0m @ coords[root]
    for step, v in enumerate(order):
        tgt = Rm @ coords[v] + tm
        placed_nb = [u for u in adj[v] if u in picks]
        cands = foam.tree.query_ball_point(tgt, 1.6 if placed_nb else 0.9)
        best_sc, best_c = None, None
        for cand in cands:
            cand = int(cand)
            if cand in used:
                continue
            s = 0.3 * float(np.sum((foam.xyz[cand] - tgt) ** 2))
            ok = True
            for u in placed_nb:
                lk, dd = foam.linked(picks[u], cand)
                if not lk:
                    ok = False
                    break
                di = float(np.linalg.norm(coords[u] - coords[v]))
                s += (dd - di) ** 2
            if not ok:
                continue
            if best_sc is None or s < best_sc:
                best_sc, best_c = s, cand
        if best_c is None:
            return None, None
        picks[v] = best_c
        used.add(best_c)
        if step >= 2 and (step % 3 == 0 or step == len(order) - 1):
            P = np.array([coords[u] for u in picks])
            Q = np.array([foam.xyz[picks[u]] for u in picks])
            try:
                Rm, tm = kabsch(P, Q)
            except np.linalg.LinAlgError:
                pass
    pl = [picks[k] for k in range(nv)]
    return pl, [(picks[a], picks[b]) for a, b in edges]


def embed_candidate(foam, cand, rng, tries=60, box_margin=4.0):
    coords, edges = cand["coords"], cand["edges"]
    if len(coords) > 16:
        tries = tries * 3
    lo, hi = box_margin, 24.0 - box_margin
    n_ok = 0
    spread = []
    best = None
    best_embed = None
    scales = (1.0, 1.05, 0.97)
    for t in range(tries):
        center = rng.uniform(lo, hi, 3)
        R = rand_rot(rng)
        sc_geom = scales[t % len(scales)]
        if t % 2 == 0 or not edges:
            picks, cell_edges = embed_once(foam, coords * sc_geom, edges,
                                           center, R)
        else:
            picks, cell_edges = embed_grow(foam, coords * sc_geom, edges,
                                           center, R, rng)
        if picks is None or len(set(picks)) < len(picks):
            continue
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
        if ds:
            spread.append(np.array(ds).std())
        w_eff, drops_map = design_basin(cand, cell_edges, ds)
        s = best_omega(foam, picks, cell_edges, w_design=w_eff,
                       heavy_ok=cand.get("heavy_ok", False),
                       span=0.04 if drops_map else 0.07, drops_map=drops_map)
        if s is not None and s.ok and (best is None or s.total < best.total):
            best = s
            best_embed = (picks, cell_edges)
    stats = {"tries": tries, "embed_ok": n_ok,
             "realizability": round(n_ok / tries, 3),
             "sigma_d_med": (round(float(np.median(spread)), 4) if spread else None)}
    return best, best_embed, stats


def paper_score(cand):
    """Score the candidate on its IDEAL geometry (r = r0 nominal): the clean
    parity/tuning number, independent of foam supply."""
    coords = np.asarray(cand["coords"], float)
    mini = Foam(coords, np.full(len(coords), R0))
    cells = list(range(len(coords)))
    ds = [float(np.linalg.norm(coords[a] - coords[b])) for a, b in cand["edges"]]
    w_eff, drops_map = design_basin(cand, cand["edges"], ds)
    s = best_omega(mini, cells, cand["edges"], w_design=w_eff,
                   heavy_ok=cand.get("heavy_ok", False),
                   span=0.04 if drops_map else 0.07, drops_map=drops_map)
    return s


# ---------------------------------------------------------------- catalog
def ring_coords(n, d, z=0.0):
    R = d / (2 * math.sin(math.pi / n))
    return np.array([[R * math.cos(2 * math.pi * k / n),
                      R * math.sin(2 * math.pi * k / n), z] for k in range(n)])


def g_ring(n, d):
    return {"coords": ring_coords(n, d),
            "edges": [(k, (k + 1) % n) for k in range(n)]}


def g_chain(n, d):
    c = np.array([[k * d, 0.15 * ((-1) ** k), 0] for k in range(n)])
    return {"coords": c, "edges": [(k, k + 1) for k in range(n - 1)]}


def g_blob(n, spacing=1.2):
    rng = np.random.default_rng(7)
    return {"coords": rng.normal(scale=spacing, size=(n, 3)), "edges": []}


def g_cube(a):
    c = np.array([[(k & 1), (k >> 1) & 1, (k >> 2) & 1] for k in range(8)],
                 float) * a
    e = [(i, i ^ (1 << b)) for i in range(8) for b in range(3)
         if i < (i ^ (1 << b))]
    return {"coords": c, "edges": e}


def g_tower(layers, a):
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
    c = np.vstack([ring_coords(ngon, d, 0.0), ring_coords(ngon, d, h)])
    e = [(k, (k + 1) % ngon) for k in range(ngon)]
    e += [(ngon + k, ngon + (k + 1) % ngon) for k in range(ngon)]
    e += [(k, ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_tube(ngon, layers, d_ring, d_ax):
    rings = [ring_coords(ngon, d_ring, z=i * d_ax) for i in range(layers)]
    c = np.vstack(rings)
    e = []
    for i in range(layers):
        e += [(i * ngon + k, i * ngon + (k + 1) % ngon) for k in range(ngon)]
    for i in range(layers - 1):
        e += [(i * ngon + k, (i + 1) * ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_trunc_oct(d):
    from itertools import permutations
    s = d / math.sqrt(2)
    pts = set()
    for p in permutations([0, 1, 2]):
        for s1 in (1, -1):
            for s2 in (1, -1):
                pts.add(tuple(np.array([p[0], p[1] * s1, p[2] * s2], float)))
    c = np.array(sorted(pts)) * s
    e = [(i, j) for i in range(len(c)) for j in range(i + 1, len(c))
         if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-6]
    return {"coords": c, "edges": e}


def g_rhombic_dodeca(d):
    s = d / math.sqrt(3)
    pts = [(sx, sy, sz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    pts += [(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)]
    c = np.array(pts, float) * s
    e = [(i, j) for i in range(len(c)) for j in range(i + 1, len(c))
         if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-6]
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
    c = np.array(v, float) * (d / 2.0)
    e = [(i, j) for i in range(12) for j in range(i + 1, 12)
         if abs(np.linalg.norm(c[i] - c[j]) - d) < 1e-4]
    return {"coords": c, "edges": e}


def g_antiprism(ngon, d):
    R = d / (2 * math.sin(math.pi / ngon))
    ch = 2 * R * math.sin(math.pi / (2 * ngon))
    h = math.sqrt(max(d * d - ch * ch, 1e-9))
    bot = ring_coords(ngon, d, 0.0)
    top = np.array([[R * math.cos(2 * math.pi * (k + 0.5) / ngon),
                     R * math.sin(2 * math.pi * (k + 0.5) / ngon), h]
                    for k in range(ngon)])
    c = np.vstack([bot, top])
    e = [(k, (k + 1) % ngon) for k in range(ngon)]
    e += [(ngon + k, ngon + (k + 1) % ngon) for k in range(ngon)]
    e += [(k, ngon + k) for k in range(ngon)]
    e += [((k + 1) % ngon, ngon + k) for k in range(ngon)]
    return {"coords": c, "edges": e}


def g_mobius(k, d):
    n = 2 * k
    w = d / 2.0
    R = d * k / (2 * math.pi) * 1.15
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
    e = set()
    for j in range(k):
        jn = (j + 1) % k
        a0, a1 = 2 * j, 2 * j + 1
        b0, b1 = (1, 0) if jn == 0 else (2 * jn, 2 * jn + 1)
        e.add(tuple(sorted((a0, b0))))
        e.add(tuple(sorted((a1, b1))))
        e.add(tuple(sorted((a0, a1))))
    return {"coords": c, "edges": sorted(e)}


def _arc_solve(na, d, D):
    """Circular arc of na chords of length d spanning straight distance D:
    returns (R, A) with span 2R sin(A/2) = D, chord 2R sin(A/(2na)) = d."""
    lo, hi = 1e-3, 2 * math.pi - 0.2
    for _ in range(80):
        A = 0.5 * (lo + hi)
        span = d * math.sin(A / 2) / math.sin(A / (2 * na))
        if span > D:
            lo = A
        else:
            hi = A
    A = 0.5 * (lo + hi)
    R = d / (2 * math.sin(A / (2 * na)))
    return R, A


def g_theta(arms, d):
    """Theta/book graph: two poles joined by len(arms) arcs (edge counts).
    Every arc spans the same pole distance; unequal arms bend more."""
    na_min = min(arms)
    D = 2 * (d / (2 * math.sin(math.pi / (2 * na_min)))) * math.sin(math.pi / 2)
    D *= 0.95
    coords = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, D])]
    edges = []
    for ai, na in enumerate(arms):
        phi = ai * 2 * math.pi / len(arms)
        R, A = _arc_solve(na, d, D)
        rot = np.array([[math.cos(phi), -math.sin(phi), 0],
                        [math.sin(phi), math.cos(phi), 0], [0, 0, 1]])
        # arc in the local xz plane from (0,0,0) to (0,0,D), bulging +x
        xc = -R * math.cos(A / 2)  # circle center x so ends at x=0
        prev = 0
        for s in range(1, na + 1):
            ang = -A / 2 + A * s / na
            local = np.array([xc + R * math.cos(ang), 0.0,
                              D / 2 + R * math.sin(ang) * (D / abs(D))])
            local = np.array([xc + R * math.cos(ang), 0.0,
                              D / 2 + R * math.sin(ang)])
            p = rot @ local
            if s == na:
                edges.append((prev, 1))
            else:
                coords.append(p)
                edges.append((prev, len(coords) - 1))
                prev = len(coords) - 1
    return {"coords": np.array(coords), "edges": edges}


def g_cube_in_cube(a):
    inner = g_cube(a)
    b = a + 2 * a / math.sqrt(3)
    outer = g_cube(b)
    c = np.vstack([inner["coords"] - a / 2, outer["coords"] - b / 2])
    e = list(inner["edges"]) + [(i + 8, j + 8) for i, j in outer["edges"]]
    e += [(k, k + 8) for k in range(8)]
    return {"coords": c, "edges": e}


def g_hopf_rings(n, d):
    R = d / (2 * math.sin(math.pi / n))
    A = ring_coords(n, d)
    rot = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], float)
    B = ring_coords(n, d) @ rot.T + np.array([R, 0, 0])
    c = np.vstack([A, B])
    e = [(k, (k + 1) % n) for k in range(n)]
    e += [(n + k, n + (k + 1) % n) for k in range(n)]
    return {"coords": c, "edges": e}


def catalog():
    """(name, class, generator, omega_design, heavy_ok, note)."""
    W15 = math.pi / 1.5    # 2.0944  (x=0.3205: the comp12-class load)
    W16 = math.pi / 1.6    # 1.9635  (x=0.3977)
    cat = []
    # anchors / calibration controls
    cat.append(("blob12", "control", g_blob(12), W15, False,
                "no intended edges; measured -0.232%/tu class"))
    cat.append(("chain12", "control", g_chain(12, 1.5), W15, False,
                "open chain; measured 1.7x worse than closed ring"))
    r = g_ring(12, 1.5); r["m_target"] = 6
    cat.append(("ring12_m6", "anchor", r, W15, False,
                "unwound strut ring m=6: died t=2221"))
    r = g_ring(12, 1.25); r["m_target"] = 5
    cat.append(("ring12_m5", "anchor", r, W15, False,
                "comp12 wound cable ring m=5: ALIVE t=5000"))
    r = g_ring(6, 1.10); r["m_target"] = 2
    cat.append(("ring6_m2", "anchor", r, 2 * math.pi / (3 * 1.10), False,
                "comp6 one-way m=2: died t=3836; check 2nd-nbr 1.905 parasites"))
    cat.append(("cube_a1.25", "anchor", g_cube(1.25), math.pi / 1.43, False,
                "H1 cube (mean-tuned 1.43): died ~1200-1800; diagonals 1.77 in-ceiling"))
    # bipartite strut shells
    cat.append(("cube_a1.5", "shell", g_cube(1.5), W15, False,
                "P2 fix: diagonals 2.12 past ceiling; edges at dbar"))
    cat.append(("tower2_a1.5", "shell", g_tower(2, 1.5), W15, False,
                "2-cube stack = capped square tube, n=12"))
    cat.append(("tower3_a1.5", "shell", g_tower(3, 1.5), W15, False,
                "3-cube stack, n=16"))
    cat.append(("hexprism_d1.5", "shell", g_prism(6, 1.5), W15, False,
                "n=12 3-regular bipartite"))
    cat.append(("hexprism_d1.6", "shell", g_prism(6, 1.6), W16, False,
                "n=12 at d=1.6: more diagonal clearance"))
    cat.append(("octprism_d1.5", "shell", g_prism(8, 1.5), W15, False,
                "n=16 3-regular bipartite"))
    cat.append(("truncoct_d1.5", "shell", g_trunc_oct(1.5), W15, False,
                "Kelvin cell n=24, bipartite 3-regular"))
    cat.append(("rhombdodec_d1.5", "shell", g_rhombic_dodeca(1.5), W15, False,
                "PREDICT frustrated: rhombus short diagonal 1.73 = 12 parasites"))
    cat.append(("cubeincube_a1.5", "shell", g_cube_in_cube(1.5), W15, False,
                "PREDICT broken: outer edges 3.23 beyond link ceiling"))
    # tubes
    cat.append(("tube6x3_d1.5", "tube", g_tube(6, 3, 1.5, 1.5), W15, False,
                "hex tube 3 rings n=18"))
    cat.append(("tube6x3_d1.6", "tube", g_tube(6, 3, 1.6, 1.6), W16, False,
                "hex tube at 1.6"))
    cat.append(("tube8x3_d1.5", "tube", g_tube(8, 3, 1.5, 1.5), W15, False,
                "oct tube 3 rings n=24"))
    cat.append(("tube6x4_d1.5", "tube", g_tube(6, 4, 1.5, 1.5), W15, False,
                "hex tube 4 rings n=24"))
    W17 = math.pi / 1.7   # 1.8480 (x = 0.4740)
    d_ring_w = (math.pi - 2 * math.pi / 12) / W17
    cat.append(("woundtube12x2_w1", "tube", g_tube(12, 2, d_ring_w, 1.7), W17, False,
                f"B1c co-rotating wound tube: d_ring={d_ring_w:.3f} (phi=5pi/6, m=5) "
                f"d_ax=1.7 (struts); diag {math.hypot(d_ring_w,1.7):.2f} clean"))
    # theta / book
    cat.append(("theta_666_d1.5", "theta", g_theta([6, 6, 6], 1.5), W15, False,
                "3 arms x 6 edges; junctions deg 3; cycles m=6"))
    cat.append(("theta_866_d1.5", "theta", g_theta([8, 6, 6], 1.5), W15, False,
                "unequal arms 8/6/6: cycle integers 7/6"))
    # negative controls
    cat.append(("octahedron_d1.5", "control", g_octahedron(1.5), W15, False,
                "NEG: odd cycles + equator 4-cycles jointly unclosable"))
    cat.append(("icosahedron_d1.5", "control", g_icosahedron(1.5), W15, False,
                "NEG: triangulated shell"))
    cat.append(("antiprism4_d1.5", "control", g_antiprism(4, 1.5), W15, False,
                "NEG: triangles; slant phi=pi/2 needs d=0.75 unreachable"))
    cat.append(("mobius6_d1.5", "control", g_mobius(6, 1.5), W15, False,
                "NEG: Mobius ladder k=6 (even) => odd cycles"))
    # wound / heavy rings (cable species)
    r = g_ring(9, 1.15); r["m_target"] = 3
    cat.append(("ring9_m3", "ring", r, 2 * math.pi / (3 * 1.15), False,
                "odd-N one-way ring phi=120deg: parity-free cables, back gate 1.5e-5"))
    r = g_ring(8, 1.30); r["m_target"] = 3
    cat.append(("ring8_m3", "ring", r, 3 * math.pi / (4 * 1.30), False,
                "phi=3pi/4 wound ring: back gate G(pi/2)=0.0039"))
    r = g_ring(10, 1.35); r["m_target"] = 4
    cat.append(("ring10_m4", "ring", r, 4 * math.pi / (5 * 1.35), False,
                "phi=4pi/5: back gate G(2pi/5)=0.0263"))
    r = g_ring(16, 1.5); r["m_target"] = 8
    cat.append(("ring16_m8", "ring", r, W15, False,
                "large strut ring n=16"))
    r = g_ring(8, 1.126); r["m_target"] = 2
    cat.append(("ring8_m2_x0.9", "ring", r, math.pi / (2 * 1.126), True,
                "UNMEASURED-HEAVY probe: phi=pi/2 at x~0.9 — the ZERO back-gate "
                "point G(pi)=0; needs short links d=1.126"))
    # odd strut ring: the cleanest parity negative control (m = 6.5 impossible)
    cat.append(("ring13_strut", "control", g_ring(13, 1.5), W15, False,
                "NEG: odd strut ring — closure needs m=6.5, the seam cannot close"))
    # Hopf link molecule
    r = g_hopf_rings(12, 1.25); r["m_target"] = 5
    cat.append(("hopf12x12_m5", "exotic", r, W15, False,
                "two comp12 rings Hopf-linked; purely topological bond"))
    return cat


# ---------------------------------------------------------------- free search
def free_search(foam, rng, omega_class, n_restarts=6, steps=1200,
                nmin=8, nmax=32, seed_states=None, log=None):
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
            cset = set(int(c) for c in seed_states[restart][0])
            edges = set((min(a, b), max(a, b)) for a, b in seed_states[restart][1])
        else:
            cset, edges = rand_seed_state()
        cur = state_score(cset, edges)
        if cur is None:
            continue
        best = cur
        best_state = (frozenset(cset), frozenset(edges))
        T0, T1 = 0.08, 0.003
        for step in range(steps):
            T = T0 * (T1 / T0) ** (step / max(1, steps - 1))
            move = rng.random()
            nc, ne = set(cset), set(edges)
            if move < 0.28 and len(nc) < nmax + 2:
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
            elif move < 0.52 and len(nc) > nmin - 2:
                degs = {}
                for a, b in ne:
                    degs[a] = degs.get(a, 0) + 1
                    degs[b] = degs.get(b, 0) + 1
                cl = sorted(nc, key=lambda c: degs.get(c, 0))
                victim = cl[rng.integers(0, min(4, len(cl)))]
                nc.discard(victim)
                ne = set(e for e in ne if victim not in e)
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
            elif move < 0.60:
                # close-an-end move: wire the lowest-degree vertex to a linked
                # set member it is not yet wired to
                degs = {c: 0 for c in nc}
                for a, b in ne:
                    degs[a] += 1
                    degs[b] += 1
                opens = sorted(nc, key=lambda c: degs[c])
                v = opens[0]
                options = [(abs(d - math.pi / omega_class), j, d)
                           for j, d in foam.adj[v] if j in nc
                           and (min(v, j), max(v, j)) not in ne]
                if not options:
                    continue
                options.sort()
                _, j, d = options[0]
                ne.add((min(v, j), max(v, j)))
            else:
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
        if log:
            log(f"    restart {restart}: leak/voice={best.per_voice:.4f} "
                f"total={best.total:.4f} n={best.n} ne={best.ne} "
                f"npar={best.npar} open={best.nopen} bridges={best.nbridge}")
    return results


# ---------------------------------------------------------------- exotics
def exotic_fifth_bridge():
    out = {}
    out["comment"] = ("3:2 is the ONLY non-unison interval whose two pitches "
                      "both fit the seedable window x in [0.25,0.9]; octave "
                      "needs ratio 2.0 > 2.9/(1.3)/ (2.9/2.08) = 1.60 max")
    w_lo = 1.45
    w_hi = 1.5 * w_lo
    out["example"] = {"w_low": w_lo, "w_high": w_hi,
                      "x_low": round(x_of_omega(w_lo), 4),
                      "x_high": round(x_of_omega(w_hi), 4)}
    Om = 2 * w_hi  # coincident partial q*w_hi = p*w_lo
    out["bridge_d_m1"] = round(math.pi / Om, 4)
    out["bridge_d_m2"] = round(2 * math.pi / Om, 4)
    gw = GAMMA_M / 6.0
    ddet_dxhi = 2 * QDET * w_hi ** 2 / W2
    ddet_dxlo = 3 * QDET * w_lo ** 2 / W2
    out["acceptance_det"] = round(gw, 5)
    out["dx_tolerance_rms"] = round(gw / math.hypot(ddet_dxhi, ddet_dxlo), 5)
    out["res_on_tooth"] = round(1.0 / 6.0, 4)
    out["verdict"] = ("bridge struts exist at d=m2 rung ~1.44 (in-window) but "
                      "lock tolerance is ~0.4% in load vs ~30% foam chaos: "
                      "molecule bonds are 6x weaker and ~100x more fragile "
                      "than unison links — paper-only until a unison particle "
                      "is certified")
    return out


def exotic_hopf_geometry(n=12, d=1.25):
    R = d / (2 * math.sin(math.pi / n))
    t = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    A = np.stack([R * np.cos(t), R * np.sin(t), 0 * t], 1)
    B = np.stack([R * np.cos(t) + R, 0 * t, R * np.sin(t)], 1)
    dmin = float(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)))
    return {"R": round(R, 4), "min_dist": round(dmin, 4),
            "worst_link_ceiling": 2.07,
            "clean": bool(dmin > 2.07),
            "verdict": "two m=5 rings interlock with min approach %.2f > ceiling: "
                       "zero channels between them — a purely topological bond" % dmin}


# ---------------------------------------------------------------- output
def edge_records(sc):
    return [{"i": int(a), "j": int(b), "label": lab,
             "psi_f": round(pf, 4), "psi_b": round(pb, 4),
             "gf": round(gf, 4), "gb": round(gb, 4)}
            for (a, b), lab, (pf, pb), (gf, gb) in sc.labels]


def write_net(path, foam, sc, cells):
    idx = {c: k for k, c in enumerate(cells)}
    with open(path, "w") as f:
        f.write(f"# morpho export: omega={sc.omega:.5f} x={sc.x:.5f} "
                f"leak/voice={sc.per_voice:.4f} npar={sc.npar} wind={sc.wind}\n")
        for c in cells:
            th = sc.theta.get(c, 0.0) % (2 * math.pi)
            p = foam.xyz[c]
            f.write(f"V {p[0]:.4f} {p[1]:.4f} {p[2]:.4f} {sc.x:.6f} {th:.6f}\n")
        for (a, b), lab, _, _ in sc.labels:
            f.write(f"E {idx[a]} {idx[b]}\n")


def row_from_score(name, klass, note, n_ideal, ne_ideal, bip, s, cells,
                   stats, paper, foam):
    labs = [l for _, l, _, _ in s.labels]
    row = {"name": name, "class": klass, "note": note, "n": n_ideal,
           "ne": ne_ideal, "bipartite": bip, "stats": stats,
           "verdict": "scored", "omega": round(s.omega, 4),
           "x": round(s.x, 4), "x_margin_skirt": round(s.x / X_SKIRT, 1),
           "mass": s.mass,
           "x_unmeasured": bool(s.x > X_MEASURED_MAX + 1e-9),
           "leak_per_voice": round(s.per_voice, 4),
           "leak_edge": round(s.leak_edge, 4),
           "leak_par": round(s.leak_par, 4),
           "leak_vac": round(s.leak_vac, 4),
           "npar": s.npar, "gpar_max": round(s.gpar_max, 4),
           "nopen": s.nopen, "nbridge": s.nbridge, "wind": s.wind,
           "flux_imbal": round(s.flux_imbal, 4),
           "n_strut": labs.count("strut"),
           "n_cable": labs.count("cable_fwd") + labs.count("cable_bwd"),
           "cycles": s.cycles,
           "cells": [int(c) for c in cells],
           "positions": [[round(float(v), 4) for v in foam.xyz[c]] for c in cells],
           "edges": edge_records(s),
           "parasites": s.par_pairs,
           "penalty": round(s.penalty, 3), "total": round(s.total, 4)}
    if paper is not None and paper.ok:
        row["paper_leak_per_voice"] = round(paper.per_voice, 4)
        row["paper_npar"] = paper.npar
    row["score_obj"] = s
    return row


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
    steps = args.steps or (250 if args.quick else 1400)

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    print(f"# loading foam {os.path.normpath(FOAM)}")
    foam = Foam.load(FOAM)
    print(f"# foam: {foam.n} cells, {len(foam.pairs)} links, "
          f"dbar={foam.pdist.mean():.4f}, deg mean={foam.deg.mean():.1f}, "
          f"link window [{foam.pdist.min():.3f},{foam.pdist.max():.3f}]")

    ledger = []
    results = {}

    print("\n# === SWEEP 1: catalog ===")
    for name, klass, cand, w_design, heavy_ok, note in catalog():
        cand["omega"] = w_design
        cand["heavy_ok"] = heavy_ok
        bip = is_bipartite(range(len(cand["coords"])), cand["edges"])
        paper = paper_score(cand) if cand["edges"] else None
        best, emb, stats = embed_candidate(foam, cand, rng, tries=tries)
        if best is None or not best.ok:
            verdict = "NOT-SUPPLIED" if stats["embed_ok"] == 0 else "NO-SCORE"
            row = {"name": name, "class": klass, "note": note,
                   "n": len(cand["coords"]), "ne": len(cand["edges"]),
                   "bipartite": bip, "stats": stats, "verdict": verdict}
            if paper is not None and paper.ok:
                row["paper_leak_per_voice"] = round(paper.per_voice, 4)
                row["paper_npar"] = paper.npar
            print(f"  {name:22s} [{klass}] n={row['n']:3d} -> {verdict} "
                  f"(embeds {stats['embed_ok']}/{stats['tries']}, paper="
                  f"{row.get('paper_leak_per_voice', 'n/a')})")
        else:
            row = row_from_score(name, klass, note, len(cand["coords"]),
                                 len(cand["edges"]), bip, best, emb[0],
                                 stats, paper, foam)
            print(f"  {name:22s} [{klass}] n={row['n']:3d} leak/voice="
                  f"{best.per_voice:.4f} (e={best.leak_edge:.3f} "
                  f"p={best.leak_par:.3f} v={best.leak_vac:.3f}) "
                  f"npar={best.npar} x={best.x:.3f} w={best.omega:.4f} "
                  f"s/c={row['n_strut']}/{row['n_cable']} wind={best.wind} "
                  f"realiz={stats['realizability']:.2f} "
                  f"paper={row.get('paper_leak_per_voice', 'n/a')}")
        ledger.append(row)
        results[name] = row

    print("\n# === SWEEP 2: free search (backtracer) ===")
    free_finds = []
    omega_classes = [2.0944, 1.9635, 2.2308, 1.8000]
    cat_sets = set()
    seed_from = []
    for nm in ("tube6x3_d1.5", "ring12_m5", "truncoct_d1.5", "woundtube12x2_w1"):
        r = results.get(nm)
        if r and r.get("verdict") == "scored":
            cells = r["cells"]
            edges = [(e["i"], e["j"]) for e in r["edges"]]
            seed_from.append((set(cells), set((min(a, b), max(a, b))
                                              for a, b in edges)))
            cat_sets.add(frozenset(cells))
    for wc in omega_classes:
        print(f"  omega class {wc:.4f} (d*={math.pi/wc:.3f}, x={x_of_omega(wc):.3f})")
        res = free_search(foam, rng, wc, n_restarts=restarts, steps=steps,
                          seed_states=seed_from if wc == 2.0944 else None,
                          log=print)
        for best, (cset, eset) in res:
            free_finds.append((best, sorted(cset), sorted(eset), wc))

    free_finds.sort(key=lambda t: t[0].total)
    polished = []
    seen_sets = set(cat_sets)
    for best, cells, edges, wc in free_finds[:14]:
        key = frozenset(cells)
        if key in seen_sets:
            continue
        seen_sets.add(key)
        s = best_omega(foam, cells, edges, w_design=wc, heavy_ok=False, span=0.10)
        if not (s and s.ok):
            continue
        # greedy closure pass: add internal links that lower the total
        cur_edges = set(edges)
        improved = True
        rounds = 0
        while improved and rounds < 6 and (s.nopen or s.nbridge):
            improved = False
            rounds += 1
            for i, j, d in foam.links_within(cells):
                key2 = (min(i, j), max(i, j))
                if key2 in cur_edges:
                    continue
                s2 = best_omega(foam, cells, sorted(cur_edges | {key2}),
                                w_design=s.omega, span=0.03)
                if s2 and s2.ok and s2.total < s.total - 1e-4:
                    s = s2
                    cur_edges.add(key2)
                    improved = True
                    break
        polished.append((s, cells, sorted(cur_edges)))
    polished.sort(key=lambda t: t[0].total)
    for k, (s, cells, edges) in enumerate(polished[:8]):
        name = f"free_{k}"
        row = row_from_score(name, "free-search", "backtracer find", s.n, s.ne,
                             is_bipartite(cells, edges), s, cells,
                             {"tries": 1, "embed_ok": 1, "realizability": 1.0,
                              "sigma_d_med": None}, None, foam)
        ledger.append(row)
        results[name] = row
        print(f"  {name}: n={s.n} ne={s.ne} leak/voice={s.per_voice:.4f} "
              f"npar={s.npar} open={s.nopen} bridges={s.nbridge} "
              f"wind={s.wind} x={s.x:.3f}")

    print("\n# === SWEEP 3: exotics (paper checks) ===")
    ex = {"fifth_bridge": exotic_fifth_bridge(),
          "hopf_geometry": exotic_hopf_geometry(),
          "beat_locked_pairs": {
              "verdict": "NOT EXPRESSIBLE as a static lock: pitch w is an "
              "instantaneous function of load (w2e = w2/(1+q*(Em+flload)/cap), "
              "cellfab.c:1715) — there is no independent frequency dynamics to "
              "entrain, so a pair with dw != 0 beats forever; conversions are "
              "atomized (quant_mode=2 credit, atoms_fire cellfab.c:1558) but "
              "the credit ledger has no phase, so nothing couples the beat "
              "phase to the atom cadence. A stroboscopic 'lock' would be a "
              "limit cycle of kappa_lock phase kicks crossing the dead-gate "
              "desert every beat — a periodic unlock/refire radiator, not a "
              "particle. The only stable dw != 0 relations are the comb's p:q "
              "locks (lp/lq, cellfab.c:1787-1798); within the seedable pitch "
              "window that means the fifth (3:2) or nothing."},
          "wound_double_tori": {
              "verdict": "true quad-mesh tori need major radius > minor+2.0 "
              "clearance => l >= ~12 rings => n >= ~72 cells at k=6: beyond "
              "the n<=32 mandate (and 2 linked tori >= ~150). The reachable "
              "linked-homology object at n<=32 is the Hopf pair of wound "
              "RINGS (see hopf_geometry: min approach 2.41 > ceiling, zero "
              "cross-channels). Sim-worthy: bond is purely topological — "
              "separation requires crossing = a parasite storm."}}
    print(f"  fifth bridge: d_m2={ex['fifth_bridge']['bridge_d_m2']} "
          f"dx_tol={ex['fifth_bridge']['dx_tolerance_rms']}")
    print(f"  hopf: min_dist={ex['hopf_geometry']['min_dist']} "
          f"clean={ex['hopf_geometry']['clean']}")

    print("\n# === SHORTLIST ===")
    scored = [r for r in ledger if r.get("verdict") == "scored"
              and r["class"] not in ("control",)
              and not r.get("x_unmeasured", False)
              and r.get("x", 0) >= X_SEED_MIN - 1e-6
              and r.get("nopen", 9) == 0 and r.get("nbridge", 9) == 0]
    scored.sort(key=lambda r: r["leak_per_voice"])
    shortlist = []
    have = set()
    for r in scored:
        if len(shortlist) >= 8:
            break
        if r["name"] in have:
            continue
        shortlist.append(r)
        have.add(r["name"])
    # diversity guarantee: >= 2 wound (wind > 0.3) candidates, best-first,
    # never importing anything grossly worse than the list it joins
    nw = sum(1 for r in shortlist if r.get("wind", 0) > 0.3)
    if nw < 2 and shortlist:
        cutoff = 4 * shortlist[-1]["leak_per_voice"] + 0.05
        wound = [r for r in scored if r.get("wind", 0) > 0.3
                 and r["name"] not in have and r["leak_per_voice"] < cutoff]
        while nw < 2 and wound:
            r = wound.pop(0)
            shortlist[-1] = r
            have.add(r["name"])
            nw += 1
    negatives = [results.get("mobius6_d1.5"), results.get("octahedron_d1.5")]

    out = {"laws": {"w2": W2, "q_detune": QDET, "cap": CAP, "p_gate": PGATE,
                    "gamma_res_m": GAMMA_M, "x_skirt": X_SKIRT,
                    "seed_window_x": [X_SEED_MIN, X_MAX],
                    "validity_window_x": [X_SEED_MIN, X_MEASURED_MAX],
                    "omega_window": [W_OMEGA_MIN, W_OMEGA_MAX]},
           "score_weights": {"W_BACK": W_BACK, "W_PAR": W_PAR, "W_VAC": W_VAC,
                             "G_RAND": G_RAND, "PEN_OPEN": PEN_OPEN,
                             "PEN_BRIDGE": PEN_BRIDGE},
           "shortlist": [], "negative_controls": [], "exotics": ex}

    def strip(r):
        return {k: v for k, v in r.items() if k != "score_obj"}

    for rank, r in enumerate(shortlist):
        rr = strip(r)
        rr["rank"] = rank + 1
        out["shortlist"].append(rr)
        print(f"  #{rank+1} {r['name']:22s} leak/voice={r['leak_per_voice']:.4f} "
              f"x={r['x']:.3f} npar={r['npar']} wind={r.get('wind', 0)} "
              f"class={r['class']}")
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

    with open(os.path.join(HERE, "ledger.tsv"), "w") as f:
        cols = ("name", "class", "n", "verdict", "bipartite", "omega", "x",
                "leak_per_voice", "leak_edge", "leak_par", "leak_vac", "npar",
                "gpar_max", "nopen", "nbridge", "wind", "n_strut", "n_cable",
                "paper_leak_per_voice")
        f.write("\t".join(cols) + "\trealizability\tnote\n")
        for r in ledger:
            f.write("\t".join(str(r.get(k, "")) for k in cols) +
                    f"\t{r['stats']['realizability'] if r.get('stats') else ''}"
                    f"\t{r.get('note','')}\n")

    with open(os.path.join(HERE, "report.json"), "w") as f:
        json.dump({"ledger": [strip(r) for r in ledger], "exotics": ex},
                  f, indent=1)
    print(f"\n# done in {time.time()-t0:.1f}s -> shortlist.json, ledger.tsv, "
          f"report.json, nets/")


if __name__ == "__main__":
    main()
