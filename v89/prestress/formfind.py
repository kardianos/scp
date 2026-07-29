#!/usr/bin/env python3
# =========================================================================
# formfind.py — v89 PRESTRESS stream B: form-finding solver on the real foam
#
# Pick cells, tune (omega, loads, theta, m), predict gates, emit kernel-ready
# .net seed files for cellfab init=net.  PURE PYTHON (numpy/scipy).  Never
# runs the simulator.  Works only from v89 material.
#
# Physics (verified against v89/cellfab.c at 2026-07-28 state):
#   gate_of  (cellfab.c ~530):  gate(psi) = ((1+cos psi)/2)^p_gate, p_gate=8
#   forward gate of link i->j:  psi_f = wrap(th_i - w_i*d/C - th_j)
#   backward:                   psi_b = wrap(th_j - w_j*d/C - th_i)
#   pitch from load:            w = w2/(1+q_detune*x);  x = (w2/w-1)/q_detune
#   seeded load:                Em = x*cap   (add + s_pull mechanics, net x*cap)
#   link rule (foam build):     linked iff d < 1.15*(r_i + r_j)
#   net reader (init=net):      V x y z xk th2 / E a b ; '#' comments ignored;
#                               xk clamped >= 0.02; th2 = fmod(th+16pi, 2pi);
#                               picks = nearest free cell, id order, strict <
#   ring seeder (init=7, mimic lines ~916-991), shell seeder (init=8, ~992+)
#
# Design algebra (LEDGER.md translation):
#   strut  = both-gate link: needs 2*w*d = 2pi*m, only m=1 reachable
#            => d = pi/w shared per component; optimal phases put
#            th_i - th_j = pi (antiphase), each gate at psi = pi - w*d.
#   cable  = forward-only link i->j: th_i - th_j = w*d exactly on a tree;
#            every independent cycle needs sum of oriented targets = 0 mod 2pi
#            (strut contributes pi, cable +-w*d)  — ring_m generalized.
#            back gate of a cable = gate(-2*w*d): the slow back-leak.
#   parasite = unintended pair of picked cells inside the link ceiling.
#
# Determinism: every stochastic step uses np.random.default_rng with fixed
# seeds derived from FSEED; identical runs give identical bytes.
# =========================================================================

import argparse
import math
import os
import sys

import numpy as np

try:
    from scipy.optimize import minimize as _scipy_minimize
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

TWO_PI = 2.0 * math.pi
FSEED = 20260727

# ---- the law table (v89/battery/laws_V2g.cfg) ---------------------------
LAW = dict(C=1.0, w1=1.65, w2=2.9, q_detune=1.2, p_gate=8, cap=2.5,
           s_pull=0.15, es_floor=0.05, e_s0=1.0, r0=0.85, rjit=0.06,
           dmin=1.0, L=24.0)
X_HARD = (0.05, 0.90)        # usable load window
X_SAFE = 0.25                # skirt margin floor (x_skirt = 0.0617)
X_SKIRT = 0.0617
BOUND_MARGIN = 2.8           # picked cells stay this far from box faces
                             # (edge sinks live within 1.6)

DEF_FOAM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "foam", "foam_s20260727.tsv")

# ---- reference numbers from the validation logs -------------------------
REF_RING = dict(Lring=16.542, closure=5.0000, n=12, d_target=1.25, m=5)
REF_SHELL = dict(abar=1.586, omega=1.9809, x=0.3866, gmin=0.001, gmean=0.597,
                 a_target=1.25, r_core=0.8)


# ------------------------------------------------------------------------
# primitives
# ------------------------------------------------------------------------

def wrap(a):
    """wrap_pi, kernel-identical semantics (vectorized)."""
    return (np.asarray(a) + math.pi) % TWO_PI - math.pi


def gate(psi):
    h = 0.5 * (1.0 + np.cos(psi))
    return h ** LAW["p_gate"]


def dgate(psi):
    p = LAW["p_gate"]
    h = 0.5 * (1.0 + np.cos(psi))
    return p * h ** (p - 1) * (-0.5 * np.sin(psi))


def d2gate(psi):
    p = LAW["p_gate"]
    h = 0.5 * (1.0 + np.cos(psi))
    hp = -0.5 * np.sin(psi)
    hpp = -0.5 * np.cos(psi)
    return p * (p - 1) * h ** (p - 2) * hp * hp + p * h ** (p - 1) * hpp


def omega_of_x(x):
    return LAW["w2"] / (1.0 + LAW["q_detune"] * x)


def x_of_omega(w):
    return (LAW["w2"] / w - 1.0) / LAW["q_detune"]


OMEGA_HARD = (omega_of_x(X_HARD[1]), omega_of_x(X_HARD[0]))   # (1.394, 2.736)
OMEGA_SAFE_HI = omega_of_x(X_SAFE)                            # 2.2308


# ------------------------------------------------------------------------
# foam
# ------------------------------------------------------------------------

class Foam:
    def __init__(self, path):
        raw = np.loadtxt(path, skiprows=1)
        self.path = path
        self.n = raw.shape[0]
        self.pos = raw[:, 1:4].copy()          # x y z
        self.r = raw[:, 4].copy()              # nominal radius cr0
        self._build_links()

    def _build_links(self):
        """Exact kernel rule: linked iff d < 1.15*(r_i+r_j).  Grid binned."""
        L = LAW["L"]
        rmax = self.r.max()
        hcut = 1.15 * 2.0 * rmax
        g = max(1, int(math.ceil(L / hcut)))
        bx = np.minimum((self.pos[:, 0] / L * g).astype(int), g - 1)
        by = np.minimum((self.pos[:, 1] / L * g).astype(int), g - 1)
        bz = np.minimum((self.pos[:, 2] / L * g).astype(int), g - 1)
        bins = {}
        for i in range(self.n):
            bins.setdefault((bx[i], by[i], bz[i]), []).append(i)
        li, lj, ld = [], [], []
        for i in range(self.n):
            b = (bx[i], by[i], bz[i])
            for ax in range(b[0] - 1, b[0] + 2):
                for ay in range(b[1] - 1, b[1] + 2):
                    for az in range(b[2] - 1, b[2] + 2):
                        for q in bins.get((ax, ay, az), ()):
                            if q <= i:
                                continue
                            dv = self.pos[q] - self.pos[i]
                            d2 = float(dv @ dv)
                            cut = 1.15 * (self.r[i] + self.r[q])
                            if d2 < cut * cut:
                                li.append(i); lj.append(q)
                                ld.append(math.sqrt(d2))
        self.li = np.array(li); self.lj = np.array(lj); self.ld = np.array(ld)
        self.nl = len(ld)
        self.adj = {}
        for a, b, d in zip(self.li, self.lj, self.ld):
            self.adj[(int(a), int(b))] = float(d)
            self.adj[(int(b), int(a))] = float(d)
        self.dbar = float(self.ld.mean())
        self.mean_degree = 2.0 * self.nl / self.n

    def linked(self, a, b):
        return (a, b) in self.adj

    def link_d(self, a, b):
        return self.adj.get((a, b))

    def cut(self, a, b):
        return 1.15 * (self.r[a] + self.r[b])

    def nearest_free(self, target, used, excl_center=None, excl_r=0.0):
        """Kernel pick loop: ascending id, strict <, skip used / core ball."""
        d2 = np.sum((self.pos - np.asarray(target)) ** 2, axis=1)
        bad = np.zeros(self.n, bool)
        if used:
            bad[list(used)] = True
        if excl_center is not None and excl_r > 0:
            c2 = np.sum((self.pos - np.asarray(excl_center)) ** 2, axis=1)
            bad |= c2 < excl_r * excl_r
        d2 = np.where(bad, np.inf, d2)
        return int(np.argmin(d2))     # first minimum == kernel strict <

    def margin_ok(self, i):
        p = self.pos[i]
        return (p.min() >= BOUND_MARGIN) and (p.max() <= LAW["L"] - BOUND_MARGIN)


# ------------------------------------------------------------------------
# kernel seeder mimics (validation)
# ------------------------------------------------------------------------

def kernel_ring_pick(foam, n=12, ring_d=1.25, center=None):
    """cellfab.c init=7 pick loop, lines ~916-934, replicated exactly."""
    if center is None:
        center = (0.5 * LAW["L"], 0.5 * LAW["L"], 0.5 * LAW["L"])
    R = ring_d / (2.0 * math.sin(math.pi / n))
    picks = []
    for k in range(n):
        ph = TWO_PI * k / n
        t = (center[0] + R * math.cos(ph), center[1] + R * math.sin(ph),
             center[2])
        picks.append(foam.nearest_free(t, set(picks)))
    P = foam.pos[picks]
    Lring = float(sum(np.linalg.norm(P[(k + 1) % n] - P[k])
                      for k in range(n)))
    return picks, Lring, R


def kernel_shell_pick(foam, a_target=1.25, center=None, r_core=0.8):
    """cellfab.c init=8 cube picks + 3 edge-uniformity refine passes."""
    if center is None:
        center = (0.5 * LAW["L"], 0.5 * LAW["L"], 0.5 * LAW["L"])
    cc = np.asarray(center)
    h = 0.5 * a_target
    corners = [cc + np.array([h if k & 1 else -h,
                              h if k & 2 else -h,
                              h if k & 4 else -h]) for k in range(8)]
    picks = []
    for k in range(8):
        picks.append(foam.nearest_free(corners[k], set(picks),
                                       excl_center=cc, excl_r=r_core))
    picks = list(picks)
    # refinement: 3 passes, sequential in k, search radius sqrt(1.44)
    core2 = r_core * r_core
    c2all = np.sum((foam.pos - cc) ** 2, axis=1)
    for _pass in range(3):
        for k in range(8):
            t = corners[k]
            d2t = np.sum((foam.pos - t) ** 2, axis=1)
            elig = d2t <= 1.44
            if r_core > 0:
                elig &= c2all >= core2
            for q in range(8):
                if q != k:
                    elig[picks[q]] = False
            best_sc, best_c = 1e30, picks[k]
            idx = np.nonzero(elig)[0]
            if len(idx):
                nb = [picks[k ^ (1 << b)] for b in range(3)]
                sc = np.zeros(len(idx))
                for pkn in nb:
                    dd = np.linalg.norm(foam.pos[idx] - foam.pos[pkn], axis=1)
                    sc += (dd - a_target) ** 2
                j = int(np.argmin(sc))
                if sc[j] < best_sc:
                    best_sc, best_c = float(sc[j]), int(idx[j])
            picks[k] = best_c
    edges = [(k, k ^ (1 << b)) for k in range(8) for b in range(3)
             if (k ^ (1 << b)) > k]
    dl = [float(np.linalg.norm(foam.pos[picks[b]] - foam.pos[picks[a]]))
          for a, b in edges]
    abar = float(np.mean(dl))
    return picks, edges, dl, abar


def kernel_shell_bfs(foam, picks, om):
    """Kernel BFS phases (actual lengths) + the kernel 12-edge gate report."""
    th = {0: 0.0}
    seen = [False] * 8
    seen[0] = True
    q = [0]
    while q:
        k = q.pop(0)
        for b in range(3):
            k2 = k ^ (1 << b)
            if seen[k2]:
                continue
            dd = float(np.linalg.norm(foam.pos[picks[k2]] - foam.pos[picks[k]]))
            th[k2] = math.fmod(th[k] - om * dd / LAW["C"] + 8.0 * TWO_PI,
                               TWO_PI)
            seen[k2] = True
            q.append(k2)
    gates = []
    for k in range(8):
        for b in range(3):
            k2 = k ^ (1 << b)
            if k2 < k:
                continue
            dd = float(np.linalg.norm(foam.pos[picks[k2]]
                                      - foam.pos[picks[k]]))
            ps = wrap(th[k] - om * dd / LAW["C"] - th[k2])
            gates.append(float(gate(ps)))
    return th, gates


# ------------------------------------------------------------------------
# designs
# ------------------------------------------------------------------------

class Design:
    """verts: (K,3) template offsets about a center; edges: (a,b,type)
    with type 'S' (strut, both gates) or 'C' (cable, forward a->b)."""

    def __init__(self, name, verts, edges, band="strut", band_cycle=None,
                 band_m=None, pool_r=1.5, notes="", pinned=None):
        self.name = name
        self.verts = np.asarray(verts, float)
        self.edges = edges
        self.band = band              # 'strut' | 'cycle' | 'free' | ('fix',w)
        self.band_cycle = band_cycle  # vertex loop for band='cycle'
        self.band_m = band_m
        self.pool_r = pool_r
        self.notes = notes
        self.pinned = pinned or []    # [(vertex_loop, m)] forced closure ints
        self.nv = len(self.verts)
        self.ne = len(edges)
        self.ext = float(np.max(np.linalg.norm(self.verts, axis=1)))

    def struts(self):
        return [e for e in self.edges if e[2] == "S"]

    def cables(self):
        return [e for e in self.edges if e[2] == "C"]


def cycle_basis(design):
    """Fundamental cycles from a BFS forest (roots in vertex order).
    Returns (cycles, parent, order, roots); cycle = [(edge_index, sign)]."""
    adj = {}
    for idx, (a, b, _t) in enumerate(design.edges):
        adj.setdefault(a, []).append((b, idx, +1))
        adj.setdefault(b, []).append((a, idx, -1))
    parent = {}
    order = []
    roots = []
    tree = set()
    for r0 in range(design.nv):
        if r0 in parent:
            continue
        parent[r0] = None
        roots.append(r0)
        order.append(r0)
        q = [r0]
        while q:
            u = q.pop(0)
            for (v, idx, sg) in adj.get(u, ()):
                if v in parent:
                    continue
                parent[v] = (u, idx, sg)
                tree.add(idx)
                order.append(v)
                q.append(v)
    def path_to_root(v):
        out = []
        while parent[v] is not None:
            u, idx, sg = parent[v]
            out.append((idx, sg, u, v))
            v = u
        return out
    cycles = []
    for idx, (a, b, _t) in enumerate(design.edges):
        if idx in tree:
            continue
        pa = path_to_root(a)
        pb = path_to_root(b)
        sa = {e[0] for e in pa}
        sb = {e[0] for e in pb}
        common = sa & sb
        cyc = [(idx, +1)]                       # a -> b via the edge
        for (i2, sg, _u, _v) in pb:             # b -> root, keep uncommon
            if i2 not in common:
                cyc.append((i2, -sg))           # traversed v->u
        for (i2, sg, _u, _v) in reversed([e for e in pa if e[0] not in common]):
            cyc.append((i2, sg))
        cycles.append(cyc)
    return cycles, parent, order, roots


def edge_targets(design, dlen, omega):
    """Oriented lock target per edge: strut pi, cable omega*d."""
    t = np.empty(design.ne)
    for i, (_a, _b, ty) in enumerate(design.edges):
        t[i] = math.pi if ty == "S" else omega * dlen[i] / LAW["C"]
    return t


def cycle_defects(design, cycles, dlen, omega):
    out = []
    for cyc in cycles:
        s = 0.0
        for idx, sg in cyc:
            ty = design.edges[idx][2]
            mu = math.pi if ty == "S" else omega * dlen[idx] / LAW["C"]
            s += sg * mu
        out.append(float(wrap(s)))
    return out


# ------------------------------------------------------------------------
# phase solver
# ------------------------------------------------------------------------

def solve_phases(design, dlen, omega, w_strut=2.0, w_cable=1.0, iters=6):
    """Tree init + wrap-iterated weighted least squares (lock Laplacian).
    Integers per cycle fixed by rounding the tree solution (wrap).
    Handles a BFS forest (disconnected intended graphs): gauge pin per
    component root."""
    nv, ne = design.nv, design.ne
    targ = edge_targets(design, dlen, omega)
    _cycles, parent, order, roots = cycle_basis(design)
    th = np.zeros(nv)
    for v in order:
        if parent[v] is None:
            continue
        u, idx, sg = parent[v]
        # sg=+1: edge stored (u,v): th_u - th_v = targ  => th_v = th_u - targ
        # sg=-1: edge stored (v,u): th_v - th_u = targ  => th_v = th_u + targ
        th[v] = th[u] - sg * targ[idx]
    W = np.array([w_strut if e[2] == "S" else w_cable for e in design.edges])
    A = np.zeros((nv, nv))
    for i, (a, b, _t) in enumerate(design.edges):
        A[a, a] += W[i]; A[b, b] += W[i]
        A[a, b] -= W[i]; A[b, a] -= W[i]
    for r0 in roots:
        A[r0, :] = 0.0; A[:, r0] = 0.0; A[r0, r0] = 1.0   # gauge pins
    for _ in range(iters):
        r = np.zeros(ne)
        for i, (a, b, _t) in enumerate(design.edges):
            r[i] = wrap(th[a] - th[b] - targ[i])
        g = np.zeros(nv)
        for i, (a, b, _t) in enumerate(design.edges):
            g[a] -= W[i] * r[i]
            g[b] += W[i] * r[i]
        for r0 in roots:
            g[r0] = 0.0
        dth = np.linalg.solve(A, g)
        if float(np.max(np.abs(dth))) < 1e-14:
            break
        th += dth
    return th


def leak_terms(design, dlen, th, omega, parasites):
    """Directional gate table + the leak accounting.

    Leak model (documented in FORMFIND.md):
      strut: both directions are wanted channels; deficit (1-g) each.
      cable: forward wanted: deficit (1-g_f); backward unwanted: openness g_b.
      parasite: both directions unwanted: openness g_f + g_b.
    leak_alt uses the task-sheet literal cable term (1 - gate(-2 phi)).
    """
    rows = []
    leak = 0.0
    leak_alt = 0.0
    back_sum = 0.0
    live = []           # gates over wanted directions (struts both, cable fwd)
    fwd = []            # kernel '# net:' comparable: forward gates as listed
    for i, (a, b, ty) in enumerate(design.edges):
        pf = float(wrap(th[a] - omega * dlen[i] / LAW["C"] - th[b]))
        pb = float(wrap(th[b] - omega * dlen[i] / LAW["C"] - th[a]))
        gf, gb = float(gate(pf)), float(gate(pb))
        rows.append((a, b, ty, dlen[i], pf, pb, gf, gb))
        fwd.append(gf)
        if ty == "S":
            leak += (1.0 - gf) + (1.0 - gb)
            leak_alt += 2.0 * (1.0 - gb)
            live += [gf, gb]
        else:
            leak += (1.0 - gf) + gb
            leak_alt += (1.0 - gf) + (1.0 - gb)
            back_sum += gb
            live.append(gf)
    prow = []
    for (a, b, d) in parasites:
        pf = float(wrap(th[a] - omega * d / LAW["C"] - th[b]))
        pb = float(wrap(th[b] - omega * d / LAW["C"] - th[a]))
        gf, gb = float(gate(pf)), float(gate(pb))
        prow.append((a, b, d, pf, pb, gf, gb))
        leak += gf + gb
        leak_alt += gf + gb
    return dict(rows=rows, prows=prow, leak=leak, leak_alt=leak_alt,
                back_sum=back_sum,
                min_live=float(min(live)) if live else 1.0,
                mean_live=float(np.mean(live)) if live else 1.0,
                min_fwd=float(min(fwd)) if fwd else 1.0,
                mean_fwd=float(np.mean(fwd)) if fwd else 1.0)


def leak_value(design, dlen, th, omega, parasites):
    return leak_terms(design, dlen, th, omega, parasites)["leak"]


def polish_phases(design, dlen, th0, omega, parasites):
    """L-BFGS on the exact leak objective (incl. parasites), gauge pinned."""
    if not HAVE_SCIPY or design.nv < 2:
        return th0
    ne = design.ne
    ea = np.array([e[0] for e in design.edges])
    eb = np.array([e[1] for e in design.edges])
    phi = np.array([omega * d / LAW["C"] for d in dlen])
    is_s = np.array([e[2] == "S" for e in design.edges])
    pa = np.array([p[0] for p in parasites], int) if parasites else np.zeros(0, int)
    pb = np.array([p[1] for p in parasites], int) if parasites else np.zeros(0, int)
    pphi = (np.array([p[2] for p in parasites]) * omega / LAW["C"]
            if parasites else np.zeros(0))

    def J_and_grad(thf):
        th = np.concatenate(([th0[0]], thf))    # gauge: th_0 fixed
        g = np.zeros(design.nv)
        psf = th[ea] - phi - th[eb]
        psb = th[eb] - phi - th[ea]
        gf, gb = gate(psf), gate(psb)
        dgf, dgb = dgate(psf), dgate(psb)
        # struts: (1-gf)+(1-gb);   cables: (1-gf) + gb
        J = float(np.sum(np.where(is_s, (1 - gf) + (1 - gb), (1 - gf) + gb)))
        cf = np.where(is_s, -1.0, -1.0)     # d/dpsi coefficient on gf
        cb = np.where(is_s, -1.0, +1.0)     # on gb
        np.add.at(g, ea, cf * dgf - cb * dgb)
        np.add.at(g, eb, -cf * dgf + cb * dgb)
        if len(pa):
            qf = th[pa] - pphi - th[pb]
            qb = th[pb] - pphi - th[pa]
            J += float(np.sum(gate(qf) + gate(qb)))
            dqf, dqb = dgate(qf), dgate(qb)
            np.add.at(g, pa, dqf - dqb)
            np.add.at(g, pb, -dqf + dqb)
        return J, g[1:]

    res = _scipy_minimize(J_and_grad, th0[1:], jac=True, method="L-BFGS-B",
                          options=dict(maxiter=300, ftol=1e-14, gtol=1e-12))
    return np.concatenate(([th0[0]], res.x))


# ------------------------------------------------------------------------
# picking (annealed) and parasites
# ------------------------------------------------------------------------

def find_parasites(foam, design, cells):
    """Unintended pairs of picked cells inside the exact link ceiling."""
    listed = set()
    for (a, b, _t) in design.edges:
        listed.add((a, b)); listed.add((b, a))
    out = []
    for a in range(design.nv):
        for b in range(a + 1, design.nv):
            if (a, b) in listed:
                continue
            d = foam.link_d(cells[a], cells[b])
            if d is not None:
                out.append((a, b, d))
    return out


def missing_links(foam, design, cells):
    return [i for i, (a, b, _t) in enumerate(design.edges)
            if not foam.linked(cells[a], cells[b])]


def edge_lengths(foam, design, cells):
    out = np.empty(design.ne)
    for i, (a, b, _t) in enumerate(design.edges):
        d = foam.link_d(cells[a], cells[b])
        if d is None:
            d = float(np.linalg.norm(foam.pos[cells[b]] - foam.pos[cells[a]]))
        out[i] = d
    return out


def omega_band(design, foam, cells, dlen):
    lo, hi = OMEGA_HARD
    lo += 1e-3; hi -= 1e-3
    if isinstance(design.band, tuple) and design.band[0] == "fix":
        w = design.band[1]
        return (w, w)
    if design.band == "strut":
        ds = [dlen[i] for i, e in enumerate(design.edges) if e[2] == "S"]
        w0 = math.pi / float(np.mean(ds)) if ds else 2.0
        return (max(lo, 0.85 * w0), min(hi, 1.15 * w0))
    if design.band == "cycle":
        Lc = 0.0
        loop = design.band_cycle
        for k in range(len(loop)):
            a, b = loop[k], loop[(k + 1) % len(loop)]
            d = foam.link_d(cells[a], cells[b])
            if d is None:
                d = float(np.linalg.norm(foam.pos[cells[b]]
                                         - foam.pos[cells[a]]))
            Lc += d
        w0 = TWO_PI * design.band_m / Lc
        return (max(lo, 0.93 * w0), min(hi, 1.07 * w0))
    return (lo, hi)


class PickCtx:
    """Precomputed structure for the pick objective: cycle basis + BFS-tree
    paths between all vertex pairs (for tree-predicted parasite gates) +
    pinned-loop edge lists (forced closure integers)."""

    def __init__(self, design):
        self.cycles, parent, order, roots = cycle_basis(design)
        edir = {}
        for idx, (a, b, _t) in enumerate(design.edges):
            edir[(a, b)] = (idx, +1)
            edir[(b, a)] = (idx, -1)
        self.pinned = []
        for loop, m in design.pinned:
            lst = []
            for k in range(len(loop)):
                a, b = loop[k], loop[(k + 1) % len(loop)]
                lst.append(edir[(a, b)])
            self.pinned.append((lst, m))
        comp = {}
        root_path = {}
        for v in order:
            if parent[v] is None:
                root_path[v] = []
                comp[v] = v
            else:
                u, idx, sg = parent[v]
                root_path[v] = root_path[u] + [(idx, sg)]
                comp[v] = comp[u]
        self.pair_path = {}
        for a in range(design.nv):
            for b in range(a + 1, design.nv):
                if comp[a] != comp[b]:
                    self.pair_path[(a, b)] = None    # cross-component
                    continue
                pa, pb = root_path[a], root_path[b]
                k = 0
                while k < len(pa) and k < len(pb) and pa[k] == pb[k]:
                    k += 1
                # path a->b: reverse(pa[k:]) then pb[k:]
                path = ([(idx, -sg) for idx, sg in reversed(pa[k:])]
                        + list(pb[k:]))
                self.pair_path[(a, b)] = path


def pick_objective(foam, design, cells, ctx, wgrid):
    """Gate-faithful pick objective, minimized over the omega grid.

    struts: exact deficit 2*(1-gate(pi - omega*d)) + 3x worst-edge term;
    cable cycles: spread proxy 2*wrap(defect)^2/n_links;
    parasites: gates at the TREE-predicted phase offset (theta_a - theta_b
    ~ pi*n_struts_on_path + omega*signed_cable_length) + 0.25 structural;
    missing intended links: +25 each."""
    dlen = edge_lengths(foam, design, cells)
    miss = missing_links(foam, design, cells)
    s_idx = [i for i, e in enumerate(design.edges) if e[2] == "S"]
    F = np.zeros(len(wgrid))
    if s_idx:
        ds = dlen[s_idx]
        psi = math.pi - wgrid[:, None] * ds[None, :]
        dfc = 1.0 - gate(psi)
        F += 2.0 * np.sum(dfc, axis=1) + 3.0 * np.max(dfc, axis=1)
    for cyc in ctx.cycles:
        ns = sum(1 for idx, _sg in cyc if design.edges[idx][2] == "S")
        Sc = sum(sg * dlen[idx] for idx, sg in cyc
                 if design.edges[idx][2] == "C")
        D = wrap(math.pi * ns + wgrid * Sc)
        F += 2.0 * D * D / max(1, len(cyc))
    for lst, m in ctx.pinned:
        ns = sum(1 for idx, sg in lst if design.edges[idx][2] == "S")
        Sc = sum(sg * dlen[idx] for idx, sg in lst
                 if design.edges[idx][2] == "C")
        D = math.pi * ns + wgrid * Sc - TWO_PI * m     # NOT wrapped
        F += 2.0 * D * D / max(1, len(lst))
    par = find_parasites(foam, design, cells)
    for (a, b, d) in par:
        path = ctx.pair_path[(a, b)]
        if path is None:
            F += 1.25          # cross-component: phases unconstrained
            continue
        ns = sum(1 for idx, _sg in path if design.edges[idx][2] == "S")
        Sc = sum(sg * dlen[idx] for idx, sg in path
                 if design.edges[idx][2] == "C")
        off = math.pi * ns + wgrid * Sc          # predicted th_a - th_b
        F += gate(off - wgrid * d) + gate(-off - wgrid * d) + 0.25
    J = float(F.min())
    J += 25.0 * len(miss)
    return J


def random_rotation(rng):
    """Uniform random rotation (QR of a gaussian matrix, det +1)."""
    M = rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(M)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 2] = -Q[:, 2]
    return Q


def pick_cells(foam, design, center, rng, restarts=8, max_sweeps=12,
               nrot=10):
    """Template-rotation search x multistart coordinate descent over
    per-vertex candidate pools.  Each sweep tries every pool candidate for
    every vertex (random order) and keeps improvements; iterate to
    convergence; randomized restarts escape local minima.  Beats the kernel
    shell seeder's 3 greedy passes.  Deterministic under the passed rng."""
    cc = np.asarray(center)
    ok_margin = np.all((foam.pos >= BOUND_MARGIN)
                       & (foam.pos <= LAW["L"] - BOUND_MARGIN), axis=1)
    ctx = PickCtx(design)
    lo, hi = OMEGA_HARD
    wgrid = np.linspace(lo + 1e-3, hi - 1e-3, 60)

    best_cells, best_J = None, math.inf
    for ir in range(max(1, nrot)):
        Rm = np.eye(3) if ir == 0 else random_rotation(rng)
        targets = design.verts @ Rm.T + cc
        pools = []
        for t in targets:
            d2 = np.sum((foam.pos - t) ** 2, axis=1)
            idx = np.nonzero((d2 <= design.pool_r ** 2) & ok_margin)[0]
            idx = idx[np.argsort(d2[idx], kind="stable")]
            pools.append(list(map(int, idx)))

        def init_cells(randomized):
            cells, used = [], set()
            for k, t in enumerate(targets):
                cand = [c for c in pools[k] if c not in used]
                if not cand:
                    got = foam.nearest_free(t, used)
                elif randomized:
                    got = cand[int(rng.integers(len(cand)))]
                else:
                    got = cand[0]
                cells.append(got); used.add(got)
            return cells

        for rs in range(restarts):
            cur = init_cells(randomized=(rs > 0))
            J = pick_objective(foam, design, cur, ctx, wgrid)
            for _sw in range(max_sweeps):
                improved = False
                for k in rng.permutation(design.nv):
                    k = int(k)
                    if not pools[k]:
                        continue
                    used = set(cur) - {cur[k]}
                    bestc, bestJ = cur[k], J
                    for c in pools[k]:
                        if c in used or c == cur[k]:
                            continue
                        old = cur[k]
                        cur[k] = c
                        J2 = pick_objective(foam, design, cur, ctx, wgrid)
                        if J2 < bestJ - 1e-12:
                            bestc, bestJ = c, J2
                        cur[k] = old
                    if bestc != cur[k]:
                        cur[k] = bestc
                        J = bestJ
                        improved = True
                if not improved:
                    break
            if J < best_J:
                best_J, best_cells = J, list(cur)
    return best_cells, best_J


# ------------------------------------------------------------------------
# tuner
# ------------------------------------------------------------------------

def tune(foam, design, cells, ngrid=240, do_polish=True):
    """Joint (omega, integers, phases): grid + golden refine on exact leak.

    Intended edges that are NOT foam links (d >= 1.15*(r_i+r_j)) are
    excluded from solving and scoring — they carry no channel; they are
    reported as MISSING (the H1 v3 cube silently had two of these)."""
    dlen_full = edge_lengths(foam, design, cells)
    miss = missing_links(foam, design, cells)
    if miss:
        mset = set(miss)
        keep = [i for i in range(design.ne) if i not in mset]
        eff = Design(design.name, design.verts,
                     [design.edges[i] for i in keep], band=design.band,
                     band_cycle=design.band_cycle, band_m=design.band_m,
                     pool_r=design.pool_r, pinned=design.pinned)
        dlen = dlen_full[keep]
    else:
        keep = list(range(design.ne))
        eff, dlen = design, dlen_full
    parasites = find_parasites(foam, design, cells)
    lo, hi = omega_band(eff, foam, cells, dlen)

    def eval_w(w):
        th = solve_phases(eff, dlen, w)
        return leak_value(eff, dlen, th, w, parasites), th

    if hi - lo < 1e-12:
        wbest = lo
        Jbest, thbest = eval_w(wbest)
    else:
        ws = np.linspace(lo, hi, ngrid)
        Js = np.empty(ngrid)
        for i, w in enumerate(ws):
            Js[i], _ = eval_w(w)
        i0 = int(np.argmin(Js))
        a = ws[max(0, i0 - 1)]; b = ws[min(ngrid - 1, i0 + 1)]
        gr = (math.sqrt(5) - 1) / 2
        c = b - gr * (b - a); d = a + gr * (b - a)
        fc, _ = eval_w(c); fd, _ = eval_w(d)
        for _ in range(28):
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - gr * (b - a)
                fc, _ = eval_w(c)
            else:
                a, c, fc = c, d, fd
                d = a + gr * (b - a)
                fd, _ = eval_w(d)
        wbest = 0.5 * (a + b)
        Jbest, thbest = eval_w(wbest)
    if do_polish and HAVE_SCIPY:
        th2 = polish_phases(eff, dlen, thbest, wbest, parasites)
        J2 = leak_value(eff, dlen, th2, wbest, parasites)
        if J2 < Jbest:
            Jbest, thbest = J2, th2
    x = x_of_omega(wbest)
    cycles, _p, _o, _r = cycle_basis(eff)
    dfs = cycle_defects(eff, cycles, dlen, wbest)
    mints = []
    for cyc in cycles:
        s = 0.0
        for idx, sg in cyc:
            ty = eff.edges[idx][2]
            mu = math.pi if ty == "S" else wbest * dlen[idx] / LAW["C"]
            s += sg * mu
        mints.append(int(round(s / TWO_PI)))
    return dict(omega=wbest, x=x, th=thbest, dlen=dlen, dlen_full=dlen_full,
                parasites=parasites, leak=Jbest, cycles=cycles, defects=dfs,
                mints=mints, miss=miss, keep=keep, eff=eff)


def stiffness_spectrum(design, dlen, th, omega):
    """Weighted lock Laplacian (weights = gate stiffness -g'' per wanted
    direction at the operating point); eigenvalues; gap above gauge zero."""
    W = np.zeros(design.ne)
    for i, (a, b, ty) in enumerate(design.edges):
        pf = wrap(th[a] - omega * dlen[i] / LAW["C"] - th[b])
        pb = wrap(th[b] - omega * dlen[i] / LAW["C"] - th[a])
        if ty == "S":
            W[i] = max(0.0, -d2gate(pf)) + max(0.0, -d2gate(pb))
        else:
            W[i] = max(0.0, -d2gate(pf))
    Lm = np.zeros((design.nv, design.nv))
    for i, (a, b, _t) in enumerate(design.edges):
        Lm[a, a] += W[i]; Lm[b, b] += W[i]
        Lm[a, b] -= W[i]; Lm[b, a] -= W[i]
    ev = np.linalg.eigvalsh(Lm)
    tol = max(1e-9, 1e-9 * max(1.0, float(ev[-1])))
    nz = ev[ev > tol]
    gap = float(nz[0]) if len(nz) else 0.0
    return ev, gap


# ------------------------------------------------------------------------
# emission
# ------------------------------------------------------------------------

def bfs_mean_phases(design, dlen, omega):
    """Control twin: BFS over intended edges with per-TYPE mean drops,
    ignoring actual per-edge lengths (the mis-tuned A/B)."""
    ds = [dlen[i] for i, e in enumerate(design.edges) if e[2] == "S"]
    dc = [dlen[i] for i, e in enumerate(design.edges) if e[2] == "C"]
    mu_s = omega * float(np.mean(ds)) / LAW["C"] if ds else math.pi
    mu_c = omega * float(np.mean(dc)) / LAW["C"] if dc else 0.0
    adj = {}
    for idx, (a, b, ty) in enumerate(design.edges):
        adj.setdefault(a, []).append((b, ty, +1))
        adj.setdefault(b, []).append((a, ty, -1))
    th = np.zeros(design.nv)
    seen = set()
    for r0 in range(design.nv):
        if r0 in seen:
            continue
        seen.add(r0)
        q = [r0]
        while q:
            u = q.pop(0)
            for (v, ty, sg) in adj.get(u, ()):
                if v in seen:
                    continue
                mu = mu_s if ty == "S" else mu_c
                th[v] = th[u] - sg * mu
                seen.add(v)
                q.append(v)
    return th


def emit_net(path, foam, design, cells, x, th, header_lines):
    with open(path, "w") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        for k in range(design.nv):
            p = foam.pos[cells[k]]
            t = math.fmod(th[k] + 8.0 * TWO_PI, TWO_PI)
            f.write("V %.4f %.4f %.4f %.6f %.6f\n"
                    % (p[0], p[1], p[2], x, t))
        for (a, b, _ty) in design.edges:
            f.write("E %d %d\n" % (a, b))


# ------------------------------------------------------------------------
# jitter robustness
# ------------------------------------------------------------------------

def jitter_study(foam, design, rng, n=20, restarts=3, max_sweeps=8, nrot=4):
    ext = design.ext + design.pool_r + 0.4
    lo = BOUND_MARGIN + ext
    hi = LAW["L"] - lo
    out = []
    for j in range(n):
        c = lo + (hi - lo) * rng.random(3)
        cells, _J = pick_cells(foam, design, c, rng, restarts=restarts,
                               max_sweeps=max_sweeps, nrot=nrot)
        tn = tune(foam, design, cells, ngrid=100, do_polish=False)
        rep = leak_terms(tn["eff"], tn["dlen"], tn["th"], tn["omega"],
                         tn["parasites"])
        out.append(dict(center=tuple(np.round(c, 2)),
                        min_gate=rep["min_live"], mean_gate=rep["mean_live"],
                        leak=rep["leak"], miss=len(tn["miss"]),
                        npar=len(tn["parasites"])))
    return out


def percentile_of(value, samples, higher_better=True):
    if not samples:
        return 100.0
    s = np.asarray(samples)
    if higher_better:
        return 100.0 * float(np.mean(s <= value))
    return 100.0 * float(np.mean(s >= value))


# ------------------------------------------------------------------------
# candidate templates
# ------------------------------------------------------------------------

def cube_design(name, a, pool_r=1.5, notes=""):
    h = 0.5 * a
    verts = [[(h if k & 1 else -h), (h if k & 2 else -h),
              (h if k & 4 else -h)] for k in range(8)]
    edges = [(k, k ^ (1 << b), "S") for k in range(8) for b in range(3)
             if (k ^ (1 << b)) > k]
    return Design(name, verts, edges, band="strut", pool_r=pool_r,
                  notes=notes)


def hexprism_design(name, a=1.5, h=1.5):
    verts = []
    for k in range(6):
        ph = TWO_PI * k / 6
        verts.append([a * math.cos(ph), a * math.sin(ph), -0.5 * h])
    for k in range(6):
        ph = TWO_PI * k / 6
        verts.append([a * math.cos(ph), a * math.sin(ph), +0.5 * h])
    edges = []
    for k in range(6):
        edges.append((k, (k + 1) % 6, "S"))
        edges.append((6 + k, 6 + (k + 1) % 6, "S"))
        edges.append((k, 6 + k, "S"))
    return Design(name, verts, edges, band="strut",
                  notes="all-strut hexagonal prism (over-constraint probe)")


def tube6_design(name, a=1.2, hfac=1.5):
    h = hfac * a
    verts = []
    for k in range(6):
        ph = TWO_PI * k / 6
        verts.append([a * math.cos(ph), a * math.sin(ph), -0.5 * h])
    for k in range(6):
        ph = TWO_PI * k / 6
        verts.append([a * math.cos(ph), a * math.sin(ph), +0.5 * h])
    edges = []
    for k in range(6):
        edges.append((k, (k + 1) % 6, "C"))            # ring 1, co-rotating
    for k in range(6):
        edges.append((6 + k, 6 + (k + 1) % 6, "C"))    # ring 2, same sense
    for k in range(6):
        edges.append((k, 6 + k, "S"))                  # axial pi-rung struts
    return Design(name, verts, edges, band="cycle",
                  band_cycle=[0, 1, 2, 3, 4, 5], band_m=2, pool_r=1.35,
                  notes="B1c co-rotating wound tube: rings m=2 (phi=2pi/3, "
                        "one-way), axial struts d=pi/omega=1.5*dbar_ring",
                  pinned=[([0, 1, 2, 3, 4, 5], 2),
                          ([6, 7, 8, 9, 10, 11], 2)])


def ring12_design(name, d_ring=1.5, m=5, as_struts=False):
    R = d_ring / (2.0 * math.sin(math.pi / 12))
    verts = [[R * math.cos(TWO_PI * k / 12), R * math.sin(TWO_PI * k / 12),
              0.0] for k in range(12)]
    ty = "S" if as_struts else "C"
    edges = [(k, (k + 1) % 12, ty) for k in range(12)]
    d = Design(name, verts, edges,
               band=("strut" if as_struts else "cycle"),
               band_cycle=list(range(12)), band_m=m,
               notes="ring12")
    return d


def ring12_chords_design(name, d_ring=1.5):
    """C4: unwound strut ring12 (m=6) + two 5-hop cable chords (omega*L=4pi,
    the chords carry winding m=2 each) + two consonant cross struts.
    Single-link chords are geometrically impossible (2R=5.8, k=2 skip
    2*d*cos15 > ceiling); on a wound m=5 ring every symmetric crossing is
    frustrated by pi/2 — worked out in FORMFIND.md."""
    R = d_ring / (2.0 * math.sin(math.pi / 12))
    verts = [[R * math.cos(TWO_PI * k / 12), R * math.sin(TWO_PI * k / 12),
              0.0] for k in range(12)]
    # chord A: v0 -> a1..a4 -> v6 along x; arch +z
    fr = [0.2, 0.4, 0.6, 0.8]
    za = [0.42, 0.65, 0.65, 0.42]
    A = [[R * (1 - 2 * f), 0.0, z] for f, z in zip(fr, za)]
    B = [[0.0, R * (1 - 2 * f), -z] for f, z in zip(fr, za)]
    verts = verts + A + B          # a1..a4 = 12..15, b1..b4 = 16..19
    edges = [(k, (k + 1) % 12, "S") for k in range(12)]
    chainA = [0, 12, 13, 14, 15, 6]
    chainB = [3, 16, 17, 18, 19, 9]
    for u, v in zip(chainA[:-1], chainA[1:]):
        edges.append((u, v, "C"))
    for u, v in zip(chainB[:-1], chainB[1:]):
        edges.append((u, v, "C"))
    edges.append((13, 17, "S"))    # a2-b2 cross strut
    edges.append((14, 18, "S"))    # a3-b3 cross strut
    return Design(name, verts, edges, band="strut", pool_r=1.4,
                  notes="unwound strut ring12 + two wound cable chords "
                        "(m=2 each) + consonant cross struts")


def torus_design(name, n1=3, n2=8, r_m=0.80, R_M=2.15):
    verts = []
    beta0 = math.pi / 3.0 if n1 == 3 else math.pi / 4.0   # inward-pointing
    for j in range(n2):
        al = TWO_PI * j / n2
        for i in range(n1):
            be = beta0 + TWO_PI * i / n1
            lat = R_M + r_m * math.cos(be)
            verts.append([lat * math.cos(al), lat * math.sin(al),
                          r_m * math.sin(be)])
    def vid(i, j):
        return (j % n2) * n1 + (i % n1)
    edges = []
    for j in range(n2):
        for i in range(n1):
            edges.append((vid(i, j), vid(i + 1, j), "C"))     # minor
    for j in range(n2):
        for i in range(n1):
            edges.append((vid(i, j), vid(i, j + 1), "C"))     # major
    d = Design(name, verts, edges, band="cycle",
               band_cycle=[vid(i, 0) for i in range(n1)], band_m=1,
               pool_r=1.6,
               notes="quad/tri-mesh torus %dx%d, homology integers (m1,m2)"
                     % (n1, n2),
               pinned=[([vid(i, j) for i in range(n1)], 1)
                       for j in range(n2)])
    return d


def truncocta_design(name, edge=1.5):
    s = edge / math.sqrt(2.0)
    base = []
    from itertools import permutations
    seen = set()
    for perm in permutations(range(3)):
        for s1 in (+1, -1):
            for s2 in (+1, -1):
                v = [0.0, 0.0, 0.0]
                v[perm[0]] = 0.0
                v[perm[1]] = s1 * s
                v[perm[2]] = s2 * 2 * s
                t = tuple(np.round(v, 9))
                if t not in seen:
                    seen.add(t)
                    base.append(list(v))
    verts = np.array(base)
    edges = []
    for a in range(len(verts)):
        for b in range(a + 1, len(verts)):
            d = np.linalg.norm(verts[b] - verts[a])
            if abs(d - edge) < 1e-6:
                edges.append((a, b, "S"))
    return Design(name, verts, edges, band="strut", pool_r=1.4,
                  notes="truncated octahedron n=24, 36 struts, bipartite")


# ------------------------------------------------------------------------
# reporting
# ------------------------------------------------------------------------

def netgate_lines(design, rep):
    out = []
    for (a, b, ty, d, pf, pb, gf, gb) in rep["rows"]:
        out.append("# NETGATE %d %d d=%.4f psi_f=%+.4f psi_b=%+.4f "
                   "gf=%.4f gb=%.4f   [%s]" % (a, b, d, pf, pb, gf, gb, ty))
    for (a, b, d, pf, pb, gf, gb) in rep["prows"]:
        out.append("# NETGATE P %d %d d=%.4f psi_f=%+.4f psi_b=%+.4f "
                   "gf=%.4f gb=%.4f" % (a, b, d, pf, pb, gf, gb))
    return out


def write_report(path, name, foam, design, cells, tn, rep, ev, gap,
                 jit, extra=()):
    with open(path, "w") as f:
        f.write("# formfind report — %s\n" % name)
        f.write("# foam %s (cells=%d links=%d dbar=%.4f)\n"
                % (os.path.basename(foam.path), foam.n, foam.nl, foam.dbar))
        for line in extra:
            f.write("# %s\n" % line)
        f.write("# omega=%.6f x=%.6f Em/voice=%.6f mass=%.4f\n"
                % (tn["omega"], tn["x"], tn["x"] * LAW["cap"],
                   design.nv * tn["x"] * LAW["cap"]))
        f.write("# gates live: min=%.4f mean=%.4f | fwd (kernel '# net:') "
                "min=%.4f mean=%.4f\n"
                % (rep["min_live"], rep["mean_live"],
                   rep["min_fwd"], rep["mean_fwd"]))
        f.write("# leak=%.5f leak_alt=%.5f back_sum=%.5f parasites=%d "
                "missing=%d\n"
                % (rep["leak"], rep["leak_alt"], rep["back_sum"],
                   len(tn["parasites"]), len(tn["miss"])))
        f.write("# cycle integers m=%s defects=%s\n"
                % (tn["mints"], ["%.4f" % d for d in tn["defects"]]))
        f.write("# stiffness spectrum: lam1=%.2e gap=%.4f lam_max=%.4f\n"
                % (ev[0], gap, ev[-1]))
        f.write("# picks (vertex -> cell id, position, radius):\n")
        for k, c in enumerate(cells):
            p = foam.pos[c]
            f.write("#   v%-3d cell %-5d (%.4f, %.4f, %.4f) r=%.4f\n"
                    % (k, c, p[0], p[1], p[2], foam.r[c]))
        for line in netgate_lines(design, rep):
            f.write(line + "\n")
        for i in tn.get("miss", []):
            a, b, ty = design.edges[i]
            dd = float(np.linalg.norm(foam.pos[cells[b]]
                                      - foam.pos[cells[a]]))
            cut = foam.cut(cells[a], cells[b])
            f.write("# MISSING %d %d d=%.4f cut=%.4f [%s] — intended edge "
                    "is NOT a foam link\n" % (a, b, dd, cut, ty))
        if jit is not None:
            mg = sorted(j["min_gate"] for j in jit)
            f.write("# jitter (%d centers): min_gate quartiles "
                    "%.4f / %.4f / %.4f  (worst %.4f best %.4f)\n"
                    % (len(jit), np.percentile(mg, 25), np.percentile(mg, 50),
                       np.percentile(mg, 75), mg[0], mg[-1]))
            for j in jit:
                f.write("#   center %s min=%.4f mean=%.4f leak=%.4f "
                        "miss=%d par=%d\n"
                        % (j["center"], j["min_gate"], j["mean_gate"],
                           j["leak"], j["miss"], j["npar"]))


# ------------------------------------------------------------------------
# candidate driver
# ------------------------------------------------------------------------

def run_candidate(foam, design, outdir, rng, jitter_n=20,
                  cells_override=None, omega_fix=None, emit=True,
                  jit_run=True, restarts=8, nrot=10):
    center = (0.5 * LAW["L"], 0.5 * LAW["L"], 0.5 * LAW["L"])
    if cells_override is not None:
        cells = list(cells_override)
    else:
        cells, _J = pick_cells(foam, design, center, rng,
                               restarts=restarts, nrot=nrot)
    if omega_fix is not None:
        design.band = ("fix", omega_fix)
    tn = tune(foam, design, cells)
    rep = leak_terms(tn["eff"], tn["dlen"], tn["th"], tn["omega"],
                     tn["parasites"])
    ev, gap = stiffness_spectrum(tn["eff"], tn["dlen"], tn["th"],
                                 tn["omega"])
    jit = jitter_study(foam, design, rng, n=jitter_n) if jit_run else None
    pct = (percentile_of(rep["min_live"], [j["min_gate"] for j in jit])
           if jit else float("nan"))
    if emit:
        os.makedirs(outdir, exist_ok=True)
        hdr = [
            "formfind v1 — candidate %s (v89 PRESTRESS stream B)" % design.name,
            "foam %s; solved deterministically (seed %d)"
            % (os.path.basename(foam.path), FSEED),
            "omega=%.6f x=%.6f Em/voice=%.4f mass=%.4f"
            % (tn["omega"], tn["x"], tn["x"] * LAW["cap"],
               design.nv * tn["x"] * LAW["cap"]),
            "gates live min=%.4f mean=%.4f; pred '# net:' mean(gf)=%.4f"
            % (rep["min_live"], rep["mean_live"], rep["mean_fwd"]),
            "leak=%.5f back_sum=%.5f parasites=%d missing=%d m=%s"
            % (rep["leak"], rep["back_sum"], len(tn["parasites"]),
               len(tn["miss"]), tn["mints"]),
            "edge types: %s"
            % " ".join("%d-%d:%s" % (a, b, t) for a, b, t in design.edges),
        ]
        emit_net(os.path.join(outdir, design.name + ".net"),
                 foam, design, cells, tn["x"], tn["th"], hdr)
        thc = bfs_mean_phases(design, tn["dlen_full"], tn["omega"])
        repc = leak_terms(tn["eff"], tn["dlen"], thc, tn["omega"],
                          tn["parasites"])
        hdrc = ["formfind v1 — CONTROL twin of %s" % design.name,
                "same cells, same loads; phases mean-tuned BFS ignoring "
                "actual lengths",
                "predicted gates live min=%.4f mean=%.4f (tuned: %.4f/%.4f)"
                % (repc["min_live"], repc["mean_live"],
                   rep["min_live"], rep["mean_live"])]
        emit_net(os.path.join(outdir, design.name + "_ctrl.net"),
                 foam, design, cells, tn["x"], thc, hdrc)
        write_report(os.path.join(outdir, design.name + "_report.txt"),
                     design.name, foam, design, cells, tn, rep, ev, gap, jit)
    return dict(design=design, cells=cells, tn=tn, rep=rep, ev=ev, gap=gap,
                jit=jit, pct=pct)


def summary_row(name, res, notes=""):
    d, tn, rep = res["design"], res["tn"], res["rep"]
    ns = len(d.struts()); nc = len(d.cables())
    jit = res["jit"]
    med = (float(np.median([j["min_gate"] for j in jit]))
           if jit else float("nan"))
    gpar_max = max([max(pr[5], pr[6]) for pr in rep["prows"]], default=0.0)
    return ("%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
            "\t%.5f\t%.5f\t%d\t%.4f\t%d\t%.4f\t%.1f\t%.4f\t%s"
            % (name, d.nv, d.ne, ns, nc, tn["omega"], tn["x"],
               tn["x"] * LAW["cap"], d.nv * tn["x"] * LAW["cap"],
               rep["min_live"], rep["mean_live"], rep["back_sum"],
               rep["leak"], rep["leak_alt"], len(tn["parasites"]),
               gpar_max, len(tn["miss"]), res["gap"], res["pct"], med, notes))


SUMMARY_HDR = ("cand\tnv\tne\tn_strut\tn_cable\tomega\tx\tem_voice\tmass"
               "\tmin_gate\tmean_gate\tback_sum\tleak\tleak_alt\tn_par"
               "\tgpar_max\tn_missing\tspec_gap\trobust_pctile"
               "\trobust_med_min\tnotes")


# ------------------------------------------------------------------------
# validations (C8 + C1 kernel mimics)
# ------------------------------------------------------------------------

def validate_ring(foam, verbose=True):
    picks, Lring, R = kernel_ring_pick(foam, n=REF_RING["n"],
                                       ring_d=REF_RING["d_target"])
    om = TWO_PI * REF_RING["m"] * LAW["C"] / Lring
    x = x_of_omega(om)
    closure = om * Lring / TWO_PI
    ok = abs(Lring - REF_RING["Lring"]) <= 0.002
    if verbose:
        print("# ring mimic: n=%d R=%.2f d_target=%.2f Lring=%.3f "
              "closure/2pi=%.4f ring_m=%d  (ref Lring=%.3f)  %s"
              % (REF_RING["n"], R, REF_RING["d_target"], Lring, closure,
                 REF_RING["m"], REF_RING["Lring"],
                 "MATCH" if ok else "MISMATCH"))
    return dict(picks=picks, Lring=Lring, omega=om, x=x, ok=ok, R=R)


def validate_shell(foam, verbose=True):
    picks, edges, dl, abar = kernel_shell_pick(
        foam, a_target=REF_SHELL["a_target"], r_core=REF_SHELL["r_core"])
    om = math.pi * LAW["C"] / abar
    x = x_of_omega(om)
    th, gates = kernel_shell_bfs(foam, picks, om)
    gmin, gmean = min(gates), float(np.mean(gates))
    ok = (abs(abar - REF_SHELL["abar"]) <= 0.002
          and abs(om - REF_SHELL["omega"]) <= 0.002
          and abs(x - REF_SHELL["x"]) <= 0.002
          and abs(gmean - REF_SHELL["gmean"]) <= 0.01
          and abs(gmin - REF_SHELL["gmin"]) <= 0.005)
    if verbose:
        print("# shell mimic: cube a_target=%.2f abar=%.3f omega=%.4f "
              "x=%.4f gates min=%.3f mean=%.3f  (ref %.3f/%.4f/%.4f/"
              "%.3f/%.3f)  %s"
              % (REF_SHELL["a_target"], abar, om, x, gmin, gmean,
                 REF_SHELL["abar"], REF_SHELL["omega"], REF_SHELL["x"],
                 REF_SHELL["gmin"], REF_SHELL["gmean"],
                 "MATCH" if ok else "MISMATCH"))
    return dict(picks=picks, edges=edges, dl=dl, abar=abar, omega=om, x=x,
                th=th, gates=gates, gmin=gmin, gmean=gmean, ok=ok)


def selftest(foam):
    print("== selftest ==")
    print("# foam: cells=%d links=%d mean_degree=%.2f dbar=%.4f "
          "(ref 9741 / 85886 / 17.63 / 1.5053)"
          % (foam.n, foam.nl, foam.mean_degree, foam.dbar))
    assert foam.n == 9741
    assert abs(foam.dbar - 1.5053) < 0.002
    assert abs(foam.nl - 85886) <= 80   # 4-dp TSV rounding wobbles the rule
    # gradient check on a small design
    d = cube_design("t", 1.5)
    rngt = np.random.default_rng(1)
    dlen = 1.5 + 0.1 * rngt.standard_normal(d.ne)
    th = solve_phases(d, dlen, 2.0)
    par = []
    if HAVE_SCIPY:
        ea = np.array([e[0] for e in d.edges]); eb = np.array([e[1] for e in d.edges])
        phi = 2.0 * dlen
        is_s = np.array([True] * d.ne)
        def J(t):
            psf = t[ea] - phi - t[eb]; psb = t[eb] - phi - t[ea]
            return float(np.sum((1 - gate(psf)) + (1 - gate(psb))))
        g0 = np.zeros(d.nv)
        psf = th[ea] - phi - th[eb]; psb = th[eb] - phi - th[ea]
        dgf, dgb = dgate(psf), dgate(psb)
        np.add.at(g0, ea, -dgf + dgb)
        np.add.at(g0, eb, dgf - dgb)
        num = np.zeros(d.nv)
        for i in range(d.nv):
            e = np.zeros(d.nv); e[i] = 1e-6
            num[i] = (J(th + e) - J(th - e)) / 2e-6
        assert np.max(np.abs(num - g0)) < 1e-5, "gradient check failed"
        print("# gradient check: PASS")
    # pair rung sanity: strut at exactly d=pi/omega must give both gates 1
    dd = math.pi / 2.0
    ps = wrap(math.pi - 2.0 * dd)
    assert abs(float(gate(ps)) - 1.0) < 1e-12
    print("# pi-rung pair: gates 1.0000/1.0000 PASS")
    r = validate_ring(foam)
    s = validate_shell(foam)
    print("== selftest %s ==" % ("PASS" if (r["ok"] and s["ok"]) else "FAIL"))
    return r["ok"] and s["ok"]


# ------------------------------------------------------------------------
# per-candidate flows
# ------------------------------------------------------------------------

def flow_c1(foam, outdir, rng, jitter_n):
    """C1 cube a=1.25 exact retune — the H1 geometry done right, plus the
    deficit decomposition (phases vs lengths vs parasites)."""
    print("\n=== C1 cube a=1.25 (H1 done right) ===")
    v = validate_shell(foam)
    # (a) kernel picks + kernel BFS phases (validated above)
    # (b) same cells + exact retune (LSQ + polish, tuned omega)
    d = cube_design("c1_cube125", 1.25,
                    notes="H1 geometry, exact retune on kernel picks")
    res = run_candidate(foam, d, outdir, rng, jitter_n=jitter_n,
                        cells_override=v["picks"], jit_run=True)
    # decomposition on the same cells at the kernel's omega (realized
    # edges only — the kernel cube had 2 phantom edges):
    eff = res["tn"]["eff"]
    dlen = res["tn"]["dlen"]
    om_k = v["omega"]
    th_k = solve_phases(eff, dlen, om_k)
    par = res["tn"]["parasites"]
    rep_k = leak_terms(eff, dlen, th_k, om_k, par)
    # (c) free annealed picks at the same template
    d2 = cube_design("c1_cube125_freepick", 1.25)
    cells2, _ = pick_cells(foam, d2, (12.0, 12.0, 12.0), rng)
    tn2 = tune(foam, d2, cells2)
    rep2 = leak_terms(tn2["eff"], tn2["dlen"], tn2["th"], tn2["omega"],
                      tn2["parasites"])
    keep = res["tn"]["keep"]
    g_real = [v["gates"][i] for i in keep]
    dec = dict(kernel_mean_fwd=v["gmean"], kernel_min_fwd=v["gmin"],
               kernel_mean_realized=float(np.mean(g_real)),
               kernel_min_realized=float(min(g_real)),
               retune_at_kernel_omega_mean=rep_k["mean_fwd"],
               retune_tuned_mean=res["rep"]["mean_fwd"],
               retune_tuned_min=res["rep"]["min_fwd"],
               freepick_mean=rep2["mean_fwd"], freepick_min=rep2["min_fwd"],
               freepick_abar=float(np.mean(tn2["dlen"])),
               kernel_abar=v["abar"],
               n_parasites_kernel_cells=len(par),
               n_parasites_freepick=len(tn2["parasites"]))
    print("# C1 decomposition (realized edges): kernel %.3f (min %.3f) -> "
          "retune %.3f (min %.3f, same cells) -> freepick %.3f (min %.3f); "
          "parasites %d (kernel cells) %d (freepick); phantom edges %d"
          % (dec["kernel_mean_realized"], dec["kernel_min_realized"],
             dec["retune_tuned_mean"], dec["retune_tuned_min"],
             dec["freepick_mean"], dec["freepick_min"],
             dec["n_parasites_kernel_cells"], dec["n_parasites_freepick"],
             len(res["tn"]["miss"])))
    print("# C1 freepick abar=%.3f (kernel abar %.3f)"
          % (dec["freepick_abar"], dec["kernel_abar"]))
    return res, dec


def flow_c8(foam, outdir, rng):
    """C8 ring12 exact re-derivation — the validation gate — and .net."""
    print("\n=== C8 ring12 exact re-derivation (validation) ===")
    v = validate_ring(foam)
    d = ring12_design("c8_ring12", d_ring=1.25, m=5)
    # kernel lock-recursion phases, th_0 = 0 (gauge; kernel used random)
    picks = v["picks"]
    n = 12
    P = foam.pos[picks]
    dlen = np.array([float(np.linalg.norm(P[(k + 1) % n] - P[k]))
                     for k in range(n)])
    om = v["omega"]
    th = np.zeros(n)
    for k in range(n - 1):
        th[k + 1] = math.fmod(th[k] - om * dlen[k] / LAW["C"] + 8.0 * TWO_PI,
                              TWO_PI)
    # wrap into the design container for reporting/emission
    d.band = ("fix", om)
    par = find_parasites(foam, d, picks)
    rep = leak_terms(d, dlen, th, om, par)
    ev, gap = stiffness_spectrum(d, dlen, th, om)
    os.makedirs(outdir, exist_ok=True)
    hdr = ["formfind v1 — candidate c8_ring12 (kernel-exact comp12 rebuild)",
           "kernel mimic: Lring=%.3f closure/2pi=%.4f ring_m=5 omega=%.6f"
           % (v["Lring"], om * v["Lring"] / TWO_PI, om),
           "x=%.6f Em/voice=%.4f mass=%.4f (v3_comp12: x=0.4392 Em=1.0979)"
           % (v["x"], v["x"] * LAW["cap"], 12 * v["x"] * LAW["cap"]),
           "back gates gate(-2*omega*d): sum=%.4f mean=%.4f (B1: 0.1001 "
           "at phi=150deg exactly; long links leak more)"
           % (rep["back_sum"], rep["back_sum"] / 12.0)]
    emit_net(os.path.join(outdir, "c8_ring12.net"), foam, d, picks,
             v["x"], th, hdr)
    thc = bfs_mean_phases(d, dlen, om)
    repc = leak_terms(d, dlen, thc, om, par)
    emit_net(os.path.join(outdir, "c8_ring12_ctrl.net"), foam, d, picks,
             v["x"], thc,
             ["formfind v1 — CONTROL twin of c8_ring12",
              "same cells/loads; phases mean-tuned BFS ignoring lengths",
              "predicted gates live min=%.4f mean=%.4f (tuned %.4f/%.4f)"
              % (repc["min_live"], repc["mean_live"],
                 rep["min_live"], rep["mean_live"])])
    tn = dict(omega=om, x=v["x"], th=th, dlen=dlen, parasites=par,
              leak=rep["leak"], cycles=[], defects=[om * v["Lring"] / TWO_PI],
              mints=[5], miss=missing_links(foam, d, picks))
    write_report(os.path.join(outdir, "c8_ring12_report.txt"), "c8_ring12",
                 foam, d, picks, tn, rep, ev, gap, None,
                 extra=["VALIDATION: Lring=%.3f vs ref %.3f -> %s"
                        % (v["Lring"], REF_RING["Lring"],
                           "MATCH" if v["ok"] else "MISMATCH")])
    return dict(design=d, cells=picks, tn=tn, rep=rep, ev=ev, gap=gap,
                jit=None, pct=float("nan")), v


def main():
    ap = argparse.ArgumentParser(description="v89 PRESTRESS form-finder")
    ap.add_argument("--foam", default=DEF_FOAM)
    ap.add_argument("--cand", default=None,
                    help="c1..c8 (or comma list)")
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "candidates"))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--jitter", type=int, default=20)
    args = ap.parse_args()

    print("# loading foam %s" % args.foam)
    foam = Foam(args.foam)
    print("# foam: cells=%d links=%d mean_degree=%.2f dbar=%.4f"
          % (foam.n, foam.nl, foam.mean_degree, foam.dbar))

    if args.selftest:
        ok = selftest(foam)
        sys.exit(0 if ok else 1)

    wanted = ([c.strip() for c in args.cand.split(",")] if args.cand
              else None)
    if args.all:
        wanted = ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]
    if not wanted:
        ap.error("give --cand or --all")

    os.makedirs(args.out, exist_ok=True)
    rows = []
    seed_ix = {c: i for i, c in enumerate(
        ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"])}

    for cand in wanted:
        rng = np.random.default_rng(FSEED + 101 * seed_ix.get(cand, 99))
        if cand == "c1":
            res, dec = flow_c1(foam, args.out, rng, args.jitter)
            rows.append(summary_row("c1_cube125", res,
                        notes="H1 cells retuned; kernel 0.597 -> %.3f"
                              % res["rep"]["mean_fwd"]))
        elif cand == "c2":
            print("\n=== C2 cube a=1.5 parasite-free ===")
            d = cube_design("c2_cube150", 1.5,
                            notes="edges at foam dbar; diagonals out of range")
            res = run_candidate(foam, d, args.out, rng, jitter_n=args.jitter)
            rows.append(summary_row("c2_cube150", res))
        elif cand == "c3":
            print("\n=== C3 hexagonal prism (all-strut over-constraint) ===")
            d = hexprism_design("c3_hexprism", a=1.5, h=1.5)
            res = run_candidate(foam, d, args.out, rng, jitter_n=args.jitter)
            rows.append(summary_row("c3_hexprism", res,
                        notes="best mixed split = c5 tube (see FORMFIND.md)"))
        elif cand == "c4":
            print("\n=== C4 ring12 + 2 consonant chords ===")
            d = ring12_chords_design("c4_ring12_chords", d_ring=1.5)
            res = run_candidate(foam, d, args.out, rng, jitter_n=args.jitter,
                                restarts=8, nrot=14)
            rows.append(summary_row("c4_ring12_chords", res,
                        notes="strut ring m=6 + wound chords m=2"))
        elif cand == "c5":
            print("\n=== C5 co-rotating wound tube (2x ring6) ===")
            d = tube6_design("c5_tube6", a=1.2)
            res = run_candidate(foam, d, args.out, rng, jitter_n=args.jitter,
                                restarts=8, nrot=14)
            rows.append(summary_row("c5_tube6", res,
                        notes="rings m=2 one-way; axial pi-struts"))
        elif cand == "c6":
            print("\n=== C6 quad-mesh torus ===")
            d38 = torus_design("c6_torus3x8", 3, 8, r_m=0.80, R_M=2.15)
            r38 = run_candidate(foam, d38, args.out, rng,
                                jitter_n=max(6, args.jitter // 3),
                                restarts=8, nrot=14)
            d46 = torus_design("c6_torus4x6", 4, 6, r_m=1.03, R_M=1.90)
            r46 = run_candidate(foam, d46, args.out, rng,
                                jitter_n=max(6, args.jitter // 3),
                                restarts=8, nrot=14, emit=True)
            rows.append(summary_row("c6_torus3x8", r38,
                        notes="row-closure Diophantine strain (see md)"))
            rows.append(summary_row("c6_torus4x6", r46,
                        notes="row-closure Diophantine strain (see md)"))
        elif cand == "c7":
            print("\n=== C7 truncated octahedron n=24 ===")
            d = truncocta_design("c7_truncocta", edge=1.5)
            res = run_candidate(foam, d, args.out, rng,
                                jitter_n=max(6, args.jitter // 2),
                                restarts=8, nrot=14)
            rows.append(summary_row("c7_truncocta", res))
        elif cand == "c8":
            res, v = flow_c8(foam, args.out, rng)
            rows.append(summary_row("c8_ring12", res,
                        notes="kernel-exact comp12; Lring %.3f %s"
                              % (v["Lring"],
                                 "MATCH" if v["ok"] else "MISMATCH")))
        else:
            print("# unknown candidate %s" % cand)

    if rows:
        spath = os.path.join(args.out, "summary.tsv")
        with open(spath, "w") as f:
            f.write(SUMMARY_HDR + "\n")
            for r in rows:
                f.write(r + "\n")
        print("\n# wrote %s (%d rows)" % (spath, len(rows)))


if __name__ == "__main__":
    main()
