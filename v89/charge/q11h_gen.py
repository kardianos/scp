#!/usr/bin/env python3
"""Q11h — holonomy current vs HARMONIC and SIZE (the charge-offset scan).

The comb vocabulary at comb_limit=6 admits exactly three addressable
rungs in the pitch band w2/(1+q) .. w2 (ratio cap 2.2):
    unison 1:1 (Tenney 1), octave 2:1 (Tenney 2), fifth 3:2 (Tenney 6).
Nets generated (pinned, dual-book slips; the existing fifth-triangle is
the 3:2 anchor):
    oct_tri   A-B 2:1, B-C 1:2, C-A 1:1   (the octave triangle)
    fifth_tri2  a SECOND fifth triangle at a different foam site
                (address independence of the 0.44 partial-lock ratio)
    fifth_hex  Z6 ring of alternating 3:2/2:3 edges (6 voices) — the
               SIZE axis: is holonomy current per edge intensive?

Pin loads from the pitch law w = w2/(1+q*x):
    fifth: xh=0.195035 (w=2.3504), xl=0.709220 (w=1.5670), ratio 3:2
    octave: xh=0.03 (w=2.7992), xl=0.89305 (w=1.3996), ratio 2:1
Phases chained so every edge but the last seeds at comb phase 0; the
last edge's residual is the cycle's closure defect (printed).
Edge convention: E a b p q encodes pitch ratio wa:wb = p:q
(coincidence q*wa = p*wb; kernel registers link-oriented).
"""
import numpy as np

FOAM = "/home/d/code/scp/v89/prestress/foam/foam_s20260727.tsv"
OUT = "/home/d/code/scp/v89/charge/nets"
W2, QD, C = 2.9, 1.2, 1.0
XH_F, XL_F = 0.195035, 0.709220          # fifth pair (P17)
XH_O, XL_O = 0.03, 0.89305               # octave pair


def w_of(x):
    return W2 / (1.0 + QD * x)


d = np.loadtxt(FOAM, skiprows=1)
P, R = d[:, 1:4], d[:, 4]
from scipy.spatial import cKDTree
tr = cKDTree(P)
pairs = tr.query_pairs(r=2.1, output_type="ndarray")
dd = np.linalg.norm(P[pairs[:, 0]] - P[pairs[:, 1]], axis=1)
ok = dd < 1.15 * (R[pairs[:, 0]] + R[pairs[:, 1]])
E, dd = pairs[ok], dd[ok]
adj = {}
for (a, b), l in zip(E, dd):
    adj.setdefault(a, {})[b] = l
    adj.setdefault(b, {})[a] = l

AVOID = {3928, 4407, 5902}   # the P17 fifth-triangle cells


def residual_of(cells, loads, edges):
    th = {cells[0]: 0.0}
    for (i, j, p, q) in edges[:-1]:
        dij = float(np.linalg.norm(P[i] - P[j]))
        wi = w_of(loads[i])
        if i in th and j not in th:
            th[j] = (q * th[i] - q * wi * dij / C) / p
        elif j in th and i not in th:
            th[i] = (p * th[j] + q * wi * dij / C) / q
    i, j, p, q = edges[-1]
    dij = float(np.linalg.norm(P[i] - P[j]))
    res = (q * th[i] - q * w_of(loads[i]) * dij / C - p * th[j])
    return (res + np.pi) % (2 * np.pi) - np.pi


def find_triangle(center, exclude, mk_loads, mk_edges, nbest=4000):
    """best-|closure-residual| linked triple with sides in band."""
    order = np.argsort(np.linalg.norm(P - center, axis=1))
    best, seen = None, 0
    for a in order:
        a = int(a)
        if a in exclude or a not in adj:
            continue
        na = [int(n) for n, l in adj[a].items()
              if 1.25 < l < 1.60 and n not in exclude]
        for bi in range(len(na)):
            for ci in range(bi + 1, len(na)):
                b, c = na[bi], na[ci]
                if c not in adj.get(b, {}) or not (1.25 < adj[b][c] < 1.60):
                    continue
                cells = [a, b, c]
                r = residual_of(cells, mk_loads(cells), mk_edges(cells))
                seen += 1
                if best is None or abs(r) < abs(best[1]):
                    best = (cells, r)
                if seen >= nbest and abs(best[1]) < 0.10:
                    return best[0]
        if seen >= 40 * nbest:
            break
    if best is None:
        raise SystemExit("no triangle found")
    return best[0]


def find_hexagon(center, exclude, mk_loads, mk_edges):
    """deterministic restart walks; best closure residual among rings."""
    order = [int(x) for x in np.argsort(np.linalg.norm(P - center, axis=1))]
    def band(u, v):
        return v in adj.get(u, {}) and 1.30 < adj[u][v] < 1.60
    best = None
    tried = 0
    for a in order[:1500]:
        if a in exclude or a not in adj:
            continue
        stack = [(a, [a])]
        expand = 0
        while stack and expand < 30000:
            u, path = stack.pop()
            expand += 1
            if len(path) == 6:
                if band(u, a) and len(set(path)) == 6:
                    cells = path
                    r = residual_of(cells, mk_loads(cells), mk_edges(cells))
                    tried += 1
                    if best is None or abs(r) < abs(best[1]):
                        best = (cells, r)
                continue
            for v in sorted(adj.get(u, {})):
                v = int(v)
                if v in path or v in exclude or not band(u, v):
                    continue
                # no chords: v must not link to earlier ring cells (except
                # the start when closing)
                if any(x in adj.get(v, {}) for x in path[1:-1]):
                    continue
                stack.append((v, path + [v]))
        if best is not None and abs(best[1]) < 0.15:
            break
        if tried > 4000:
            break
    if best is None:
        raise SystemExit("no hexagon found")
    return best[0]


def chain_phases(cells, loads, edges):
    """Set th[0]=0; zero each edge's comb phase in order; last = residual."""
    th = {cells[0]: 0.0}
    for (i, j, p, q) in edges[:-1]:
        dij = float(np.linalg.norm(P[i] - P[j]))
        wi = w_of(loads[i])
        if i in th and j not in th:
            th[j] = (q * th[i] - q * wi * dij / C) / p
        elif j in th and i not in th:
            th[i] = (p * th[j] + q * wi * dij / C) / q
    i, j, p, q = edges[-1]
    dij = float(np.linalg.norm(P[i] - P[j]))
    res = (q * th[i] - q * w_of(loads[i]) * dij / C - p * th[j])
    res = (res + np.pi) % (2 * np.pi) - np.pi
    return th, res


def emit(name, cells, loads, edges, note):
    idx = {c: k for k, c in enumerate(cells)}
    th, res = chain_phases(cells, loads, edges)
    lines = [f"# Q11h {name}: {note} closure_residual={res:+.4f} rad"]
    for c in cells:
        lines.append(f"V {P[c][0]:.4f} {P[c][1]:.4f} {P[c][2]:.4f} "
                     f"{loads[c]:.6f} {th[c] % (2*np.pi):.6f}")
    for (i, j, p, q) in edges:
        lines.append(f"E {idx[i]} {idx[j]} {p} {q}")
    path = f"{OUT}/q11h_{name}.net"
    open(path, "w").write("\n".join(lines) + "\n")
    sides = [f"{np.linalg.norm(P[i]-P[j]):.3f}" for i, j, _, _ in edges]
    print(f"{name}: cells={cells} sides={sides} residual={res:+.4f} -> {path}")


cen = P.mean(axis=0)


def tri_loads_oct(cells):
    a, b, c = cells
    return {a: XH_O, b: XL_O, c: XH_O}
def tri_edges_oct(cells):
    a, b, c = cells
    return [(a, b, 2, 1), (b, c, 1, 2), (c, a, 1, 1)]
def tri_loads_f(cells):
    a, b, c = cells
    return {a: XH_F, b: XL_F, c: XH_F}
def tri_edges_f(cells):
    a, b, c = cells
    return [(a, b, 3, 2), (b, c, 2, 3), (c, a, 1, 1)]
def hex_loads(cells):
    return {u: (XH_F if k % 2 == 0 else XL_F) for k, u in enumerate(cells)}
def hex_edges(cells):
    ee = []
    for k in range(6):
        u, v = cells[k], cells[(k + 1) % 6]
        ee.append((u, v, 3, 2) if k % 2 == 0 else (u, v, 2, 3))
    return ee


tri_o = find_triangle(cen, AVOID, tri_loads_oct, tri_edges_oct)
emit("oct_tri", tri_o, tri_loads_oct(tri_o), tri_edges_oct(tri_o),
     "A-B octave 2:1, B-C 1:2, C-A unison")

tri_f2 = find_triangle(cen + np.array([6.0, 5.0, -4.0]),
                       AVOID | set(tri_o), tri_loads_f, tri_edges_f)
emit("fifth_tri2", tri_f2, tri_loads_f(tri_f2), tri_edges_f(tri_f2),
     "second-site fifth triangle (address independence)")

hexc = find_hexagon(cen - np.array([5.0, -4.0, 3.0]),
                    AVOID | set(tri_o) | set(tri_f2), hex_loads, hex_edges)
emit("fifth_hex", hexc, hex_loads(hexc), hex_edges(hexc),
     "Z6 alternating 3:2/2:3 ring")
