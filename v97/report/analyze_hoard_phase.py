#!/usr/bin/env python3
"""v93 R3 PHASE AUDIT (RECONCILIATION R3 — the hinge). Input: fcsdump -mode cells.
The natural Face-C condensate is seeded with UNIFORM phase (seed_mw=0); B.4
measured the SEEDED-winding retention and found it dead, but no one has measured
the NATURAL hoard's internal phase. This identifies hoards (spatial clusters of
high-Em cells via union-find at r_cut) and, per hoard, measures:
  R2_m0   |<e^{i th2}>|              over hoard cells (m=0 coherence)
  best-m  argmax_m |<e^{i(th2-m*phi)}>|  about the hoard centroid (m in -4..4)
  R2_best the peak coherence
  W       discrete azimuthal winding of th2 about the centroid (sum wrapped
          increments / 2pi, sorted by phi) — topological charge
Plus the BATH noise floor: R2_m0 over an equal-count random sample of low-Em
cells.  VERDICT: if hoard R2 >> bath floor (and/or W locks to an integer), the
natural hoard carries hidden phase order (Stage 1 achievable in-kernel); if at
floor, q_detune is confirmed phase-blind."""
import sys, math, numpy as np

THRESH = 0.5      # Em threshold for "in a hoard"
R_CUT = 1.7       # clustering distance (dbar ~1.4; just over a cell gap)
MINHOARD = 5      # min cells to score a hoard


def union_find_clusters(xy, rcut):
    n = len(xy)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    # grid for O(N) neighbor search
    cells = {}
    for i, (x, y) in enumerate(xy):
        cells.setdefault((int(x / rcut), int(y / rcut)), []).append(i)
    for (gx, gy), idxs in cells.items():
        nb = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                nb += cells.get((gx + dx, gy + dy), [])
        for i in idxs:
            for j in nb:
                if j <= i:
                    continue
                if (xy[i][0] - xy[j][0])**2 + (xy[i][1] - xy[j][1])**2 < rcut * rcut:
                    parent[find(i)] = find(j)
    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)
    return list(clusters.values())


def main():
    f = open(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "-" else sys.stdin
    f.readline()
    frames = {}
    for line in f:
        p = line.split()
        if len(p) < 14:
            continue
        # t i x y z r es em ee xload tag fa1 fa2 th2
        frames.setdefault(float(p[0]), []).append(
            (float(p[2]), float(p[3]), float(p[7]), float(p[13])))
    tmax = max(frames)
    print(f"(last frame t={tmax})")
    a = np.array(frames[tmax])
    x, y, em, th = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    L = max(x.max() - x.min(), y.max() - y.min()) + 1.0
    wrap = lambda v: (v + L / 2) % L - L / 2

    hmask = em > THRESH
    nh = int(hmask.sum())
    print(f"hoard cells(Em>{THRESH}): {nh}   bath cells: {len(em) - nh}")
    # bath noise floor: R2_m0 over equal-count random low-Em cells
    bath = np.where(~hmask)[0]
    if len(bath) > 100:
        rs = np.random.choice(bath, size=min(nh, len(bath)), replace=False)
        r2bath = abs(np.exp(1j * th[rs]).mean())
        print(f"bath R2_m0 (noise floor, {len(rs)} cells): {r2bath:.4f}")
    if nh < MINHOARD:
        print("too few hoard cells"); return

    xy = list(zip(x[hmask], y[hmask]))
    cl = union_find_clusters(xy, R_CUT)
    cl = [c for c in cl if len(c) >= MINHOARD]
    hx, hy, he, ht = x[hmask], y[hmask], em[hmask], th[hmask]
    print(f"hoards(>={MINHOARD} cells): {len(cl)}")
    rows = []
    for c in cl:
        cx, cyc = hx[c].mean(), hy[c].mean()
        phi = np.arctan2(wrap(hy[c] - cyc), wrap(hx[c] - cx))
        tt = ht[c]
        r2m0 = abs(np.exp(1j * tt).mean())
        bestm, bestr2 = 0, r2m0
        for m in range(-4, 5):
            r = abs(np.exp(1j * (tt - m * phi)).mean())
            if r > bestr2:
                bestr2, bestm = r, m
        # winding
        o = np.argsort(phi)
        tr = tt[o]
        d = np.diff(np.concatenate([tr, tr[:1]]))
        d = (d + np.pi) % (2 * np.pi) - np.pi
        W = d.sum() / (2 * np.pi)
        rows.append((len(c), he[c].sum(), he[c].max(), r2m0, bestm, bestr2, W))
    rows.sort(reverse=True)
    print(f"{'n':>4} {'Em_tot':>7} {'Em_max':>6} {'R2_m0':>6} {'bestm':>5} {'R2_best':>7} {'W':>6}")
    for r in rows[:25]:
        print(f"{r[0]:4d} {r[1]:7.3f} {r[2]:6.3f} {r[3]:6.3f} {r[4]:5d} {r[5]:7.3f} {r[6]:6.2f}")
    import numpy as _np
    arr = _np.array([r[3] for r in rows])
    arrb = _np.array([r[5] for r in rows])
    print(f"R2_m0  median={_np.median(arr):.3f} p90={_np.percentile(arr,90):.3f} max={arr.max():.3f}")
    print(f"R2_best median={_np.median(arrb):.3f} p90={_np.percentile(arrb,90):.3f} max={arrb.max():.3f}")


if __name__ == "__main__":
    main()
