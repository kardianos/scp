#!/usr/bin/env python3
"""Progressive-cell meters (v93/theory/PROGRESSIVE_CELL.md).

Input: fcsdump -mode cells text (file or stdin).
Rebuilds a contact graph from positions+radii (spatial bins), reconstructs
dense psi = sqrt(Em)*exp(i*th2), reports pattern-vs-capsule meters.

Usage:
  fcsdump -mode cells run.fcs | python3 analyze_progressive.py [L] [tau]
  python3 analyze_progressive.py dump.txt [L] [tau]

Optional env:
  PROG_STRIDE=k   process every k-th frame only (default 1)
"""
from __future__ import annotations

import math
import os
import sys
from collections import defaultdict

import numpy as np

L = 16.0
TAU = 0.15


def parse_args():
    global L, TAU
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    path = None
    nums = []
    for a in args:
        try:
            nums.append(float(a))
        except ValueError:
            path = a
    if len(nums) >= 1:
        L = nums[0]
    if len(nums) >= 2:
        TAU = nums[1]
    return path


def wr(dx, box):
    return (dx + 0.5 * box) % box - 0.5 * box


def load_frames(f):
    f.readline()
    frames = defaultdict(list)
    for line in f:
        p = line.split()
        if len(p) < 14:
            continue
        t = float(p[0])
        frames[t].append(
            (
                int(p[1]),
                float(p[2]),
                float(p[3]),
                float(p[4]),
                float(p[5]),
                float(p[7]),
                float(p[8]),
                float(p[10]),
                float(p[11]),
                float(p[12]),
                float(p[13]),
            )
        )
    return frames


def build_contacts(x, y, z, rad, box):
    """Return list of (i,j,d,dx,dy,dz) contact pairs via spatial bins."""
    n = len(x)
    max_r = float(rad.max()) if n else 1.0
    bin_sz = max(2.1 * max_r, 0.5)
    g = max(int(math.floor(box / bin_sz)), 1)
    bin_sz = box / g
    heads = defaultdict(list)
    for i in range(n):
        bx = int(x[i] / bin_sz) % g
        by = int(y[i] / bin_sz) % g
        bz = int(z[i] / bin_sz) % g
        heads[(bx, by, bz)].append(i)

    pairs = []
    neighbors = [[] for _ in range(n)]
    for (bx, by, bz), members in heads.items():
        for dxb in (-1, 0, 1):
            for dyb in (-1, 0, 1):
                for dzb in (-1, 0, 1):
                    key = ((bx + dxb) % g, (by + dyb) % g, (bz + dzb) % g)
                    other = heads.get(key)
                    if not other:
                        continue
                    for i in members:
                        for j in other:
                            if j <= i:
                                continue
                            ddx = wr(x[j] - x[i], box)
                            ddy = wr(y[j] - y[i], box)
                            ddz = wr(z[j] - z[i], box)
                            d2 = ddx * ddx + ddy * ddy + ddz * ddz
                            contact = rad[i] + rad[j]
                            if d2 >= contact * contact or d2 <= 1e-20:
                                continue
                            d = math.sqrt(d2)
                            pairs.append((i, j, d, ddx, ddy, ddz))
                            neighbors[i].append(j)
                            neighbors[j].append(i)
    return pairs, neighbors


def connected_components(neighbors, mask):
    n = len(mask)
    seen = [False] * n
    comps = []
    for s in range(n):
        if seen[s] or not mask[s]:
            continue
        stack = [s]
        seen[s] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in neighbors[u]:
                if not seen[v] and mask[v]:
                    seen[v] = True
                    stack.append(v)
        comps.append(comp)
    return comps


def analyze_frame(rows, box, tau, prev_support):
    n = len(rows)
    idx = np.array([r[0] for r in rows], dtype=np.int32)
    x = np.array([r[1] for r in rows])
    y = np.array([r[2] for r in rows])
    z = np.array([r[3] for r in rows])
    rad = np.array([r[4] for r in rows])
    em = np.array([r[5] for r in rows])
    ee = np.array([r[6] for r in rows])
    tag = np.array([r[7] for r in rows])
    fa1 = np.array([r[8] for r in rows])
    fa2 = np.array([r[9] for r in rows])
    th2 = np.array([r[10] for r in rows])

    em_tot = float(em.sum())
    ee_tot = float(ee.sum())
    em_max = float(em.max()) if n else 0.0
    pr = (em_tot * em_tot / float((em * em).sum())) if em_tot > 0 else 0.0
    f_tag = float(em[tag > 0.5].sum()) / em_tot if em_tot > 0 else 0.0

    pairs, neighbors = build_contacts(x, y, z, rad, box)
    J_acc = gphi_acc = w_acc = 0.0
    Jf_acc = wf_acc = 0.0
    for i, j, d, *_ in pairs:
        ri, rj = em[i], em[j]
        if ri > 1e-12 and rj > 1e-12:
            dphi = math.atan2(math.sin(th2[j] - th2[i]), math.cos(th2[j] - th2[i]))
            J = 2.0 * math.sqrt(ri * rj) * math.sin(dphi)
            w = math.sqrt(ri * rj)
            J_acc += abs(J) * w
            gphi_acc += abs(dphi) / d * w
            w_acc += w
        ei = fa1[i] * fa1[i] + fa2[i] * fa2[i]
        ej = fa1[j] * fa1[j] + fa2[j] * fa2[j]
        if ei > 1e-14 and ej > 1e-14:
            Jf = 2.0 * (fa1[i] * fa2[j] - fa2[i] * fa1[j])
            wf = math.sqrt(ei * ej)
            Jf_acc += abs(Jf) * wf
            wf_acc += wf

    J_mean = J_acc / w_acc if w_acc > 0 else 0.0
    gphi_mean = gphi_acc / w_acc if w_acc > 0 else 0.0
    Jf_mean = Jf_acc / wf_acc if wf_acc > 0 else 0.0

    mask = em >= tau
    comps = connected_components(neighbors, mask)
    n_comp = len(comps)
    em_cmax = r2_cmax = 0.0
    size_cmax = 0
    if comps:
        best = max(comps, key=lambda c: float(em[c].sum()))
        size_cmax = len(best)
        em_c = em[best]
        em_cmax = float(em_c.sum())
        s = em_cmax
        if s > 0:
            cx = float(np.sum(em_c * np.cos(th2[best])))
            sx = float(np.sum(em_c * np.sin(th2[best])))
            r2_cmax = math.hypot(cx, sx) / s

    support = set(int(idx[i]) for i in range(n) if mask[i])
    if prev_support is None or not support or not prev_support:
        jaccard = 1.0 if support else 0.0
    else:
        inter = len(support & prev_support)
        union = len(support | prev_support)
        jaccard = inter / union if union else 0.0

    c = 0.5 * box

    def com(weights, use_tag=False):
        w = weights * (tag > 0.5) if use_tag else weights
        tw = float(w.sum())
        if tw <= 0:
            return (float("nan"),) * 3
        dx = np.array([wr(xx - c, box) for xx in x])
        dy = np.array([wr(yy - c, box) for yy in y])
        dz = np.array([wr(zz - c, box) for zz in z])
        return (
            c + float(np.average(dx, weights=w)),
            c + float(np.average(dy, weights=w)),
            c + float(np.average(dz, weights=w)),
        )

    com_all = com(em, False)
    com_tag = com(em, True)
    dcom = float("nan")
    if not math.isnan(com_all[0]) and not math.isnan(com_tag[0]):
        dcom = math.sqrt(
            wr(com_all[0] - com_tag[0], box) ** 2
            + wr(com_all[1] - com_tag[1], box) ** 2
            + wr(com_all[2] - com_tag[2], box) ** 2
        )

    ee_mask = ee > 1e-8
    if ee_mask.any() and float(ee[ee_mask].sum()) > 0:
        w = ee[ee_mask]
        ph = np.arctan2(fa2[ee_mask], fa1[ee_mask])
        r2f = math.hypot(
            float(np.sum(w * np.cos(ph))), float(np.sum(w * np.sin(ph)))
        ) / float(w.sum())
    else:
        r2f = 0.0

    return {
        "em_tot": em_tot,
        "ee_tot": ee_tot,
        "em_max": em_max,
        "pr": pr,
        "f_tag": f_tag,
        "n_comp": n_comp,
        "em_cmax": em_cmax,
        "size_cmax": size_cmax,
        "r2_cmax": r2_cmax,
        "J_mean": J_mean,
        "gphi_mean": gphi_mean,
        "Jf_mean": Jf_mean,
        "r2f": r2f,
        "jaccard": jaccard,
        "n_links": len(pairs),
        "com_all": com_all,
        "com_tag": com_tag,
        "dcom": dcom,
        "support": support,
    }


def main():
    path = parse_args()
    stride = max(int(os.environ.get("PROG_STRIDE", "1")), 1)
    f = open(path) if path else sys.stdin
    frames = load_frames(f)
    if path:
        f.close()

    ts = sorted(frames)
    if stride > 1:
        ts = ts[::stride]

    print(f"# progressive meters L={L} tau={TAU} n_frames={len(ts)} stride={stride}")
    print(
        "t Em_tot Em_max PR F_tag n_comp Em_cmax size_cmax R2_cmax "
        "J_mean gphi Jf_mean R2f Jaccard n_links "
        "comx_all comx_tag dcom"
    )

    prev_support = None
    rows_out = []
    for t in ts:
        m = analyze_frame(frames[t], L, TAU, prev_support)
        prev_support = m["support"]
        ca, ct = m["com_all"], m["com_tag"]
        print(
            f"{t:7.2f} {m['em_tot']:9.4f} {m['em_max']:7.4f} {m['pr']:7.1f} "
            f"{m['f_tag']:6.3f} {m['n_comp']:4d} {m['em_cmax']:8.3f} "
            f"{m['size_cmax']:4d} {m['r2_cmax']:6.3f} "
            f"{m['J_mean']:7.4f} {m['gphi_mean']:7.4f} {m['Jf_mean']:7.4f} "
            f"{m['r2f']:6.3f} {m['jaccard']:6.3f} {m['n_links']:6d} "
            f"{ca[0]:7.3f} {ct[0]:7.3f} {m['dcom']:6.3f}"
        )
        rows_out.append((t, m))

    if len(rows_out) >= 2:
        t0, m0 = rows_out[0]
        t1, m1 = rows_out[-1]
        # mean jaccard excluding first
        jacs = [m["jaccard"] for _, m in rows_out[1:]]
        mean_j = sum(jacs) / len(jacs) if jacs else 0.0
        print(
            f"# SUMMARY t0={t0} t1={t1} "
            f"F_tag {m0['f_tag']:.3f}->{m1['f_tag']:.3f} "
            f"PR {m0['pr']:.1f}->{m1['pr']:.1f} "
            f"Em_max {m0['em_max']:.3f}->{m1['em_max']:.3f} "
            f"R2 {m0['r2_cmax']:.3f}->{m1['r2_cmax']:.3f} "
            f"mean_Jaccard={mean_j:.3f} "
            f"J_mean {m0['J_mean']:.4f}->{m1['J_mean']:.4f} "
            f"n_comp {m0['n_comp']}->{m1['n_comp']}"
        )


if __name__ == "__main__":
    main()
