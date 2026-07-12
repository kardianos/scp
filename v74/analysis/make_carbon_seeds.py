#!/usr/bin/env python3
"""Stage 2A: build liquid-drop carbon multi-ball seeds via gen_qball_multi.

Usage:
  python3 make_carbon_seeds.py [--N 192] [--L 36] [--which c6,c12]
"""
from __future__ import annotations

import argparse
import itertools
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V74 = ROOT / "v74"
GEN = ROOT / "bin" / "gen_qball_multi"
PROF_LIGHT = V74 / "profiles" / "f_w146_g005.txt"
W = 1.46
DELTA = 0.0
D = 10.0  # nearest-neighbour edge length


def octahedron_centers(edge: float) -> list[tuple[float, float, float]]:
    """Regular octahedron with edge length `edge` (axial distance edge/sqrt(2))."""
    a = edge / math.sqrt(2.0)
    return [
        (a, 0.0, 0.0),
        (-a, 0.0, 0.0),
        (0.0, a, 0.0),
        (0.0, -a, 0.0),
        (0.0, 0.0, a),
        (0.0, 0.0, -a),
    ]


def icosahedron_centers(edge: float) -> list[tuple[float, float, float]]:
    """Regular icosahedron vertices with edge length `edge`."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    verts: list[tuple[float, float, float]] = []
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            verts.append((0.0, s1, s2 * phi))
            verts.append((s1, s2 * phi, 0.0))
            verts.append((s2 * phi, 0.0, s1))
    md = min(math.dist(u, v) for u, v in itertools.combinations(verts, 2))
    s = edge / md
    return [(s * x, s * y, s * z) for x, y, z in verts]


def build(name: str, centers: list[tuple[float, float, float]], N: int, L: float) -> Path:
    out = V74 / "seeds" / f"{name}_N{N}.sfa"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(GEN), str(N), str(L), str(out)]
    for x, y, z in centers:
        cmd += [str(PROF_LIGHT), str(W), str(DELTA), f"{x:.10f}", f"{y:.10f}", f"{z:.10f}"]
    print("RUN:", " ".join(cmd[:6]), f"... ({len(centers)} balls)")
    subprocess.check_call(cmd)
    print(f"  -> {out}  ({out.stat().st_size / 1e6:.1f} MB)")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=192)
    ap.add_argument("--L", type=float, default=36.0)
    ap.add_argument("--which", default="c6,c12", help="comma list: c6,c12")
    args = ap.parse_args()
    if not GEN.is_file():
        print(f"FATAL: {GEN} missing — make -C sfa install", file=sys.stderr)
        return 1
    if not PROF_LIGHT.is_file():
        print(f"FATAL: profile {PROF_LIGHT} missing", file=sys.stderr)
        return 1

    which = {w.strip() for w in args.which.split(",") if w.strip()}
    if "c6" in which:
        c = octahedron_centers(D)
        print(f"c6 octahedron: {len(c)} balls, axial={D/math.sqrt(2):.4f}, edge D={D}")
        build("c6_light", c, args.N, args.L)
    if "c12" in which:
        c = icosahedron_centers(D)
        print(f"c12 icosahedron: {len(c)} balls, edge D={D}")
        build("c12_light", c, args.N, args.L)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
