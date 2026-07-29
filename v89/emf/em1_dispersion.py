#!/usr/bin/env python3
"""EM1 — ω(k) / group-velocity scan (field pulse, apparatus only).

Runs a small set of carrier kx values under laws_V2g + e2-like pulse
apparatus, parses front_speed, writes emf/runs/em1_dispersion.tsv.

  python3 v89/emf/em1_dispersion.py [--threads 2]
"""
import argparse
import os
import re
import subprocess
import sys

ROOT = "/home/d/code/scp"
V89 = os.path.join(ROOT, "v89")
BAT = os.path.join(V89, "battery")
OUT = os.path.join(V89, "emf", "runs")
BIN = os.path.join(V89, "cellfab")

# kx values spanning the battery e2 calibration (kx=1.1) and neighbors
KX = [0.4, 0.7, 1.0, 1.1, 1.4, 1.7, 2.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    laws = open(os.path.join(BAT, "laws_V2g.cfg")).read()
    rows = []
    for kx in KX:
        name = f"em1_kx{kx:.2f}".replace(".", "p")
        cfg = os.path.join(OUT, name + ".cfg")
        log = os.path.join(OUT, name + ".log")
        with open(cfg, "w") as f:
            f.write(laws)
            f.write(f"""
# --- EM1 apparatus (pulse; kx sweep) ---
mode=field
init=pulse
seed=20260727
L=36
T=24
dt=0.02
amp=0.25
sigma=2.5
px=6
kx={kx}
prealign=1
diag_every=50
""")
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
        with open(log, "w") as lf:
            rc = subprocess.call([BIN, cfg], stdout=lf, stderr=subprocess.STDOUT,
                                 env=env)
        text = open(log, errors="replace").read()
        m = re.search(
            r"# RESULT front_speed v=([\d.eE+-]+) v_over_C=([\d.eE+-]+)", text)
        v = float(m.group(1)) if m else None
        vc = float(m.group(2)) if m else None
        rows.append((kx, v, vc, rc))
        print(f"kx={kx:.2f}  v={v}  v/C={vc}  rc={rc}")
    tsv = os.path.join(OUT, "em1_dispersion.tsv")
    with open(tsv, "w") as f:
        f.write("kx\tv\tv_over_C\trc\n")
        for kx, v, vc, rc in rows:
            f.write(f"{kx}\t{v}\t{vc}\t{rc}\n")
    print(f"wrote {tsv}")
    # crude note: massless would be v≈const; quadratic/massive-like: v rises with k
    ok = [r for r in rows if r[1] is not None]
    if len(ok) >= 3:
        print(f"# EM1 note: v range [{min(r[1] for r in ok):.4f}, "
              f"{max(r[1] for r in ok):.4f}] across kx — "
              f"{'k-dependent (massive-like branch)' if max(r[1] for r in ok)-min(r[1] for r in ok) > 0.05 else 'roughly flat'}")


if __name__ == "__main__":
    main()
