#!/usr/bin/env python3
# PRESTRESS run harness — one candidate .net through cellfab, with renders.
# Composes laws_V2g + standard apparatus, runs, renders last snapshot
# (slab + 3d; standing visual practice), and appends a score row.
#
#   run_net.py --name p2_cube15 --net candidates/p2_cube15.net [--T 3000]
#              [--seed 20260727] [--threads 4] [--extra k=v ...] [--laws PATH]

import argparse
import os
import subprocess
import sys

ROOT = "/home/d/code/scp"
RUNS = os.path.join(ROOT, "v89/prestress/runs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--net", required=True)
    ap.add_argument("--T", type=float, default=3000)
    ap.add_argument("--seed", type=int, default=20260727)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--laws", default=os.path.join(ROOT, "v89/battery/laws_V2g.cfg"))
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--snap_every", type=int, default=25000)
    ap.add_argument("--no-render", action="store_true")
    a = ap.parse_args()

    os.makedirs(RUNS, exist_ok=True)
    snapdir = os.path.join(RUNS, a.name + "_snaps")
    os.makedirs(snapdir, exist_ok=True)
    cfg = os.path.join(RUNS, a.name + ".cfg")
    with open(a.laws) as f:
        laws = f.read()
    with open(cfg, "w") as f:
        f.write(laws)
        f.write(f"""
L=24
T={a.T}
dt=0.02
seed={a.seed}
init=net
net_file={a.net}
edge_sink=1.6
lump_diag=1
diag_every=50
snap_every={a.snap_every}
snap_dir={snapdir}
""")
        for kv in a.extra:
            f.write(kv + "\n")

    log = os.path.join(RUNS, a.name + ".log")
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
    with open(log, "w") as lf:
        rc = subprocess.call([os.path.join(ROOT, "v89/cellfab"), cfg],
                             stdout=lf, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        print(f"FAILED rc={rc} {a.name} — see {log}")
        sys.exit(rc)

    if not a.no_render:
        snaps = sorted(p for p in os.listdir(snapdir) if p.startswith("cells_"))
        if snaps:
            last = os.path.join(snapdir, snaps[-1])
            for mode, out in ((["--mode3d"], a.name + ".3d.png"),
                              ([], a.name + ".slab.png")):
                subprocess.call([sys.executable,
                                 os.path.join(ROOT, "v89/viz/render_slice.py"),
                                 last, "--out", os.path.join(RUNS, out)] + mode,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    subprocess.call([sys.executable,
                     os.path.join(ROOT, "v89/prestress/score_net.py"), log])
    print(f"done {a.name}")


if __name__ == "__main__":
    main()
