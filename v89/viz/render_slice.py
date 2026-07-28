#!/usr/bin/env python3
# v89 field renderer — the fresh eye (VIZ-1, MASS tech tree D1).
# Renders kernel snapshot TSVs (snap_every / snap_dir) as slab scatter
# panels: dense Em, field Ee, space Es, and dense phase th2 (shown only
# where dense lives). v89-native: no SFA, no volview, no lattice.
#
#   render_slice.py cells_000500.tsv [--z Z0] [--dz 1.5] [--out out.png]
#   render_slice.py --series snapdir/ [--panel em] [--out series.png]
#
# Standing practice (user, 2026-07-28): every MASS/EMF campaign produces
# frames alongside analytics — visuals have caught what numbers missed.

import argparse
import glob
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    d = np.genfromtxt(path, names=True)
    return d


def slab(d, z0, dz):
    if z0 is None:
        z0 = 0.5 * (d["z"].min() + d["z"].max())
    m = np.abs(d["z"] - z0) <= dz
    return d[m], z0


def panel(ax, d, vals, title, cmap, vmin=None, vmax=None, phase=False):
    s = 22.0 * d["r"] ** 2
    if phase:
        live = d["Em"] > 0.02
        ax.scatter(d["x"][~live], d["y"][~live], s=s[~live] * 0.5,
                   c="#dddddd", linewidths=0)
        sc = ax.scatter(d["x"][live], d["y"][live], s=s[live],
                        c=vals[live], cmap=cmap, vmin=-math.pi,
                        vmax=math.pi, linewidths=0)
    else:
        sc = ax.scatter(d["x"], d["y"], s=s, c=vals, cmap=cmap,
                        vmin=vmin, vmax=vmax, linewidths=0)
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    return sc


def render_one(path, z0, dz, out):
    d = load(path)
    d, z0 = slab(d, z0, dz)
    step = os.path.basename(path).split("_")[-1].split(".")[0]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    p = panel(axes[0], d, d["Em"], f"dense Em  (step {step}, z≈{z0:.1f})",
              "inferno", vmin=0)
    plt.colorbar(p, ax=axes[0], fraction=0.046)
    p = panel(axes[1], d, d["Ee"], "field Ee", "viridis", vmin=0)
    plt.colorbar(p, ax=axes[1], fraction=0.046)
    p = panel(axes[2], d, d["Es"], "space Es", "coolwarm",
              vmin=0.4, vmax=1.2)
    plt.colorbar(p, ax=axes[2], fraction=0.046)
    ph = np.mod(d["th2"] + math.pi, 2 * math.pi) - math.pi
    p = panel(axes[3], d, ph, "dense phase th2 (where Em>0.02)",
              "twilight", phase=True)
    plt.colorbar(p, ax=axes[3], fraction=0.046)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"wrote {out}  ({len(d)} cells in slab)")


def render_series(snapdir, which, z0, dz, out):
    paths = sorted(glob.glob(os.path.join(snapdir, "cells_*.tsv")))
    if not paths:
        sys.exit(f"no snapshots in {snapdir}")
    n = len(paths)
    cols = min(n, 6)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.0 * cols, 3.2 * rows),
                             squeeze=False)
    key = {"em": ("Em", "inferno"), "ee": ("Ee", "viridis"),
           "es": ("Es", "coolwarm")}[which]
    vmax = None
    for k, path in enumerate(paths):
        d = load(path)
        d, _ = slab(d, z0, dz)
        if vmax is None:
            vmax = max(1e-9, np.percentile(d[key[0]], 99.9))
        ax = axes[k // cols][k % cols]
        step = os.path.basename(path).split("_")[-1].split(".")[0]
        panel(ax, d, d[key[0]], f"step {step}", key[1], vmin=0 if which != "es" else 0.4,
              vmax=vmax if which != "es" else 1.2)
    for k in range(n, rows * cols):
        axes[k // cols][k % cols].axis("off")
    fig.suptitle(f"{key[0]} series — {snapdir}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"wrote {out}  ({n} frames)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", help="snapshot TSV, or snapshot dir with --series")
    ap.add_argument("--series", action="store_true")
    ap.add_argument("--panel", default="em", choices=["em", "ee", "es"])
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--dz", type=float, default=1.5)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or (a.target.rstrip("/") + (".series.png" if a.series else ".png"))
    if a.series:
        render_series(a.target, a.panel, a.z, a.dz, out)
    else:
        render_one(a.target, a.z, a.dz, out)


if __name__ == "__main__":
    main()
