#!/usr/bin/env python3
"""render_slices.py — render sfa_slice output to a PNG panel grid.

Rows = derived quantities, columns = frames. Per-quantity normalization is
shared across frames so time evolution is visually comparable.

  rho2, s, th2 : log10 color scale (these span decades)
  rhoQ         : symmetric diverging scale about 0 (charge sign visible)
  Emag         : linear

Usage:
  render_slices.py prefix [prefix2 ...] [--out out.png] [--max-frames K]

`prefix` is the --out value given to sfa_slice (expects prefix.meta.txt and
prefix.f32). With multiple prefixes, each gets its own figure saved as
<prefix>.png (the --out option only applies to a single prefix).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm


def load(prefix):
    meta = {}
    with open(prefix + ".meta.txt") as f:
        for line in f:
            k, _, v = line.strip().partition("=")
            meta[k] = v
    N = int(meta["N"]); L = float(meta["L"]); nq = int(meta["nq"])
    quants = meta["quantities"].split(",")
    nf = int(meta["nframes"])
    times = [float(t) for t in meta["times"].split(",")]
    data = np.fromfile(prefix + ".f32", dtype=np.float32)
    data = data.reshape(nf, nq, N, N)
    return meta, N, L, quants, nf, times, data


def render(prefix, out=None, max_frames=8):
    meta, N, L, quants, nf, times, data = load(prefix)
    axis = meta["axis"]
    plane_axes = {"z": ("x", "y"), "y": ("x", "z"), "x": ("y", "z")}[axis]

    if nf > max_frames:
        sel = np.unique(np.linspace(0, nf - 1, max_frames).astype(int))
    else:
        sel = np.arange(nf)
    nfs = len(sel)
    nq = len(quants)

    fig, axes = plt.subplots(nq, nfs, figsize=(2.3 * nfs + 1.2, 2.3 * nq + 0.8),
                             squeeze=False, constrained_layout=True)
    extent = [-L, L, -L, L]

    for qi, qname in enumerate(quants):
        qdata = data[:, qi]
        finite = qdata[np.isfinite(qdata)]
        if qname.startswith("rho2") or qname in ("s", "th2"):
            vmax = float(np.percentile(finite, 99.9)) if finite.size else 1.0
            vmax = max(vmax, 1e-12)
            vmin = vmax * 1e-6
            norm = LogNorm(vmin=vmin, vmax=vmax)
            cmap = "inferno"
        elif qname == "rhoQ":
            vamp = float(np.percentile(np.abs(finite), 99.9)) if finite.size else 1.0
            vamp = max(vamp, 1e-12)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-vamp, vmax=vamp)
            cmap = "RdBu_r"
        else:  # Emag
            vmax = float(np.percentile(finite, 99.9)) if finite.size else 1.0
            norm = plt.Normalize(vmin=0.0, vmax=max(vmax, 1e-12))
            cmap = "viridis"

        im = None
        for ci, fi in enumerate(sel):
            ax = axes[qi][ci]
            img = qdata[fi].copy()
            if qname.startswith("rho2") or qname in ("s", "th2"):
                img = np.clip(img, norm.vmin, None)
            img = np.nan_to_num(img, nan=0.0, posinf=0.0, neginf=0.0)
            if qname.startswith("rho2") or qname in ("s", "th2"):
                img[img <= 0] = norm.vmin
            im = ax.imshow(img, origin="lower", extent=extent, norm=norm,
                           cmap=cmap, interpolation="nearest")
            if qi == 0:
                ax.set_title(f"t={times[fi]:.0f}", fontsize=9)
            if ci == 0:
                ax.set_ylabel(f"{qname}\n{plane_axes[1]}", fontsize=9)
            else:
                ax.set_yticklabels([])
            if qi == nq - 1:
                ax.set_xlabel(plane_axes[0], fontsize=8)
            else:
                ax.set_xticklabels([])
            ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=[axes[qi][c] for c in range(nfs)],
                     shrink=0.85, aspect=12, pad=0.01)

    fig.suptitle(f"{prefix.split('/')[-1]}  (slice {axis}={float(meta['slice_pos']):.2f},"
                 f" N={N}, L={L:g})", fontsize=11)
    out_path = out or (prefix + ".png")
    fig.savefig(out_path, dpi=105)
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    out = None
    max_frames = 8
    prefixes = []
    i = 0
    while i < len(args):
        if args[i] == "--out":
            out = args[i + 1]; i += 2
        elif args[i] == "--max-frames":
            max_frames = int(args[i + 1]); i += 2
        else:
            prefixes.append(args[i]); i += 1
    if not prefixes:
        print(__doc__); sys.exit(1)
    for p in prefixes:
        render(p, out=out if len(prefixes) == 1 else None, max_frames=max_frames)
