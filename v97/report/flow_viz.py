#!/usr/bin/env python3
"""FLOW campaign analysis + figures (FLOW.md; deliverable = exhibit).
Reads arm logs (# BED / # HOLE / # PAR / # BEDMAP lines) and fcs cell
dumps (via fcsdump text piped to files), renders PNGs.
Palette: the dataviz reference instance — diverging blue#2a78d6 /
gray#f0efec / red#e34948 for bed weight (neutral = 1); sequential
blues for well depth; categorical order blue/orange/aqua for series.
Usage: flow_viz.py <outdir>  (expects runs/flow/*.log and
pre-dumped   runs/flow/<arm>.cells  text files)"""
import sys, re, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection

BLUE, ORANGE, AQUA, RED, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#f0efec"
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
DIV = LinearSegmentedColormap.from_list("bed", [BLUE, GRAY, RED])
plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.dpi": 130})

RD = "runs/flow"

def read_cells(arm):
    frames = {}
    with open(f"{RD}/{arm}.cells") as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) < 14:
                continue
            frames.setdefault(float(p[0]), []).append(
                (float(p[2]), float(p[3]), float(p[6]), float(p[7]),
                 float(p[8])))          # x y es em ee
    return {t: np.array(v) for t, v in frames.items()}

def read_bedmap(arm):
    beds = {}
    for line in open(f"{RD}/{arm}.log"):
        m = re.match(r"# BEDMAP t=([\d.]+) (\d+) (\d+) ([\d.]+)", line)
        if m:
            beds.setdefault(float(m.group(1)), []).append(
                (int(m.group(2)), int(m.group(3)), float(m.group(4))))
    return beds

def read_series(arm, tag, fields):
    out = {k: [] for k in ["t"] + fields}
    for line in open(f"{RD}/{arm}.log"):
        if not line.startswith(f"# {tag} "):
            continue
        m = dict(re.findall(r"(\w+)=\[?([\d.eE+-]+)", line))
        q = re.search(r"q=\[([\d.]+) ([\d.]+) ([\d.]+)\]", line)
        out["t"].append(float(m["t"]))
        for k in fields:
            if k in ("q25", "q50", "q75") and q:
                out[k].append(float(q.group(("q25", "q50", "q75").index(k) + 1)))
            elif k in m:
                out[k].append(float(m[k]))
            else:
                out[k].append(np.nan)
    return {k: np.array(v) for k, v in out.items()}

def bed_map_panel(ax, arm, t, title, hole_r=0, hole_sep=0, L=64):
    fr = read_cells(arm)
    tk = min(fr.keys(), key=lambda x: abs(x - t))
    a = fr[tk]
    x, y, es, em = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    ax.scatter(x, y, s=6, c=es, cmap="Blues_r", vmin=0, vmax=1.2,
               linewidths=0, alpha=0.75, rasterized=True)
    live = em >= 0.05
    ax.scatter(x[live], y[live], s=10, c=ORANGE, linewidths=0, alpha=0.8,
               rasterized=True)
    beds = read_bedmap(arm)
    if beds:
        bk = min(beds.keys(), key=lambda k: abs(k - tk))
        segs, cols, lws = [], [], []
        for i, j, w in beds[bk]:
            if i >= len(x) or j >= len(x):
                continue
            dx = x[j] - x[i]; dy = y[j] - y[i]
            if abs(dx) > L / 2 or abs(dy) > L / 2:
                continue                      # skip wrapped segments
            segs.append([(x[i], y[i]), (x[j], y[j])])
            cols.append(DIV((np.clip(w, 0.2, 3.0) - 1) / 4 + 0.5))
            lws.append(0.6 + 2.4 * min(abs(w - 1), 2))
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws))
        n = len(segs)
    else:
        n = 0
    cx = L / 2
    for sx in ([cx] if hole_sep == 0 else [cx - hole_sep / 2, cx + hole_sep / 2]):
        if hole_r > 0:
            ax.add_patch(plt.Circle((sx, cx), hole_r, fill=False,
                                    color=INK, lw=1.2, ls="--"))
    ax.set_xlim(8, L - 8); ax.set_ylim(8, L - 8)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{title}  (t={tk:.0f}, {n} beds shown)", fontsize=10)

def radial_profile(arm, t, L=64, what="pi"):
    fr = read_cells(arm)
    tk = min(fr.keys(), key=lambda x: abs(x - t))
    a = fr[tk]
    x, y, es, em, ee = (a[:, i] for i in range(5))
    r = np.hypot(x - L / 2, y - L / 2)
    val = es + 0.3 * (em + ee)
    bins = np.arange(0, 26, 2)
    prof = [val[(r >= b) & (r < b + 2)].mean() if ((r >= b) & (r < b + 2)).sum() > 3
            else np.nan for b in bins[:-1]]
    return bins[:-1] + 1, np.array(prof)

def bed_anisotropy(arm, t, L=64):
    """mean sbed of radial-ish vs tangential-ish links by shell"""
    fr = read_cells(arm)
    tk = min(fr.keys(), key=lambda x: abs(x - t))
    a = fr[tk]; x, y = a[:, 0], a[:, 1]
    beds = read_bedmap(arm)
    if not beds:
        return None
    bk = min(beds.keys(), key=lambda k: abs(k - tk))
    shells = np.arange(4, 22, 3)
    rad = {s: [] for s in shells}; tan = {s: [] for s in shells}
    for i, j, w in beds[bk]:
        if i >= len(x) or j >= len(x):
            continue
        mx, my = (x[i] + x[j]) / 2 - L / 2, (y[i] + y[j]) / 2 - L / 2
        rr = np.hypot(mx, my)
        lx, ly = x[j] - x[i], y[j] - y[i]
        if abs(lx) > L / 2 or abs(ly) > L / 2 or (lx == 0 and ly == 0):
            continue
        cosr = abs((mx * lx + my * ly) / (np.hypot(lx, ly) * (rr + 1e-9)))
        for s in shells:
            if s <= rr < s + 3:
                (rad if cosr > 0.7 else tan if cosr < 0.3 else {s: []})[s] \
                    .append(w) if cosr > 0.7 or cosr < 0.3 else None
    out = []
    for s in shells:
        out.append((s + 1.5,
                    np.mean(rad[s]) if rad[s] else np.nan,
                    np.mean(tan[s]) if tan[s] else np.nan))
    return np.array(out)

def main(outdir):
    os.makedirs(outdir, exist_ok=True)

    # FIG 1 — the eater digs its river (FH-1 epochs)
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.4))
    for ax, t in zip(axs, [100, 600, 1500]):
        bed_map_panel(ax, "fh1", t, "eater + beds", hole_r=4)
    fig.suptitle("FH-1: a feeding hole digs its bed — space-channel weights "
                 "(red = grown, blue = starved), matter orange, Es shading", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig1_eater_beds.png"); plt.close(fig)

    # FIG 2 — two eaters: the corridor (FH-2 vs control)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 5.0))
    bed_map_panel(axs[0], "fh2", 1500, "two eaters + beds", hole_r=3, hole_sep=16)
    bed_map_panel(axs[1], "fh2c", 1500, "control (bed off)", hole_r=3, hole_sep=16)
    fig.suptitle("FH-2: do two eaters grow a channel between them?", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig2_double.png"); plt.close(fig)

    # FIG 3 — quench with beds: lawful matter stays bed-free
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.6))
    bed_map_panel(axs[0], "fq1", 1500, "quench cloud + beds (P-F3)")
    ax = axs[1]
    for arm, col, lab in [("fq1", BLUE, "bed on")]:
        s = read_series(arm, "PAR", ["gids"])
        # PAR line: gids n=...
        ts, ns = [], []
        for line in open(f"{RD}/{arm}.log"):
            m = re.match(r"# PAR t=([\d.]+) gids n=(\d+)", line)
            if m:
                ts.append(float(m.group(1))); ns.append(int(m.group(2)))
        ax.plot(ts, ns, color=col, lw=2, label=lab)
    ts, ns = [], []
    for line in open("runs/quench/q2_beam8.log"):
        m = re.match(r"# PAR t=([\d.]+) gids n=(\d+)", line)
        if m:
            ts.append(float(m.group(1))); ns.append(int(m.group(2)))
    ax.plot(ts, ns, color=ORANGE, lw=2, ls="--", label="bed off (archived Q2)")
    ax.set_xlabel("t"); ax.set_ylabel("live episodes")
    ax.legend(frameon=False); ax.set_title("cloud census, bed on vs off", fontsize=10)
    fig.suptitle("FQ-1: the quench cloud under the channel law", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig3_quench.png"); plt.close(fig)

    # FIG 4 — growth curves: bed feedback + anti-ignition
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 3.9))
    ax = axs[0]
    for arm, col, lab in [("fh1", BLUE, "FH-1 bed on")]:
        s = read_series(arm, "HOLE", ["Eh"])
        ax.plot(s["t"], s["Eh"], color=col, lw=2, label=lab)
    s0 = read_series("hz0_ref", "HOLE", ["Eh"]) if os.path.exists(f"{RD}/hz0_ref.log") else None
    ts, es_ = [], []
    for line in open("runs/horizon/hz0.log"):
        m = re.match(r"# HOLE t=([\d.]+) Eh=([\d.]+)", line)
        if m:
            ts.append(float(m.group(1))); es_.append(float(m.group(2)))
    ax.plot(ts, es_, color=ORANGE, lw=2, ls="--", label="HZ-0 bed off (archived)")
    ax.set_xlabel("t"); ax.set_ylabel("Eh swallowed")
    ax.legend(frameon=False); ax.set_title("P-F4: does the bed feed the river?", fontsize=10)
    ax = axs[1]
    for arm, col, lab in [("fh1", BLUE, "FH-1 eater"), ("fq1", AQUA, "FQ-1 quench"),
                          ("fwb", ORANGE, "FW-B bath")]:
        s = read_series(arm, "BED", ["q50", "q75", "max"])
        if len(s["t"]):
            ax.plot(s["t"], s["max"], color=col, lw=2, label=lab)
    ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("t"); ax.set_ylabel("max sbed")
    ax.legend(frameon=False); ax.set_title("bed maxima: eater vs cloud vs bath", fontsize=10)
    ax = axs[2]
    for arm, col, lab in [("fh1", BLUE, "FH-1"), ("fwb", ORANGE, "FW-B bath")]:
        s = read_series(arm, "BED", ["q25", "q75"])
        if len(s["t"]):
            ax.fill_between(s["t"], s["q25"], s["q75"], color=col, alpha=0.25,
                            linewidth=0)
            ax.plot(s["t"], s["q75"], color=col, lw=1.5, label=f"{lab} q25–q75")
    ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("t"); ax.set_ylabel("sbed quartile band")
    ax.legend(frameon=False)
    ax.set_title("P-F2 anti-ignition: the bath must hug 1", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig4_growth.png"); plt.close(fig)

    # FIG 5 — the well + anisotropy
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 3.9))
    ax = axs[0]
    r, p = radial_profile("fh1", 1500)
    ax.plot(r, p, color=BLUE, lw=2, label="FH-1 bed on, t=1500")
    r0, p0 = None, None
    try:
        rr, pp = [], []
        # archived hz0 profile from its fcs is not re-dumped here; skip
    except Exception:
        pass
    ax.set_xlabel("r"); ax.set_ylabel("π(r)")
    ax.set_title("the standing well under the channel law", fontsize=10)
    ax.legend(frameon=False)
    ax = axs[1]
    an = bed_anisotropy("fh1", 1500)
    if an is not None:
        ax.plot(an[:, 0], an[:, 1], color=RED, lw=2, marker="o", ms=5,
                label="radial links")
        ax.plot(an[:, 0], an[:, 2], color=BLUE, lw=2, marker="o", ms=5,
                label="tangential links")
        ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("shell r"); ax.set_ylabel("mean sbed")
    ax.set_title("P-F4: radial anisotropy of the dug beds", fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{outdir}/fig5_well.png"); plt.close(fig)
    print("figures written to", outdir)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs/flow/figs")
