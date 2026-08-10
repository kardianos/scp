#!/usr/bin/env python3
"""fig_substrate: (a) one cell with its mode fiber e-chi-m + Es; (b) the prism graph."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.6),
                               gridspec_kw={"width_ratios": [1, 1.6]})
for ax in (ax1, ax2):
    ax.set_aspect("equal"); ax.axis("off")

# ---------------- panel (a): one cell ----------------
ax1.set_xlim(-1.6, 1.8); ax1.set_ylim(-1.9, 2.1)
ax1.add_patch(Circle((0, 0), 1.45, fc="#f2f2f2", ec=GRAY, lw=0.8, zorder=0))
ax1.text(0, -1.72, "one cell (roster site)", ha="center", fontsize=8)
nodes = {"e": (0, 1.0), "chi": (0, 0.1), "m": (0, -0.8)}
lbl = {"e": r"$\psi_e$  field", "chi": r"$\chi$  register", "m": r"$\psi_m$  matter"}
col = {"e": BLUE, "chi": GRAY, "m": ORANGE}
for k, (x, y) in nodes.items():
    ax1.add_patch(Circle((x, y), 0.16, fc="white", ec=col[k], lw=1.6, zorder=3))
    ax1.text(x + 0.26, y, lbl[k], va="center", fontsize=8)
# fiber edges
ax1.annotate("", xy=(0, 0.27), xytext=(0, 0.83),
             arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
ax1.text(-0.12, 0.55, "continuous\n(pass P)", ha="right", fontsize=7)
ax1.annotate("", xy=(0, -0.63), xytext=(0, -0.07),
             arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2,
                             linestyle=(0, (2, 2))))
ax1.text(-0.12, -0.35, "click, in\natoms $\\varepsilon$", ha="right", fontsize=7)
ax1.text(0.95, -1.25, "$E_s$ (space)", fontsize=8, color=INK)

# ---------------- panel (b): prism graph over 4 cells ----------------
ax2.set_xlim(-0.6, 4.4); ax2.set_ylim(-1.6, 1.6)
xs = [0, 1.3, 2.6, 3.9]
ye, yc, ym = 1.0, 0.15, -0.7
# spatial edges
for i in range(3):
    ax2.plot([xs[i] + 0.14, xs[i + 1] - 0.14], [ye, ye], color=BLUE, lw=1.6)
    ax2.plot([xs[i] + 0.14, xs[i + 1] - 0.14], [ym, ym], color=ORANGE, lw=1.6)
for x in xs:
    ax2.plot([x, x], [yc + 0.12, ye - 0.12], color=INK, lw=1.0)
    ax2.plot([x, x], [ym + 0.12, yc - 0.12], color=INK, lw=1.0,
             linestyle=(0, (2, 2)))
    for y, c in ((ye, BLUE), (yc, GRAY), (ym, ORANGE)):
        ax2.add_patch(Circle((x, y), 0.11, fc="white", ec=c, lw=1.5, zorder=3))
ax2.text(1.95, 1.32, "spatial field links (pass F)", color=BLUE,
         ha="center", fontsize=8)
ax2.text(1.95, -1.05, "spatial matter links (pass U)", color=ORANGE,
         ha="center", fontsize=8)
ax2.text(4.05, ye, r"$\psi_e$", fontsize=9, va="center")
ax2.text(4.05, yc, r"$\chi$", fontsize=9, va="center")
ax2.text(4.05, ym, r"$\psi_m$", fontsize=9, va="center")
ax2.text(-0.25, 0.56, "fiber\n(mode)\nlinks", fontsize=7, ha="center")
ax2.text(1.95, -1.45, "cells  (the roster is dynamical: cells are born and die)",
         ha="center", fontsize=7.5, color=GRAY)

fig.tight_layout()
fig.savefig("fig_substrate.pdf")
print("wrote fig_substrate.pdf")
