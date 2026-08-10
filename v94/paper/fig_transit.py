#!/usr/bin/env python3
"""fig_transit: energy in transit — flight parcel on a link (space-transit)
vs register inside a cell (mode-transit)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.4, 2.9))
for ax in (ax1, ax2):
    ax.set_aspect("equal"); ax.axis("off")

# (a) flight: parcel mid-link between two cells
ax1.set_xlim(-0.4, 6.4); ax1.set_ylim(-1.15, 1.3)
for cx, name in ((0.8, "cell A"), (5.2, "cell B")):
    ax1.add_patch(Circle((cx, 0), 0.62, fc="#f2f2f2", ec=GRAY, lw=1.0))
    ax1.text(cx, -0.95, name, ha="center", fontsize=7.5)
ax1.plot([1.42, 4.58], [0, 0], color=GRAY, lw=1.0)
ax1.add_patch(Circle((2.9, 0), 0.17, fc=ORANGE, ec="none", zorder=4))
ax1.annotate("", xy=(3.7, 0.0), xytext=(3.15, 0.0),
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.3))
ax1.text(2.9, 0.42, "flight parcel:\nenergy crossing SPACE,\nbelongs to neither cell yet",
         ha="center", fontsize=7, color=INK)
ax1.text(-0.35, 1.12, "(a) space-transit (a link)", fontsize=8)

# (b) register: energy crossing between modes inside one cell
ax2.set_xlim(-0.4, 6.4); ax2.set_ylim(-1.5, 1.35)
ax2.add_patch(Circle((2.9, 0), 1.28, fc="#f2f2f2", ec=GRAY, lw=1.0))
ax2.text(2.9, -1.42, "one cell", ha="center", fontsize=7.5)
ax2.add_patch(Circle((1.95, 0.45), 0.20, fc="white", ec=BLUE, lw=1.6, zorder=3))
ax2.text(1.13, 0.72, "field", color=BLUE, fontsize=7.5)
ax2.add_patch(Circle((3.85, -0.45), 0.20, fc="white", ec=INK, lw=1.6, zorder=3))
ax2.text(4.18, -0.55, "matter", color=INK, fontsize=7.5)
ax2.add_patch(Circle((2.9, 0), 0.17, fc=ORANGE, ec="none", zorder=4))
ax2.annotate("", xy=(3.55, -0.30), xytext=(2.22, 0.36),
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
ax2.text(4.55, 0.55, "register $\\chi$:\nenergy crossing\nbetween MODES,\nnot settled yet",
         ha="left", fontsize=7)
ax2.text(-0.35, 1.18, "(b) mode-transit (the register)", fontsize=8)

fig.tight_layout()
fig.savefig("fig_transit.pdf")
print("wrote fig_transit.pdf")
