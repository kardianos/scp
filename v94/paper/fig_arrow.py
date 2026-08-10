#!/usr/bin/env python3
"""fig_arrow: complex amplitude as an arrow; multiplying by i = quarter-turn; energy = radius^2."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY})

fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.15))
for ax in axes:
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th), np.sin(th), color=GRAY, lw=0.6, alpha=0.6)
    ax.annotate("", xy=(1.4, 0), xytext=(-1.4, 0),
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.6))
    ax.annotate("", xy=(0, 1.4), xytext=(0, -1.4),
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.6))

# (a) the arrow: length + clock angle
a = np.deg2rad(30)
axes[0].annotate("", xy=(np.cos(a), np.sin(a)), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.0))
axes[0].plot([0.35 * np.cos(t) for t in np.linspace(0, a, 40)],
             [0.35 * np.sin(t) for t in np.linspace(0, a, 40)],
             color=INK, lw=0.9)
axes[0].text(0.42, 0.10, "phase\n(clock angle)", fontsize=7)
axes[0].text(0.28, 0.75, "length $r$", color=BLUE, fontsize=8, rotation=30)
axes[0].text(0, -1.32, "(a) an amplitude is an arrow:\na length and a clock angle",
             ha="center", fontsize=7.5)

# (b) multiply by i = quarter-turn
b = np.deg2rad(30)
axes[1].annotate("", xy=(np.cos(b), np.sin(b)), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.0))
axes[1].annotate("", xy=(np.cos(b + np.pi/2), np.sin(b + np.pi/2)), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.0))
arc = np.linspace(b + 0.10, b + np.pi/2 - 0.10, 50)
axes[1].plot(0.72 * np.cos(arc), 0.72 * np.sin(arc), color=INK, lw=0.9)
axes[1].annotate("", xy=(0.72 * np.cos(arc[-1]), 0.72 * np.sin(arc[-1])),
                 xytext=(0.72 * np.cos(arc[-3]), 0.72 * np.sin(arc[-3])),
                 arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
axes[1].text(0.86, 0.48, r"$\psi$", color=BLUE, fontsize=9)
axes[1].text(-0.72, 0.86, r"$i\,\psi$", color=ORANGE, fontsize=9)
axes[1].text(0, -1.32, "(b) multiplying by $i$ turns the\narrow a quarter-turn ($90^\\circ$)",
             ha="center", fontsize=7.5)

# (c) energy = radius squared
c = np.deg2rad(55)
axes[2].annotate("", xy=(np.cos(c), np.sin(c)), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.0))
axes[2].plot([np.cos(c), np.cos(c)], [0, np.sin(c)], color=GRAY, lw=0.9,
             linestyle=(0, (2, 2)))
axes[2].plot([0, np.cos(c)], [np.sin(c), np.sin(c)], color=GRAY, lw=0.9,
             linestyle=(0, (2, 2)))
axes[2].text(np.cos(c) + 0.06, 0.24, r"$f_{a1}$", fontsize=8)
axes[2].text(0.10, np.sin(c) + 0.08, r"$f_{a2}$", fontsize=8)
axes[2].text(0, -1.32, "(c) energy = (length)$^2$\n$= f_{a1}^2 + f_{a2}^2$ — never negative",
             ha="center", fontsize=7.5)

fig.tight_layout()
fig.savefig("fig_arrow.pdf")
print("wrote fig_arrow.pdf")
