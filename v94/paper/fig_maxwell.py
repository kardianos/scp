#!/usr/bin/env python3
"""fig_maxwell: (a) the two quadratures inducing each other; (b) a travelling packet."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.4),
                               gridspec_kw={"width_ratios": [1, 1.8]})

# panel (a): rotation in the (fa1, fa2) plane
ax1.set_aspect("equal"); ax1.set_xlim(-1.5, 1.5); ax1.set_ylim(-1.5, 1.5)
th = np.linspace(0, 2 * np.pi, 200)
ax1.plot(np.cos(th), np.sin(th), color=GRAY, lw=0.8)
ax1.annotate("", xy=(1.45, 0), xytext=(-1.45, 0),
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.7))
ax1.annotate("", xy=(0, 1.45), xytext=(0, -1.45),
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.7))
a = np.deg2rad(35)
ax1.annotate("", xy=(np.cos(a), np.sin(a)), xytext=(0, 0),
             arrowprops=dict(arrowstyle="->", color=INK, lw=1.6))
ax1.annotate("", xy=(np.cos(a + 0.65), np.sin(a + 0.65)),
             xytext=(np.cos(a + 0.08), np.sin(a + 0.08)),
             arrowprops=dict(arrowstyle="->", color=INK, lw=1.0,
                             connectionstyle="arc3,rad=0.25"))
ax1.text(1.12, -0.18, r"$f_{a1}$", color=BLUE, fontsize=9)
ax1.text(0.08, 1.2, r"$f_{a2}$", color=ORANGE, fontsize=9)
ax1.text(0.30, 0.85, r"$\Omega_1$", fontsize=9)
ax1.text(0, -1.72, "(a)", ha="center", fontsize=9)
ax1.axis("off")

# panel (b): travelling packet, quarter-wave offset quadratures
x = np.linspace(0, 12, 600)
env = np.exp(-((x - 4.5) ** 2) / 4.0)
k = 3.0
f1 = env * np.cos(k * x)
f2 = env * np.sin(k * x)
ax2.plot(x, f1, color=BLUE, lw=1.4)
ax2.plot(x, f2, color=ORANGE, lw=1.4)
ax2.plot(x, env, color=GRAY, lw=0.8, linestyle=(0, (2, 2)))
ax2.plot(x, -env, color=GRAY, lw=0.8, linestyle=(0, (2, 2)))
ax2.annotate("", xy=(9.4, 1.02), xytext=(7.0, 1.02),
             arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
ax2.text(8.2, 1.12, r"moves at $c_{\mathrm{sub}}$", ha="center", fontsize=8)
ax2.text(1.15, 0.86, r"$f_{a1}$", color=BLUE, fontsize=9)
ax2.text(2.15, -0.9, r"$f_{a2}$", color=ORANGE, fontsize=9)
ax2.set_ylim(-1.35, 1.35); ax2.set_xlim(0, 12)
ax2.set_yticks([]); ax2.set_xticks([])
ax2.set_xlabel("position along a chain of cells", fontsize=8)
ax2.text(11.7, -1.2, "(b)", ha="right", fontsize=9)
for s in ("top", "right", "left"):
    ax2.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig("fig_maxwell.pdf")
print("wrote fig_maxwell.pdf")
