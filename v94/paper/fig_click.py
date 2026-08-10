#!/usr/bin/env python3
"""fig_click: (a) the beat phase winding + click events; (b) continuous register vs quantized matter."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.4, 3.4), sharex=True)

t = np.linspace(0, 30, 3000)
beat = 0.9                                    # relative pitch w1e - w2e
phase = (beat * t) % (2 * np.pi)
clicks = t[np.where(np.diff(phase) < 0)]      # wrap instants

# (a) winding accumulator
ax1.plot(t, phase, color=INK, lw=1.2)
ax1.axhline(2 * np.pi, color=GRAY, lw=0.7, linestyle=(0, (2, 2)))
ax1.text(0.4, 2 * np.pi + 0.25, r"$2\pi$ — one full winding", fontsize=7.5)
for tc in clicks:
    ax1.axvline(tc, color=ORANGE, lw=0.8, alpha=0.85)
ax1.text(clicks[0] + 0.3, 1.15, "click", color=ORANGE, fontsize=8, rotation=90)
ax1.set_ylabel(r"beat phase $\varphi_D$")
ax1.set_yticks([0, 2 * np.pi]); ax1.set_yticklabels(["0", r"$2\pi$"])
ax1.set_ylim(0, 7.4)

# (b) register (continuous) and matter (staircase)
rng = np.random.default_rng(3)
chi2 = 0.25 + 0.22 * np.sin(0.9 * t) ** 2 + 0.12 * np.sin(2.7 * t + 1.0)
chi2 = np.clip(chi2 + 0.35 * (1 - np.exp(-t / 8)), 0, None)
eps = 0.37
Em = np.zeros_like(t)
level = 0.0
for i, ti in enumerate(t):
    fired = clicks[clicks <= ti]
    n = 0
    for tc in fired:
        j = np.searchsorted(t, tc)
        if chi2[j] > eps:
            n += 1
    Em[i] = n * eps
ax2.plot(t, chi2, color=GRAY, lw=1.1)
ax2.plot(t, Em, color=ORANGE, lw=1.6)
ax2.axhline(eps, color=BLUE, lw=0.7, linestyle=(0, (2, 2)))
ax2.text(0.4, eps + 0.05, r"one atom  $\varepsilon = A_0\omega/2\pi$",
         color=BLUE, fontsize=7.5)
ax2.text(20.5, 0.62, r"register $|\chi|^2$ (continuous,", color=GRAY, fontsize=7.5)
ax2.text(20.5, 0.50, "reversible)", color=GRAY, fontsize=7.5)
ax2.text(13.2, 1.75, r"matter $E_m$ (whole atoms only)", color=ORANGE, fontsize=7.5)
ax2.set_xlabel("time"); ax2.set_ylabel("energy")
ax2.set_xlim(0, 30); ax2.set_ylim(0, 2.4)
for ax in (ax1, ax2):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("fig_click.pdf")
print("wrote fig_click.pdf")
