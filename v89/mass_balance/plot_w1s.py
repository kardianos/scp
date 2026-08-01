#!/usr/bin/env python3
"""M-B5 figure: free-object drain trails vs the M-B4 pinned-pair curve.

Each candidate walks left as it sheds; the trail (x(t), drain_pv(t))
overlays the measured balance landscape. Output:
runs/w1s_drain_vs_curve.png
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import score_w1s as S   # parsers + curve (module-level scoring prints run)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "w1s_drain_vs_curve.png")

fig, ax = plt.subplots(figsize=(8.5, 5.5))

xs = [x / 1000 for x in range(80, 840, 4)]
ax.plot(xs, [S.B_pv(x) for x in xs], "k-", lw=1.2,
        label="M-B4 pinned-pair curve B_pv(x)")
ax.plot([p[0] for p in S.MB4], [p[1] / 2 for p in S.MB4], "ks", ms=4)
ax.plot([0.4156], [4.009e-3 / 2], "k^", ms=7, label="fifth tongue spike")

COL = {"free1": "tab:blue", "ring12_m5": "tab:green",
       "ring12_m5_ctrl": "tab:olive", "ring12_m6": "tab:red",
       "ring8_m3": "tab:purple"}
for name, nv, _ in S.RUNS:
    log = os.path.join(HERE, "runs", f"w1s_{name}.log")
    if not os.path.exists(log):
        continue
    wd = S.windowed_drain(S.lump_series(log, nv), nv)
    pts = [(x, d) for t, x, d in wd if t > 200 and d > 0]
    if not pts:
        continue
    ax.plot([p[0] for p in pts], [p[1] for p in pts], ".",
            ms=2.5, alpha=0.45, color=COL[name])
    ax.plot([], [], "o", ms=5, color=COL[name], label=name)

for xv, lab in [(S.X_SKIRT, "skirt"), (S.X_FIFTH, "3:2"), (0.8333, "2:1")]:
    ax.axvline(xv, color="gray", ls=":", lw=0.8)
    ax.text(xv, 2.5e-3, lab, ha="center", fontsize=8, color="gray")

ax.set_yscale("log")
ax.set_xlim(0.05, 0.85)
ax.set_ylim(2e-6, 4e-3)
ax.set_xlabel("address x = Em/(cap) per voice")
ax.set_ylabel("per-voice drain  [Em/t.u.]")
ax.set_title("M-B5: free desert riders on V2s vs the pinned balance landscape")
ax.legend(fontsize=8, loc="upper center", ncol=2)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUT, dpi=140)
print(f"wrote {OUT}")
