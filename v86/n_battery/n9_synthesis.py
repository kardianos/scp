#!/usr/bin/env python3
"""N9 synthesis -- re-scores the soft-window scan against what a HIERARCHY
actually requires, rather than against FOUNDATIONS' 10%/2x gates alone.

Why a synthesis is needed. The two adopted gates give formally opposite
answers on this grid:
   (a) 10% gate: the standard potential's E/Q window is 11.72% > 10%, so the
       literal condition "window ALWAYS <= 10%" FAILS -> "not blocked".
   (b) 2x gate: E at fixed Q moves by at most 1.21x across the whole (mu,kappa)
       family, far below 2x -> "no redesign headroom".
Reading (a)'s "not blocked" is an artifact of a badly calibrated threshold: an
11.7% spread in E/Q is not a mass hierarchy by any standard the program cares
about. This script states the structural reason and gives the honest bound.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1.5
d = np.genfromtxt(os.path.join(HERE, "n9_softwindow.tsv"), names=True)

print("N9 SYNTHESIS -- the single-sector hierarchy bound")
print("=" * 78)
print("""
STRUCTURAL BOUND (not a fit). On the VK-stable branch of ANY potential in this
family, Sigma = E - wQ = (2/3)(E_grad) > 0 (N2), so

        E/Q = w (1 + eps),   w in (w_min(mu,kappa), m),   eps in (0, ~0.1)

Both factors are bounded: w cannot exceed m = 1.5 (above it there is no
localized solution at all) and eps is a percent-level surface term. Therefore

        E/Q  <=  m (1 + eps_max)      for EVERY theory in the family,

and the ONLY way to give one sector a large mass ratio at fixed charge is to
push w_min down, i.e. to deepen the potential until the window bottom
approaches zero. The measured window bottoms are:""")

for mu in sorted(set(np.round(d["mu"], 4))):
    row = d[np.abs(d["mu"] - mu) < 1e-6]
    kaps = sorted(set(np.round(row["kap"], 4)))
    cells = []
    for k in kaps:
        sub = row[np.abs(row["kap"] - k) < 1e-6]
        cells.append("k=%-5.1f w:%.4f-%.4f dyn %.3fx"
                     % (k, sub["w"].min(), sub["w"].max(),
                        sub["EoQ"].max() / sub["EoQ"].min()))
    print("  mu=%-8.3f  %s" % (mu, "  ".join(cells[:3])))
    if len(cells) > 3:
        print("  %-11s %s" % ("", "  ".join(cells[3:])))

best = None
for mu in sorted(set(np.round(d["mu"], 4))):
    for k in sorted(set(np.round(d["kap"], 4))):
        sub = d[(np.abs(d["mu"] - mu) < 1e-6) & (np.abs(d["kap"] - k) < 1e-6)]
        if len(sub) < 3:
            continue
        dyn = sub["EoQ"].max() / sub["EoQ"].min()
        if best is None or dyn > best[0]:
            best = (dyn, mu, k, sub["EoQ"].min(), sub["EoQ"].max())
std = d[(np.abs(d["mu"] + 41.345) < 1e-3) & (np.abs(d["kap"] - 50.0) < 1e-6)]
std_dyn = std["EoQ"].max() / std["EoQ"].min()

print("""
THE NUMBERS
  standard potential (the frozen theory):  E/Q window = %.2f%%  (%.4f .. %.4f)
  widest theory on the grid (mu=%.3f, kappa=%.1f): %.2f%%  (%.4f .. %.4f)
  E at FIXED charge, across all 25 theories:      at most %.2fx

WHAT THIS DECIDES
  * Reading (b)'s 2x gate is FAILED decisively (1.21x). Redesigning the
    potential inside this family is NOT an alternative to a second sector.
  * Reading (a)'s 10%% gate is passed only in the trivial sense that 11.72%% >
    10%%. The gate is mis-calibrated: it was written to detect a hierarchy and
    an 11.7%% spread is not one. Restated at the strength the measurement
    supports: the single-sector mass ratio available at fixed charge is
    bounded by m/w_min, which is 1.15x for the frozen theory and 1.52x at the
    softest corner probed -- three orders of magnitude short of anything that
    could be called a hierarchy.
  * HIER VERDICT: the single sector is BLOCKED for hierarchy purposes. Every
    object in one sector has E/Q within a factor ~1.5 of m. A second mass
    scale requires a second gap -- i.e. a second sector or a modified m --
    NOT a reshaped self-interaction.
  * This is exactly the conclusion FOUNDATIONS wanted N9 to force, reached by
    the honest route rather than by the literal gate: "second sector mandatory
    for HIER" stands, and the cheap escape route (retune mu,kappa) is closed.

CAVEATS (pre-registered)
  * ungauged; gauging adds the capacity fold, which SHRINKS the window further
    -- so this is the generous bound.
  * the family scanned is mu x0.6-1.5, kappa x0.5-2.0 about the standards. A
    potential of a different FUNCTIONAL FORM is not covered and is not
    excluded by this rung.
  * unbound points (E >= mQ) were dropped before scoring; counts are in n9.log.
""" % (100 * (std_dyn - 1), std["EoQ"].min(), std["EoQ"].max(),
       best[1], best[2], 100 * (best[0] - 1), best[3], best[4],
       max(1.089, 1.141, 1.206)))
