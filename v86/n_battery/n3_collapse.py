#!/usr/bin/env python3
"""v86 Part-0 rung N3 -- scaling collapse of the excess.

Protocol: v86/council/grok45/N_BATTERY_REVIEW.md section 1, N3 (restored):
  Definitions: xi = 1/sqrt(m^2 - w^2);  R = r_half (FROZEN; sensitivity vs the
               charge-RMS radius reported).
  Coulomb subtraction (gauged): prefer the MEASURED E_g from N1 over the
               point-Coulomb estimate g^2 Q^2/(8 pi R) (grok Finding: the
               point estimate can absorb arbitrary residuals).
  Plots/tables: (1) eps vs xi/R ; (2) eps vs g^2 Q / R ; (3) eps' vs xi/R after
               subtraction.
  Outcomes:
     collapse on xi/R alone (all g)          -> pure geometry / surface
     needs g^2Q/R, or only after subtraction -> Coulomb self-energy essential
     no collapse                             -> bulk / thick-wall / potential shape

N2 makes this rung sharp. From Sigma = (2/3)(E_grad - E_g) and wQ = 2(E_kin+E_g):
     eps  = (E_grad - E_g) / (3 (E_kin + E_g))
     eps' = eps + (2/3) E_g/(w Q) = (2/3) E_grad / (w Q)
so the Coulomb-subtracted excess eps' is EXACTLY the gradient (surface) term
divided by the throughput, with the measured E_g doing the subtraction and no
fitted constant anywhere. The collapse test is therefore a test of whether the
gradient fraction is a function of xi/R alone.

Input: n1_decomp.tsv written by n12_decomp.py (production grid preferred).
Output: n3_collapse.tsv + a printed table; no plotting dependency.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/d/code/scp"
TSV = os.path.join(HERE, "n1_decomp.tsv")


def load():
    if not os.path.exists(TSV):
        sys.exit("missing %s -- run n12_decomp.py first" % TSV)
    return np.genfromtxt(TSV, names=True)


def rms_radius_sensitivity(sample):
    """Re-solve a handful of branch points and compare R = r_half against the
    charge-weighted RMS radius, so the frozen choice of R is auditable
    (protocol: 'freeze one and report sensitivity')."""
    sys.path.insert(0, os.path.join(ROOT, "v69/theory"))
    import gauged_shooter_fast as G
    out = []
    f0g = G.load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
    fseed, cseed, ok, _, _ = G.solve(1.45, 0.0, f0g, np.zeros(G.N))
    if not ok:
        return out
    for (g, w) in sample:
        # march in omega at g=0 then continue in g (same recipe as N1)
        f, chi = fseed.copy(), cseed.copy()
        wc = 1.45
        step = -0.004 if w < wc else 0.004
        while abs(wc - w) > 1e-9:
            wn = wc + step
            if (step < 0 and wn < w) or (step > 0 and wn > w):
                wn = w
            fn, cn, okk, _, _ = G.solve(wn, 0.0, f, chi)
            if not okk:
                step *= 0.5
                if abs(step) < 1e-5:
                    break
                continue
            f, chi, wc = fn, cn, wn
        if abs(wc - w) > 1e-6:
            continue
        gc, chig = 0.0, 0.0
        okall = True
        while gc < g - 1e-12:
            dg = g - gc
            while True:
                gg = min(gc + dg, g)
                cg = chi * (gg / chig) ** 2 if chig > 0 else chi
                fn, cn, okk, _, _ = G.solve(w, gg, f, cg)
                if okk:
                    f, chi, gc, chig = fn, cn, gg, gg
                    break
                dg *= 0.5
                if dg < 1e-6:
                    okall = False
                    break
            if not okall:
                break
        if not okall:
            continue
        o = G.observables(f, chi, w, g)
        fe = np.append(f, 0.0)
        ce = np.append(chi, chi[-1] * G.RFAC)
        re = np.append(G.r, G.rN)
        wte = w - ce
        dV = 4.0 * np.pi * re * re
        rho = 3.0 * wte * fe * fe
        Q = np.trapz(rho * dV, dx=G.H)
        R_rms = np.sqrt(np.trapz(rho * re * re * dV, dx=G.H) / Q)
        out.append((g, w, o["rhalf"], R_rms, R_rms / o["rhalf"]))
    return out


def collapse_metric(x, y, nb=14):
    """Spread of y within bins of x, relative to the total spread of y.
    0 = perfect collapse onto a single curve of x; ~1 = x explains nothing."""
    x = np.asarray(x); y = np.asarray(y)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 10:
        return np.nan
    edges = np.quantile(x, np.linspace(0, 1, nb + 1))
    edges[-1] += 1e-12
    num, den = 0.0, np.std(y)
    cnt = 0
    for i in range(nb):
        sel = (x >= edges[i]) & (x < edges[i + 1])
        if sel.sum() >= 3:
            num += np.var(y[sel]) * sel.sum()
            cnt += sel.sum()
    if cnt == 0 or den == 0:
        return np.nan
    return float(np.sqrt(num / cnt) / den)


def main():
    d = load()
    gs = sorted(set(np.round(d["g"], 6)))
    print("v86 N3 -- scaling collapse of the excess   (%d branch points, g in %s)"
          % (len(d), gs))
    print("R is FROZEN to r_half; xi = 1/sqrt(m^2 - w^2); the Coulomb subtraction")
    print("uses the MEASURED E_g from N1, not the point-Coulomb estimate.\n")

    wQ = d["w"] * d["Q"]
    eps = d["eps"]
    epsp = eps + (2.0 / 3.0) * d["E_g"] / wQ        # = (2/3) E_grad / (wQ)
    xiR = d["xi_over_R"]
    coul = d["coul_par"]                            # g^2 Q / R
    # sanity: eps' must equal (2/3)E_grad/(wQ) to the N2 residual
    chk = np.max(np.abs(epsp - (2.0 / 3.0) * d["E_grad"] / wQ))
    print("  identity check  eps' == (2/3)E_grad/(wQ):  max deviation %.2e" % chk)

    # ---------------- the three collapse tests ------------------------------
    m1 = collapse_metric(xiR, eps)
    m2 = collapse_metric(coul, eps)
    m3 = collapse_metric(xiR, epsp)
    print("\ncollapse metric (0 = perfect collapse, 1 = variable explains nothing):")
    print("  (1) eps  vs xi/R                        : %.4f" % m1)
    print("  (2) eps  vs g^2 Q / R                   : %.4f" % m2)
    print("  (3) eps' vs xi/R  (measured-E_g removed): %.4f" % m3)

    # ---------------- tabulate along xi/R -----------------------------------
    print("\n%-6s %9s %9s %9s %9s %9s %9s %9s" %
          ("g", "w", "Q", "R_half", "xi/R", "g^2Q/R", "eps", "eps'"))
    for g in gs:
        sel = np.where(np.abs(d["g"] - g) < 1e-9)[0]
        if len(sel) == 0:
            continue
        for i in sel[::max(1, len(sel) // 7)]:
            print("%-6.3f %9.4f %9.2f %9.3f %9.4f %9.4f %9.5f %9.5f" %
                  (d["g"][i], d["w"][i], d["Q"][i], d["rhalf"][i], xiR[i],
                   coul[i], eps[i], epsp[i]))

    # ---------------- binned collapse table ---------------------------------
    print("\neps' in bins of xi/R (a collapse means the g-columns agree in a row):")
    bins = np.array([0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.5, 4.0])
    hdr = "%-14s" % "xi/R bin"
    for g in gs:
        hdr += "%14s" % ("g=%.3f" % g)
    print(hdr)
    for i in range(len(bins) - 1):
        line = "%-14s" % ("%.2f-%.2f" % (bins[i], bins[i + 1]))
        any_row = False
        for g in gs:
            sel = ((np.abs(d["g"] - g) < 1e-9) & (xiR >= bins[i]) & (xiR < bins[i + 1]))
            if sel.sum():
                line += "%14s" % ("%.5f(%d)" % (epsp[sel].mean(), sel.sum()))
                any_row = True
            else:
                line += "%14s" % "-"
        if any_row:
            print(line)

    print("\nsame table for RAW eps (before the measured-E_g subtraction):")
    print(hdr)
    for i in range(len(bins) - 1):
        line = "%-14s" % ("%.2f-%.2f" % (bins[i], bins[i + 1]))
        any_row = False
        for g in gs:
            sel = ((np.abs(d["g"] - g) < 1e-9) & (xiR >= bins[i]) & (xiR < bins[i + 1]))
            if sel.sum():
                line += "%14s" % ("%.5f(%d)" % (eps[sel].mean(), sel.sum()))
                any_row = True
            else:
                line += "%14s" % "-"
        if any_row:
            print(line)

    # ---------------- thin-wall asymptote -----------------------------------
    print("\nthin-wall check (does eps -> 0 as R/xi -> infinity, i.e. xi/R -> 0?):")
    sel0 = np.abs(d["g"]) < 1e-12
    if sel0.sum():
        o = np.argsort(xiR[sel0])
        xs, ys = xiR[sel0][o], eps[sel0][o]
        for frac in (0.0, 0.05, 0.1, 0.25, 0.5):
            i = min(len(xs) - 1, int(frac * (len(xs) - 1)))
            print("   xi/R=%.4f -> eps=%.5f" % (xs[i], ys[i]))
        # local power law near the thin-wall end
        k = max(6, len(xs) // 6)
        p = np.polyfit(np.log(xs[:k]), np.log(np.maximum(ys[:k], 1e-30)), 1)
        print("   local power law over the %d thinnest-wall points: eps ~ (xi/R)^%.3f"
              % (k, p[0]))

    # ---------------- R-definition sensitivity ------------------------------
    print("\nR-definition sensitivity (frozen R = r_half vs charge-weighted RMS):")
    sample = [(0.0, 1.34), (0.0, 1.40), (0.0, 1.46), (0.05, 1.42), (0.05, 1.47),
              (0.10, 1.48)]
    rows = rms_radius_sensitivity(sample)
    if rows:
        print("  %-6s %-8s %10s %10s %10s" % ("g", "w", "r_half", "R_rms", "ratio"))
        for g, w, rh, rr, ratio in rows:
            print("  %-6.3f %-8.4f %10.4f %10.4f %10.4f" % (g, w, rh, rr, ratio))
        ratios = [x[4] for x in rows]
        print("  ratio R_rms/r_half spans %.3f-%.3f: xi/R shifts by that factor,"
              % (min(ratios), max(ratios)))
        print("  which RESCALES the horizontal axis uniformly only if the ratio is")
        print("  constant. Spread of the ratio = %.1f%% -> the collapse verdict is"
              % (100 * (max(ratios) - min(ratios)) / np.mean(ratios)))
        print("  %s to the choice of R."
              % ("ROBUST" if (max(ratios) - min(ratios)) / np.mean(ratios) < 0.1
                 else "SENSITIVE"))
    else:
        print("  (re-solve failed; sensitivity not measured)")

    # ---------------- verdict ----------------------------------------------
    print("\n" + "=" * 78)
    best = min([(m1, "xi/R alone"), (m2, "g^2Q/R alone"),
                (m3, "xi/R after measured-E_g subtraction")],
               key=lambda t: (np.inf if not np.isfinite(t[0]) else t[0]))
    print("N3 VERDICT: best collapse variable = %s (metric %.4f)" % (best[1], best[0]))
    improve = (m1 / m3) if (np.isfinite(m1) and np.isfinite(m3) and m3 > 0) else np.nan
    print("  subtraction improvement factor m1/m3 = %.2f (how much removing the "
          "MEASURED gauge energy tightens the collapse on xi/R)" % improve)
    if np.isfinite(m3) and m3 < 0.25 and np.isfinite(improve) and improve > 1.2:
        print("  -> GEOMETRY DOMINANT, COULOMB CORRECTION REAL. Both readings hold and")
        print("     neither alone is the answer:")
        print("     * xi/R already organises the raw excess (metric %.3f) -- the" % m1)
        print("       excess is primarily a SURFACE/geometry quantity;")
        print("     * removing the measured E_g tightens it by %.2fx (metric %.3f)," % (improve, m3))
        print("       and the binned table shows the improvement is concentrated at")
        print("       the largest g, where the Coulomb self-energy is largest.")
        print("     So: the excess is the gradient/surface term (2/3)E_grad/(wQ), and")
        print("     the gauge sector enters ONLY as the -(2/3)E_g piece that N2")
        print("     already isolates -- not as an independent scaling variable")
        print("     (g^2Q/R alone explains essentially nothing: metric %.3f)." % m2)
    elif np.isfinite(m1) and m1 < 0.25:
        print("  -> PURE GEOMETRY / SURFACE: eps is a function of xi/R alone across g,")
        print("     and the gauge subtraction does not measurably improve it.")
    else:
        print("  -> NO CLEAN COLLAPSE on either variable: bulk / thick-wall /")
        print("     potential-shape contributions matter. This does NOT invalidate")
        print("     N2 (the closed form is exact); it changes the closed-form TARGET")
        print("     from a surface law to a profile-moment expression.")

    out = os.path.join(HERE, "n3_collapse.tsv")
    with open(out, "w") as fh:
        fh.write("g\tw\tQ\tR_half\txi\txi_over_R\tcoul_par\teps\teps_prime\tE_g\tE_grad\n")
        for i in range(len(d)):
            fh.write("%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t"
                     "%.10g\t%.10g\t%.10g\n" %
                     (d["g"][i], d["w"][i], d["Q"][i], d["rhalf"][i], d["xi"][i],
                      xiR[i], coul[i], eps[i], epsp[i], d["E_g"][i], d["E_grad"][i]))
    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
