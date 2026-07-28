#!/usr/bin/env python3
"""v86 HC-3-volume -- the partition-space scan of the charge-response signature.

WHY THIS RUNG EXISTS
  HC-3 (n4_hc3_flavored.py) computed n(D) for the signature of
  D_ab = dQ_a/dw_b along ONE detuning ray: w0 down, w1 = w2 up. It found
  n(D) = 1 at every point and flagged the limitation itself -- "this is ONE
  detuning ray at ONE total charge ... a full (Q0,Q1,Q2) volume scan is HC-3's
  follow-on". NEXT_PROGRAM.md promotes that follow-on to the gate on HC-6:
  either the volume contains a signature change (HC-6 has targets) or it does
  not (HC-6 is unrunnable, which is itself the result).

THE GEOMETRY, AND WHAT THE OLD SCAN ACTUALLY COVERED
  Fixed total charge means moving in the TRACELESS plane of (w0,w1,w2). Every
  traceless detuning can be written
      delta_a(theta, rho) = rho * cos(theta - 120 deg * a),   a = 0,1,2
  which satisfies sum_a delta_a = 0 identically. Two facts follow:

    * theta = 180 deg gives delta = rho * (-1, +0.5, +0.5) -- EXACTLY the old
      ray. So HC-3 sampled a single value of theta.
    * The 3 component permutations act as rotations by 120 deg plus a
      reflection, so the fundamental domain is theta in [0, 60] deg. Everything
      outside is a relabelling of something inside.

  theta = 180 deg is equivalent to theta = 60 deg, i.e. the old ray sits on ONE
  EDGE of the fundamental domain. The opposite edge, theta = 0 deg -- one
  component detuned UP and two DOWN -- has never been sampled. That is the
  gap this rung closes, and it is not a small one: the two edges are physically
  different (one flavour enriched vs two flavours enriched).

WHAT IS AND IS NOT FIXED
  The old rung Newton-corrected a uniform w shift to pin Q_tot, because N4 was
  comparing Sigma at fixed Q. HC-3 asks a different question -- the SIGNATURE of
  D at a point of the branch -- which is a property of the point, not of a
  constraint. So no Q correction is applied here; Q_tot is reported for every
  point, and sweeping the base frequency wbar sweeps the total charge. The
  (Q0,Q1,Q2) coverage is therefore an output, printed as a coverage table.

Reuses the solver, observables and continuation from n4_hc3_flavored.

Usage: hc3_volume.py [--quick] [--out prefix]
"""
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import n4_hc3_flavored as H3           # noqa: E402  (solver machinery)

DEG = np.pi / 180.0


def delta(theta, rho):
    """Traceless detuning: sum_a delta_a == 0 to machine precision."""
    return rho * np.cos(theta - 120.0 * DEG * np.arange(3))


def charge_response(F, w, hw=1e-3):
    """D_ab = dQ_a/dw_b by central differences, warm-started from (F,w).

    Returns (eigenvalues, n_neg, n_pos, n_zero, asym) or None if any of the six
    perturbed solves fails."""
    D = np.zeros((3, 3))
    for b in range(3):
        Qpm = []
        for sgn in (+1, -1):
            wp = np.array(w, float)
            wp[b] += sgn * hw
            Fp, wpp = H3.walk_to(F, w, wp, nstep=1)
            if Fp is None:
                return None
            Qpm.append(H3.observables(Fp, wpp)["Q"])
        D[:, b] = (Qpm[0] - Qpm[1]) / (2 * hw)
    Dsym = 0.5 * (D + D.T)
    ev = np.linalg.eigvalsh(Dsym)
    tol = 1e-6 * max(1e-30, np.max(np.abs(ev)))
    nneg = int((ev < -tol).sum())
    npos = int((ev > tol).sum())
    asym = np.max(np.abs(D - D.T)) / max(1e-30, np.max(np.abs(D)))
    return ev, nneg, npos, 3 - nneg - npos, asym


def main():
    quick = "--quick" in sys.argv
    if "--from-tsv" in sys.argv:
        import csv
        rows = []
        with open(os.path.join(HERE, "hc3_volume.tsv")) as fh:
            for d_ in csv.DictReader(fh, delimiter="\t"):
                rows.append(dict(wbar=float(d_["wbar"]), theta=float(d_["theta"]),
                                 rho=float(d_["rho"]),
                                 w=(float(d_["w0"]), float(d_["w1"]), float(d_["w2"])),
                                 Q=(float(d_["Q0"]), float(d_["Q1"]), float(d_["Q2"])),
                                 Qtot=float(d_["Qtot"]), E=float(d_["E"]),
                                 eps=float(d_["eps"]),
                                 ev=(float(d_["ev0"]), float(d_["ev1"]), float(d_["ev2"])),
                                 nneg=int(d_["nneg"]), npos=int(d_["npos"]),
                                 nzero=int(d_["nzero"]), asym=float(d_["asym"])))
        globals()["_PRELOADED"] = rows
    pfx = "hc3_volume"
    if "--out" in sys.argv:
        pfx = sys.argv[sys.argv.index("--out") + 1]

    if quick:
        wbars = [1.42]
        thetas = [0.0, 60.0]
        rho_max, rho_step = 0.03, 0.015
    else:
        # spans the ungauged Q-ball window (1.3087, 1.5); wbar sweeps Q_tot
        wbars = [1.36, 1.39, 1.42, 1.45, 1.48]
        thetas = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        rho_max, rho_step = 0.10, 0.005

    print("v86 HC-3-volume -- charge-response signature over the partition plane")
    print("  solver: H=%g RMAX=%g N=%d   (ungauged, 3 components)"
          % (H3.H, H3.RMAX, H3.N))
    print("  wbar   : %s" % wbars)
    print("  theta  : %s deg   (fundamental domain [0,60]; old HC-3 ray = 60)"
          % thetas)
    print("  rho    : 0 .. %g step %g" % (rho_max, rho_step))
    print("  NOTE: total charge is NOT constrained here -- the signature of D is")
    print("        a property of the branch point. Q_tot is reported per point.")
    print()

    t0 = time.time()
    rows = []
    if "_PRELOADED" in globals():
        rows = globals()["_PRELOADED"]
        wbars = []
        print("  (re-deriving the verdict from hc3_volume.tsv, %d rows)" % len(rows))
    for wbar in wbars:
        try:
            F0 = H3.symmetric_seed(wbar)
        except Exception as e:
            print("wbar=%.3f: symmetric seed FAILED (%s) -- branch edge" % (wbar, e))
            continue
        w0 = np.array([wbar, wbar, wbar])
        o0 = H3.observables(F0, w0)
        print("=" * 96)
        print("wbar = %.3f   symmetric point: Q_tot = %.3f  E = %.3f  eps = %.5f"
              % (wbar, o0["Qtot"], o0["E"], o0["eps"]))
        print("%-7s %-7s %10s %10s %10s %31s %5s %5s  %s"
              % ("theta", "rho", "Q0", "Q1", "Q2", "eig(dQ_a/dw_b)",
                 "n(D)", "p(D)", "verdict"))
        for th in thetas:
            F, w = F0.copy(), w0.copy()
            rho = 0.0
            while rho <= rho_max + 1e-12:
                wt = w0 + delta(th * DEG, rho)
                Fn, wn = H3.walk_to(F, w, wt, nstep=1 if rho > 0 else 1)
                if Fn is None:
                    print("%-7.1f %-7.3f   branch end (no convergence)" % (th, rho))
                    break
                F, w = Fn, wn
                res = charge_response(F, w)
                if res is None:
                    print("%-7.1f %-7.3f   FD failed" % (th, rho))
                    break
                ev, nneg, npos, nzero, asym = res
                o = H3.observables(F, w)
                verdict = ("n(D)=1" if nneg == 1
                           else "*** n(D)=%d -- HC-6 TARGET ***" % nneg)
                if nzero:
                    verdict += "  [%d near-zero]" % nzero
                print("%-7.1f %-7.3f %10.3f %10.3f %10.3f  [%9.1f %9.1f %9.1f] "
                      "%5d %5d  %s"
                      % (th, rho, o["Q"][0], o["Q"][1], o["Q"][2],
                         ev[0], ev[1], ev[2], nneg, npos, verdict))
                rows.append(dict(wbar=wbar, theta=th, rho=rho, w=tuple(w),
                                 Q=tuple(o["Q"]), Qtot=o["Qtot"], E=o["E"],
                                 eps=o["eps"], ev=tuple(ev), nneg=nneg,
                                 npos=npos, nzero=nzero, asym=asym))
                rho += rho_step
        print()

    # --------------------------------------------------- window filter first
    # A converged Newton solution is not automatically a physical one: the
    # continuation can walk a component past the Q-ball window edge omega=1.5,
    # where the "solution" is a solver artifact (Q jumps by an order of
    # magnitude, one eigenvalue blows up to ~1e6). Those points are dropped
    # BEFORE any n(D) claim is made.
    WLO, WHI = 1.30870, 1.5
    kept = [r for r in rows if all(WLO < wv < WHI for wv in r["w"])]
    drop = [r for r in rows if r not in kept]
    print("=" * 96)
    print("WINDOW FILTER: %d converged, %d inside (%.4f, %.4f), %d dropped"
          % (len(rows), len(kept), WLO, WHI, len(drop)))
    if drop:
        dn = sorted(set(r["nneg"] for r in drop))
        print("  dropped points carried n(D) in %s -- including any n(D)=2, which"
              % dn)
        print("  is a continuation artifact past the window edge, not a target.")
    rows = kept
    print()

    # ------------------------------------------------------------- synthesis
    print("=" * 96)
    print("HC-3-VOLUME SYNTHESIS   (%d converged partition points, %.1f s)"
          % (len(rows), time.time() - t0))
    if not rows:
        print("  no points converged -- nothing to conclude")
        return
    seen = sorted(set(r["nneg"] for r in rows))
    masym = max(r["asym"] for r in rows)
    evmin = min(min(abs(e) for e in r["ev"]) for r in rows)
    Qt = np.array([r["Qtot"] for r in rows])
    Qa = np.array([r["Q"] for r in rows])
    frac = Qa / Qt[:, None]
    print("  total charge covered      : Q_tot in [%.2f, %.2f]" % (Qt.min(), Qt.max()))
    print("  partition fractions covered: Q_a/Q_tot in [%.4f, %.4f]"
          % (frac.min(), frac.max()))
    print("      (symmetric point is 1/3 = 0.3333; the scan reaches %.4f from it)"
          % max(abs(frac.min() - 1.0 / 3), abs(frac.max() - 1.0 / 3)))
    print("  D symmetry max|D-D^T|/max|D| : %.2e%s"
          % (masym, "   *** LARGE -- n(D) NOT TRUSTWORTHY ***" if masym > 1e-2 else ""))
    print("  smallest |eigenvalue| seen   : %.3e" % evmin)
    print("  n(D) values seen             : %s" % seen)
    print()
    print("  n(D) by base frequency:")
    for wb in sorted(set(r["wbar"] for r in rows)):
        sub = [r for r in rows if r["wbar"] == wb]
        print("    wbar=%.3f : %4d pts, n(D) values %s"
              % (wb, len(sub), sorted(set(r["nneg"] for r in sub))))
    print()
    if seen == [1]:
        print("  RESULT: n(D) = 1 at EVERY converged point of the accessible")
        print("  partition volume, including the theta = 0 edge (one flavour up,")
        print("  two down) that the original ray never sampled, and across a")
        print("  Q_tot range of %.2f .. %.2f." % (Qt.min(), Qt.max()))
        print()
        print("  CONSEQUENCE FOR HC-6: no GSS-unstable partition exists to seed.")
        print("  HC-6 as designed is UNRUNNABLE on the flavoured branch of this")
        print("  potential -- not because the test is bad but because the target")
        print("  set is empty. That converts HC-3 from a local check into a")
        print("  global statement about this sector, which is what the redesign")
        print("  was for.")
        print()
        print("  This is CONSISTENT with, and now numerically supports, HC-1's")
        print("  analytic result L_x^flav = L_0 - A with A = 2 P_0 < 0: the")
        print("  flavour channels cannot supply a negative direction in this")
        print("  potential, so no partition can carry n(D) > 1.")
    else:
        print("  RESULT: n(D) CHANGES across the partition volume: %s" % seen)
        print("  The index boundary lies inside the accessible region. Points")
        print("  with n(D) != 1 are HC-6's decay targets -- listed below.")
        for r in rows:
            if r["nneg"] != 1:
                print("    wbar=%.3f theta=%.1f rho=%.3f  w=(%.4f,%.4f,%.4f)  "
                      "Q=(%.2f,%.2f,%.2f)  n(D)=%d"
                      % (r["wbar"], r["theta"], r["rho"], r["w"][0], r["w"][1],
                         r["w"][2], r["Q"][0], r["Q"][1], r["Q"][2], r["nneg"]))
    print()
    print("  CLAIM LIMIT (carried from HC-3): this is an n(D) map. The GSS")
    print("  criterion is n(H_omega) = n(D), and n(H_omega) is HC-1's")
    print("  deliverable, assumed = 1 for this one-hump family. UNGAUGED.")

    path = os.path.join(HERE, pfx + ".tsv")
    with open(path, "w") as fh:
        fh.write("wbar\ttheta\trho\tw0\tw1\tw2\tQ0\tQ1\tQ2\tQtot\tE\teps"
                 "\tev0\tev1\tev2\tnneg\tnpos\tnzero\tasym\n")
        for r in rows:
            fh.write("%.6g\t%.6g\t%.6g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g"
                     "\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%d\t%d\t%d\t%.3e\n"
                     % (r["wbar"], r["theta"], r["rho"], r["w"][0], r["w"][1],
                        r["w"][2], r["Q"][0], r["Q"][1], r["Q"][2], r["Qtot"],
                        r["E"], r["eps"], r["ev"][0], r["ev"][1], r["ev"][2],
                        r["nneg"], r["npos"], r["nzero"], r["asym"]))
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
