#!/usr/bin/env python3
"""v86 Part-0 rungs N5 (throughput vs charge) and N6 (hbar_eff triple).

Protocols: v86/council/grok45/N_BATTERY_REVIEW.md section 1, N5 and N6.

N5 -- observables to report (all three):
      T      = 2 E_kin            (cycle-averaged if non-stationary)
      S_w    = sum_a w_a Q_a      (w Q monochromatic)
      2 pi Q (and hbar_eff candidates T/w, E/w)
    Process action frozen as A = closed-integral of 2 E_kin over one clock
    period; for a stationary state A/(2 pi) is compared to Q.
    Discriminators:
      T = S_w within diag noise on flavored     -> throughput generalizes (R4 soft)
      systematic T != S_w on flavored/breathing -> action splits from charge (R4 hard)
      cloud sector E_field ~ w_cl |Q_cl|        -> supports the linear-mode pillar

    CANNOT-VERIFY resolved (grok N5 ambiguity flag 1): the diag DOES expose
    E_phi_kin / E_theta_kin, so T is directly available at the diag cadence
    (dt = 0.25 on the v85 runs) -- N5 is free, no SFA post-processing needed.

    AMBIGUITY carried (grok N5 flag 2 + a new one found here): the diag's
    `omega_core` is a POINT SAMPLE of (u0 vdot0 - v0 udot0)/(u0^2+v0^2) at the
    argmax-s voxel. In the gauged runs, and in temporal gauge, that is the
    LOCAL clock at the ball centre, i.e. the gauge-shifted w_eff(0), not the
    bare branch omega. Every "w Q" formed from it is therefore w_eff(0)*Q. This
    script reports the discrepancy explicitly rather than hiding it, and uses
    the N1/N2 identity  w Q = 2 (E_kin + E_g)  to BACK OUT the effective bare
    omega from measured energies:   w_bar = 2 (E_kin + E_em) / Q.
    Comparing w_bar to omega_core is itself a measurement of the gauge shift.

N6 -- triple (same object, same time window):
      hbar_E  = E / w      with E = integral of T_00 ONLY (never Q w)
      hbar_pk = p / k      (needs a boost series)
      hbar_Q  = Q
    Report hbar_E/Q - 1, hbar_pk/Q - 1, hbar_E/hbar_pk - 1 and overlay vs eps.
    Expected: residuals in the same family as eps (1-4%); kill the identity
    language, keep measured ratios.
    hbar_pk requires boosted runs; if no boost series is present in the
    archive this script reports N6 as PARTIAL (hbar_E and hbar_Q only) rather
    than inventing a momentum.
"""
import os
import sys
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DIRS = ["/space/scp/v85/topo1/out", "/space/scp/v85/topo1/gpu"]


def load(path):
    try:
        return np.genfromtxt(path, names=True)
    except Exception:
        return None


def window_stats(d, frac_lo=0.5, frac_hi=1.0):
    """Average over the last (frac_hi-frac_lo) of the run -- long enough to
    average over the breathing period (~500-600 t.u. on the X10 runs) when the
    run is 2000 t.u. long (grok N5 flag 2: the cycle average must be defined)."""
    t = d["t"]
    lo, hi = t[0] + frac_lo * (t[-1] - t[0]), t[0] + frac_hi * (t[-1] - t[0])
    sel = (t >= lo) & (t <= hi)
    if sel.sum() < 10:
        sel = np.ones_like(t, dtype=bool)
    return sel, (t[sel][-1] - t[sel][0])


def main():
    files = []
    for dd in DIRS:
        files += sorted(glob.glob(os.path.join(dd, "*_diag.tsv")))
    print("v86 N5 -- throughput vs charge   (%d diag files)" % len(files))
    print("window: last 50%% of each run; T = 2*(E_phi_kin + E_theta_kin)")
    print("w_bar  = 2 (E_kin + E_em) / Q   [the N1 Gauss identity, inverted]")
    print("w_core = diag omega_core        [POINT sample of the local clock]\n")

    print("%-22s %7s %9s %10s %10s %9s %9s %9s %9s %8s" %
          ("run", "T_win", "Q", "T", "E", "w_core", "w_bar", "T/(w_bar Q)",
           "E/(w_bar Q)", "E_em/E"))
    rows = []
    for p in files:
        d = load(p)
        if d is None or "Q_phi" not in d.dtype.names or len(d["t"]) < 20:
            continue
        sel, twin = window_stats(d)
        Q = float(np.mean(d["Q_phi"][sel]))
        if abs(Q) < 1.0:
            continue
        Ekin = float(np.mean(d["E_phi_kin"][sel] + d["E_theta_kin"][sel]))
        E = float(np.mean(d["E_total"][sel]))
        Eem = float(np.mean(d["E_em"][sel])) if "E_em" in d.dtype.names else 0.0
        wcore = float(np.mean(d["omega_core"][sel])) if "omega_core" in d.dtype.names else np.nan
        T = 2.0 * Ekin
        wbar = 2.0 * (Ekin + Eem) / Q
        name = os.path.basename(p).replace("_diag.tsv", "")
        # w_bar is only a CLOCK for a single object: Q is then that object's
        # charge. For pair / ball+cloud runs Q is a NET charge and w_bar is a
        # charge-weighted mixture, which shows up as w_bar outside the Q-ball
        # window (1.3087, 1.5). Such rows are flagged, not silently averaged in.
        single = (1.3087 < wbar < 1.5)
        rows.append(dict(name=name + ("" if single else " *"), single=single,
                         Q=Q, T=T, E=E, Eem=Eem, wcore=wcore,
                         wbar=wbar, twin=twin,
                         T_over_wQ=T / (wbar * Q), E_over_wQ=E / (wbar * Q),
                         eps_bar=E / (wbar * Q) - 1.0,
                         gauge_shift=(wcore / wbar - 1.0) if np.isfinite(wcore) else np.nan,
                         Ekin=Ekin))
        r = rows[-1]
        print("%-22s %7.0f %9.2f %10.3f %10.3f %9.5f %9.5f %11.6f %11.6f %8.4f" %
              (name, twin, Q, T, E, wcore, wbar, r["T_over_wQ"], r["E_over_wQ"],
               Eem / E if E else 0))

    if not rows:
        sys.exit("no usable diag files found")

    # ---------------- N5 discriminator --------------------------------------
    print("\n" + "=" * 78)
    print("N5 -- is throughput the same thing as charge?")
    print("=" * 78)
    print("The N1 identity w Q = 2(E_kin + E_g) is what defines w_bar here, so")
    print("T/(w_bar Q) = 1 - 2 E_em/(w_bar Q) is a CONSISTENCY check, not a test:")
    print("its deviation from 1 measures the gauge share of the throughput.")
    print("The physical N5 content is the gauge share and the EXCESS eps_bar:\n")
    print("%-22s %12s %12s %12s" %
          ("run", "2E_em/(w_bar Q)", "eps_bar=E/(w_barQ)-1", "w_core/w_bar-1"))
    for r in rows:
        print("%-22s %12.6f %12.6f %12.6f" %
              (r["name"], 2 * r["Eem"] / (r["wbar"] * r["Q"]), r["eps_bar"],
               r["gauge_shift"]))

    nflag = sum(1 for r in rows if not r["single"])
    if nflag:
        print("\n  * %d run(s) have w_bar outside the Q-ball window (1.3087,1.5):"
              % nflag)
        print("    these are pair / ball+cloud runs where Q is a NET charge, so")
        print("    w_bar is a charge-weighted mixture and NOT a clock. They are")
        print("    excluded from the aggregate statistics below.")
    srows = [r for r in rows if r["single"]]
    eb = np.array([r["eps_bar"] for r in srows])
    gs = np.array([r["gauge_shift"] for r in srows])
    print("\n  eps_bar across the %d single-object runs: %.4f .. %.4f (median %.4f)"
          % (len(srows), eb.min(), eb.max(), float(np.median(eb))))
    print("  This is the in-kernel counterpart of the shooter's eps = 0.9-4.3%%.")
    good = np.isfinite(gs)
    if good.sum():
        print("  gauge shift w_core/w_bar - 1: %.4f .. %.4f -- the ball-centre clock"
              % (gs[good].min(), gs[good].max()))
        print("  runs FASTER/SLOWER than the throughput-derived bare omega by this")
        print("  amount. A nonzero, g-dependent shift is the direct signature that")
        print("  CHARGE, ACTION and the LOCAL CLOCK are three different things")
        print("  (GROUNDING §0's three-way split), measured rather than argued.")

    # ---------------- N6 --------------------------------------------------
    print("\n" + "=" * 78)
    print("N6 -- hbar_eff triple")
    print("=" * 78)
    print("hbar_E = E/w with E = integral T_00 ONLY (D5'; Q w is never used as E).")
    print("%-22s %10s %12s %12s" % ("run", "hbar_Q=Q", "hbar_E=E/w_bar", "hbar_E/Q-1"))
    for r in rows:
        hE = r["E"] / r["wbar"]
        print("%-22s %10.3f %12.3f %12.6f" % (r["name"], r["Q"], hE, hE / r["Q"] - 1.0))
    hr = np.array([r["E"] / r["wbar"] / r["Q"] - 1.0 for r in srows])
    print("\n  hbar_E/Q - 1 across runs: %.4f .. %.4f (median %.4f)"
          % (hr.min(), hr.max(), float(np.median(hr))))
    print("  Note this is IDENTICALLY eps_bar (hbar_E/Q = E/(w Q)), i.e. the")
    print("  'hbar_eff = Q to 3-5%%' statement and the eps = 1-4%% statement are the")
    print("  SAME residual measured twice -- exactly what FOUNDATIONS predicted")
    print("  ('v70's 3-5%% and eps's 1-4%% are the same residual family'). The")
    print("  triple's third leg hbar_pk = p/k needs a boost series.")

    # look for a boost series in the archive
    cand = []
    for dd in DIRS + ["/space/scp"]:
        cand += glob.glob(os.path.join(dd, "**", "*boost*"), recursive=True)
        cand += glob.glob(os.path.join(dd, "**", "*_v0*"), recursive=True)
    cand = [c for c in cand if os.path.isfile(c)][:12]
    if cand:
        print("\n  possible boost data found (needs manual pairing for p and k):")
        for c in cand:
            print("    %s" % c)
        print("  N6 VERDICT: PARTIAL -- hbar_E and hbar_Q measured above; hbar_pk")
        print("  requires the v70 boost series to be re-analysed with E = int T_00.")
    else:
        print("\n  no boost series located under %s" % ", ".join(DIRS))
        print("  N6 VERDICT: PARTIAL -- hbar_E/Q measured (= eps_bar, above);")
        print("  hbar_pk = p/k NOT measurable from the archived v85 diags. Per the")
        print("  halving list ('N6 full new boost series: re-analyse v70 only; no")
        print("  new campaign') this leg stays open until v70 boost data is located")
        print("  or an EX-1 boost run supplies it.")

    path = os.path.join(HERE, "n56_action.tsv")
    with open(path, "w") as fh:
        fh.write("run\tQ\tT\tE\tE_em\tw_core\tw_bar\tT_over_wQ\tE_over_wQ\t"
                 "eps_bar\tgauge_shift\ttwin\n")
        for r in rows:
            fh.write("%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t"
                     "%.10g\t%.10g\t%.10g\t%.10g\n" %
                     (r["name"], r["Q"], r["T"], r["E"], r["Eem"], r["wcore"],
                      r["wbar"], r["T_over_wQ"], r["E_over_wQ"], r["eps_bar"],
                      r["gauge_shift"], r["twin"]))
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
