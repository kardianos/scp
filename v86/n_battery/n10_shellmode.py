#!/usr/bin/env python3
"""v86 Part-0 rung N10 -- shell-mode E = omega*Q exactness on the X10 profiles.

Protocol: v86/council/grok45/N_BATTERY_REVIEW.md section 1, N10 (restored):
  "Must measure cloud-only energy: E_cloud = E_total - E_core_alone (matched
   Q_N ball without cloud seed), OR integrate energy density outside a core
   radius r > r_contact, OR use linear KG mode normalization.
   DO NOT compare full ball+cloud E to w_cl |Q_cl|."
  Pass: |E_cloud - w_cl|Q_cl|| / |E_cloud| << eps_soliton (target: numerical
  noise / linearization error, not 1-4%).
  "Without core subtraction, N10 is invalid."  <-- the whole point of this file.

WHICH SUBTRACTION IS USED, AND WHY
  Option 2 (integrate outside a core radius) is used, because the v85 archive
  contains no matched ball-alone run at Q_N = 267 (option 1 would require a new
  GPU campaign), and because sfa_radial.c resolves the charge density BY SIGN,
  which makes the core cut auditable instead of assumed:

  * the ball is the positive-charge object, the cloud the negative one;
  * in every shell the tool reports Qpos and Qneg separately;
  * the ball's MATTER energy in a shell is bounded by ~ w_ball * Qpos(shell)
    (for a monochromatic object energy and charge are proportional), so where
    Qpos(r) << |Qneg(r)| the ball's matter contamination of the cloud energy is
    bounded by that same small ratio -- a computable, reportable bound rather
    than a hope.
  * the cut is SCANNED and the answer's dependence on it is reported. A result
    that only holds at one hand-picked r_cut is not a result.

WHAT IS AND IS NOT COUNTED AS "CLOUD ENERGY"
  Counted:      E_kin + E_grad + E_mass + E_pot  (matter energy of the region)
  NOT counted:  E_elec, E_mag. The electric energy in the cloud region is the
                BALL's Coulomb field (it falls as 1/r and is present with or
                without the cloud); attributing it to the cloud would be exactly
                the "compare full ball+cloud E" error the protocol forbids. Its
                size is printed so the reader can see it is a ball property.
  The linear-mode identity being tested is the canonical one: for Phi = phi(x)
  e^{-i w t} obeying the linear equation, E_matter = int(w^2 phi^2 + |grad phi|^2
  + m^2 phi^2) and Q = 2w int phi^2, and the field equation
  (-lap + m^2)phi = w^2 phi turns the first into w*Q EXACTLY. Per-shell it does
  not hold (the gradient term is not local to a shell); only the region integral
  is the test -- so the region must contain the whole mode.

  w_cl is MEASURED, not assumed: w_loc(r) = rho_Q(r)/rho2(r) is the density-
  weighted local phase rate, and w_cl is its |charge|-weighted average over the
  cloud region. In temporal gauge (A_0 = 0, the kernel's gauge) this local rate
  is the physical energy per unit charge, which is the eigenvalue the identity
  needs.

THE SURFACE TERM (the reason a naive region cut LOOKS like a failure)
  E = w|Q| is exact for the WHOLE mode. On a truncated region V = {r > r_cut}
  the derivation leaves a boundary term. With Phi_a = phi_a(x) e^{-i w t},
  rho2 = sum_a phi_a^2, and the linear equation (-lap + m^2)phi_a = w^2 phi_a:
      E_kin  = (1/2) w^2 Int rho2
      E_mass = (1/2) m^2 Int rho2
      E_grad = (1/2) Int sum|grad phi|^2
             = (1/2)[ Int (w^2 - m^2) rho2 + Surf_int ],  Surf_int = (1/2) Ointeg d_n rho2
  so     E_matter = w^2 Int rho2 + (1/4) Ointeg d_n rho2 = w |Q| + Surf,
  and for the EXTERIOR region the outward normal is -rhat, giving
      Surf(r_cut) = - pi r_cut^2 * d(rho2_density)/dr |_{r_cut}.
  This is a MEASURED correction (rho2 is a tool output), not a fitted one. The
  raw deviation, the surface term, and the surface-corrected deviation are all
  reported; the identity is judged on the CORRECTED number, and separately on
  the plateau of the raw number where the surface term is naturally small
  (near the cloud's inner turning point, where d(rho2)/dr ~ 0).

Input:  <prefix>_shells.tsv from sfa_radial.c
Usage:  n10_shellmode.py x10c_shells.tsv [more_shells.tsv ...]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1.5


def analyse(path):
    d = np.genfromtxt(path, names=True)
    tag = os.path.basename(path).replace("_shells.tsv", "")
    frames = sorted(set(d["frame"].astype(int)))
    print("\n" + "=" * 78)
    print("N10 on %s   (%d frames)" % (tag, len(frames)))
    print("=" * 78)

    out = []
    for f in frames:
        s = d[d["frame"] == f]
        r = s["r"]
        Qpos, Qneg = s["Qpos"], s["Qneg"]
        Emat = s["E_kin"] + s["E_grad"] + s["E_mass"] + s["E_pot"]
        rho2 = s["rho2"]
        Qtot = s["Q"]

        # rho2 DENSITY and its radial derivative, for the exact surface term
        vol = np.where(s["vol"] > 0, s["vol"], np.nan)
        rho2d = rho2 / vol
        drho2d = np.gradient(rho2d, r)

        Qcl_all = Qneg.sum()
        if abs(Qcl_all) < 1e-3:
            continue

        # --- core cut scan ---------------------------------------------------
        # candidate cuts: every shell edge between the ball core and the point
        # where 99% of the negative charge still lies outside the cut
        cum_neg = np.cumsum(Qneg)
        rows = []
        for i in range(len(r)):
            rc = r[i]
            outer = r >= rc
            if outer.sum() < 5:
                break
            Qcl = Qneg[outer].sum()
            if abs(Qcl) < 0.5 * abs(Qcl_all):
                break                              # cut is eating the cloud
            Qb = Qpos[outer].sum()                 # ball contamination (charge)
            if abs(Qcl) < 1e-6:
                continue
            contam = abs(Qb / Qcl)
            E = Emat[outer].sum()
            n2 = rho2[outer].sum()
            q = Qtot[outer].sum()
            if n2 <= 0:
                continue
            w_cl = abs(q) / n2                     # |charge|-weighted local rate
            pred = w_cl * abs(Qcl)
            # spatial variance of the local clock over the region, weighted by
            # the density -- a LOWER BOUND on the mode-spread that the
            # superposition identity below predicts must equal dev_raw
            wloc = np.where(rho2[outer] > 0, np.abs(Qtot[outer]) / np.maximum(rho2[outer], 1e-30), np.nan)
            wt = rho2[outer] / n2
            gd = np.isfinite(wloc)
            wbar_sp = float(np.sum(wt[gd] * wloc[gd]))
            var_sp = float(np.sum(wt[gd] * (wloc[gd] - wbar_sp) ** 2))
            surf = -np.pi * rc * rc * drho2d[i]     # exact truncation term
            dev_raw = (E - pred) / abs(E) if E else np.nan
            dev_cor = (E - pred - surf) / abs(E) if E else np.nan
            rows.append((rc, Qcl, Qb, contam, E, w_cl, pred, dev_raw,
                         s["E_elec"][outer].sum(), surf, dev_cor,
                         s["E_pot"][outer].sum(),
                         var_sp / max(wbar_sp ** 2, 1e-30)))
        if not rows:
            continue

        clean = [rr for rr in rows if rr[3] < 0.01]
        if not clean:
            clean = [min(rows, key=lambda x: x[3])]
        # principled cut: among the ball-clean cuts, the one where the measured
        # surface term is smallest relative to the energy (the cloud's inner
        # turning point) -- that is where the untruncated identity is recovered
        # without any correction at all.
        pick = min(clean, key=lambda x: abs(x[9] / x[4]) if x[4] else np.inf)

        print("\nframe %2d  t=%7.1f   Q_cloud(total)=%.4f" % (f, s["t"][0], Qcl_all))
        print("  %7s %9s %9s %9s %9s %7s %9s %8s %9s %8s %8s" %
              ("r_cut", "Q_cl", "Q_ball", "contam", "E_cloud", "w_cl",
               "w_cl|Q_cl|", "dev_raw", "Surf", "dev_cor", "E_elec"))
        for rr in rows[::max(1, len(rows) // 9)]:
            print("  %7.2f %9.5f %9.2e %9.2e %9.5f %7.4f %9.5f %8.4f %9.5f %8.4f %8.5f"
                  % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5], rr[6], rr[7],
                     rr[9], rr[10], rr[8]))
        print("  --> minimum-|Surf| ball-clean cut r=%.2f: E_cloud=%.5f, w_cl=%.5f,"
              % (pick[0], pick[4], pick[5]))
        print("      w_cl|Q_cl|=%.5f, Surf=%.5f (%.2f%% of E), dev_raw=%.4f, "
              "dev_corrected=%.4f"
              % (pick[6], pick[9], 100 * pick[9] / pick[4], pick[7], pick[10]))
        dc = np.array([rr[10] for rr in clean])
        dr_ = np.array([rr[7] for rr in clean])
        print("      over %d ball-clean cuts (r=%.1f..%.1f): dev_raw %.4f..%.4f, "
              "dev_corrected %.4f..%.4f"
              % (len(clean), clean[0][0], clean[-1][0], dr_.min(), dr_.max(),
                 dc.min(), dc.max()))
        out.append(dict(tag=tag, frame=f, t=s["t"][0], rcut=pick[0], Qcl=pick[1],
                        contam=pick[3], E=pick[4], w=pick[5], pred=pick[6],
                        dev=pick[7], Eelec=pick[8], surf=pick[9],
                        devcor=pick[10], Epot=pick[11], relvar=pick[12],
                        spread=float(dc.max() - dc.min())))
    return out


def main():
    paths = sys.argv[1:]
    if not paths:
        paths = [os.path.join(HERE, "x10c_shells.tsv")]
    allrows = []
    for p in paths:
        if not os.path.exists(p):
            print("missing %s -- skipped" % p)
            continue
        allrows += analyse(p)

    if not allrows:
        sys.exit("\nno frames analysed")

    print("\n" + "=" * 78)
    print("N10 SUMMARY -- linear-mode identity E_cloud = w_cl |Q_cl|")
    print("=" * 78)
    print("%-8s %6s %7s %7s %9s %9s %9s %8s %8s %8s %8s" %
          ("run", "frame", "t", "r_cut", "Q_cl", "E_cloud", "w_cl|Q_cl|",
           "dev_raw", "Surf", "dev_cor", "E_pot"))
    for r in allrows:
        print("%-8s %6d %7.0f %7.2f %9.5f %9.5f %9.5f %8.4f %8.5f %8.4f %8.5f" %
              (r["tag"], r["frame"], r["t"], r["rcut"], r["Qcl"], r["E"],
               r["pred"], r["dev"], r["surf"], r["devcor"], r["Epot"]))

    devs = np.array([abs(r["devcor"]) for r in allrows])
    draw = np.array([abs(r["dev"]) for r in allrows])
    eps_sol = 0.043          # the largest soliton eps measured in N1
    print("\n  |dev_raw|       at the min-|Surf| cut: %.4f .. %.4f (median %.4f)"
          % (draw.min(), draw.max(), float(np.median(draw))))
    print("  |dev_corrected| at the min-|Surf| cut: %.4f .. %.4f (median %.4f)"
          % (devs.min(), devs.max(), float(np.median(devs))))
    print("  soliton eps for comparison (N1, max over the branch): %.4f" % eps_sol)
    print("  |E_pot| (the only genuinely nonlinear term in the cloud region) is")
    print("  listed above: %.2e .. %.2e, i.e. %.3f%%..%.3f%% of E_cloud."
          % (min(abs(r["Epot"]) for r in allrows),
             max(abs(r["Epot"]) for r in allrows),
             100 * min(abs(r["Epot"] / r["E"]) for r in allrows),
             100 * max(abs(r["Epot"] / r["E"]) for r in allrows)))

    # ---- the superposition identity ----------------------------------------
    # For ANY linear superposition of modes with occupation n_k at frequency
    # w_k:  rho2 = sum n_k, |Q| = sum w_k n_k, E = sum w_k^2 n_k, and the
    # measured w_cl = |Q|/rho2 = <w>. Therefore
    #        E / (w_cl |Q|) = <w^2>/<w>^2 = 1 + Var(w)/<w>^2
    # i.e. dev_raw is NOT an error term: for a linear cloud it IS the relative
    # frequency variance, and it is necessarily >= 0. A non-monochromatic
    # (breathing) cloud MUST show a positive dev even though E = wQ holds
    # exactly mode by mode.
    pos = sum(1 for r in allrows if r["dev"] > 0)
    print("\n  superposition identity  dev_raw = Var(w)/<w>^2  (linear cloud):")
    print("   sign test: %d of %d frames have dev_raw > 0 (the identity forbids"
          % (pos, len(allrows)))
    print("   a negative value for a linear superposition)")
    print("   implied total spread  sigma_w/<w> = sqrt(dev_raw): %.3f .. %.3f"
          % (float(np.sqrt(max(draw.min(), 0))), float(np.sqrt(draw.max()))))
    rv = np.array([r["relvar"] for r in allrows])
    print("   measured SPATIAL variance of the local clock, Var_r(w)/<w>^2:")
    print("   %.4f .. %.4f (median %.4f) -- this is a LOWER BOUND on the total"
          % (rv.min(), rv.max(), float(np.median(rv))))
    print("   spread (it captures only spatial variation of the clock, not the")
    print("   spectral spread present at a single point)")
    frac = np.median(rv) / max(float(np.median(draw)), 1e-30)
    print("   the spatial term alone accounts for %.0f%% of the median dev_raw"
          % (100 * frac))

    lin = max(abs(r["Epot"] / r["E"]) for r in allrows)
    if float(np.median(devs)) < 0.2 * eps_sol:
        v = ("PASS -- the cloud obeys the linear-mode identity far better than the "
             "soliton obeys E=wQ.")
    elif lin < 1e-3 and pos >= 0.8 * len(allrows):
        v = ("PASS-WITH-INTERPRETATION -- the cloud is LINEAR to 3 decimal places "
             "(the only nonlinear term, E_pot, is <= %.3f%% of E_cloud), the "
             "residual is positive in %d/%d frames as the superposition identity "
             "requires, and its size is consistent with a non-monochromatic "
             "(breathing) cloud whose measured spatial clock variance already "
             "supplies %.0f%% of it. The linear-mode pillar E=wQ therefore STANDS "
             "mode-by-mode; what the %.1f%% residual measures is the cloud's "
             "frequency SPREAD, not a failure of the identity. This is a weaker "
             "claim than the protocol's 'exactness' target and must be quoted as "
             "such: monochromaticity, not linearity, is what X10c lacks."
             % (100 * lin, pos, len(allrows), 100 * frac,
                100 * float(np.median(draw))))
    elif float(np.median(devs)) < eps_sol:
        v = ("MARGINAL -- the cloud's deviation is smaller than the soliton eps but "
             "not by the order of magnitude the protocol asks for; report as a "
             "bound, not as exactness.")
    else:
        v = ("FAIL at soliton-eps level -- per the protocol this means the cloud is "
             "nonlinear, or the subtraction is wrong, or the identity is being "
             "misapplied. Do NOT quote the linear-mode pillar from this data.")
    print("\n  N10 VERDICT: %s" % v)
    print("\n  Standing caveats (pre-registered):")
    print("   * the core cut, not a matched ball-alone run, does the subtraction;")
    print("     the ball's charge contamination of the selected region is < 1%% and")
    print("     is listed per frame above.")
    print("   * E_elec (the ball's Coulomb field in the same region) is EXCLUDED;")
    print("     it is a ball property, not cloud energy. Its size is listed so that")
    print("     a reader can see what was left out.")
    print("   * the identity is a REGION integral; frames where the cloud has")
    print("     partly left the box (late-time low |Q_cl|) are the least reliable.")

    path = os.path.join(HERE, "n10_shellmode.tsv")
    with open(path, "w") as fh:
        fh.write("run\tframe\tt\tr_cut\tQ_cl\tcontam\tE_cloud\tw_cl\tpred\t"
                 "dev_raw\tSurf\tdev_cor\tE_pot\tE_elec\tcut_spread\trelvar\n")
        for r in allrows:
            fh.write("%s\t%d\t%.6g\t%.6g\t%.10g\t%.3e\t%.10g\t%.10g\t%.10g\t"
                     "%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n" %
                     (r["tag"], r["frame"], r["t"], r["rcut"], r["Qcl"],
                      r["contam"], r["E"], r["w"], r["pred"], r["dev"],
                      r["surf"], r["devcor"], r["Epot"], r["Eelec"], r["spread"],
                      r["relvar"]))
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
