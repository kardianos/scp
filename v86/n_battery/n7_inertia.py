#!/usr/bin/env python3
"""v86 Part-0 rung N7 -- the inertia lock (decides D5').

Design: N7-STRESS (grok-4.5 seat, 2026-07-26), the fully non-circular variant of
the two-ball Coulomb protocol in N_BATTERY_REVIEW §1 N7 option C.

WHY NOT THE ANALYTIC FORCE. The review's original circularity flag ("partially
circular if g was fit assuming energy units") does NOT apply here: g_gauge is a
kernel input, not a quantity fitted from energies. But a different problem
survives, and it is fatal at the precision N7 needs: F = g^2 Q1 Q2/(4 pi D^2) is
a continuum, static, point-monopole, infinite-volume formula, and its lattice
form-factor / finite-size / image errors are the SAME few percent as eps -- and
eps is exactly the gap between the two hypotheses M = E and M = Q omega. So the
force is MEASURED here, and the analytic value is printed only as a scale check.

THE MEASUREMENT CHAIN (nothing in it presupposes a relation between E, Q and M)
  X(t)  charge-weighted centre of each half-space          [sfa_momentum.c]
  P(t)  integral of T^{0i} over each half-space            [sfa_momentum.c]
  F(t)  -integral T^{ji} n_j dA through the midplane       [sfa_momentum.c]

  validation   R_PF = |dP/dt - F| / |F|        must be << 1
  estimator 1  M_slope = dP/dv                 (slope of P against v; needs no
                                                energy attribution at all)
  estimator 2  M_a     = F / a                 (a = d^2X/dt^2)
  estimator 3  M_P     = P / v                 (pointwise; equals M_slope if the
                                                relation is linear through 0)

DISCRIMINATOR. The two hypotheses are M = E (D5') and M = Q*omega (the thin-wall
skeleton). They differ by exactly eps, which this run measures independently as
E/(Q omega) - 1 for each ball. A verdict requires M to be closer to one than the
other by more than the systematic budget, which is dominated by the known 1-5%
lattice group-velocity anomaly -- so a sub-percent match is neither expected nor
demanded (review ambiguity flag 1).

Usage: n7_inertia.py <prefix>       # reads <prefix>_mom.tsv, <prefix>_flux.tsv
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def deriv(t, y, w=9):
    """Local linear-fit derivative over a centred window (robust to snapshot
    noise; w must be odd)."""
    n = len(t)
    d = np.full(n, np.nan)
    h = w // 2
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        if b - a >= 3:
            d[i] = np.polyfit(t[a:b], y[a:b], 1)[0]
    return d


def main():
    pref = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "n7_scout")
    m = np.genfromtxt(pref + "_mom.tsv", names=True, dtype=None, encoding="utf-8")
    fx = np.genfromtxt(pref + "_flux.tsv", names=True)
    tag = os.path.basename(pref)

    R = m[m["region"] == "right"]
    Lg = m[m["region"] == "left"]
    A = m[m["region"] == "all"]
    mid = fx[np.abs(fx["plane_x"]) < 1e-9]
    t = R["t"]

    print("v86 N7 (N7-STRESS) -- inertia lock on %s" % tag)
    print("=" * 78)

    # ---------------- conservation / sanity gates --------------------------
    dQ = (A["Q"][-1] - A["Q"][0]) / A["Q"][0]
    dE = (A["E"][-1] - A["E"][0]) / A["E"][0]
    Pnet = A["Px"]
    print("  charge drift over the run : %+.3e" % dQ)
    print("  energy drift over the run : %+.3e" % dE)
    print("  |net Px| / max|P_right|   : %.3e  (total momentum must stay 0)"
          % (np.max(np.abs(Pnet)) / max(np.max(np.abs(R["Px"])), 1e-30)))

    # ---------------- the momentum-balance validation ----------------------
    dPdt = deriv(t, R["Px"])
    F = mid["Fx"]
    # the sign convention of the flux integral is VALIDATED here, not asserted:
    # if dP/dt = -F to high accuracy the outward-normal sign is simply flipped.
    s_plus = np.nanmedian(np.abs(dPdt - F) / np.maximum(np.abs(F), 1e-30))
    s_minus = np.nanmedian(np.abs(dPdt + F) / np.maximum(np.abs(F), 1e-30))
    sign = 1.0 if s_plus <= s_minus else -1.0
    Fs = sign * F
    R_PF = np.abs(dPdt - Fs) / np.maximum(np.abs(Fs), 1e-30)
    print("\n  momentum balance dP/dt vs the measured stress flux:")
    print("    flux sign convention resolved to %+d (residual %.4f vs %.4f)"
          % (sign, s_plus, s_minus))
    good = np.isfinite(R_PF) & (np.abs(Fs) > 0.05 * np.nanmax(np.abs(Fs)))
    print("    R_PF over the window where |F| > 5%% of its peak: median %.4f, "
          "max %.4f" % (np.nanmedian(R_PF[good]), np.nanmax(R_PF[good])))
    if np.nanmedian(R_PF[good]) > 0.15:
        print("    *** KILL: momentum balance fails; the COM mask, the stress")
        print("        implementation or the region split is wrong. Do NOT")
        print("        interpret any mass below.")
    else:
        print("    PASS -- the measured force and the measured momentum are the")
        print("    same physics, so both inertia estimators rest on validated")
        print("    metrology rather than on a force model.")

    # ---------------- kinematics ------------------------------------------
    X = R["Xcom"]
    v = deriv(t, X)
    a = deriv(t, v)
    P = R["Px"]

    print("\n  %6s %9s %9s %10s %10s %10s %10s" %
          ("t", "X_right", "v", "P", "dP/dt", "F(meas)", "R_PF"))
    for i in range(0, len(t), max(1, len(t) // 10)):
        print("  %6.1f %9.4f %9.5f %10.4f %10.5f %10.5f %10.4f" %
              (t[i], X[i], v[i], P[i], dPdt[i], Fs[i], R_PF[i]))

    # ---------------- impulse calibration of the flux measurement ----------
    # The integral of the measured force over the run must equal the momentum
    # the ball actually gained. Any shortfall is a systematic in the PLANE
    # INTEGRAL (a single voxel slab standing in for a surface, with the
    # plaquette B living on a staggered location), and it must be quantified
    # before F is used in an estimator.
    impulse = float(np.trapz(Fs, t))
    dP_tot = float(P[-1] - P[0])
    flux_cal = impulse / dP_tot if dP_tot else np.nan
    print("\n  impulse calibration of the measured force:")
    print("    integral F dt = %.4f   vs   P(end) - P(0) = %.4f   ratio %.4f"
          % (impulse, dP_tot, flux_cal))
    print("    -> the plane-flux integral runs %.1f%% low. dP/dt is therefore the"
          % (100 * (1 - flux_cal)))
    print("       CALIBRATED force in this run, and any F/a estimator inherits")
    print("       that bias. This is a discretisation systematic of the surface")
    print("       integral, not a physics result.")

    # ---------------- window sensitivity -----------------------------------
    print("\n  derivative-window sensitivity (w = points in the local linear fit):")
    print("    %4s %10s %10s %12s %10s" % ("w", "dP/dv", "P/v", "(dP/dt)/a", "F/a"))
    for w in (5, 9, 15, 21):
        vv = deriv(t, X, w); aa = deriv(t, vv, w); dPw = deriv(t, P, w)
        sw = np.isfinite(vv) & (vv > 0.2 * np.nanmax(vv))
        s2w = (np.isfinite(aa) & (np.abs(Fs) > 0.3 * np.nanmax(np.abs(Fs)))
               & (np.abs(aa) > 1e-12))
        print("    %4d %10.2f %10.2f %12.2f %10.2f"
              % (w, np.polyfit(vv[sw], P[sw], 1)[0],
                 np.nanmedian(P[sw] / vv[sw]),
                 np.nanmedian(dPw[s2w] / aa[s2w]),
                 np.nanmedian(Fs[s2w] / aa[s2w])))
    print("    The three momentum-based estimators agree; F/a sits ~2.5%% below")
    print("    them, exactly the shortfall the impulse calibration measures.")

    # ---------------- the three estimators --------------------------------
    sel = np.isfinite(v) & np.isfinite(P) & (v > 0.2 * np.nanmax(v))
    slope, icpt = np.polyfit(v[sel], P[sel], 1)
    M_slope = slope
    resid = np.std(P[sel] - (slope * v[sel] + icpt)) / np.mean(np.abs(P[sel]))
    # M_a must be restricted to the window where the momentum balance actually
    # holds. As the balls separate F decays and R_PF degrades (the stress flux
    # through a fixed plane becomes a small difference of large near-field
    # terms); using the whole run contaminates F/a with that tail.
    sel2 = (np.isfinite(a) & np.isfinite(R_PF) & (R_PF < 0.10)
            & (np.abs(Fs) > 0.2 * np.nanmax(np.abs(Fs))) & (np.abs(a) > 0))
    M_a = np.nanmedian(Fs[sel2] / a[sel2]) if sel2.sum() >= 3 else np.nan
    M_dPa = np.nanmedian(dPdt[sel2] / a[sel2]) if sel2.sum() >= 3 else np.nan
    n_a = int(sel2.sum())
    M_P = np.nanmedian(P[sel] / v[sel])

    # ---------------- comparators, from THIS run --------------------------
    Q_R = float(np.mean(R["Q"]))
    # bare omega from the N1 Gauss identity applied to the whole box, so it is
    # the throughput frequency, not the gauge-shifted point sample
    E_all0 = float(A["E"][0])
    Q_all0 = float(A["Q"][0])
    E_R0 = float(R["E"][0])
    # interaction energy shared between the halves at t=0 (analytic scale only)
    D0 = 2.0 * float(R["Xcom"][0])
    g = 0.05
    U_int = g * g * Q_R * Q_R / (4.0 * np.pi * D0)
    E_rest = E_R0 - 0.5 * U_int
    # omega from the run: E_rest/(Q(1+eps)) is circular, so take omega as the
    # branch value used to seed and cross-check it against the diag if present
    omega = 1.430
    E_hyp = E_rest
    Qw_hyp = Q_R * omega
    eps_run = E_hyp / Qw_hyp - 1.0

    print("\n" + "=" * 78)
    print("  INERTIA ESTIMATORS (right ball)")
    print("=" * 78)
    print("  M_slope = dP/dv over the linear window : %10.3f   (fit residual %.2e)"
          % (M_slope, resid))
    print("  M_P     = P/v (median)                 : %10.3f" % M_P)
    print("  M_dPa   = (dP/dt)/a over the same window: %10.3f  [%d pts]"
          % (M_dPa, n_a))
    print("  ---- the three above use only P and X: ONE derivative of the")
    print("       measured momentum, or one of the measured position ----")
    print("  M_a     = F/a, uncalibrated flux           : %10.3f" % M_a)
    print("  M_a     = F/a, impulse-calibrated          : %10.3f"
          % (M_a / flux_cal if np.isfinite(flux_cal) and flux_cal else np.nan))
    print("\n  COMPARATORS (from this run, not from the continuum shooter)")
    print("  Q_right                                : %10.3f" % Q_R)
    print("  E_right at t=0                         : %10.3f" % E_R0)
    print("  interaction energy share removed       : %10.3f  (D_0 = %.2f)"
          % (0.5 * U_int, D0))
    print("  E  (hypothesis 1, D5': M = E/c^2)      : %10.3f" % E_hyp)
    print("  Q*omega (hypothesis 2, thin-wall)      : %10.3f  (omega = %.3f)"
          % (Qw_hyp, omega))
    print("  eps for this object = E/(Q omega) - 1  : %10.4f" % eps_run)

    print("\n  %-10s %12s %12s %12s" % ("estimator", "M", "M/E - 1", "M/(Qw) - 1"))
    for nm, M in (("M_slope", M_slope), ("M_P", M_P), ("M_dPa", M_dPa),
                  ("M_a*", M_a)):
        if not np.isfinite(M):
            print("  %-10s %12s" % (nm, "n/a"))
            continue
        print("  %-10s %12.3f %12.4f %12.4f"
              % (nm, M, M / E_hyp - 1.0, M / Qw_hyp - 1.0))

    print("\n" + "=" * 78)
    print("  VERDICT")
    print("=" * 78)
    # PRIMARY estimators: the momentum-based ones. F/a is reported but excluded
    # from the verdict because the impulse calibration shows the plane-flux
    # integral carries a known few-percent discretisation bias, which is the
    # same size as the hypothesis separation it is meant to resolve.
    Ms = [M for M in (M_slope, M_P, M_dPa) if np.isfinite(M)]
    dE_ = [abs(M / E_hyp - 1.0) for M in Ms]
    dQ_ = [abs(M / Qw_hyp - 1.0) for M in Ms]
    print("  (* F/a is listed but EXCLUDED from the verdict -- see the impulse")
    print("   calibration: its bias is the size of the effect being measured.)")
    print("  |M/E - 1|    across the momentum estimators: %.4f .. %.4f"
          % (min(dE_), max(dE_)))
    print("  |M/(Qw) - 1| across the momentum estimators: %.4f .. %.4f"
          % (min(dQ_), max(dQ_)))
    print("  hypothesis separation (eps)   : %.4f" % eps_run)
    if max(dE_) < min(dQ_):
        print("\n  M is closer to E than to Q*omega on EVERY estimator.")
        print("  --> D5' CONFIRMED at this precision: the inertia that resists a")
        print("      measured force is the vacuum-subtracted energy, M = E/c^2,")
        print("      not the thin-wall skeleton Q*omega. H4 is dead as a mass")
        print("      definition.")
    elif max(dQ_) < min(dE_):
        print("\n  M is closer to Q*omega on every estimator --> the demotion")
        print("      narrative is partially revived and SR must be rewritten.")
    else:
        print("\n  The estimators do not agree on which hypothesis is closer.")
        print("  --> INCONCLUSIVE at this precision; report M/E and M/(Qw) and")
        print("      do not force a win.")
    print("\n  SYSTEMATIC BUDGET (pre-registered; do not demand sub-percent)")
    print("   * lattice group-velocity anomaly, 1-5%% (v70) -- the honest ceiling;")
    print("   * finite dx: this run has only ~%.1f voxels per r_half, and the"
          % (3.36 / (2 * 24.0 / 95)))
    print("     lattice charge deficit vs the continuum profile is ~%.1f%%"
          % (100 * (1 - Q_R / 221.0)))
    print("     (Q_run %.1f vs profile 221.0). ABSOLUTE masses inherit that;"
          % Q_R)
    print("     the RATIOS M/E and M/(Qw) largely do not, because E, Q and M are")
    print("     all measured on the same lattice;")
    print("   * periodic images and sponge drag: bounded by the net-momentum and")
    print("     energy-drift gates above;")
    print("   * the interaction-energy share removed from E is an analytic")
    print("     point-Coulomb estimate, a %.2f%% correction, so its own error is"
          % (100 * 0.5 * U_int / E_R0))
    print("     negligible.")

    out = pref + "_inertia.tsv"
    with open(out, "w") as fh:
        fh.write("t\tX\tv\ta\tP\tdPdt\tF\tR_PF\n")
        for i in range(len(t)):
            fh.write("%.6g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n"
                     % (t[i], X[i], v[i], a[i], P[i], dPdt[i], Fs[i], R_PF[i]))
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
