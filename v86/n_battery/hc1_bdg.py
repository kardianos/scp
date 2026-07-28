#!/usr/bin/env python3
"""v86 census rung HC-1 -- ball-background linear-response (BdG) spectrum.

Deliverables (PROPOSAL I.4 + Amendment 1):
  * the MODE CATALOG the whole census scores against, with MULTIPOLE content
    (HC-2 classifies multipole-FIRST, so l must be a catalog column);
  * the CONTINUUM EDGES (GROUNDING §2: gap arithmetic runs against the matter
    band, and near-edge verdicts need the true edge, not m=1.5 folklore);
  * n(H_omega) -- the constrained negative index of the linearized energy,
    which the corrected GSS criterion (GROUNDING §1) needs and which HC-3
    currently assumes to be 1.

SCOPE: ungauged (g=0). This is deliberate and matches HC-3: GROUNDING §1
caveat (i) says the gauged system is Coulomb-phase and gauged GSS is a
heuristic with the UNGAUGED theorem as the anchor. HC-1 supplies that anchor.
A gauged BdG needs the A_0 perturbation solved self-consistently and is a
separate build.

--------------------------------------------------------------------------
DERIVATION (three components, product potential, symmetric background)
--------------------------------------------------------------------------
Energy density  e = sum_a (1/2)(|Phidot_a|^2 + |grad Phi_a|^2 + m^2|Phi_a|^2)
                    + Vt(s),  s = prod_a |Phi_a|^2,  Vt(s) = (mu/2)s/(1+kappa s)
EOM             Phiddot_a = lap Phi_a - m^2 Phi_a - P_a Phi_a,
                P_a = 2 Vt'(s) s/|Phi_a|^2 ;  background P0 = 2Vt'(f^6) f^4.

Write Phi_a = e^{i w t}(f(r) + eta_a), eta_a = x_a + i y_a. To first order
    rho_a = f^2 + 2 f x_a,   s = f^6 + 2 f^5 sum_b x_b,
    f * delta P_a = B sum_b x_b + A (sum_b x_b - x_a),
    A = 4 f^4 Vt'(f^6) = 2 P0,      B = 4 f^10 Vt''(f^6) = -4 kappa mu f^10/(1+kappa f^6)^3
giving the coupled linear system
    xddot_a = -(L_x x)_a + 2 w ydot_a
    yddot_a = -(L0  y)_a - 2 w xdot_a
with
    L0        = -lap + m^2 - w^2 + P0
    (L_x x)_a = L0 x_a - A x_a + (A+B) sum_b x_b .
In the flavor basis (symmetric mode s = sum_a/sqrt3, two orthogonal modes):
    L_x^sym  = L0 + 2A + 3B = -lap + m^2 - w^2 + 5 P0 + 3B
    L_x^flav = L0 -  A      = -lap + m^2 - w^2 -   P0
    L_y      = L0            (all three channels)

CHECKS THIS SCRIPT RUNS (each is a real falsification opportunity):
  1. L0 must annihilate u = r f at l=0 -- that is the U(1) phase Goldstone.
     Three such modes exist (one per global U(1)); the code verifies the l=0
     eigenvalue of L0 is zero to the discretization floor.
  2. L0 >= 0 (its zero mode is nodeless, hence the ground state) -> the y block
     contributes NO negative directions.
  3. L_x^flav = L0 - A with A = 2P0 and P0 = mu f^4/(1+kappa f^6)^2 < 0 for
     mu < 0, so -A > 0 and L_x^flav >= L0 >= 0 -> the flavor channels
     contribute NO negative directions either.
  4. Therefore n(H_omega) is the negative count of L_x^sym alone, expected to
     be exactly 1 (the l=0 dilational direction) -- the VK/GSS anchor.
  Any of 1-4 failing is a real result, not a bug to be papered over; the script
  reports what it finds.

Second variation: d(w) = E - w Q has  delta^2 d = <x, L_x x> + <y, L0 y>
(the cross terms are total derivatives), so n(H_omega) = n(L_x) + n(L0)
counted over ALL l with the (2l+1) degeneracy noted but not multiplied in
(the index is per-mode; the degeneracy is reported separately).

Mode frequencies: with (x,y) ~ e^{lam t},
    (L_x + lam^2) x = -4 w^2 lam^2 (L0 + lam^2)^{-1} x
which this script solves as a quadratic eigenvalue problem projected onto the
lowest eigenvectors of L0 and L_x (companion linearization), giving Omega with
lam = -i Omega. Real Omega = oscillation, Re lam > 0 = instability.

THE CONTINUUM EDGE IN Omega (this is the number the census needs, and it is
NOT sqrt(m^2 - w^2)). At large r the channel potentials vanish and the coupled
system reduces to, for a plane wave of wavenumber k and G = m^2 - w^2 + k^2,
    (Omega^2 - G)^2 - 4 w^2 Omega^2 = 0   ->   Omega = -w +/- sqrt(m^2 + k^2),
so the lowest radiating co-rotating frequency is
    Omega_c = m - w                        <-- the matter continuum edge
which is exactly the lab-frame statement "a perturbation at co-rotating Omega
appears at lab frequency w + Omega and radiates iff w + Omega >= m" (GROUNDING
§2's band [m, w_max]). Using sqrt(m^2 - w^2) here would overstate the gap by a
factor of ~6 at w = 1.42 and would wrongly promote box-discretized continuum
states to "bound internal modes". The lattice band TOP is w_max^2 = m^2 +
12/dx^2, i.e. Omega_max = w_max - w, and per the council warning any near-edge
verdict is refinement-gated.

BOX STATES vs BOUND STATES. In a finite box the continuum is discretized with
spacing set by Rmax; near threshold Omega_n ~ sqrt(m^2 + (n pi/Rmax)^2) - w.
This script therefore runs EVERY spectrum at two box sizes and reports which
sub-threshold frequencies MOVE with Rmax (box artifacts) and which do not
(genuine bound internal modes). A mode catalog that does not do this cannot
tell an internal mode from a standing wave.

Radial discretization: u(r) = r phi(r) turns the radial operator into the
standard 1D Schrodinger form -u'' + [l(l+1)/r^2 + V(r)] u = lambda u on
(0, Rmax] with u(0) = u(Rmax) = 0 -- symmetric tridiagonal, so eigenvalues are
real by construction and the negative count is trustworthy.
"""
import os
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal, eig

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/d/code/scp"
sys.path.insert(0, os.path.join(ROOT, "v69/theory"))
import gauged_shooter_fast as G

M2, MU, KAP, M = G.M2, G.MU, G.KAP, G.M


def background(w, hr=0.01, Rmax=40.0):
    """Ungauged radial profile f(r) on a uniform grid r = h..Rmax."""
    f0g = G.load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
    f, chi, ok, _, rn = G.solve(1.45, 0.0, f0g, np.zeros(G.N))
    if not ok:
        raise RuntimeError("seed failed")
    wc, step = 1.45, (0.004 if w > 1.45 else -0.004)
    while abs(wc - w) > 1e-9:
        wn = wc + step
        if (step < 0 and wn < w) or (step > 0 and wn > w):
            wn = w
        fn, cn, okk, _, _ = G.solve(wn, 0.0, f, chi)
        if not okk:
            step *= 0.5
            if abs(step) < 1e-6:
                raise RuntimeError("continuation stalled at w=%.5f" % wc)
            continue
        f, chi, wc = fn, cn, wn
    n = int(Rmax / hr)
    r = (np.arange(n) + 1) * hr
    fr = np.interp(r, G.r, f, right=0.0)
    return r, fr, hr, rn


def potentials(f):
    """P0, A, B and the three channel potentials (without the l-barrier)."""
    d = 1.0 + KAP * f ** 6
    P0 = MU * f ** 4 / d ** 2                 # = 2 Vt'(f^6) f^4
    A = 2.0 * P0
    B = -4.0 * KAP * MU * f ** 10 / d ** 3
    return P0, A, B


def channel_V(w, f, name):
    P0, A, B = potentials(f)
    base = M2 - w * w
    if name == "L0":
        return base + P0
    if name == "Lx_sym":
        return base + 5.0 * P0 + 3.0 * B
    if name == "Lx_flav":
        return base - P0
    raise ValueError(name)


def eigs(V, r, h, l, k=14):
    """Lowest k eigenpairs of -u'' + [l(l+1)/r^2 + V] u on the u = r phi grid."""
    diag = 2.0 / h ** 2 + l * (l + 1) / r ** 2 + V
    off = -np.ones(len(r) - 1) / h ** 2
    vals, vecs = eigh_tridiagonal(diag, off, select="i",
                                  select_range=(0, min(k, len(r)) - 1))
    return vals, vecs


def main():
    ws = [float(x) for x in (sys.argv[1:] or ["1.34", "1.39", "1.42", "1.45", "1.47"])]
    LMAX = 4
    print("v86 HC-1 -- BdG mode catalog on the UNGAUGED Q-ball background")
    print("grid: u=r*phi, h=0.01, ADAPTIVE Rmax = max(40, 5 r_rms), l=0..%d\n" % LMAX)

    catalog = []
    boundlist = []
    rejected = []
    for w in ws:
        # ADAPTIVE BOX. A fixed Rmax truncates the huge thin-wall balls, which
        # makes the Goldstone check fail and then silently reports n(H_omega)=0
        # (the review caught exactly this at w=1.315, where Q = 3.7e5). Size the
        # box to the object from a scouting solve.
        r0, f0, h0, _ = background(w, Rmax=60.0)
        rr0 = float(np.trapz(f0 * f0 * r0 ** 4, dx=h0)
                    / max(np.trapz(f0 * f0 * r0 * r0, dx=h0), 1e-30)) ** 0.5
        Rmax = max(40.0, 5.0 * rr0)
        r, f, h, rn = background(w, Rmax=Rmax)
        gap = M2 - w * w                    # operator-level edge: P0,A,B -> 0
        # so every channel potential tends to m^2 - w^2 and lambda < gap is a
        # bound eigenvalue OF THAT OPERATOR. The physical radiation threshold of
        # the COUPLED system is Omega_c = m - w (see the header derivation).
        Om_c = M - w
        xi = 1.0 / np.sqrt(gap)
        Q = np.trapz(3.0 * w * f * f * 4 * np.pi * r * r, dx=h)
        print("=" * 78)
        print("omega = %.4f   Q = %.2f   xi = %.3f   (shooter resid %.1e)"
              % (w, Q, xi, rn))
        print("operator edge lambda_c = m^2 - w^2 = %.5f   |   PHYSICAL radiation "
              "threshold Omega_c = m - w = %.5f" % (gap, Om_c))
        print("=" * 78)

        # ---- check 1: the U(1) Goldstone --------------------------------
        v0, u0 = eigs(channel_V(w, f, "L0"), r, h, 0, k=6)
        gold = v0[0]
        print("  check 1  L0 lowest eigenvalue at l=0 (the U(1) Goldstone, must be 0):"
              "  %.3e   [floor ~ h^2 = %.1e]" % (gold, h * h))
        ov = np.abs(np.dot(u0[:, 0] / np.linalg.norm(u0[:, 0]),
                           (r * f) / np.linalg.norm(r * f)))
        print("           overlap of that eigenvector with u = r*f: %.6f" % ov)
        print("  box: Rmax=%.1f = %.1f x r_rms (r_rms=%.2f)" % (Rmax, Rmax / rr0, rr0))
        # HARD GATE. A failed Goldstone means the box does not contain the
        # object, and nothing measured at such an omega may be quoted -- least
        # of all a spurious n(H_omega) = 0.
        if abs(gold) > 100.0 * h * h or ov < 0.99:
            print("  *** OMEGA REJECTED: Goldstone check FAILED "
                  "(|lambda_0| = %.2e vs floor %.1e, overlap %.4f < 0.99)."
                  % (abs(gold), 100.0 * h * h, ov))
            print("      The box does not resolve this object. NOTHING from this")
            print("      omega -- including n(H_omega) -- is reported or quoted.")
            rejected.append((w, Q, gold, ov, Rmax, rr0))
            continue

        nneg_tot = 0
        rows = []
        for name in ("L0", "Lx_sym", "Lx_flav"):
            V = channel_V(w, f, name)
            for l in range(LMAX + 1):
                vals, _ = eigs(V, r, h, l, k=14)
                nneg = int((vals < -10.0 * h * h).sum())
                bound = vals[vals < gap - 1e-9]
                mult = 3 if name == "L0" else (1 if name == "Lx_sym" else 2)
                nneg_tot += nneg * mult
                rows.append((name, l, nneg, mult, vals[:6], bound))
        print("\n  %-9s %3s %6s %6s %s" %
              ("channel", "l", "n_neg", "mult", "lowest eigenvalues "
                                                "(bound ones are < lambda_c)"))
        for name, l, nneg, mult, vals, bound in rows:
            mark = "".join("*" if v < gap - 1e-9 else " " for v in vals[:6])
            print("  %-9s %3d %6d %6d  %s   [%s]" %
                  (name, l, nneg, mult,
                   " ".join("%9.5f" % v for v in vals[:6]), mark))
        print("\n  (* = below the continuum edge lambda_c = %.5f, i.e. a BOUND "
              "radial mode)" % gap)
        print("  n(H_omega) = sum over channels/l of n_neg x multiplicity = %d"
              % nneg_tot)

        # ---- checks 2,3 -------------------------------------------------
        n_L0 = sum(nn * mu for nm, l, nn, mu, _, _ in rows if nm == "L0")
        n_fl = sum(nn * mu for nm, l, nn, mu, _, _ in rows if nm == "Lx_flav")
        n_sy = sum(nn * mu for nm, l, nn, mu, _, _ in rows if nm == "Lx_sym")
        print("  check 2  L0 contributes %d negative directions (expected 0)" % n_L0)
        print("  check 3  Lx_flav contributes %d (expected 0: A = 2 P0 < 0 so "
              "Lx_flav = L0 - A >= L0)" % n_fl)
        print("  check 4  Lx_sym contributes %d (expected 1: the l=0 dilational "
              "direction)" % n_sy)
        verdict = ("n(H_omega) = 1 -- the GSS/VK anchor CONFIRMED, and the "
                   "provisional n(H)=1 used by HC-3 is now MEASURED"
                   if nneg_tot == 1 else
                   "n(H_omega) = %d -- NOT 1. HC-3's provisional assumption is "
                   "WRONG at this omega and the GSS matching must be redone."
                   % nneg_tot)
        print("  --> %s" % verdict)

        # ---- the mode catalog with multipole labels ----------------------
        print("\n  MODE CATALOG (bound radial modes, by multipole) -- HC-2 input")
        print("  %-9s %3s %4s %11s %11s %10s %s" %
              ("channel", "l", "node", "lambda", "Omega_est", "radiates?",
               "classification"))
        cat_rows = []
        for name, l, nneg, mult, vals, bound in rows:
            for k, lam in enumerate(bound):
                # Frequency estimate from the COUPLED dispersion, not sqrt(lam):
                # the asymptotic relation Omega = -w + sqrt(m^2 + k^2) with
                # lam = m^2 - w^2 + k^2 gives k^2 = lam - (m^2 - w^2), hence
                #     Omega(lam) = sqrt(w^2 + lam) - w
                # which maps lam = 0 -> Omega = 0 (Goldstone) and
                # lam = lam_c = m^2 - w^2 -> Omega = m - w = Omega_c exactly.
                # Using sqrt(lam) here (the naive single-operator guess) inflates
                # every frequency by roughly a factor 2w/Omega and would put
                # sub-threshold modes above the continuum edge.
                Om = np.sqrt(w * w + lam) - w if lam > -w * w else np.nan
                # MULTIPOLE FIRST (GROUNDING §2 / council fix #3): any mode with
                # a time-varying l>=1 moment couples to the massless gauge
                # channel, open at ALL frequencies -> golden rule, not arithmetic
                # exact zero modes are symmetries, not excitations
                zero = abs(lam) < 50.0 * h * h
                if zero and name == "L0" and l == 0:
                    rad, cls = "n/a", "U(1) PHASE GOLDSTONE (x3, one per global U(1))"
                elif zero and name == "Lx_sym" and l == 1:
                    rad, cls = "n/a", "TRANSLATION GOLDSTONE (x3: l=1, m=-1,0,+1)"
                elif l >= 1:
                    rad = "YES (l>=1)"
                    cls = "golden-rule width vs the massless A channel (gauged case)"
                else:
                    rad = "no (l=0)"
                    if lam < -50.0 * h * h:
                        cls = "NEGATIVE direction (GSS index contribution)"
                    else:
                        cls = ("monopole: gap arithmetic vs Omega_c = m - w = %.4f"
                               % Om_c)
                print("  %-9s %3d %4d %11.5f %11s %10s %s" %
                      (name, l, k, lam,
                       ("%.5f" % Om) if np.isfinite(Om) else "imag",
                       rad, cls))
                cat_rows.append((w, name, l, k, lam, Om, mult))
        catalog += cat_rows

        # ---- coupled quadratic pencil for the l=0 symmetric sector -------
        # (x in Lx_sym, y in L0) -- this is where the real internal-mode
        # frequency of the breathing/dilational sector lives.
        print("\n  coupled (x,y) l=0 symmetric sector -- BOX-SIZE TEST")
        print("  (a sub-threshold Omega that MOVES with Rmax is a discretized")
        print("   continuum state, not an internal mode)")
        res = {}
        for RB in (Rmax, 1.5 * Rmax):
            rb, fb, hb, _ = background(w, Rmax=RB)
            Vx = channel_V(w, fb, "Lx_sym")
            Vy = channel_V(w, fb, "L0")
            nb = 40
            ax, Ux = eigs(Vx, rb, hb, 0, k=nb)
            ay, Uy = eigs(Vy, rb, hb, 0, k=nb)
            S = Ux.T @ Uy                              # basis overlap
            Lxm = np.diag(ax)
            L0m = S @ np.diag(ay) @ S.T
            Z = np.zeros((nb, nb))
            I = np.eye(nb)
            Amat = np.block([
                [Z, Z, I, Z],
                [Z, Z, Z, I],
                [-Lxm, Z, Z, 2 * w * I],
                [Z, -L0m, -2 * w * I, Z],
            ])
            ev = eig(Amat, right=False)
            res[RB] = (np.sort(np.abs(ev.imag[ev.imag > 1e-8])),
                       ev.real[ev.real > 1e-6])
        Om_a, gr_a = res[Rmax]
        Om_b, _ = res[1.5 * Rmax]
        print("    Rmax=%-5g Omega: %s" % (Rmax, " ".join("%.5f" % o for o in Om_a[:7])))
        print("    Rmax=%-5g Omega: %s" % (1.5 * Rmax,
                                           " ".join("%.5f" % o for o in Om_b[:7])))
        print("    predicted box-mode ladder sqrt(m^2+(n pi/Rmax)^2) - w:")
        for RB in (Rmax, 1.5 * Rmax):
            lad = [np.sqrt(M2 + (n * np.pi / RB) ** 2) - w for n in range(1, 6)]
            print("      Rmax=%-5g %s" % (RB, " ".join("%.5f" % x for x in lad)))
        sub_a = Om_a[(Om_a > 1e-6) & (Om_a < Om_c - 1e-6)]
        sub_b = Om_b[(Om_b > 1e-6) & (Om_b < Om_c - 1e-6)]
        # Review Finding 10: a flat 2% match is too weak near threshold, where
        # the box ladder is dense and accidental alignment is likely. Require
        # the shift between boxes to be small compared with the LOCAL LEVEL
        # SPACING of the larger box -- the scale an accidental match would have.
        spacing_b = (float(np.min(np.diff(sub_b))) if len(sub_b) > 1
                     else (Om_c if len(sub_b) else np.inf))
        tol = min(0.02 * max(float(np.max(sub_a)) if len(sub_a) else 1.0, 1e-9),
                  0.25 * spacing_b)
        stable = [o for o in sub_a
                  if len(sub_b) and np.min(np.abs(sub_b - o)) < tol]
        print("    box-stability tolerance %.2e = min(2%% of the highest "
              "sub-threshold mode, 0.25 x larger-box level spacing %.2e)"
              % (tol, spacing_b))
        print("    sub-threshold (Omega < Omega_c = %.5f) frequencies: %d at "
              "Rmax=%g, %d at Rmax=%g; %d survive the box change"
              % (Om_c, len(sub_a), Rmax, len(sub_b), 1.5 * Rmax, len(stable)))
        if stable:
            print("    --> GENUINE BOUND INTERNAL MODE(S): %s"
                  % " ".join("%.5f" % o for o in stable))
        else:
            print("    --> NO bound internal mode in this sector: every")
            print("        sub-threshold frequency moves with the box, i.e. is a")
            print("        discretized continuum state. Every genuine excitation of")
            print("        the ball therefore sits AT OR ABOVE the radiation")
            print("        threshold and CARRIES A WIDTH -- the census claim C1's")
            print("        'every other localized excitation carries a width'")
            print("        holds here by absence of bound structure, not by")
            print("        arithmetic on a mode that does not exist.")
        rh = float(np.trapz(f * f * r ** 4, dx=h) / max(np.trapz(f * f * r * r, dx=h), 1e-30)) ** 0.5
        if Rmax < 2.0 * rh:
            print("    *** BOX WARNING: Rmax=%g is only %.2f x the ball's rms radius"
                  " (%.2f); the bound-mode COUNT is not converged in box size and"
                  " must be re-run larger before being quoted." % (Rmax, Rmax / rh, rh))
        boundlist.append((w, Q, rh, Om_c, list(stable), len(sub_a),
                          len(sub_b), nneg_tot, Rmax))
        print("    growth rates Re(lambda) > 0 at Rmax=%g: %s" %
              (Rmax, " ".join("%.5f" % gr for gr in np.sort(gr_a)[::-1][:4])
               if len(gr_a) else "none (no linear instability in this sector)"))

    path = os.path.join(HERE, "hc1_catalog.tsv")
    with open(path, "w") as fh:
        fh.write("w\tchannel\tl\tnode\tlambda\tOmega\tmultiplicity\n")
        for row in catalog:
            fh.write("%.6g\t%s\t%d\t%d\t%.10g\t%.10g\t%d\n" %
                     (row[0], row[1], row[2], row[3], row[4],
                      row[5] if np.isfinite(row[5]) else float("nan"), row[6]))
    print("\nwrote %s (%d catalog entries)" % (path, len(catalog)))

    bp = os.path.join(HERE, "hc1_bound_modes.tsv")
    with open(bp, "w") as fh:
        fh.write("w\tQ\tr_rms\tRmax\tOmega_c\tn_boxstable\tOmegas\t"
                 "n_sub_small\tn_sub_large\tn_H_omega\n")
        for (w, Q, rh, Oc, st, na, nb, nh1, RB) in boundlist:
            fh.write("%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%s\t%d\t%d\t%d\n" %
                     (w, Q, rh, RB, Oc, len(st),
                      ",".join("%.6g" % o for o in st) if st else "-", na, nb, nh1))
    print("wrote %s" % bp)
    print("\n" + "=" * 78)
    print("HC-1 SUMMARY -- bound internal modes along the branch")
    print("=" * 78)
    print("%-8s %11s %8s %7s %9s %8s %6s %s" %
          ("w", "Q", "r_rms", "Rmax", "Omega_c", "n(H_w)", "count",
           "box-stable Omega"))
    for (w, Q, rh, Oc, st, na, nb, nh1, RB) in boundlist:
        print("%-8.4f %11.1f %8.2f %7.1f %9.5f %8d %6d %s" %
              (w, Q, rh, RB, Oc, nh1, len(st),
               " ".join("%.5f" % o for o in st) if st else "(none)"))
    if rejected:
        print("\nREJECTED omegas (Goldstone check failed -- box too small):")
        for (w, Q, gold, ov, RB, rr) in rejected:
            print("  w=%.4f Q=%.1f lambda_0=%.2e overlap=%.4f Rmax=%.1f r_rms=%.2f"
                  % (w, Q, gold, ov, RB, rr))
        print("  NOTHING from these omegas is quoted anywhere in the results.")
    nh = sorted(set(b[7] for b in boundlist))
    print("\nn(H_omega) over the %d ACCEPTED omegas: %s" % (len(boundlist), nh))
    if nh == [1]:
        print("The GSS/VK anchor n(H_omega) = 1 is a MEASURED input to HC-3 over")
        print("the accepted range [%.4f, %.4f] -- and ONLY over that range. It is"
              % (min(b[0] for b in boundlist), max(b[0] for b in boundlist)))
        print("not established outside it, and this run is UNGAUGED, so the")
        print("Coulomb-phase (production, g=0.05) case remains heuristic.")
    else:
        print("n(H_omega) is NOT uniformly 1 over the accepted range -- HC-3's")
        print("provisional assumption fails somewhere and the GSS matching must")
        print("be redone omega by omega.")


if __name__ == "__main__":
    main()
