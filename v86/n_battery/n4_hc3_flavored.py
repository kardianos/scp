#!/usr/bin/env python3
"""v86 rungs N4 (multi-flavor Legendre excess) and HC-3 (corrected GSS signature).

Protocols:
  N4  -- v86/council/grok45/N_BATTERY_REVIEW.md section 1, N4 (restored):
         Sigma(w_vec) = E - sum_a w_a Q_a at fixed total Q = sum Q_a, varying the
         partition. Discriminator table:
           Sigma ~ constant vs partition at fixed Q -> excess is geometric
           Sigma strongly partition-dependent       -> V(s) interference is a
                                                       leading eps channel (gate
                                                       for HC-3)
           Sigma < 0 anywhere                       -> sign convention or
                                                       non-minimizing branch
  HC-3 -- v86/GROUNDING.md section 1 (council-corrected GSS):
         D_ab = dQ_a/dw_b by finite differences on the flavored branch; record
         the FULL signature; orbital stability requires
             n(H_omega) = n(D)        [NEGATIVE index of the charge-response matrix]
         Single-charge reduction: dQ/dw < 0 -> n(D)=1 = n(H_omega) = 1 (VK).
         Caveat carried in the output: this run is UNGAUGED, which is exactly
         the regime where the GSS theorem is an anchor rather than a heuristic.
         n(H_omega) itself comes from HC-1 and is NOT computed here; HC-3
         delivers n(D) plus the stability MAP under the standing assumption
         n(H_omega) = 1 for the one-hump family (flagged, not assumed silently).

Physics (ungauged, 3 components, phase-blind product potential), following
v71/analysis/flavored_qball.py but with a corrected energy functional:
    Phi_a = f_a(r) e^{i w_a t},  s = prod_a f_a^2
    f_a'' + (2/r) f_a' = (m^2 - w_a^2) f_a + 2 Vt'(s) f_a prod_{b!=a} f_b^2
    Q_a = int w_a f_a^2 dV
    E   = int [ sum_a (1/2)(w_a^2 f_a^2 + f_a'^2 + m^2 f_a^2) + Vt(s) ] dV
(the symmetric limit f_a=f, w_a=w reproduces the v66/v69 convention exactly:
 Q = int 3 w f^2 dV, e = (3/2)(w^2f^2+f'^2+m^2f^2) + Vt(f^6) -- asserted below.)

Derived prediction carried over from N2 (this file also tests it):
    Sigma = E - sum_a w_a Q_a = E_grad + E_m + E_V - E_kin,  and Derrick gives
    R_vir = E_grad + 3(E_m + E_V - E_kin) = 0, hence  Sigma = (2/3) E_grad
for ANY partition. So N4's "is Sigma partition-dependent?" has a sharp
prediction: Sigma tracks E_grad exactly, whatever the flavor split does.
"""
import os
import sys
import itertools
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/d/code/scp"

M2, MU, KAP = 2.25, -41.345, 50.0
M = np.sqrt(M2)
H, RMAX = 0.01, 60.0
N = int(RMAX / H)
r = (np.arange(N) + 1) * H          # r = h .. RMAX ; ghost at r=0 enforces f'(0)=0
dV = 4.0 * np.pi * r * r * H


def Vt(s):
    return 0.5 * MU * s / (1.0 + KAP * s)


def Vp(s):
    return 0.5 * MU / (1.0 + KAP * s) ** 2


def Vpp(s):
    return -KAP * MU / (1.0 + KAP * s) ** 3


def _lap(f):
    fp = np.empty(N); fm = np.empty(N)
    fp[:-1] = f[1:]; fp[-1] = 0.0            # f(RMAX) = 0
    fm[1:] = f[:-1]; fm[0] = f[1]            # ghost f(-h) = f(h) -> f'(0) = 0
    return (fp - 2 * f + fm) / H ** 2 + (2.0 / r) * (fp - fm) / (2 * H)


def residual(F, w2):
    res = np.zeros_like(F)
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    vp = Vp(s)
    for a in range(3):
        prest = np.ones(N)
        for b in range(3):
            if b != a:
                prest *= F[b] ** 2
        res[a] = _lap(F[a]) - (M2 - w2[a]) * F[a] - 2 * vp * F[a] * prest
    return res


def jacobian(F, w2):
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    vp, vpp = Vp(s), Vpp(s)
    J = lil_matrix((3 * N, 3 * N))
    cm = 1.0 / H ** 2 - 1.0 / (r * H)
    cc = -2.0 / H ** 2
    cp = 1.0 / H ** 2 + 1.0 / (r * H)
    for a in range(3):
        prest = np.ones(N)
        for b in range(3):
            if b != a:
                prest *= F[b] ** 2
        # diagonal block: laplacian + d/df_a of the algebraic terms
        diag = cc - (M2 - w2[a]) - 2 * vp * prest \
            - 2 * vpp * F[a] * prest * (2 * F[a] * prest)
        idx = np.arange(N)
        J[a * N + idx, a * N + idx] = diag
        J[a * N + idx[1:], a * N + idx[1:] - 1] = cm[1:]
        J[a * N + idx[:-1], a * N + idx[:-1] + 1] = cp[:-1]
        J[a * N + 0, a * N + 1] += cm[0]                 # ghost folds back
        # off-diagonal blocks
        for b in range(3):
            if b == a:
                continue
            dprest = 2 * F[b]
            for c in range(3):
                if c != a and c != b:
                    dprest = dprest * F[c] ** 2
            ds_db = 2 * F[b] * (s / np.maximum(F[b] ** 2, 1e-300))
            off = -2 * (vpp * ds_db * F[a] * prest + vp * F[a] * dprest)
            J[a * N + idx, b * N + idx] = off
    return csr_matrix(J)


def solve(F, w, tol=1e-10, itmax=60):
    w2 = np.asarray(w) ** 2
    F = F.copy()
    for it in range(itmax):
        res = residual(F, w2)
        rn = np.max(np.abs(res))
        if rn < tol:
            return F, True, it, rn
        try:
            dF = spsolve(jacobian(F, w2), -res.reshape(-1)).reshape(3, N)
        except Exception:
            return F, False, it, rn
        lam = 1.0
        ok = False
        for _ in range(12):
            Ft = F + lam * dF
            if np.all(np.isfinite(Ft)) and np.max(np.abs(residual(Ft, w2))) < rn:
                F, ok = Ft, True
                break
            lam *= 0.5
        if not ok:
            return F, False, it, rn
    res = residual(F, np.asarray(w) ** 2)
    rn = np.max(np.abs(res))
    return F, rn < 1e-8, itmax, rn


def observables(F, w):
    """Energy/charge with the v66-consistent normalization (asserted in main)."""
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    Q = np.array([np.sum(w[a] * F[a] ** 2 * dV) for a in range(3)])
    E_kin = sum(0.5 * np.sum(w[a] ** 2 * F[a] ** 2 * dV) for a in range(3))
    E_grad = sum(0.5 * np.sum(np.gradient(F[a], H, edge_order=2) ** 2 * dV)
                 for a in range(3))
    E_mass = sum(0.5 * np.sum(M2 * F[a] ** 2 * dV) for a in range(3))
    E_V = np.sum(Vt(s) * dV)
    E = E_kin + E_grad + E_mass + E_V
    Sigma = E - float(np.dot(w, Q))
    R_vir = E_grad + 3.0 * (E_mass + E_V - E_kin)
    Sigma_pred = (2.0 / 3.0) * E_grad
    Qtot = float(Q.sum())
    return dict(w=tuple(w), Q=Q, Qtot=Qtot, E=E, E_kin=E_kin, E_grad=E_grad,
                E_mass=E_mass, E_V=E_V, Sigma=Sigma, Sigma_pred=Sigma_pred,
                eps=Sigma / float(np.dot(w, Q)),
                R_vir_rel=R_vir / E, I4_rel=(Sigma - Sigma_pred) / E,
                f0=tuple(F[a][0] for a in range(3)))


# ---------------------------------------------------------------- solving grid
def symmetric_seed(wb):
    src = np.loadtxt(os.path.join(ROOT, "v69/theory/gprofile_w142_g005.txt"))
    f0 = np.interp(r, src[:, 0], src[:, 1], right=0.0)
    F = np.vstack([f0, f0, f0])
    F, ok, it, rn = solve(F, np.array([wb, wb, wb]))
    if not ok:
        raise RuntimeError("symmetric seed failed at w=%.4f (resid %.2e)" % (wb, rn))
    return F


def walk_to(F, w_from, w_to, nstep=8):
    """Linear continuation in omega-space from w_from to w_to."""
    F = F.copy()
    w_from = np.asarray(w_from, float)
    w_to = np.asarray(w_to, float)
    for k in range(1, nstep + 1):
        w = w_from + (w_to - w_from) * (k / nstep)
        Fn, ok, it, rn = solve(F, w)
        if not ok:
            return None, None
        F = Fn
    return F, w_to


def main():
    wb = 1.42
    print("v86 N4 + HC-3 -- flavored (multi-clock) branch  [UNGAUGED, H=%g R=%g]"
          % (H, RMAX))

    # ---- normalization assert against the symmetric v66/v69 convention -----
    F = symmetric_seed(wb)
    o = observables(F, np.array([wb, wb, wb]))
    fsym = F[0]
    Q_v66 = np.sum(3.0 * wb * fsym ** 2 * dV)
    e_v66 = (1.5 * (wb ** 2 * fsym ** 2
                    + np.gradient(fsym, H, edge_order=2) ** 2
                    + M2 * fsym ** 2) + Vt(fsym ** 6))
    E_v66 = np.sum(e_v66 * dV)
    assert abs(o["Qtot"] - Q_v66) / Q_v66 < 1e-12, "Q normalization mismatch"
    assert abs(o["E"] - E_v66) / E_v66 < 1e-12, "E normalization mismatch"
    print("  normalization assert vs v66 symmetric convention: PASS "
          "(Q=%.4f, E=%.4f)" % (o["Qtot"], o["E"]))
    print("  symmetric point: eps=%.5f  R_vir/E=%.2e  I4/E=%.2e"
          % (o["eps"], o["R_vir_rel"], o["I4_rel"]))

    # ================= N4: partition scan at (nearly) fixed Q ================
    # Detune w_0 down and w_1,w_2 up so that sum_a w_a stays ~constant, which
    # keeps Q_tot close; then correct onto a Q target by a scalar shift in w.
    print("\n" + "=" * 78)
    print("N4 -- Legendre excess vs flavor partition")
    print("=" * 78)
    Qtarget = o["Qtot"]
    print("  Q target (symmetric point at w=%.3f): %.3f" % (wb, Qtarget))
    print("%-26s %9s %9s %9s %9s %9s %10s %10s %10s" %
          ("partition (w0,w1,w2)", "Qtot", "Q0/Qtot", "E", "Sigma", "eps",
           "2/3E_grad", "R_vir/E", "I4/E"))

    rows = []
    Fcur, wcur = F.copy(), np.array([wb, wb, wb])
    for dw in [0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050]:
        wt = np.array([wb - dw, wb + 0.5 * dw, wb + 0.5 * dw])
        Fn, wn = walk_to(Fcur, wcur, wt, nstep=max(2, int(dw / 0.005) + 2))
        if Fn is None:
            print("  dw=%.3f: NO CONVERGENCE -- flavored branch end" % dw)
            break
        # Newton-correct a uniform shift in w to hit Qtarget exactly
        for _ in range(12):
            oo = observables(Fn, wn)
            err = oo["Qtot"] - Qtarget
            if abs(err) / Qtarget < 1e-6:
                break
            h = 2e-4
            Fp, wp = walk_to(Fn, wn, wn + h, nstep=1)
            if Fp is None:
                break
            dQ = (observables(Fp, wp)["Qtot"] - oo["Qtot"]) / h
            shift = -err / dQ
            shift = float(np.clip(shift, -0.01, 0.01))
            Fs, ws = walk_to(Fn, wn, wn + shift, nstep=2)
            if Fs is None:
                break
            Fn, wn = Fs, ws
        oo = observables(Fn, wn)
        # fix 1 (grok Finding 10): the fixed-Q constraint must be reported, and
        # a row that failed to land on Q* is not comparable to the others.
        oo["Qerr"] = abs(oo["Qtot"] - Qtarget) / Qtarget
        if oo["Qerr"] > 1e-4:
            print("%-26s %9.3f  ROW REJECTED: |Q-Q*|/Q* = %.2e > 1e-4 "
                  "(fixed-Q control failed)"
                  % ("(%.4f,%.4f,%.4f)" % tuple(wn), oo["Qtot"], oo["Qerr"]))
            Fcur, wcur = Fn, wn
            continue
        rows.append(oo)
        print("%-26s %9.3f %9.4f %9.3f %9.4f %9.5f %10.4f %10.2e %10.2e  "
              "|dQ|/Q*=%.1e" %
              ("(%.4f,%.4f,%.4f)" % tuple(wn), oo["Qtot"],
               oo["Q"][0] / oo["Qtot"], oo["E"], oo["Sigma"], oo["eps"],
               oo["Sigma_pred"], oo["R_vir_rel"], oo["I4_rel"], oo["Qerr"]))
        Fcur, wcur = Fn, wn

    if len(rows) >= 2:
        S = np.array([d["Sigma"] for d in rows])
        Qr = np.array([d["Qtot"] for d in rows])
        frac = (S.max() - S.min()) / S.mean()
        print("\n  Sigma at fixed Q=%.2f (+/- %.1e): %.4f .. %.4f "
              "-> spread %.2f%% of the mean"
              % (Qr.mean(), (Qr.max() - Qr.min()) / Qr.mean(),
                 S.min(), S.max(), 100 * frac))
        maxI4 = max(abs(d["I4_rel"]) for d in rows)
        Eg = np.array([d["E_grad"] for d in rows])
        gfrac = (Eg.max() - Eg.min()) / Eg.mean()
        # fix 2 (grok Finding 10/13-2): if I4 still holds while Sigma moves, the
        # movement is carried by E_grad -- that is GEOMETRY (surface/profile
        # shape responding to the flavor split), not automatically V(s)
        # interference. Only a Sigma move that BREAKS the I4 identity would
        # implicate the potential coupling directly.
        if S.min() < 0:
            v = "SIGN ALARM -- Sigma < 0 somewhere; investigate before GSS"
        elif frac < 0.02:
            v = ("GEOMETRIC (flat) -- excess is a total-charge/geometry property; "
                 "spread %.2f%% < 2%%" % (100 * frac))
        elif maxI4 < 0.01 * frac:
            v = ("GEOMETRY-DRIVEN PARTITION DEPENDENCE (%.1f%%) -- Sigma moves, but "
                 "the I4 identity Sigma=(2/3)E_grad still holds to %.1e, and "
                 "E_grad moves by %.1f%%: the excess tracks the gradient/surface "
                 "term, so the flavor split acts through PROFILE SHAPE, not through "
                 "a separate V(s) interference channel. HC-3 gate: single-charge eps "
                 "formulas may not be reused numerically, but the closed form is "
                 "unchanged." % (100 * frac, maxI4, 100 * gfrac))
        else:
            v = ("PARTITION-DEPENDENT (%.1f%%) WITH I4 BREAKDOWN (%.1e) -- V(s) "
                 "interference is a leading eps channel; gate for HC-3: do NOT "
                 "reuse single-charge eps formulas naively" % (100 * frac, maxI4))
        print("  N4 VERDICT: %s" % v)
        print("  cross-check: max |Sigma - (2/3)E_grad|/E over the scan = %.2e"
              % maxI4)
        print("  E_grad spread over the scan: %.2f%% (Sigma spread %.2f%%)"
              % (100 * gfrac, 100 * frac))
        print("  PATH LIMITATION (pre-registered, grok Finding 13-1): this is ONE "
              "detuning ray (w0 down, w1=w2 up) at ONE total charge Q*=%.2f. It "
              "does not gate flavored physics in general -- only 'on this path'. A "
              "full (Q0,Q1,Q2) volume scan is HC-3's follow-on." % Qtarget)

    # ================= HC-3: charge-response matrix D_ab ====================
    print("\n" + "=" * 78)
    print("HC-3 -- corrected GSS: signature of D_ab = dQ_a/dw_b")
    print("   criterion  n(H_omega) = n(D)   [NEGATIVE indices must match]")
    print("   UNGAUGED run: theorem regime (the gauged case is heuristic).")
    print("   n(H_omega) is HC-1's deliverable; the single-hump family is")
    print("   ASSUMED n(H_omega)=1 here -- flagged, to be replaced by HC-1.")
    print("=" * 78)
    print("%-26s %9s %28s %8s %6s %s" %
          ("partition (w0,w1,w2)", "Qtot", "eigenvalues of dQ_a/dw_b",
           "n(D)", "p(D)", "GSS verdict (assuming n(H)=1)"))

    hw = 1e-3
    hc3 = []
    for base in rows:
        w0 = np.array(base["w"], float)
        # re-solve the base point (walk from the symmetric seed each time keeps
        # the finite differences honest and independent of scan history)
        Fb, wb2 = walk_to(F, np.array([wb, wb, wb]), w0,
                          nstep=max(2, int(np.max(np.abs(w0 - wb)) / 0.005) + 2))
        if Fb is None:
            continue
        D = np.zeros((3, 3))
        good = True
        for b in range(3):
            Qpm = []
            for sgn in (+1, -1):
                wpert = w0.copy()
                wpert[b] += sgn * hw
                Fp, wp = walk_to(Fb, w0, wpert, nstep=1)
                if Fp is None:
                    good = False
                    break
                Qpm.append(observables(Fp, wp)["Q"])
            if not good:
                break
            D[:, b] = (Qpm[0] - Qpm[1]) / (2 * hw)
        if not good:
            print("%-26s   finite differences failed" % ("(%.4f,%.4f,%.4f)" % tuple(w0)))
            continue
        # D = -Hess d is symmetric on a smooth branch (mixed partials of d);
        # symmetrizing removes finite-difference noise, and asym is reported so a
        # genuine asymmetry (= FD breakdown) is visible rather than hidden.
        Dsym = 0.5 * (D + D.T)
        ev = np.linalg.eigvalsh(Dsym)
        # fix 5 (grok Finding 12): count negatives with a RELATIVE threshold --
        # a near-zero eigenvalue must not flip n(D) on FD noise.
        tol = 1e-6 * max(1e-30, np.max(np.abs(ev)))
        nneg = int((ev < -tol).sum())
        npos = int((ev > tol).sum())
        nzero = 3 - nneg - npos
        asym = np.max(np.abs(D - D.T)) / max(1e-30, np.max(np.abs(D)))
        # fix 4: no unqualified "GSS-stable" -- n(H_omega) is HC-1's deliverable.
        verdict = ("n(D)=1 region (provisional n(H)=1)" if nneg == 1
                   else "n(D)=%d != provisional n(H)=1 -> HC-6 decay target" % nneg)
        if nzero:
            verdict += "  [%d near-zero eigenvalue(s)]" % nzero
        hc3.append(dict(w=tuple(w0), Q=base["Qtot"], ev=ev, nneg=nneg,
                        npos=npos, nzero=nzero, asym=asym,
                        evmin=float(np.min(np.abs(ev)))))
        print("%-26s %9.3f %28s %8d %6d %s" %
              ("(%.4f,%.4f,%.4f)" % tuple(w0), base["Qtot"],
               "[%9.2f %9.2f %9.2f]" % tuple(ev), nneg, npos, verdict))

    if hc3:
        masym = max(h["asym"] for h in hc3)
        print("\n  symmetry check max|D-D^T|/max|D| (should be ~FD noise): %.2e%s"
              % (masym, "   *** ASYMMETRY LARGE -- n(D) NOT TRUSTWORTHY ***"
                 if masym > 1e-2 else ""))
        print("  smallest |eigenvalue| seen: %.3e (near-zero modes can flip n(D))"
              % min(h["evmin"] for h in hc3))
        allneg = set(h["nneg"] for h in hc3)
        print("  n(D) values seen across the scanned partitions: %s" % sorted(allneg))
        if allneg == {1}:
            print("  HC-3 RESULT: every scanned flavored partition carries n(D)=1 --")
            print("  a contiguous n(D)=1 region, not a single VK point. Under the")
            print("  PROVISIONAL n(H_omega)=1 for this one-hump family the GSS")
            print("  indices match, i.e. the region is stability-CONSISTENT.")
        else:
            print("  HC-3 RESULT: n(D) CHANGES across the partition scan %s -- the")
            print("  index boundary lies inside the scanned region; partitions with")
            print("  n(D) != n(H_omega) are HC-6's decay targets." % sorted(allneg))
        print("\n  CLAIM LIMIT (grok Finding 12, mandatory): this is an n(D) MAP under")
        print("  a provisional n(H_omega)=1, NOT a full GSS verdict. Multi-charge")
        print("  ground states need not have n(H)=1 everywhere (the decoupled limit")
        print("  gives n(H)=3). Full GSS is deferred to HC-1's n(H_omega); this run")
        print("  does not close C1 and must not be quoted as doing so.")

    # dump
    path = os.path.join(HERE, "n4_flavored.tsv")
    with open(path, "w") as fh:
        fh.write("w0\tw1\tw2\tQ0\tQ1\tQ2\tQtot\tE\tSigma\teps\tSigma_pred\t"
                 "E_grad\tE_kin\tE_mass\tE_V\tR_vir_rel\tI4_rel\n")
        for d in rows:
            fh.write("\t".join("%.10g" % x for x in
                               list(d["w"]) + list(d["Q"]) +
                               [d["Qtot"], d["E"], d["Sigma"], d["eps"],
                                d["Sigma_pred"], d["E_grad"], d["E_kin"],
                                d["E_mass"], d["E_V"], d["R_vir_rel"],
                                d["I4_rel"]]) + "\n")
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
