#!/usr/bin/env python3
# v69 gauged radial shooter -- GAUGE_DESIGN.md section 6 "pre-kernel theory task".
#
# Static gauged Q-ball in the Coulomb/A_0 picture (GAUGE_DESIGN section 2 conventions:
# D_mu = d_mu + i g A_mu on charge-(+1) fields, E_i = d_i A_0 - d_t A_i, div E = +g rho_Q).
#
# Ansatz: Phi_a = f(r) e^{i omega t} (a=0,1,2), Theta = 0, A_0 = a0(r), A_i = 0.
#   D_t Phi = i(omega + g a0) f e^{i omega t}  =>  local frequency wt(r) = omega + g a0(r).
#   rho_Q = 3 wt f^2 ;   Gauss: laplacian(a0) = g rho_Q  (E_r = + d_r a0).
#
# We solve in the positive-definite variable  chi(r) := -g a0(r) >= 0  (chi = g * a0
# in the conventional-electrostatic-sign convention used in the task statement):
#
#   f''   + (2/r) f'   = (m^2 - (omega - chi)^2) f + mu f^5/(1+kappa f^6)^2
#   chi'' + (2/r) chi' = -3 g^2 (omega - chi) f^2
#   f'(0)=0, f(inf)=0, chi'(0)=0, chi(inf)=0   (chi ~ g^2 Q/(4 pi r))
#
# This is EXACTLY GAUGE_DESIGN section 6's system with wt = omega + g a0 = omega - chi;
# the task-statement form "omega - g*a0" corresponds to a0_task = chi/g = -a0_GD.
# G12 limit (chi -> 0): E_r = -chi'/g -> (3 g omega/r^2) int f^2 s^2 ds  (checked below).
#
# Gauge-invariant observables:
#   Q        = int 3 (omega-chi) f^2 dV
#   E_matter = int [ (3/2)((omega-chi)^2 f^2 + f'^2 + m^2 f^2) + Vt(f^6) ] dV
#   E_field  = int (1/2) E_r^2 dV = int chi'^2/(2 g^2) dV
#   a0(0)    = -chi(0)/g   (GAUGE_DESIGN sign; negative for Q>0)
#
# Method: damped Newton on the coupled finite-difference BVP (banded Jacobian,
# interleaved unknowns [f_0,chi_0,f_1,chi_1,...]), continuation in g and omega.
# g=0 must reproduce v66/results/scan.tsv (cross-checked below).
#
# Outputs: v69/theory/gscan.tsv, gprofile_w139_g005.txt, gprofile_w142_g005.txt,
#          stdout (redirect to gauged_shooter.out).

import os
import numpy as np
from scipy.linalg import solve_banded

ROOT = "/home/d/code/scp"
OUT  = os.path.join(ROOT, "v69/theory")

M2, MU, KAP = 2.25, -41.345, 50.0
M = np.sqrt(M2)
WMIN0 = np.sqrt(M2 + MU/9.0*(2.0/KAP)**(2.0/3.0))   # ungauged window bottom 1.308700

H    = 0.01
RMAX = 100.0
N    = int(round(RMAX/H))          # unknowns at nodes 0..N-1; node N is the boundary
r    = np.arange(N)*H
rN   = RMAX
RFAC = r[N-1]/rN                   # Coulomb-tail Robin factor: chi_N = chi_{N-1}*RFAC

H2 = H*H

def Wpot(f):
    f6 = f**6; d = 1.0 + KAP*f6
    return MU*f**5/d**2            # = 2 Vt'(f^6) f^5

def dWpot(f):
    f6 = f**6; d = 1.0 + KAP*f6
    return MU*f**4*(5.0 - 7.0*KAP*f6)/d**3

def Vt(s):
    return 0.5*MU*s/(1.0 + KAP*s)

# radial Laplacian coefficients (i>=1):  cm f_{i-1} + cc f_i + cp f_{i+1}
i_int = np.arange(1, N)
CM = 1.0/H2 - 1.0/(r[1:]*H)
CP = 1.0/H2 + 1.0/(r[1:]*H)
CC = -2.0/H2

def lap(y, ybnd):
    """Radial Laplacian of y (length N) with y(node N)=ybnd; regularity at r=0."""
    out = np.empty(N)
    out[0] = 6.0*(y[1]-y[0])/H2
    yp = np.empty(N-1); yp[:-1] = y[2:]; yp[-1] = ybnd
    out[1:] = CM*y[0:N-1] + CC*y[1:N] + CP*yp
    return out

def residual(f, chi, w, g):
    wt = w - chi
    Rf = lap(f, 0.0) - (M2 - wt*wt)*f - Wpot(f)
    Rc = lap(chi, chi[-1]*RFAC) + 3.0*g*g*wt*f*f
    # Robin tail for f? f decays exponentially; Dirichlet 0 at RMAX is exact to e^-60.
    return Rf, Rc

def newton(f, chi, w, g, maxit=60, tol=1e-9):
    nun = 2*N
    for it in range(maxit):
        Rf, Rc = residual(f, chi, w, g)
        R = np.empty(nun); R[0::2] = Rf; R[1::2] = Rc
        rn = np.max(np.abs(R))
        if rn < tol:
            return f, chi, True, it, rn
        wt = w - chi
        # banded Jacobian, (l,u)=(2,2): ab[2+i-j, j]
        ab = np.zeros((5, nun))
        idxf = 2*np.arange(N); idxc = idxf + 1
        # f rows
        Jff = np.empty(N)
        Jff[0]  = -6.0/H2
        Jff[1:] = CC
        Jff -= (M2 - wt*wt) + dWpot(f)
        ab[2, idxf] = Jff
        ab[1, idxc] = -2.0*wt*f                      # dRf_i/dchi_i
        ab[0, idxf[1]]   = 6.0/H2                    # dRf_0/df_1
        if N > 2:
            ab[0, idxf[2:]] = CP[:-1]                # dRf_i/df_{i+1}, i=1..N-2
        ab[4, idxf[:-1]] = CM                        # dRf_i/df_{i-1}, i=1..N-1
        # chi rows
        Jcc = np.empty(N)
        Jcc[0]  = -6.0/H2
        Jcc[1:] = CC
        Jcc -= 3.0*g*g*f*f
        Jcc[-1] += CP[-1]*RFAC                       # Robin tail closure
        ab[2, idxc] = Jcc
        ab[3, idxf] = 6.0*g*g*wt*f                   # dRc_i/df_i
        ab[0, idxc[1]]   = 6.0/H2
        if N > 2:
            ab[0, idxc[2:]] = CP[:-1]
        ab[4, idxc[:-1]] = CM
        try:
            dlt = solve_banded((2, 2), ab, -R)
        except Exception:
            return f, chi, False, it, rn
        df, dc = dlt[0::2], dlt[1::2]
        lam, ok = 1.0, False
        while lam > 1e-4:
            fn, cn = f + lam*df, chi + lam*dc
            Rfn, Rcn = residual(fn, cn, w, g)
            rnn = max(np.max(np.abs(Rfn)), np.max(np.abs(Rcn)))
            if rnn < (1.0 - 0.25*lam)*rn or rnn < tol:
                ok = True; break
            lam *= 0.5
        if not ok:
            return f, chi, False, it, rn
        f, chi = fn, cn
    Rf, Rc = residual(f, chi, w, g)
    rn = max(np.max(np.abs(Rf)), np.max(np.abs(Rc)))
    return f, chi, rn < 1e-7, maxit, rn

def solve(w, g, f0, chi0):
    f, chi, ok, its, rn = newton(f0.copy(), chi0.copy(), w, g)
    if ok and (f[0] < 0.2 or np.max(f) > 5.0 or np.any(np.isnan(f))):
        ok = False                                    # trivial / runaway
    return f, chi, ok, its, rn

def observables(f, chi, w, g):
    fe  = np.append(f, 0.0)
    ce  = np.append(chi, chi[-1]*RFAC)
    re  = np.append(r, rN)
    wte = w - ce
    fpe = np.gradient(fe, H, edge_order=2)
    cpe = np.gradient(ce, H, edge_order=2)
    dV  = 4.0*np.pi*re*re
    Q   = np.trapz(3.0*wte*fe*fe*dV, dx=H)
    em  = 1.5*(wte*wte*fe*fe + fpe*fpe + M2*fe*fe) + Vt(fe**6)
    Em  = np.trapz(em*dV, dx=H)
    Ef  = np.trapz(cpe*cpe*dV, dx=H)/(2.0*g*g) if g > 0 else 0.0
    # Coulomb-tail remainder of E_field beyond RMAX: int_R^inf (C/r^2)^2/2 * 4 pi r^2
    if g > 0:
        Cc = ce[-1]*rN                 # chi ~ Cc/r
        Ef += 4.0*np.pi*Cc*Cc/(2.0*g*g*rN)
    h0 = 0.5*f[0]
    above = np.where(fe >= h0)[0]
    rhalf = 0.0
    if len(above):
        i = above[-1]
        if i < N:
            rhalf = re[i] + H*(fe[i]-h0)/(fe[i]-fe[i+1])
    return dict(w=w, g=g, f0=f[0], chi0=chi[0],
                a00=(-chi[0]/g if g > 0 else 0.0),
                Q=Q, Em=Em, Ef=Ef, Et=Em+Ef,
                EomQ=(Em+Ef)/(M*Q), rhalf=rhalf, weff0=w-chi[0])

def load_v66_profile(path):
    rr, ff = np.loadtxt(path, comments="#", unpack=True)
    f = np.interp(r, rr, ff, right=0.0)
    return f

# ---------------------------------------------------------------- sweeps
def sweep(g, w_start, sol_start, w_stop, step0, store):
    """Continue from (w_start, sol_start) toward w_stop in steps of step0 (signed).
    Adaptive halving on failure; returns (last_good_w, last_good_sol, terminated)."""
    w_cur = w_start
    f, chi = sol_start
    step = step0
    terminated = False
    while True:
        if (step0 > 0 and w_cur >= w_stop - 1e-12) or (step0 < 0 and w_cur <= w_stop + 1e-12):
            break
        w_next = w_cur + step
        if (step0 > 0 and w_next > w_stop): w_next = w_stop
        if (step0 < 0 and w_next < w_stop): w_next = w_stop
        fn, cn, ok, its, rn = solve(w_next, g, f, chi)
        if ok:
            w_cur, f, chi = w_next, fn, cn
            store[w_cur] = (f.copy(), chi.copy())
            step = step0
        else:
            step *= 0.5
            if abs(step) < 4e-6:
                terminated = True
                break
    return w_cur, (f, chi), terminated

def main():
    print("v69 gauged radial shooter  (h=%g, Rmax=%g, N=%d nodes, %d unknowns)"
          % (H, RMAX, N, 2*N))
    print("conventions: D=d+igA, E_i=d_iA_0-d_tA_i, divE=+g rho_Q;  wt(r)=omega+g a0(r)=omega-chi(r)")
    print("ungauged window: omega in (%.6f, %.6f)\n" % (WMIN0, M))

    gs = [0.0, 0.02, 0.05, 0.1]
    branches = {}          # g -> {w: (f,chi)}
    terminfo = {}          # g -> (w_min_reached, terminated?, w_max_reached)

    g_prev = None
    for g in gs:
        # ---- obtain a seed solution at some omega for this g
        seed = None
        if g == 0.0:
            wseed = 1.45
            f0g = load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
            f, chi, ok, its, rn = solve(wseed, 0.0, f0g, np.zeros(N))
            if ok: seed = (wseed, f, chi)
        else:
            for wtry in (1.45, 1.46, 1.47, 1.48, 1.49):
                wseed = min(branches[g_prev], key=lambda w: abs(w - wtry))
                fp_, cp_ = branches[g_prev][wseed]
                # adaptive g-continuation from g_prev to g
                gcur, chig = g_prev, g_prev
                f, chi = fp_.copy(), cp_.copy()
                dg = g - g_prev
                failed = False
                while gcur < g - 1e-12:
                    gg = min(gcur + dg, g)
                    cguess = chi*(gg/chig)**2 if chig > 0 else chi
                    fn, cn, ok, its, rn = solve(wseed, gg, f, cguess)
                    if ok:
                        f, chi, gcur, chig = fn, cn, gg, gg
                        dg = min(dg*2.0, g - gcur) if g - gcur > 0 else dg
                    else:
                        dg *= 0.5
                        if dg < 0.0015: failed = True; break
                if not failed:
                    seed = (wseed, f, chi); break
                print("  (g=%.3f: seed at w=%.4f failed, trying higher omega)" % (g, wseed))
        if seed is None:
            print("FATAL: could not seed branch at g=%.3f" % g); return
        wseed, fS, cS = seed
        store = {wseed: (fS.copy(), cS.copy())}
        wtop, _, term_up = sweep(g, wseed, (fS, cS), 1.4975, +0.0025, store)
        floor = 1.3125 if g == 0.0 else 1.3090
        wbot, _, term_dn = sweep(g, wseed, (fS, cS), floor, -0.0025, store)
        branches[g] = store
        terminfo[g] = (wbot, term_dn, wtop, term_up)
        g_prev = g
        print("g=%.3f: branch solved at %d omegas in [%.6f, %.6f]  (low-end %s, top %s)"
              % (g, len(store), wbot, wtop,
                 "TERMINATED (no solution below)" if term_dn else "floor reached (not terminated)",
                 "terminated" if term_up else "cap 1.4975 reached"))

    # ---------------- tables
    rows = []
    for g in gs:
        for w in sorted(branches[g]):
            f, chi = branches[g][w]
            rows.append(observables(f, chi, w, g))
    # dQ/domega sign per (g) branch
    for g in gs:
        sub = [o for o in rows if o["g"] == g]
        sub.sort(key=lambda o: o["w"])
        for k, o in enumerate(sub):
            if k == 0:               dQ = sub[1]["Q"] - sub[0]["Q"]
            elif k == len(sub)-1:    dQ = sub[-1]["Q"] - sub[-2]["Q"]
            else:                    dQ = sub[k+1]["Q"] - sub[k-1]["Q"]
            o["sgn"] = -1 if dQ < 0 else 1

    tsv = os.path.join(OUT, "gscan.tsv")
    with open(tsv, "w") as fp:
        fp.write("g\tomega\tf0\tchi0\ta0_0\tQ\tE_matter\tE_field\tE_total\tE_over_mQ\tdQdw_sign\tr_half\tweff0\n")
        for o in rows:
            fp.write("%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.3f\t%.2f\t%.5f\t%d\t%.3f\t%.6f\n"
                     % (o["g"], o["w"], o["f0"], o["chi0"], o["a00"], o["Q"],
                        o["Em"], o["Ef"], o["Et"], o["EomQ"], o["sgn"], o["rhalf"], o["weff0"]))
    print("\nwrote %s (%d rows)" % (tsv, len(rows)))

    # ---------------- g=0 verification vs v66 scan.tsv
    print("\n=== g=0 verification vs v66/results/scan.tsv ===")
    scan = {}
    with open(os.path.join(ROOT, "v66/results/scan.tsv")) as fp:
        next(fp)
        for line in fp:
            c = line.split()
            scan[round(float(c[0]), 4)] = (float(c[1]), float(c[2]), float(c[3]))  # f0,E,Q
    worst = 0.0; worst_what = ""
    print("omega    f0(here)  f0(v66)   dE%%      dQ%%      df0%%")
    for o in [x for x in rows if x["g"] == 0.0]:
        key = round(o["w"], 4)
        if key in scan:
            f0v, Ev, Qv = scan[key]
            d0 = abs(o["f0"]-f0v)/f0v*100
            dE = abs(o["Em"]-Ev)/Ev*100
            dQ = abs(o["Q"]-Qv)/Qv*100
            for d, what in ((d0, "f0"), (dE, "E"), (dQ, "Q")):
                if d > worst: worst, worst_what = d, "%s@w=%.4f" % (what, key)
            if key in (1.315, 1.33, 1.36, 1.39, 1.42, 1.45, 1.47, 1.49):
                print("%.4f   %.6f  %.6f  %+.4f  %+.4f  %+.4f" % (key, o["f0"], f0v, dE, dQ, d0))
    print("worst relative deviation over all matched rows: %.4f%%  (%s)" % (worst, worst_what))

    # ---------------- dE/dw = w dQ/dw identity
    print("\n=== gauged Q-ball identity dE/dw = w dQ/dw (centered FD on branch) ===")
    for g in gs:
        sub = sorted([x for x in rows if x["g"] == g], key=lambda o: o["w"])
        devs = []
        for k in range(1, len(sub)-1):
            dE = sub[k+1]["Et"] - sub[k-1]["Et"]
            dQ = sub[k+1]["Q"]  - sub[k-1]["Q"]
            if abs(dQ) > 1e-6:
                devs.append(abs(dE/dQ - sub[k]["w"])/sub[k]["w"])
        print("g=%.3f: median |dE/dQ - w|/w = %.2e, max = %.2e" %
              (g, np.median(devs), np.max(devs)))

    # ---------------- headline numbers
    print("\n=== window erosion / Q_max / delta-omega_self ===")
    head = {}
    for g in gs:
        sub = sorted([x for x in rows if x["g"] == g], key=lambda o: o["w"])
        wmin, term, wmax, _ = terminfo[g]
        Qmax = max(o["Q"] for o in sub)
        oQmax = max(sub, key=lambda o: o["Q"])
        # E=mQ crossing (absolute stability edge)
        wabs = None
        for k in range(len(sub)-1):
            a, b = sub[k], sub[k+1]
            if (a["EomQ"]-1.0)*(b["EomQ"]-1.0) < 0:
                t = (1.0-a["EomQ"])/(b["EomQ"]-a["EomQ"])
                wabs = a["w"] + t*(b["w"]-a["w"])
        # omega at fixed Q=482 (interpolate on the dQ/dw<0 part)
        wq = None
        mono = [o for o in sub if o["sgn"] < 0]
        for k in range(len(mono)-1):
            a, b = mono[k], mono[k+1]
            if (a["Q"]-482.0)*(b["Q"]-482.0) <= 0:
                t = (482.0-a["Q"])/(b["Q"]-a["Q"])
                wq = a["w"] + t*(b["w"]-a["w"])
        head[g] = dict(wmin=wmin, term=term, Qmax=Qmax, wQmax=oQmax["w"],
                       weff0_Qmax=oQmax["weff0"], wabs=wabs, wq=wq,
                       chi0_Qmax=oQmax["chi0"])
        print("g=%.3f: omega_min=%.5f%s  Q_max=%.0f (at w=%.4f, weff0=%.4f)  "
              "omega_abs(E=mQ)=%s  omega(Q=482)=%s" %
              (g, wmin, " [TERMINATED]" if term else " [floor]",
               Qmax, oQmax["w"], oQmax["weff0"],
               "%.4f" % wabs if wabs else "none-in-scan",
               "%.5f" % wq if wq else "Q=482 NOT on branch"))
    if head[0.0]["wq"]:
        print("\ndelta-omega_self at fixed Q=482 (vs g=0 omega=%.5f):" % head[0.0]["wq"])
        for g in gs[1:]:
            if head[g]["wq"]:
                d = head[g]["wq"] - head[0.0]["wq"]
                print("  g=%.3f: delta_w = %+.5f  (= %+.1f g^2)   [design est. +13.1 g^2 = %+.4f]"
                      % (g, d, d/(g*g), 13.1*g*g))
            else:
                print("  g=%.3f: Q=482 not on branch (Q_max=%.0f)" % (g, head[g]["Qmax"]))

    # ---------------- G12 limit + seed profiles
    print("\n=== kernel seed profiles (temporal gauge) ===")
    for wtar, tag in ((1.39, "w139"), (1.42, "w142")):
        g = 0.05
        wkey = min(branches[g], key=lambda w: abs(w - wtar))
        f, chi = branches[g][wkey]
        o = observables(f, chi, wkey, g)
        wt = wkey - chi
        # Gauss integral of the SOLVED rho_Q = 3 wt f^2:
        integ = np.concatenate(([0.0], np.cumsum(0.5*H*(3*wt*f*f*r*r)[:-1] + 0.5*H*(3*wt*f*f*r*r)[1:])))
        Er = np.zeros(N); Er[1:] = g*integ[1:]/(r[1:]**2)
        # cross-check vs -chi'/g
        chip = np.gradient(np.append(chi, chi[-1]*RFAC), H, edge_order=2)[:N]
        xdev = np.max(np.abs(Er + chip/g))
        # G12 (bare-omega) comparison
        integ0 = np.concatenate(([0.0], np.cumsum(0.5*H*(3*wkey*f*f*r*r)[:-1] + 0.5*H*(3*wkey*f*f*r*r)[1:])))
        Er0 = np.zeros(N); Er0[1:] = g*integ0[1:]/(r[1:]**2)
        print("w=%.4f g=%.2f: f0=%.6f Q=%.2f E_m=%.2f E_f=%.3f E_t=%.2f chi0=%.5f a0(0)=%.5f "
              "weff0=%.5f r_half=%.3f" %
              (wkey, g, o["f0"], o["Q"], o["Em"], o["Ef"], o["Et"], o["chi0"], o["a00"],
               o["weff0"], o["rhalf"]))
        print("   max|Er(Gauss of solved rho) - (-chi'/g)| = %.2e   "
              "max|Er - Er_bare-omega(G12)| = %.2e (the O(g^2) chi correction)" % (xdev, np.max(np.abs(Er-Er0))))
        path = os.path.join(OUT, "gprofile_%s_g005.txt" % tag)
        step = int(round(0.02/H)); nmax = int(round(60.0/0.02))
        with open(path, "w") as fp:
            fp.write("# v69 gauged_shooter profile: omega=%.6f g=%.6f m2=%.6f mu=%.6f kappa=%.6f\n"
                     % (wkey, g, M2, MU, KAP))
            fp.write("# f0=%.6f Q=%.2f E_matter=%.2f E_field=%.3f E_total=%.2f chi0=%.6f a0_GD(0)=%.6f weff0=%.6f\n"
                     % (o["f0"], o["Q"], o["Em"], o["Ef"], o["Et"], o["chi0"], o["a00"], o["weff0"]))
            fp.write("# TEMPORAL-GAUGE SEED (exact gauge transform of this static solution):\n")
            fp.write("#   u_a=f(r), v_a=0, udot_a=0, vdot_a=weff(r)*f(r), theta links=0, E_i=Er*rhat_i\n")
            fp.write("#   weff(r)=omega-chi(r)=omega+g*a0_GD(r); rho_Q(t=0)=3*weff*f^2 satisfies divE=g*rho_Q EXACTLY\n")
            fp.write("# r f Er weff\n")
            for k in range(nmax+1):
                idx = k*step
                fp.write("%.6f %.9f %.9e %.9f\n" % (r[idx], max(f[idx], 0.0), Er[idx], wt[idx]))
        print("   wrote %s (rmax=60, dr=0.02, columns: r f Er weff)" % path)

    # ---------------- positronium refinement
    print("\n=== positronium (opposite-charge pair) refinement, g=0.05 ===")
    g = 0.05
    for wtar in (1.39,):
        wkey = min(branches[g], key=lambda w: abs(w - wtar))
        f, chi = branches[g][wkey]
        o = observables(f, chi, wkey, g)
        Mb, Qb = o["Et"], o["Q"]
        D1 = 8.0*np.pi*Mb/(wkey**2 * g**2 * Qb**2)
        vrel = 2.0/(wkey*D1)
        Torb = np.pi*wkey*D1*D1
        print("solved ball at bare w=%.4f: M=E_total=%.2f, Q=%.2f" % (wkey, Mb, Qb))
        print("  D_1 = 8 pi M/(w^2 g^2 Q^2) = %.2f   v_rel = 2/(w D_1) = %.4f   T_orb = pi w D_1^2 = %.0f"
              % (D1, vrel, Torb))
        print("  [design estimate used M=692, Q=482 ungauged: D_1=15.5, T_orb=1050]")
    if head[g]["wq"]:
        wq = head[g]["wq"]
        wkey = min(branches[g], key=lambda w: abs(w - wq))
        f, chi = branches[g][wkey]
        o = observables(f, chi, wkey, g)
        D1 = 8.0*np.pi*o["Et"]/(wkey**2 * g**2 * o["Q"]**2)
        print("fixed-Q variant (w=%.4f, Q=%.1f, M=%.1f): D_1=%.2f, v_rel=%.4f, T_orb=%.0f"
              % (wkey, o["Q"], o["Et"], D1, 2.0/(wkey*D1), np.pi*wkey*D1*D1))

    # ---------------- discretization sanity at w=1.39, g=0
    print("\ndone.")

if __name__ == "__main__":
    main()
