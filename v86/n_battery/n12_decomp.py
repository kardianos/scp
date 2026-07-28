#!/usr/bin/env python3
"""v86 Part-0 rungs N1 (component decomposition of E/Sigma) and N2 (virial identity).

Protocol: v86/council/grok45/N_BATTERY_REVIEW.md section 1, N1 + N2.
Reviewed by the grok-4.5 implementation seat (verdict RUN WITH FIXES; all seven
required fixes applied here -- see CONVENTION banner, Phase A, per-row asserts,
Newton-residual gate, I2-as-duplicate marking, N12_FAST=0 default, outcome hooks).

--------------------------------------------------------------------------
DERIVATION (what this script tests; results in v86/n_battery/N1_N2_RESULTS.md)
--------------------------------------------------------------------------
Shooter conventions (v69/theory/gauged_shooter.py):
    Phi_a = f(r) e^{i w t}, a=0,1,2 ; chi(r) = -g a_0(r) ; wt(r) = w - chi(r)
    f''  + (2/r) f'  = (m^2 - wt^2) f + W(f),   W(f) = mu f^5/(1+kappa f^6)^2
    chi''+ (2/r) chi'= -3 g^2 wt f^2
Energy pieces (dV = 4 pi r^2 dr):
    E_kin = int (3/2) wt^2 f^2 dV        E_grad = int (3/2) f'^2   dV
    E_m   = int (3/2) m^2 f^2 dV         E_V    = int Vt(f^6)      dV
    E_g   = int chi'^2/(2 g^2) dV        E      = E_kin+E_grad+E_m+E_V+E_g
    Q     = int 3 wt f^2 dV

I1 (Gauss identity).  rho_Q = 3 wt f^2 = -lap(chi)/g^2, so
        int chi rho_Q dV = int |grad chi|^2 dV/g^2 = 2 E_g,
    and w Q = int (wt + chi) rho_Q dV = 2 E_kin + 2 E_g.        <-- exact identity
    (Truncated-domain IBP leaves a surface term at RMAX equal to 2*E_g_tail for
     chi ~ C/r; the same Coulomb tail correction that observables() adds to E_g
     is therefore exactly what I1 needs.)

I2 (Sigma = the reduced functional).  The free-extremum functional reproducing
    BOTH radial ODEs (gauge energy enters with a MINUS: electrostatics is a saddle
    point in a_0) is
        F[f,chi] = E_grad + E_m + E_V - E_kin - E_g.
    Using I1,  Sigma = E - w Q = F.
    NOTE: on the Gauss constraint I2 is algebraically -I1, i.e. a rewrite check,
    NOT an independent physical constraint. Reported as a duplicate for
    transparency only (grok Finding 3/13).

I3 (Derrick/virial).  F is stationary under the free variation
    f -> f(r/L), chi -> chi(r/L) at fixed w.  In 3D E_grad,E_g ~ L and
    E_m,E_V,E_kin ~ L^3, so dF/dL|_1 = 0 gives
        R_vir = (E_grad - E_g) + 3 (E_m + E_V - E_kin) = 0.      <-- N2 identity
    (The scaled pair need not obey Gauss; Derrick only needs a variation of F.)

I4 (closed form for the excess).  I2 + I3 eliminate the bulk terms:
        Sigma = (2/3) (E_grad - E_g),
        eps   = Sigma/(w Q) = (E_grad - E_g) / (3 (E_kin + E_g)).
    Ungauged: Sigma = (2/3) E_grad > 0 -- derives the measured positivity of eps.

--------------------------------------------------------------------------
CONVENTION BANNER (grok Findings 6,7,18 -- read before comparing to any doc)
--------------------------------------------------------------------------
  * E_kin here is the TRUE integrated kinetic energy int (3/2) wt^2 f^2 dV.
    The "throughput" quantity E_kin^tp = w Q / 2 is a DIFFERENT number when
    gauged: by I1, w Q / 2 = E_kin + E_g.
  * Correct gauged excess:  Sigma = E_grad + E_m + E_V - E_kin - E_g
                                  = E_grad + E_m + E_V - E_kin^tp.
  * FOUNDATIONS R2 / GROUNDING §0 / THEORY_v86 A4 write the excess as
    "gradient + potential + GAUGE - kinetic". That +E_g is WRONG once gauged
    (it is off by 2*E_g if combined with the true E_kin). This instrument's
    derivation supersedes that prose; the docs must be corrected.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/d/code/scp"
sys.path.insert(0, os.path.join(ROOT, "v69/theory"))

# fix 6: production grid is the default; N12_FAST=1 is an explicit scout opt-in.
FAST = os.environ.get("N12_FAST", "0") == "1"
if FAST:
    import gauged_shooter_fast as G
else:
    import gauged_shooter as G

M2, MU, KAP, M = G.M2, G.MU, G.KAP, G.M
H, N, rN, RFAC = G.H, G.N, G.rN, G.RFAC
r = G.r

GSCAN = os.path.join(ROOT, "v69/theory/gscan.tsv")


def decompose(f, chi, w, g, newton_res=np.nan):
    """Full energy split on a converged shooter profile. Mirrors G.observables()
    exactly for Q, Em, Ef, rhalf, and asserts agreement (fix 2)."""
    fe = np.append(f, 0.0)
    ce = np.append(chi, chi[-1] * RFAC)
    re = np.append(r, rN)
    wte = w - ce
    fpe = np.gradient(fe, H, edge_order=2)
    cpe = np.gradient(ce, H, edge_order=2)
    dV = 4.0 * np.pi * re * re

    Q = np.trapz(3.0 * wte * fe * fe * dV, dx=H)
    E_kin = np.trapz(1.5 * wte * wte * fe * fe * dV, dx=H)
    E_grad = np.trapz(1.5 * fpe * fpe * dV, dx=H)
    E_mass = np.trapz(1.5 * M2 * fe * fe * dV, dx=H)
    E_V = np.trapz(G.Vt(fe ** 6) * dV, dx=H)
    if g > 0:
        E_g = np.trapz(cpe * cpe * dV, dx=H) / (2.0 * g * g)
        Cc = ce[-1] * rN                      # chi ~ Cc/r tail beyond RMAX
        E_g += 4.0 * np.pi * Cc * Cc / (2.0 * g * g * rN)
    else:
        E_g = 0.0

    E_matter = E_kin + E_grad + E_mass + E_V
    E = E_matter + E_g
    Sigma = E - w * Q
    eps = Sigma / (w * Q)

    # --- fix 2: hard cross-check against the frozen shooter observables --------
    o = G.observables(f, chi, w, g)
    dE = abs(E - o["Et"]) / abs(o["Et"])
    dQ = abs(Q - o["Q"]) / abs(o["Q"])
    dEm = abs(E_matter - o["Em"]) / abs(o["Em"])
    if not (dE < 1e-10 and dQ < 1e-10 and dEm < 1e-10):
        raise AssertionError(
            "N1 KILL: energy pieces do not rebuild observables() at (g=%g,w=%g): "
            "dE=%.3e dQ=%.3e dEm=%.3e" % (g, w, dE, dQ, dEm))

    # --- identity residuals ---------------------------------------------------
    I1 = w * Q - 2.0 * (E_kin + E_g)                       # Gauss identity
    I2 = Sigma - (E_grad + E_mass + E_V - E_kin - E_g)     # duplicate of -I1
    R_vir = (E_grad - E_g) + 3.0 * (E_mass + E_V - E_kin)  # N2 Derrick identity
    Sigma_pred = (2.0 / 3.0) * (E_grad - E_g)              # N2 closed form
    I4 = Sigma - Sigma_pred

    h0 = 0.5 * f[0]
    above = np.where(fe >= h0)[0]
    rhalf = 0.0
    if len(above):
        i = above[-1]
        if i < N:
            rhalf = re[i] + H * (fe[i] - h0) / (fe[i] - fe[i + 1])
    xi = 1.0 / np.sqrt(max(M2 - w * w, 1e-30))

    # --- fix 7: N1 outcome-table hooks ---------------------------------------
    # surface hypothesis: Sigma ~ 4 pi R^2 sigma_eff  -> sigma_eff should be flat
    sigma_eff = Sigma / (4.0 * np.pi * rhalf ** 2) if rhalf > 0 else np.nan
    # point-Coulomb estimate of the gauge self-energy
    Eg_point = g * g * Q * Q / (8.0 * np.pi * rhalf) if (g > 0 and rhalf > 0) else 0.0
    coul_par = g * g * Q / rhalf if rhalf > 0 else np.nan   # N3 second axis

    return dict(g=g, w=w, Q=Q, E=E, E_kin=E_kin, E_grad=E_grad, E_mass=E_mass,
                E_V=E_V, E_g=E_g, E_matter=E_matter, Sigma=Sigma, eps=eps,
                Sigma_pred=Sigma_pred, eps_pred=Sigma_pred / (w * Q),
                I1=I1, I2=I2, R_vir=R_vir, I4=I4,
                I1_rel=I1 / E, I2_rel=I2 / E, R_vir_rel=R_vir / E, I4_rel=I4 / E,
                rhalf=rhalf, xi=xi, xi_over_R=(xi / rhalf if rhalf > 0 else np.nan),
                sigma_eff=sigma_eff, Eg_point=Eg_point, coul_par=coul_par,
                f0=f[0], chi0=chi[0], newton_res=newton_res)


def _sweep_w(g, w0, f0, chi0, wlo, whi, dw0, store, dw_min=2.5e-5):
    """March in omega at fixed g in both directions from w0 with adaptive step
    halving (mirrors the v69 shooter's sweep()), storing (f, chi, newton_res)
    on the requested dw0 lattice only."""
    for direction in (+1, -1):
        f, chi, w = f0.copy(), chi0.copy(), w0
        dw = dw0
        w_last_stored = w0
        while True:
            wn = round(w + direction * dw, 8)
            if wn < wlo - 1e-12 or wn > whi + 1e-12:
                break
            fn, cn, ok, _, rn = G.solve(wn, g, f, chi)
            if ok:
                f, chi, w = fn, cn, wn
                # thin the stored set to ~dw0 spacing (adaptive steps can be finer)
                if abs(wn - w_last_stored) >= 0.999 * dw0 or abs(dw) < dw0:
                    store[round(wn, 8)] = (fn.copy(), cn.copy(), rn)
                    w_last_stored = wn
                dw = dw0
            else:
                dw *= 0.5
                if abs(dw) < dw_min:
                    break


def _continue_in_g(w, f, chi, g_from, g_to):
    """Adaptive continuation in the gauge coupling at fixed omega."""
    f, chi = f.copy(), chi.copy()
    gc, chig, rn = g_from, g_from, np.nan
    dg = g_to - gc
    while gc < g_to - 1e-12:
        gg = min(gc + dg, g_to)
        cguess = chi * (gg / chig) ** 2 if chig > 0 else chi
        fn, cn, ok, _, rnn = G.solve(w, gg, f, cguess)
        if ok:
            f, chi, gc, chig, rn = fn, cn, gg, gg, rnn
            dg = g_to - gc
        else:
            dg *= 0.5
            if dg < 1e-6:
                return None
    return f, chi, rn


def build_branches(gs, wlo, whi, dw):
    """v69 continuation strategy: seed at w=1.45, g=0 from the v66 profile; sweep
    omega at g=0; then for each successive g continue in the coupling at the first
    anchor omega that works (window rises with g) and march omega at fixed g."""
    f0g = G.load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
    f, chi, ok, _, rn = G.solve(1.45, 0.0, f0g, np.zeros(N))
    if not ok:
        raise RuntimeError("seed solve failed")

    out = {}
    store0 = {1.45: (f.copy(), chi.copy(), rn)}
    _sweep_w(0.0, 1.45, f, chi, wlo, whi, dw, store0)
    out[0.0] = store0
    print("  branch g=0.000: %d omega points  [%.4f .. %.4f]"
          % (len(store0), min(store0), max(store0)), flush=True)

    g_prev = 0.0
    for g in [x for x in gs if x > 0.0]:
        prev = out[g_prev]
        anchor = None
        # window bottom rises with g: try progressively higher omega anchors
        for wtry in sorted(prev, reverse=True):
            if wtry < 1.42:
                break
            fp, cp, _ = prev[wtry]
            res = _continue_in_g(wtry, fp, cp, g_prev, g)
            if res is not None:
                anchor = (wtry, res)
                break
        if anchor is None:
            print("  branch g=%.3f: CONTINUATION FAILED (no viable anchor)" % g,
                  flush=True)
            continue
        wa, (fa, ca, ra) = anchor
        store = {wa: (fa.copy(), ca.copy(), ra)}
        _sweep_w(g, wa, fa, ca, wlo, whi, dw, store)
        out[g] = store
        print("  branch g=%.3f: %d omega points  [%.4f .. %.4f]  (anchor w=%.4f)"
              % (g, len(store), min(store), max(store), wa), flush=True)
        g_prev = g
    return out


def phase_a(rows):
    """Fix 3 -- N1 Phase A: cross-check against the frozen v69 gscan.tsv artifact
    and report the eps(g=0) vs eps(g>0) gap attributed to E_g."""
    print("\n" + "=" * 78)
    print("N1 PHASE A -- cross-check vs frozen v69/theory/gscan.tsv")
    print("=" * 78)
    if not os.path.exists(GSCAN):
        print("  gscan.tsv NOT FOUND -- Phase A skipped")
        return
    dat = np.genfromtxt(GSCAN, names=True)
    print("%-6s %-8s %12s %12s %12s %12s %12s" %
          ("g", "w", "Q(gscan)", "dQ/Q", "E(gscan)", "dE/E", "dEfield/E"))
    nbad = 0
    for g in sorted(set(d["g"] for d in rows)):
        sub = [d for d in rows if d["g"] == g]
        sel = dat[np.abs(dat["g"] - g) < 1e-9]
        if len(sel) == 0:
            continue
        for d in sub[::max(1, len(sub) // 5)]:
            j = int(np.argmin(np.abs(sel["omega"] - d["w"])))
            if abs(sel["omega"][j] - d["w"]) > 1e-4:
                continue
            dQ = (d["Q"] - sel["Q"][j]) / sel["Q"][j]
            dE = (d["E"] - sel["E_total"][j]) / sel["E_total"][j]
            dF = ((d["E_g"] - sel["E_field"][j]) / d["E"]) if g > 0 else 0.0
            print("%-6.3f %-8.4f %12.4f %12.2e %12.4f %12.2e %12.2e" %
                  (g, d["w"], sel["Q"][j], dQ, sel["E_total"][j], dE, dF))
            if abs(dQ) > 2e-3 or abs(dE) > 2e-3:
                nbad += 1
    print("  Phase A verdict: %s (%d rows off by >0.2%%; grid differences of "
          "O(1e-4) are expected between H=%.3f/RMAX=%g and gscan's H=0.004/RMAX=150)"
          % ("CONSISTENT" if nbad == 0 else "DISCREPANT", nbad, H, rN))

    # eps(g=0) vs eps(g>0) gap attributed to E_g
    print("\n  eps gap vs E_g (is Coulomb the g-piece?):")
    print("  %-8s %10s %10s %12s %12s %12s" %
          ("w", "eps(g=0)", "eps(g)", "d_eps", "-2Eg/3wQ", "ratio"))
    base = {d["w"]: d for d in rows if d["g"] == 0.0}
    for g in sorted(set(d["g"] for d in rows)):
        if g == 0.0:
            continue
        sub = [d for d in rows if d["g"] == g]
        for d in sub[::max(1, len(sub) // 4)]:
            b = base.get(d["w"])
            if b is None:
                continue
            deps = d["eps"] - b["eps"]
            pred = -(2.0 / 3.0) * d["E_g"] / (d["w"] * d["Q"])
            print("  g=%-6.3f w=%-6.3f %10.5f %10.5f %12.5f %12.5f %12.3f" %
                  (g, d["w"], b["eps"], d["eps"], deps, pred,
                   deps / pred if pred != 0 else np.nan))


def main():
    gs = [0.0, 0.02, 0.05, 0.10]
    print("v86 N1/N2 instrument -- shooter grid: H=%g RMAX=%g N=%d  (FAST=%s)"
          % (H, rN, N, FAST))
    print(__doc__[__doc__.index("CONVENTION BANNER"):])

    branches = build_branches(gs, 1.312, 1.496, 0.004)

    rows = []
    for g in gs:
        for w in sorted(branches[g]):
            f, chi, rn = branches[g][w]
            rows.append(decompose(f, chi, w, g, rn))

    cols = ["g", "w", "Q", "E", "E_kin", "E_grad", "E_mass", "E_V", "E_g",
            "Sigma", "eps", "Sigma_pred", "eps_pred",
            "I1_rel", "I2_rel", "R_vir_rel", "I4_rel",
            "rhalf", "xi", "xi_over_R", "sigma_eff", "Eg_point", "coul_par",
            "f0", "chi0", "newton_res"]
    path = os.path.join(HERE, "n1_decomp.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for d in rows:
            fh.write("\t".join("%.10g" % d[c] for c in cols) + "\n")
    print("\nwrote %s  (%d rows)" % (path, len(rows)))

    # ---- N1 decomposition table --------------------------------------------
    print("\n" + "=" * 78)
    print("N1 -- component decomposition (E_kin = TRUE integrated kinetic, not wQ/2)")
    print("=" * 78)
    print("%-5s %-6s %9s %10s %9s %9s %9s %9s %9s %9s %9s" %
          ("g", "w", "Q", "E", "E_kin", "E_grad", "E_mass", "E_V", "E_g",
           "eps", "eps_pred"))
    for g in gs:
        sub = [d for d in rows if d["g"] == g]
        for d in sub[::max(1, len(sub) // 6)]:
            print("%-5.3f %-6.3f %9.2f %10.2f %9.2f %9.2f %9.2f %9.2f %9.3f %9.5f %9.5f" %
                  (d["g"], d["w"], d["Q"], d["E"], d["E_kin"], d["E_grad"],
                   d["E_mass"], d["E_V"], d["E_g"], d["eps"], d["eps_pred"]))

    # ---- N2 identity gate ---------------------------------------------------
    print("\n" + "=" * 78)
    print("N2 -- identity residuals (relative to E), max over branch")
    print("   I2 is algebraically -I1 (duplicate, shown for transparency only)")
    print("=" * 78)
    print("%-6s %5s %11s %11s %11s %11s %11s %11s" %
          ("g", "n", "I1(Gauss)", "I2(dup)", "R_vir(N2)", "I4(closed)",
           "newton", "eps range"))
    verdict_ok = True
    for g in gs:
        sub = [d for d in rows if d["g"] == g]
        if not sub:
            continue
        mI1 = max(abs(d["I1_rel"]) for d in sub)
        mI2 = max(abs(d["I2_rel"]) for d in sub)
        mRv = max(abs(d["R_vir_rel"]) for d in sub)
        mI4 = max(abs(d["I4_rel"]) for d in sub)
        mnr = max(d["newton_res"] for d in sub)
        emin = min(d["eps"] for d in sub)
        emax = max(d["eps"] for d in sub)
        print("%-6.3f %5d %11.3e %11.3e %11.3e %11.3e %11.3e %s" %
              (g, len(sub), mI1, mI2, mRv, mI4, mnr,
               "%.4f-%.4f" % (emin, emax)))
        # fix 4: pass gate -- residuals must be far below the physical eps signal
        if max(mRv, mI4, mI1) > 0.01 * abs(emin):
            verdict_ok = False
    print("\nN2 GATE: residuals must satisfy max(|I1|,|R_vir|,|I4|)/E < 1%% of min(eps).")
    print("N2 VERDICT: %s" % ("PASS" if verdict_ok else "FAIL -- see falsification list"))

    # ---- N1 outcome table (fix 7) ------------------------------------------
    print("\n" + "=" * 78)
    print("N1 outcome hooks -- which origin does Sigma have?")
    print("=" * 78)
    print("%-5s %-6s %9s %9s %9s %9s %9s %9s %9s" %
          ("g", "w", "Sigma", "R_half", "sigma_eff", "xi/R", "E_g", "Eg_point",
           "Eg/Egpt"))
    for g in gs:
        sub = [d for d in rows if d["g"] == g]
        for d in sub[::max(1, len(sub) // 5)]:
            print("%-5.3f %-6.3f %9.3f %9.3f %9.4f %9.4f %9.3f %9.3f %9.3f" %
                  (g, d["w"], d["Sigma"], d["rhalf"], d["sigma_eff"], d["xi_over_R"],
                   d["E_g"], d["Eg_point"],
                   d["E_g"] / d["Eg_point"] if d["Eg_point"] > 0 else np.nan))

    print("\nSigma composition after the I4 reduction (Sigma = 2/3 (E_grad - E_g)):")
    print("%-5s %-6s %10s %10s %10s %10s" %
          ("g", "w", "Sigma", "2/3 E_grad", "-2/3 E_g", "frac gauge"))
    for g in gs:
        sub = [d for d in rows if d["g"] == g]
        for d in sub[::max(1, len(sub) // 4)]:
            a = (2.0 / 3.0) * d["E_grad"]
            b = -(2.0 / 3.0) * d["E_g"]
            print("%-5.3f %-6.3f %10.4f %10.4f %10.4f %10.4f" %
                  (g, d["w"], d["Sigma"], a, b, -b / a if a else np.nan))

    phase_a(rows)


if __name__ == "__main__":
    main()
