#!/usr/bin/env python3
"""v86 Part-0 rung N9 -- soft-window (mu, kappa) shooter scan.

Protocol: v86/council/grok45/N_BATTERY_REVIEW.md section 1, N9 (restored).
  "2D shooter grid about the standards (-41.345, 50). Map existence window
   width, eps range, and E/Q dynamic range at fixed Q (or fixed mid-branch
   omega)."
  Discriminator:
    E at fixed Q varies by >= 2x inside VK-stable region -> single-sector HIER
        softens; second sector less mandatory for hierarchy alone.
    Window always <= 10% in E/Q -> single-sector hierarchy blocked; the second
        sector is mandatory for HIER.

AMBIGUITY FLAG (pre-registered, and THRESHOLD ASSIGNMENT corrected by the
grok-4.5 review seat, Findings 1-2-5). The adopted one-liner does not say
whether "how much can M/Q move in one sector" means
  (a) the dynamic range of E/Q ALONG the VK-stable branch of ONE fixed theory
      (mu,kappa) -- the "one theory, one sector" window; or
  (b) the range of E at FIXED Q ACROSS the (mu,kappa) family -- POTENTIAL-
      REDESIGN freedom (different (mu,kappa) are different theories/vacua, NOT
      different states of one conserved sector; (b) answers "could a potential
      redesign replace a second sector?", never "one sector hosts a hierarchy").
FOUNDATIONS' two numeric gates map onto these as follows:
  * the 2x gate belongs to (b): E_max/E_min at fixed Q* across theories.
  * the 10% gate belongs to (a): the along-branch E/Q window, and "window
    ALWAYS <= 10%" means EVERY grid theory, not the best one.
PRE-REGISTERED, before running: the along-branch E/Q range is bounded above by
roughly m/omega_min (the skeleton E/Q ~ omega(1+eps)); on the standard
potential omega_min=1.3087 gives ~1.15x, and the softest corner of this grid
gives ~1.5x. A 2x along-branch range would need omega_min ~ 0.75, which is
OUTSIDE this grid. Applying the 2x gate to (a) would therefore be a
pre-doomed test, and is not done. For the frozen program theory (mu,kappa
fixed at the standards) only reading (a) at the STANDARD point binds HIER;
the grid-best (a) is reported separately as redesign headroom.

Scope: UNGAUGED (g=0). Rationale: the gauge sector adds the Coulomb capacity
fold, which SHRINKS the accessible window (v86 A3); the ungauged branch is
therefore the generous bound on single-sector freedom. Reported as such.

Cost: free CPU. Grid: shooter continuation in omega at each (mu,kappa).
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/d/code/scp"
sys.path.insert(0, os.path.join(ROOT, "v69/theory"))
import gauged_shooter_fast as G          # H=0.01, RMAX=100 (ungauged: plenty)

M2 = G.M2
M = G.M
N, H, rN = G.N, G.H, G.rN
r = G.r

MU0, KAP0 = -41.345, 50.0


def set_potential(mu, kap):
    """The shooter reads MU/KAP as module globals inside Wpot/dWpot/Vt, so the
    potential family can be swapped in place."""
    G.MU, G.KAP = mu, kap


def window_bottom(mu, kap):
    """Ungauged Q-ball window bottom: omega_min^2 = m^2 + (mu/9)(2/kappa)^(2/3).
    (Same formula the shooter uses for the standard potential; it is the
    condition that the effective potential first touches zero.)"""
    arg = M2 + (mu / 9.0) * (2.0 / kap) ** (2.0 / 3.0)
    return np.sqrt(arg) if arg > 0 else np.nan


def observables_ungauged(f, w):
    fe = np.append(f, 0.0)
    re = np.append(r, rN)
    fpe = np.gradient(fe, H, edge_order=2)
    dV = 4.0 * np.pi * re * re
    Q = np.trapz(3.0 * w * fe * fe * dV, dx=H)
    E_kin = np.trapz(1.5 * w * w * fe * fe * dV, dx=H)
    E_grad = np.trapz(1.5 * fpe * fpe * dV, dx=H)
    E_mass = np.trapz(1.5 * M2 * fe * fe * dV, dx=H)
    E_V = np.trapz(G.Vt(fe ** 6) * dV, dx=H)
    E = E_kin + E_grad + E_mass + E_V
    Sigma = E - w * Q
    return dict(w=w, Q=Q, E=E, E_kin=E_kin, E_grad=E_grad, E_mass=E_mass,
                E_V=E_V, Sigma=Sigma, eps=Sigma / (w * Q), EoQ=E / Q,
                # N2 closed form carries over unchanged in the ungauged case
                I4_rel=(Sigma - (2.0 / 3.0) * E_grad) / E,
                # evaporation cut (grok Finding 6): E >= m*Q means the object is
                # heavier than Q free quanta and is unbound against dispersal;
                # such points must not be allowed to inflate the E/Q extremes.
                bound=bool(E < M * Q), EoverMQ=E / (M * Q),
                f0=f[0])


W_ANCHOR = 1.45         # omega at which the (mu,kappa) continuation is carried


def continue_potential(f_from, mu_from, kap_from, mu_to, kap_to, w=W_ANCHOR):
    """Adaptive continuation of the anchor solution in the potential parameters.
    Re-seeding every theory from the standard v66 profile is what made the first
    version of this rung intractable (each miss costs a full 60-iteration Newton
    with line search); walking (mu,kappa) instead keeps every solve warm."""
    f = f_from.copy()
    mu, kap = mu_from, kap_from
    for target, which in ((mu_to, "mu"), (kap_to, "kap")):
        cur = mu if which == "mu" else kap
        step = target - cur
        while abs(cur - target) > 1e-12:
            nxt = cur + step
            if (step > 0 and nxt > target) or (step < 0 and nxt < target):
                nxt = target
            set_potential(nxt if which == "mu" else mu,
                          kap if which == "mu" else nxt)
            fn, cn, ok, _, rn = G.solve(w, 0.0, f, np.zeros(N))
            if ok:
                f, cur = fn, nxt
                if which == "mu":
                    mu = cur
                else:
                    kap = cur
                step = target - cur
            else:
                step *= 0.5
                if abs(step) < 1e-6 * max(1.0, abs(target)):
                    return None
    set_potential(mu_to, kap_to)
    return f


def branch(mu, kap, seed, dw=0.006, dw_min=1e-4):
    """Sweep the ungauged branch from the anchor omega in both directions,
    returning a list of observables dicts ordered in omega. `seed` is a
    converged profile for THIS (mu,kappa) at W_ANCHOR."""
    set_potential(mu, kap)
    wb = window_bottom(mu, kap)
    if not np.isfinite(wb) or wb >= M - 1e-3:
        return [], wb
    if seed is None:
        return [], wb
    start = (W_ANCHOR, seed)
    wlo, whi = wb, M
    sols = {}
    for direction in (-1, +1):
        w, f = start[0], start[1].copy()
        sols[round(w, 8)] = f.copy()
        step = dw
        wlast = w
        w_edge, f_edge = None, None
        while True:
            wn = w + direction * step
            if wn <= wb + 1e-6 or wn >= M - 1e-6:
                break
            fn, cn, ok, _, rn = G.solve(wn, 0.0, f, np.zeros(N))
            if ok:
                f, w = fn, wn
                w_edge, f_edge = w, f.copy()      # fix 3: never lose the edge
                if abs(w - wlast) >= 0.999 * dw:
                    sols[round(w, 8)] = f.copy()
                    wlast = w
                step = dw
            else:
                step *= 0.5
                if step < dw_min:
                    break
        # fix 3 (grok Finding 6): force-store the last converged edge state, so
        # adaptive sub-dw steps near the branch bottom cannot silently truncate
        # Qmax / the E/Q extremes that N9 actually scores.
        if w_edge is not None:
            sols[round(w_edge, 8)] = f_edge
    rows = [observables_ungauged(sols[w], w) for w in sorted(sols)]
    return rows, wb


def vk_stable(rows):
    """VK criterion: dQ/dw < 0 marks the stable branch."""
    out = []
    for i, d in enumerate(rows):
        j0 = max(0, i - 1)
        j1 = min(len(rows) - 1, i + 1)
        if j1 == j0:
            out.append(False)
            continue
        dQdw = (rows[j1]["Q"] - rows[j0]["Q"]) / (rows[j1]["w"] - rows[j0]["w"])
        out.append(dQdw < 0.0)
    return out


def main():
    mus = [-41.345 * s for s in (0.6, 0.8, 1.0, 1.25, 1.5)]
    kaps = [50.0 * s for s in (0.5, 0.75, 1.0, 1.5, 2.0)]

    print("v86 N9 -- soft-window (mu,kappa) scan   (UNGAUGED, H=%g RMAX=%g)"
          % (H, rN))
    print("standards: mu=%.3f kappa=%.1f\n" % (MU0, KAP0))

    # anchor solution for the STANDARD potential, then walk the grid by
    # continuation (never re-seed from the v66 profile: see continue_potential)
    set_potential(MU0, KAP0)
    f0g = G.load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
    f_std, _, ok, _, rn = G.solve(W_ANCHOR, 0.0, f0g, np.zeros(N))
    if not ok:
        sys.exit("standard-potential anchor solve failed (resid %.2e)" % rn)
    print("anchor: standard potential at w=%.3f solved (resid %.2e)\n" % (W_ANCHOR, rn))

    # row anchors: walk mu at the standard kappa
    row_anchor = {}
    for mu in mus:
        fa = continue_potential(f_std, MU0, KAP0, mu, KAP0)
        row_anchor[mu] = fa
        if fa is None:
            print("  (mu=%.3f, kappa=%.1f) anchor continuation FAILED" % (mu, KAP0))

    table = []
    print("%-9s %-7s %8s %8s %9s %10s %10s %10s %8s %10s" %
          ("mu", "kappa", "w_bot", "w_span", "n_VK", "Q range", "E/Q min",
           "E/Q max", "EoQ dyn", "eps range"))
    for mu in mus:
        if row_anchor[mu] is None:
            continue
        for kap in kaps:
            seed = continue_potential(row_anchor[mu], mu, KAP0, mu, kap)
            rows, wb = branch(mu, kap, seed)
            if not rows:
                print("%-9.3f %-7.1f %8s  -- no branch (anchor %s) --"
                      % (mu, kap, "%.4f" % wb,
                         "failed" if seed is None else "ok, sweep empty"))
                continue
            st = vk_stable(rows)
            sr_all = [d for d, s in zip(rows, st) if s]
            sr = [d for d in sr_all if d["bound"]]
            ndrop = len(sr_all) - len(sr)
            if not sr:
                print("%-9.3f %-7.1f %8.4f  -- no bound VK-stable segment "
                      "(%d VK points all unbound, E>=mQ) --" % (mu, kap, wb, ndrop))
                continue
            EoQ = [d["EoQ"] for d in sr]
            Qs = [d["Q"] for d in sr]
            eps = [d["eps"] for d in sr]
            rec = dict(mu=mu, kap=kap, wb=wb,
                       wspan=max(d["w"] for d in sr) - min(d["w"] for d in sr),
                       n=len(sr), Qmin=min(Qs), Qmax=max(Qs),
                       EoQmin=min(EoQ), EoQmax=max(EoQ), dyn=max(EoQ) / min(EoQ),
                       epsmin=min(eps), epsmax=max(eps), rows=sr, ndrop=ndrop)
            table.append(rec)
            print("%-9.3f %-7.1f %8.4f %8.4f %9d %10s %10.5f %10.5f %8.4f %10s%s" %
                  (mu, kap, wb, rec["wspan"], rec["n"],
                   "%.0f-%.0f" % (rec["Qmin"], rec["Qmax"]),
                   rec["EoQmin"], rec["EoQmax"], rec["dyn"],
                   "%.3f-%.3f" % (rec["epsmin"], rec["epsmax"]),
                   ("  [%d unbound dropped]" % ndrop) if ndrop else ""))

    if not table:
        print("\nNo branches found -- N9 BLOCKED.")
        return

    # ---------------- reading (a): within one theory ------------------------
    print("\n" + "=" * 78)
    print("N9 reading (a) -- E/Q window ALONG the VK-stable branch of ONE theory")
    print("   gate: the 10%% test. 'window ALWAYS <= 10%%' means EVERY theory on")
    print("   the grid, not the best one. The 2x gate does NOT belong here.")
    print("=" * 78)
    best = max(table, key=lambda t: t["dyn"])
    std = [t for t in table if abs(t["mu"] - MU0) < 1e-6 and abs(t["kap"] - KAP0) < 1e-6]
    if std:
        s = std[0]
        print("STANDARD potential (mu=%.3f, kappa=%.1f) -- the frozen program theory,"
              % (s["mu"], s["kap"]))
        print("  the only point that binds HIER for the theory as it stands:")
        print("  E/Q in [%.5f, %.5f] -> window %.2f%% over Q=%.0f..%.0f"
              % (s["EoQmin"], s["EoQmax"], 100 * (s["dyn"] - 1.0),
                 s["Qmin"], s["Qmax"]))
    print("grid-best (redesign headroom, NOT the frozen theory): "
          "mu=%.3f kappa=%.1f -> window %.2f%%"
          % (best["mu"], best["kap"], 100 * (best["dyn"] - 1.0)))
    n_le10 = sum(1 for t in table if t["dyn"] <= 1.10)
    print("theories with window <= 10%%: %d of %d" % (n_le10, len(table)))
    if n_le10 == len(table):
        va = ("BLOCKED -- window <= 10%% for EVERY theory on the grid; a single "
              "sector cannot host a mass hierarchy, second sector mandatory for HIER")
    elif std and std[0]["dyn"] <= 1.10:
        va = ("BLOCKED FOR THE FROZEN THEORY -- the standard potential's window is "
              "%.2f%% <= 10%%; %d of %d grid theories exceed 10%%, so a potential "
              "REDESIGN (not a second sector) could widen it"
              % (100 * (std[0]["dyn"] - 1.0), len(table) - n_le10, len(table)))
    else:
        va = ("NOT BLOCKED -- the standard potential's own window is %.2f%% > 10%%; "
              "FOUNDATIONS' 'always <= 10%%' condition fails"
              % (100 * (std[0]["dyn"] - 1.0) if std else np.nan))
    print("VERDICT (a) [10%% gate]: %s" % va)

    # ---------------- reading (b): across the potential family --------------
    print("\n" + "=" * 78)
    print("N9 reading (b) -- E at FIXED Q across the (mu,kappa) family")
    print("   = POTENTIAL-REDESIGN freedom. Different (mu,kappa) are different")
    print("   theories/vacua, NOT different states of one conserved sector, so a")
    print("   large (b) never means 'one sector hosts a hierarchy'; it means a")
    print("   redesigned potential is a candidate alternative to a second sector.")
    print("   gate: the 2x test.")
    print("=" * 78)
    Qtargets = [120.0, 200.0, 400.0]
    print("%-8s %-9s %-7s %10s %10s %9s %9s" %
          ("Q*", "mu", "kappa", "w", "E", "E/Q", "|dQ|/Q*"))
    dyn_b = {}
    for Qt in Qtargets:
        vals = []
        for t in table:
            rows = t["rows"]
            Qs = np.array([d["Q"] for d in rows])
            Es = np.array([d["E"] for d in rows])
            ws = np.array([d["w"] for d in rows])
            if Qt < Qs.min() or Qt > Qs.max():
                continue
            # fix 6: interpolate to the exact Q* rather than nearest-neighbour
            order = np.argsort(Qs)
            Ei = float(np.interp(Qt, Qs[order], Es[order]))
            wi = float(np.interp(Qt, Qs[order], ws[order]))
            qerr = float(np.min(np.abs(Qs - Qt)) / Qt)
            vals.append((t["mu"], t["kap"], wi, Ei, Ei / Qt, qerr))
        if len(vals) < 2:
            print("Q*=%-6.0f  (fewer than 2 theories reach this charge)" % Qt)
            continue
        Es = [v[3] for v in vals]
        lo = vals[int(np.argmin(Es))]
        hi = vals[int(np.argmax(Es))]
        for v in (lo, hi):
            print("%-8.0f %-9.3f %-7.1f %10.4f %10.2f %9.5f %9.2e" %
                  (Qt, v[0], v[1], v[2], v[3], v[4], v[5]))
        dyn_b[Qt] = max(Es) / min(Es)
        print("        -> E dynamic range at Q=%.0f across %d theories: %.3fx"
              % (Qt, len(vals), dyn_b[Qt]))
    if dyn_b:
        bb = max(dyn_b.values())
        vb = ("REDESIGN HEADROOM CONFIRMED (>=2x at fixed Q across the family) -- "
              "a potential redesign is a live alternative to a second sector for "
              "HIER; it is NOT single-sector hierarchy" if bb >= 2.0 else
              "NO REDESIGN HEADROOM ON THIS GRID (%.2fx < 2x at fixed Q)" % bb)
        print("VERDICT (b) [2x gate]: %s" % vb)

    # ---------------- eps behaviour (ties N1-N3) ----------------------------
    print("\n" + "=" * 78)
    print("N9 -- eps and the N2 closed form across the family")
    print("=" * 78)
    print("max |Sigma - (2/3)E_grad| / E over ALL (mu,kappa,omega) points: %.3e"
          % max(abs(d["I4_rel"]) for t in table for d in t["rows"]))
    print("(the N2 identity is potential-shape independent by construction -- this")
    print(" is a check that the Derrick derivation did not smuggle in mu/kappa)")

    # ---------------- dump --------------------------------------------------
    path = os.path.join(HERE, "n9_softwindow.tsv")
    cols = ["mu", "kap", "w", "Q", "E", "EoQ", "Sigma", "eps", "E_grad",
            "E_kin", "E_mass", "E_V", "I4_rel", "f0"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for t in table:
            for d in t["rows"]:
                d2 = dict(d, mu=t["mu"], kap=t["kap"])
                fh.write("\t".join("%.10g" % d2[c] for c in cols) + "\n")
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
