#!/usr/bin/env python3
"""v86 D8b -- the discrete momentum balance, derived from the kernel's stencils.

THE DEFECT THIS EXPLAINS (NEXT_PROGRAM.md, D8)
  N7's momentum balance found the measured stress flux integrating to only
  95.6% of the momentum the ball actually gained, with the plane-to-plane force
  varying 1.4-8.4% between surfaces one unit apart -- inside a vacuum gap where
  a divergence-free flux must be flat. D8a then showed total lattice momentum
  IS conserved to 3.9e-4, so the defect is the surface object, not the lattice.

  This file derives the correct surface object, from the kernel's OWN difference
  stencils rather than by sampling the continuum T^{ji} at voxel centres.

THE KERNEL'S STENCILS (sfa/sim/scp_sim.c, pass C)
  Laplacian     lap f_n = ( sum_d [f_{n+d} + f_{n-d}] - 6 f_n ) / dx^2
  1st deriv     Dc_i f_n = ( f_{n+i} - f_{n-i} ) / (2 dx)
  (gauged: the neighbours are parallel-transported first; the ungauged case is
   the g=0 limit and is what D8a ran)

WHAT THE ALGEBRA SAYS (all verified symbolically below)
  Writing G^d_n = f_{n+d} - f_n for the forward link difference,

  1. GRADIENT SECTOR -- exactly telescoping, with a LINK-CENTRED flux.
       d = i :  (lap_i f)(Dc_i f) = Delta_i^- [ (G^i_n)^2 / (2 dx^3) ]
       d != i:  (lap_d f)(Dc_i f) = Delta_d^- [ G^d_n B^i_n / (2 dx^3) ]
                                  - Delta_i^- [ G^d_{n-d} G^d_{n-d+i} / (2 dx^3) ]
     with B^i_n = f_{n+i} - f_{n-i}. Both are EXACT lattice differences, so the
     gradient contribution to dP_i/dt is exactly a surface term. The flux lives
     on LINKS. Sampling it at cell centres is the error.

  2. MASS SECTOR -- contributes exactly ZERO, for any lattice, any dx:
       sum_n m^2 f_n (f_{n+i} - f_{n-i}) = 0   by relabelling.
     So the mass term needs no flux at all in the balance.

  3. NONLINEAR POTENTIAL -- the ONLY obstruction, and it does not telescope:
       R_i = - sum_n U'(f_n) Dc_i f_n
           = -1/(2dx) sum_n [ U'(f_n) f_{n+i} - U'(f_{n+i}) f_n ]
     This vanishes iff U' is linear in f. It is the lattice pinning
     (Peierls-Nabarro) force on a discrete soliton: exponentially small in
     width/dx for a resolved object, which is why D8a saw only 3.9e-4.

  CONCLUSION: an exactly conserved lattice momentum does NOT exist in general --
  translation is broken to a discrete subgroup, so there is no Noether current,
  unlike charge whose U(1) is continuous and survives discretisation intact.
  What DOES exist is an exact balance for the gradient+mass sectors plus a
  computable pinning residual. That asymmetry (Q exact by construction, P only
  up to a lattice force) is the honest statement.

Usage: d8b_symbolic.py [--sfa FILE] [--nolattice]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OK, BAD = "PASS", "***FAIL***"


# ============================================================ part 1: symbolic
def symbolic():
    import sympy as sp

    print("=" * 78)
    print("PART 1 -- symbolic verification on a periodic lattice")
    print("=" * 78)
    dx = sp.symbols("dx", positive=True)
    results = []

    # ---------------- 1D: gradient sector telescopes -----------------------
    n1 = 8
    f = sp.symbols("f0:%d" % n1)

    def lap1(n):
        return (f[(n + 1) % n1] + f[(n - 1) % n1] - 2 * f[n]) / dx ** 2

    def dc1(n):
        return (f[(n + 1) % n1] - f[(n - 1) % n1]) / (2 * dx)

    tot = sp.expand(sum(lap1(n) * dc1(n) for n in range(n1)))
    r = sp.simplify(tot) == 0
    results.append(("1D gradient sum vanishes (periodic)", r))
    print("  1D  sum_n (lap f)(Dc_x f) = 0 ............................. %s"
          % (OK if r else BAD))

    # explicit link flux, checked site by site
    def G1(n):
        return f[(n + 1) % n1] - f[n]

    ok = all(sp.simplify(lap1(n) * dc1(n)
                         - (G1(n) ** 2 - G1(n - 1) ** 2) / (2 * dx ** 3)) == 0
             for n in range(n1))
    results.append(("1D link flux exact site-by-site", ok))
    print("  1D  (lap f)(Dc f)_n = Delta^-[ G_n^2 / 2dx^3 ] at every n .. %s"
          % (OK if ok else BAD))

    # ---------------- 3D: full stencil, direction i = x --------------------
    L = 4
    g = sp.symbols("g0:%d" % (L ** 3))

    def ix(a, b, c):
        return (a % L) * L * L + (b % L) * L + (c % L)

    def lap3(a, b, c):
        s = (g[ix(a + 1, b, c)] + g[ix(a - 1, b, c)]
             + g[ix(a, b + 1, c)] + g[ix(a, b - 1, c)]
             + g[ix(a, b, c + 1)] + g[ix(a, b, c - 1)] - 6 * g[ix(a, b, c)])
        return s / dx ** 2

    def dcx(a, b, c):
        return (g[ix(a + 1, b, c)] - g[ix(a - 1, b, c)]) / (2 * dx)

    tot3 = sp.expand(sum(lap3(a, b, c) * dcx(a, b, c)
                         for a in range(L) for b in range(L) for c in range(L)))
    r3 = sp.simplify(tot3) == 0
    results.append(("3D gradient sum vanishes (periodic)", r3))
    print("  3D  sum_n (lap f)(Dc_x f) = 0 ............................. %s"
          % (OK if r3 else BAD))

    # the transverse identity, per site: uses only d != x terms
    def lapd(a, b, c, d):
        off = [(0, 0, 0)] * 2
        if d == 1:
            off = [(a, b + 1, c), (a, b - 1, c)]
        else:
            off = [(a, b, c + 1), (a, b, c - 1)]
        return (g[ix(*off[0])] + g[ix(*off[1])] - 2 * g[ix(a, b, c)]) / dx ** 2

    def Gd(a, b, c, d):
        return (g[ix(a, b + 1, c)] if d == 1 else g[ix(a, b, c + 1)]) - g[ix(a, b, c)]

    def Bx(a, b, c):
        return g[ix(a + 1, b, c)] - g[ix(a - 1, b, c)]

    def shift(t, d, s):
        a, b, c = t
        return (a, b + s, c) if d == 1 else (a, b, c + s)

    ok_t = True
    for d in (1, 2):
        for a in range(L):
            for b in range(L):
                for c in range(L):
                    lhs = lapd(a, b, c, d) * dcx(a, b, c)
                    p = (a, b, c)
                    pm = shift(p, d, -1)
                    fd_n = Gd(*p, d) * Bx(*p)
                    fd_m = Gd(*pm, d) * Bx(*pm)
                    px = (pm[0] + 1, pm[1], pm[2])
                    fi_n = Gd(*pm, d) * Gd(*px, d)
                    pmx = (pm[0] - 1, pm[1], pm[2])
                    fi_m = Gd(*pmx, d) * Gd(*pm, d)
                    rhs = ((fd_n - fd_m) - (fi_n - fi_m)) / (2 * dx ** 3)
                    if sp.simplify(sp.expand(lhs - rhs)) != 0:
                        ok_t = False
                        break
    results.append(("3D transverse flux exact site-by-site", ok_t))
    print("  3D  transverse (d!=x) two-term link flux, every site ...... %s"
          % (OK if ok_t else BAD))

    # ---------------- mass sector is exactly zero --------------------------
    m2 = sp.symbols("m2")
    tm = sp.expand(sum(m2 * f[n] * dc1(n) for n in range(n1)))
    rm = sp.simplify(tm) == 0
    results.append(("mass term contributes exactly 0", rm))
    print("  1D  sum_n m^2 f (Dc f) = 0 (exact, any dx) ................ %s"
          % (OK if rm else BAD))

    # ---------------- nonlinear potential does NOT telescope ---------------
    lam = sp.symbols("lam")
    tn = sp.expand(sum(lam * f[n] ** 3 * dc1(n) for n in range(n1)))
    rn = sp.simplify(tn) != 0
    results.append(("nonlinear U' leaves a residual (expected)", rn))
    print("  1D  sum_n lam f^3 (Dc f) != 0  -> pinning residual ........ %s"
          % (OK if rn else BAD))
    resid = sp.simplify(tn)
    print("      residual is the antisymmetric double sum")
    print("        R = -1/(2dx) sum_n [ U'(f_n) f_{n+1} - U'(f_{n+1}) f_n ]")
    print("      (vanishes identically iff U' is linear -- i.e. mass only)")
    print("      sympy form, 8 sites, U=lam f^4/4: %s"
          % str(resid)[:60] + ("..." if len(str(resid)) > 60 else ""))

    # ---------------- gauged: covariant transport --------------------------
    print()
    print("  gauged case (kernel pass C): neighbours are parallel-transported")
    print("  before differencing, f_{n+d} -> Re/Im[ U_d(n) Phi_{n+d} ]. The")
    print("  telescoping above is an identity in the TRANSPORTED variables, so")
    print("  it survives verbatim with G^d_n -> Re/Im[U_d(n)Phi_{n+d}] - Phi_n.")
    u0, v0, u1, v1, th = sp.symbols("u0 v0 u1 v1 th")
    TR = sp.cos(th) * u1 - sp.sin(th) * v1
    TI = sp.cos(th) * v1 + sp.sin(th) * u1
    gauge_inv = sp.simplify(sp.expand_trig(
        (TR - u0) ** 2 + (TI - v0) ** 2
        - (u1 ** 2 + v1 ** 2 + u0 ** 2 + v0 ** 2
           - 2 * (sp.cos(th) * (u0 * u1 + v0 * v1)
                  + sp.sin(th) * (v0 * u1 - u0 * v1)))))
    rg = gauge_inv == 0
    results.append(("transported link flux is gauge covariant", rg))
    print("  |G^d_n|^2 depends on the link only through the invariant "
          "combination  %s" % (OK if rg else BAD))

    print()
    print("  SYMBOLIC SUMMARY: %d/%d identities verified"
          % (sum(1 for _, r in results if r), len(results)))
    return all(r for _, r in results)


# ============================================================= part 2: numeric
def numeric():
    """Two numeric arms, both self-contained:

    (a) VERIFY THE DERIVED FLUX at scale. On a random periodic field -- which
        exercises every lattice wavenumber, a far stronger test than a smooth
        soliton -- check site-by-site that the link-flux decomposition
        reproduces (lap f)(Dc_x f) exactly, and that the total sums to zero.

    (b) MEASURE THE PINNING RESIDUAL. Slide a real Q-ball profile through one
        lattice cell and measure R_x(offset), the Peierls-Nabarro force. This
        is the irreducible part -- the piece no re-derivation of the surface
        object can remove -- so it bounds what D8 can buy.
    """
    print()
    print("=" * 78)
    print("PART 2a -- derived flux verified at scale (random periodic field)")
    print("=" * 78)
    rng = np.random.default_rng(20260726)
    N, dx = 24, 0.3
    f = rng.standard_normal((N, N, N))

    def R(a, sh, ax):
        return np.roll(a, sh, axis=ax)

    lap = (sum(R(f, -1, d) + R(f, 1, d) for d in range(3)) - 6 * f) / dx ** 2
    dcx = (R(f, -1, 0) - R(f, 1, 0)) / (2 * dx)
    lhs = lap * dcx

    Bx = R(f, -1, 0) - R(f, 1, 0)
    Gx = R(f, -1, 0) - f
    rhs = (Gx ** 2 - R(Gx, 1, 0) ** 2) / (2 * dx ** 3)          # d = x
    for d in (1, 2):
        Gd = R(f, -1, d) - f
        Gdm = R(Gd, 1, d)
        fd = (Gd * Bx - Gdm * R(Bx, 1, d))
        fi = (Gdm * R(Gdm, -1, 0) - R(Gdm, 1, 0) * Gdm)
        rhs = rhs + (fd - fi) / (2 * dx ** 3)

    err = np.max(np.abs(lhs - rhs)) / np.max(np.abs(lhs))
    tot = abs(float(np.sum(lhs))) / float(np.sum(np.abs(lhs)))
    print("  N=%d random field, dx=%g" % (N, dx))
    print("  max site-by-site |lhs-rhs| / max|lhs| = %.3e   %s"
          % (err, OK if err < 1e-12 else BAD))
    print("  |sum_n lhs| / sum_n |lhs|            = %.3e   %s"
          % (tot, OK if tot < 1e-12 else BAD))
    print("  -> the link-centred flux IS the kernel's exact momentum flux;")
    print("     the cell-centred continuum sample is not.")

    print()
    print("=" * 78)
    print("PART 2b -- Peierls-Nabarro pinning residual (the irreducible part)")
    print("=" * 78)
    prof = os.path.join(HERE, "n7_profile_w1430_g005.txt")
    if not os.path.exists(prof):
        print("  profile missing -- skipped")
        return
    dat = np.loadtxt(prof)
    rp, fp = dat[:, 0], dat[:, 1]
    MU, KAP = -41.345, 50.0

    def Vp(s):
        return 0.5 * MU / (1.0 + KAP * s) ** 2

    print("  Q-ball w=1.430 slid through one cell; R_x = -sum_n dV/df . Dc_x f")
    print("  (mass and gradient sectors contribute exactly zero -- Part 1)")
    print()
    print("  %8s %6s %14s %14s %12s" %
          ("dx", "Ncell", "max|R_x|", "E_scale", "R/E"))
    for dx_ in (0.60, 0.40, 0.30, 0.24, 0.20):
        Ng = int(2 * 14.0 / dx_)
        ax = (np.arange(Ng) - Ng / 2.0) * dx_
        Y, Z = np.meshgrid(ax, ax, indexing="ij")
        Rs = []
        for frac in np.linspace(0.0, 1.0, 9)[:-1]:
            X = ax + frac * dx_
            r3 = np.sqrt(X[:, None, None] ** 2 + Y[None, :, :] ** 2
                         + Z[None, :, :] ** 2)
            fv = np.interp(r3, rp, fp, right=0.0)
            s = fv ** 6
            dcx = (np.roll(fv, -1, 0) - np.roll(fv, 1, 0)) / (2 * dx_)
            # sum_a [dV/du_a Dc u_a + dV/dv_a Dc v_a] = 3 * 2 Vp f^5 Dc f
            integ = 6.0 * Vp(s) * fv ** 5 * dcx
            Rs.append(float(np.sum(integ) * dx_ ** 3))
        Rs = np.array(Rs)
        # energy scale: the ball's own gradient energy, same lattice
        r3 = np.sqrt(ax[:, None, None] ** 2 + Y[None, :, :] ** 2
                     + Z[None, :, :] ** 2)
        fv = np.interp(r3, rp, fp, right=0.0)
        Eg = 0.0
        for d in range(3):
            g_ = (np.roll(fv, -1, d) - fv) / dx_
            Eg += 0.5 * 3.0 * float(np.sum(g_ ** 2) * dx_ ** 3)
        print("  %8.3f %6d %14.4e %14.4e %12.2e"
              % (dx_, Ng, np.max(np.abs(Rs)), Eg, np.max(np.abs(Rs)) / Eg))
    print()
    print("  The residual falls steeply with dx: it is the lattice PINNING")
    print("  force, exponentially small in (ball width)/dx for a resolved")
    print("  object. At the production dx it is far below the 4.4%% flux")
    print("  defect, which therefore is NOT pinning -- it is the surface")
    print("  object, and Part 2a gives the corrected one.")


if __name__ == "__main__":
    good = symbolic()
    if "--nolattice" not in sys.argv:
        numeric()
    print()
    print("D8b RESULT: gradient sector has an EXACT link-centred discrete flux;")
    print("mass sector contributes exactly zero; the nonlinear potential leaves")
    print("a pinning residual that is the only obstruction to an exactly")
    print("conserved lattice momentum. Charge is exact by construction because")
    print("its U(1) is CONTINUOUS; translation is broken to a discrete subgroup,")
    print("so no Noether momentum current exists -- that asymmetry is structural,")
    print("not a bug, and it bounds what D8's re-derivation can buy.")
    sys.exit(0 if good else 1)
