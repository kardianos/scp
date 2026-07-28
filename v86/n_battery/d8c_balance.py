#!/usr/bin/env python3
"""v86 D8c -- the discrete momentum balance CLOSED on real field data.

WHY THIS EXISTS: the grok-4.5 review, Finding 3.5 (MAJOR), was right --

    "D8b did not re-integrate the link flux on the N7 run and show the residual
     collapse. 'Defect is the surface object and Part 2a gives the corrected
     one' is a derivation claim, NOT a closed re-measurement."

This file closes it. D8b derived the flux; this measures it against dP/dt on the
D8a run (one boosted ball, PERIODIC box, ungauged, no sponge, N=64, L=16) --
deliberately the clean case, with no sponge and no gauge sector to muddy the
comparison.

THE FULL EXACT BALANCE (extending D8b, which only did the gradient sector)

Every term of dP_i/dt telescopes except one. With pi = phidot,
G^d_n = f_{n+d} - f_n, B^i_n = f_{n+i} - f_{n-i}:

  d/dt P_i = sum_n [ (lap f)(Dc_i f) + pi (Dc_i pi) - m^2 f (Dc_i f)
                     - U'_nl(f) (Dc_i f) ]

  term            telescopes?   link flux
  --------------  -----------   ---------------------------------------
  gradient        YES           D8b Part 1 (two pieces for d != i)
  kinetic         YES           pi_n pi_{n+i} / (2 dx)
  mass            YES           -m^2 f_n f_{n+i} / (2 dx)
  nonlinear pot.  NO            Peierls-Nabarro residual, D8b Part 2b

So the exact statement is

    d/dt P_i(region) = [link flux through the boundary] + [PN residual inside]

and BOTH sides are computed here from the same frames. The comparison is against
the naive object -- continuum T^{ji} sampled at cell centres -- which is what
sfa_momentum.c currently integrates and what produced N7's 4.4% shortfall.

Usage: d8c_balance.py [--sfa FILE] [--frames a,b] [--out prefix]
"""
import os
import sys
import subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M2, MU, KAP = 2.25, -41.345, 50.0
NF = 3


def Vt(s):
    return 0.5 * MU * s / (1.0 + KAP * s)


def Vp(s):
    return 0.5 * MU / (1.0 + KAP * s) ** 2


def build_dumper():
    exe = os.path.join(HERE, "_d8b_dump")
    src = exe + ".c"
    if os.path.exists(exe):
        return exe
    with open(src, "w") as fh:
        fh.write(r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
int main(int c,char**v){SFA*s=sfa_open(v[1]);if(!s)return 1;
 size_t nb=0;for(uint32_t i=0;i<s->n_columns;i++)nb+=s->N_total*sfa_dtype_size[s->columns[i].dtype];
 void*b=malloc(nb); int fr=atoi(v[2]); if(sfa_read_frame(s,fr,b)!=0)return 2;
 fprintf(stderr,"%u %.10g %u\n",s->Nx,s->Lx,s->total_frames);
 const char*want[12]={"phi_x","phi_y","phi_z","phiim_x","phiim_y","phiim_z",
   "phi_vx","phi_vy","phi_vz","phiim_vx","phiim_vy","phiim_vz"};
 for(int w=0;w<12;w++){uint64_t off=0;int found=-1;
  for(uint32_t i=0;i<s->n_columns;i++){if(!strncmp(s->columns[i].name,want[w],32)){found=i;break;}
   off+=s->N_total*sfa_dtype_size[s->columns[i].dtype];}
  if(found<0){fprintf(stderr,"missing %s\n",want[w]);return 3;}
  fwrite((char*)b+off,4,s->N_total,stdout);} return 0;}
''')
    r = subprocess.run(["gcc", "-O2", "-o", exe, src, "-lzstd", "-lm"],
                       cwd=HERE, capture_output=True)
    if r.returncode:
        print(r.stderr.decode()[:400])
        sys.exit(1)
    return exe


def load(path, frame):
    exe = build_dumper()
    p = subprocess.run([exe, path, str(frame)], capture_output=True)
    if p.returncode:
        print("dump failed:", p.stderr.decode()[:200])
        sys.exit(1)
    tok = p.stderr.decode().split()
    N, L, nfr = int(tok[0]), float(tok[1]), int(tok[2])
    a = np.frombuffer(p.stdout, dtype=np.float32).astype(np.float64)
    a = a.reshape(12, N, N, N)
    return (list(a[0:3]), list(a[3:6]), list(a[6:9]), list(a[9:12])), N, L, nfr


def R(a, sh, ax):
    return np.roll(a, sh, axis=ax)


def momentum_x(u, v, ud, vd, dx):
    """P_x density with the SAME centred difference the D8b identity uses."""
    p = np.zeros_like(u[0])
    for a in range(NF):
        p += ud[a] * (R(u[a], -1, 0) - R(u[a], 1, 0)) / (2 * dx)
        p += vd[a] * (R(v[a], -1, 0) - R(v[a], 1, 0)) / (2 * dx)
    return p


def link_flux_x(u, v, ud, vd, dx):
    """The EXACT discrete x-flux of x-momentum, per D8b + the kinetic and mass
    terms derived in this file's header. Lives on the link (n, n+x)."""
    F = np.zeros_like(u[0])
    for fld, fdot in ((u, ud), (v, vd)):
        for a in range(NF):
            f, pi = fld[a], fdot[a]
            Gx = R(f, -1, 0) - f
            Bx = R(f, -1, 0) - R(f, 1, 0)
            # NOTE ON dx BOOKKEEPING (this cost a wrong answer once): the
            # telescoping identity gives A_n = F_n - F_{n-1} as a DENSITY, so
            # summing A over a slab with weight dx^3 leaves (F_R - F_L)*dx^3.
            # The per-AREA flux is therefore F*dx, which is what is returned.
            # With that factor the terms are exactly the continuum stress with
            # two-point link averaging:
            #   1/2 (d_x phi)^2 - 1/2 (d_y phi)^2 - 1/2 (d_z phi)^2
            #   + 1/2 phidot^2 - 1/2 m^2 phi^2
            # and NO potential term -- Vt(s) does not telescope, so on the
            # lattice it lives in the BULK residual, not on the surface. That
            # is a structural difference from the continuum T^{xx}, not an
            # approximation.
            # gradient, d = x
            F += Gx ** 2 / (2 * dx ** 2)
            # gradient, d != x : the two-term piece that carries x-flux
            for d in (1, 2):
                Gd = R(f, -1, d) - f
                Gdm = R(Gd, 1, d)
                F -= Gdm * R(Gdm, -1, 0) / (2 * dx ** 2)
            # kinetic
            F += pi * R(pi, -1, 0) / 2.0
            # mass
            F -= M2 * f * R(f, -1, 0) / 2.0
    return F


def naive_flux_x(u, v, ud, vd, dx):
    """What sfa_momentum.c integrates: continuum T^{xx} at cell centres."""
    def cd(f, ax):
        return (R(f, -1, ax) - R(f, 1, ax)) / (2 * dx)
    s = np.ones_like(u[0])
    for a in range(NF):
        s = s * (u[a] ** 2 + v[a] ** 2)
    gradsq = np.zeros_like(u[0])
    dxx = np.zeros_like(u[0])
    kin = np.zeros_like(u[0])
    mass = np.zeros_like(u[0])
    for a in range(NF):
        for d in range(3):
            gradsq += cd(u[a], d) ** 2 + cd(v[a], d) ** 2
        dxx += cd(u[a], 0) ** 2 + cd(v[a], 0) ** 2
        kin += ud[a] ** 2 + vd[a] ** 2
        mass += u[a] ** 2 + v[a] ** 2
    return dxx + 0.5 * kin - 0.5 * gradsq - 0.5 * M2 * mass - Vt(s)


def pn_residual_x(u, v, dx):
    """The non-telescoping nonlinear-potential term, summed over a region."""
    s = np.ones_like(u[0])
    for a in range(NF):
        s = s * (u[a] ** 2 + v[a] ** 2)
    vp = Vp(s)
    res = np.zeros_like(u[0])
    for fld in (u, v):
        for a in range(NF):
            f = fld[a]
            prest = np.ones_like(f)
            for b in range(NF):
                if b != a:
                    prest = prest * (u[b] ** 2 + v[b] ** 2)
            dUdf = 2.0 * vp * f * prest
            res -= dUdf * (R(f, -1, 0) - R(f, 1, 0)) / (2 * dx)
    return res


def main():
    path = "/space/scp/v86/n7/d8a.sfa"
    if "--sfa" in sys.argv:
        path = sys.argv[sys.argv.index("--sfa") + 1]
    fa, fb = 2, 3
    if "--frames" in sys.argv:
        fa, fb = [int(x) for x in sys.argv[sys.argv.index("--frames") + 1].split(",")]

    print("=" * 92)
    print("v86 D8c -- discrete momentum balance closed on real data")
    print("=" * 92)
    (uA, vA, udA, vdA), N, L, nfr = load(path, fa)
    (uB, vB, udB, vdB), _, _, _ = load(path, fb)
    dx = 2.0 * L / (N - 1)
    print("  %s   N=%d L=%g dx=%.6f frames=%d  (using %d,%d)"
          % (os.path.basename(path), N, L, dx, nfr, fa, fb))

    # snapshot cadence from the diag file if present
    dt_snap = None
    dg = os.path.join(os.path.dirname(path), "d8a_diag.tsv")
    if os.path.exists(dg):
        t = np.loadtxt(dg, skiprows=1, usecols=0)
        dt_snap = (t[-1] - t[0]) / max(nfr - 1, 1)
    if dt_snap is None or dt_snap <= 0:
        dt_snap = 1.0
    print("  snapshot spacing dt = %.6f (from diag total / frames)" % dt_snap)

    vol = dx ** 3
    dA = dx ** 2

    PA = momentum_x(uA, vA, udA, vdA, dx)
    PB = momentum_x(uB, vB, udB, vdB, dx)
    print("  total P_x: frame %d = %.10e   frame %d = %.10e"
          % (fa, PA.sum() * vol, fb, PB.sum() * vol))

    # mid-frame fields for the flux (centred in time, matching the centred dP/dt)
    uM = [0.5 * (uA[a] + uB[a]) for a in range(NF)]
    vM = [0.5 * (vA[a] + vB[a]) for a in range(NF)]
    udM = [0.5 * (udA[a] + udB[a]) for a in range(NF)]
    vdM = [0.5 * (vdA[a] + vdB[a]) for a in range(NF)]

    Flink = link_flux_x(uM, vM, udM, vdM, dx)
    Fnaive = naive_flux_x(uM, vM, udM, vdM, dx)
    PN = pn_residual_x(uM, vM, dx)

    dPdt = (PB - PA) / dt_snap                       # per-cell rate

    xs = -L + np.arange(N) * dx
    # ---------------------------------------------------------------- TEST 1
    # PURE ALGEBRA, no time derivative, no integrator, no EOM. The claim is
    #   sum_slab [ (lap f)(Dc f) + pi(Dc pi) - m^2 f(Dc f) ]  ==  F(right)-F(left)
    # If this fails the flux expression is wrong. If it passes, any discrepancy
    # with dP/dt is about the TIME integrator, which is a different question.
    def lap(f):
        return (sum(R(f, -1, d) + R(f, 1, d) for d in range(3)) - 6 * f) / dx ** 2

    def dcx(f):
        return (R(f, -1, 0) - R(f, 1, 0)) / (2 * dx)

    A = np.zeros_like(uM[0])
    for fld, fdot in ((uM, udM), (vM, vdM)):
        for a in range(NF):
            f, pi = fld[a], fdot[a]
            A += lap(f) * dcx(f) + pi * dcx(pi) - M2 * f * dcx(f)
    print()
    print("  TEST 1 -- PURE ALGEBRA (no dt, no EOM): does the link flux reproduce")
    print("  the summand exactly?   sum_slab A_n  vs  F(right) - F(left)")
    print("  %8s %18s %18s %12s" % ("X", "sum A (slab)", "F(R)-F(L)", "rel err"))
    e1 = []
    for i in range(4, N - 4, max(1, N // 12)):
        lhs = A[:i, :, :].sum() * vol
        rhs = (Flink[i - 1, :, :].sum() - Flink[-1, :, :].sum()) * dA
        sc = max(abs(lhs), abs(rhs), 1e-30)
        e = abs(rhs - lhs) / sc
        e1.append(e)
        print("  %8.3f %18.10e %18.10e %12.2e" % (xs[i], lhs, rhs, e))
    e1 = np.array(e1)
    print("  max relative error = %.3e   %s" % (e1.max(),
          "PASS -- the flux expression is exact on real data"
          if e1.max() < 1e-8 else "*** FAIL -- flux expression is WRONG ***"))

    print()
    print("  SLAB BALANCE   d/dt P_x(0 <= x < X)  vs  [F(X) - F(wrap)] + PN(slab)")
    print("  NOTE: the box is PERIODIC, so the region x < X has TWO boundaries --")
    print("  the plane at X and the wrap at x = -L. Both must be counted; using")
    print("  only the plane is an error of the same size as the answer.")
    print("  A correct discrete flux makes the residual vanish at EVERY X.")
    print()
    print("  %8s %16s %16s %10s %16s %10s" %
          ("X", "dP/dt", "link -F+PN", "rel err", "naive -F", "rel err"))
    rows = []
    for i in range(4, N - 4, max(1, N // 24)):
        lhs = dPdt[:i, :, :].sum() * vol
        rhs_link = ((Flink[i - 1, :, :].sum() - Flink[-1, :, :].sum()) * dA
                    + PN[:i, :, :].sum() * vol)
        rhs_naive = (Fnaive[i - 1, :, :].sum() - Fnaive[-1, :, :].sum()) * dA
        sc = max(abs(lhs), 1e-30)
        e_l = abs(rhs_link - lhs) / sc
        e_n = abs(rhs_naive - lhs) / sc
        rows.append((xs[i], lhs, rhs_link, e_l, rhs_naive, e_n))
        print("  %8.3f %16.6e %16.6e %9.2e %16.6e %9.2e"
              % (xs[i], lhs, rhs_link, e_l, rhs_naive, e_n))

    el = np.array([r[3] for r in rows])
    en = np.array([r[5] for r in rows])
    print()
    print("  median relative error   link  = %.3e" % np.median(el))
    print("                          naive = %.3e" % np.median(en))
    print("  max    relative error   link  = %.3e" % el.max())
    print("                          naive = %.3e" % en.max())
    if np.median(el) < np.median(en):
        print("  -> the link flux is better by %.1fx (median)"
              % (np.median(en) / max(np.median(el), 1e-30)))
    else:
        print("  -> NO improvement. The derivation does not close the balance on")
        print("     this data, and D8b's attribution of the 4.4%% defect FAILS.")
    print()
    print("  PN residual over the whole box: %.6e  (should be ~0 by D8b Part 2b)"
          % (PN.sum() * vol))
    print("  total dP_x/dt over the box    : %.6e  (periodic -> should match PN)"
          % (dPdt.sum() * vol))


if __name__ == "__main__":
    main()


def n7_flatness():
    """TEST 3 -- the N7 defect, addressed directly and with no time derivative.

    N7 reported the plane-to-plane force varying 1.4-8.4% between surfaces one
    unit apart in the gap between two balls, where a divergence-free flux must
    be flat. That is a statement about the SURFACE OBJECT alone, so it can be
    tested on a single frame. Compares the naive cell-centred T^xx (what
    sfa_momentum.c integrates) against the exact link flux, MATTER CHANNEL ONLY
    -- N7 is gauged and the Maxwell part is reported separately as Fx_gauge.
    """
    path = "/space/scp/v86/n7/n7_N112.sfa"
    if not os.path.exists(path):
        print("  N7 sfa missing; skipped")
        return
    (u, v, ud, vd), N, L, nfr = load(path, 0)
    dx = 2.0 * L / (N - 1)
    dA = dx * dx
    print()
    print("=" * 92)
    print("TEST 3 -- N7 plane-to-plane flatness (matter channel), frame 0")
    print("=" * 92)
    print("  %s  N=%d L=%g dx=%.6f" % (os.path.basename(path), N, L, dx))
    Fl = link_flux_x(u, v, ud, vd, dx)
    Fn = naive_flux_x(u, v, ud, vd, dx)
    PN = pn_residual_x(u, v, dx)
    vol = dx ** 3
    xs = -L + np.arange(N) * dx
    sel = np.where(np.abs(xs) <= 2.01)[0]
    print("  The premise 'a divergence-free flux must be flat in the gap' is NOT")
    print("  exactly right on a lattice: the potential term does not telescope, so")
    print("  its momentum transfer is a BULK residual, not a surface flux. The")
    print("  X-independent object is therefore F(X) + sum_{x<X} PN, not F(X).")
    print()
    print("  %8s %16s %16s %18s" % ("x", "naive centre", "link", "link + bulk PN"))
    a, b, c = [], [], []
    for i in sel:
        va = Fn[i, :, :].sum() * dA
        vb = Fl[i, :, :].sum() * dA
        vc = vb + PN[:i, :, :].sum() * vol
        a.append(va); b.append(vb); c.append(vc)
        print("  %8.4f %16.8e %16.8e %18.8e" % (xs[i], va, vb, vc))
    a, b, c = np.array(a), np.array(b), np.array(c)

    def spread(z):
        m = np.mean(np.abs(z))
        return 100.0 * (z.max() - z.min()) / m if m > 0 else float("nan")
    print()
    print("  plane-to-plane spread over |x| <= 2 (0 required for a")
    print("  divergence-free flux in the gap):")
    print("    naive cell-centred      : %8.3f %%" % spread(a))
    print("    exact link flux         : %8.3f %%" % spread(b))
    print("    link flux + bulk PN     : %8.3f %%" % spread(c))
    print()
    print("    link vs naive      : %.2fx flatter" % (spread(a) / max(spread(b), 1e-12)))
    print("    link+PN vs naive   : %.2fx flatter" % (spread(a) / max(spread(c), 1e-12)))


if "--n7" in sys.argv:
    n7_flatness()
