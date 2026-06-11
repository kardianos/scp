#!/usr/bin/env python3
"""phase_tilt.py — test the v67 DEBROGLIE "hbar_eff = Q" fingerprint on
banked boost SFAs (v68 relativity ladder, bs03/05/07).

For a boosted Q-ball Phi_a = f(r') e^{i gamma omega (t - v x)} the in-core
phase tilt is k = gamma*omega*v. The v67 claim is that the conserved charge
plays the role of Planck's constant: with p = gamma*M*v (M = E0 the rest
energy), hbar_eff = p/k = M/omega = Q(1+eps), eps = (E - omega Q)/(omega Q)
from the Pohozaev identity. The discriminating alternative (a "true quantum
composite" with hbar=1) would be k = p itself — a factor ~Q larger,
sub-voxel, unrepresentable; so the test is quantitative agreement of the
measured k with gamma*omega*v using the MEASURED v, and the resulting
hbar_eff = p/k compared against Q.

Inputs: <prefix>.line.tsv files produced by
  sfa_slice bsXX.sfa --axis z --line x --out <prefix>
The lineout follows the per-frame rho2-max voxel, so each frame's line
passes through the moving core.

Ball constants (v66 THEORY table, omega=1.39 branch): E0=691.9, Q=482.2.
"""
import sys
import numpy as np

E0 = 691.9          # rest energy of the omega=1.39 ball (v66 THEORY section 4)
Q0 = 482.2          # its Noether charge
OMEGA = 1.39        # internal rotation frequency (rest frame)
WINDOW_FRAC = 0.2   # fit phase only where rho2 > frac * max (the core)


def analyze(prefix):
    d = np.genfromtxt(prefix + ".line.tsv", names=True)
    frames = np.unique(d["frame"]).astype(int)
    L = -d["coord"].min()

    # ---- measure v from core position vs t (unwrap periodic jumps) ----
    tlist, xlist = [], []
    for fi in frames:
        m = d["frame"] == fi
        rho2 = d["rho2"][m]
        coord = d["coord"][m]
        # rho2-weighted centroid near max, periodic-safe via angle average
        ang = (coord + L) / (2 * L) * 2 * np.pi
        w = np.where(rho2 > WINDOW_FRAC * rho2.max(), rho2, 0.0)
        cx = np.angle(np.sum(w * np.exp(1j * ang))) / (2 * np.pi) * 2 * L - L
        tlist.append(d["t"][m][0])
        xlist.append(cx)
    t = np.array(tlist)
    x = np.unwrap(np.array(xlist), discont=L, period=2 * L)
    vfit = np.polyfit(t, x, 1)
    v_meas = vfit[0]
    gamma_meas = 1.0 / np.sqrt(1.0 - v_meas ** 2)

    # ---- per-frame phase-tilt fit ----
    rows = []
    for fi in frames:
        m = d["frame"] == fi
        coord = d["coord"][m]
        rho2 = d["rho2"][m]
        # contiguous core window around the max
        imax = int(np.argmax(rho2))
        thr = WINDOW_FRAC * rho2[imax]
        lo = imax
        while lo > 0 and rho2[lo - 1] > thr:
            lo -= 1
        hi = imax
        while hi < len(rho2) - 1 and rho2[hi + 1] > thr:
            hi += 1
        if hi - lo < 5:
            continue
        ks = []
        for a in "xyz":
            u = d[f"phi_{a}"][m][lo:hi + 1]
            v = d[f"phiim_{a}"][m][lo:hi + 1]
            amp2 = u * u + v * v
            ok = amp2 > 1e-4
            if ok.sum() < 5:
                continue
            th = np.unwrap(np.arctan2(v[ok], u[ok]))
            cc = coord[lo:hi + 1][ok]
            k = np.polyfit(cc, th, 1)[0]
            ks.append(k)
        if ks:
            rows.append((fi, d["t"][m][0], np.mean(ks), np.std(ks)))

    return v_meas, gamma_meas, rows


def main(prefixes):
    print(f"# ball constants: E0={E0}, Q={Q0}, omega={OMEGA} (v66 THEORY)")
    print(f"# hbar_eff(theory) = E0/omega = {E0/OMEGA:.2f} = Q*(1+eps), "
          f"eps = {(E0/OMEGA - Q0)/Q0:.4f}")
    hdr = (f"{'run':8s} {'v_meas':>8s} {'gamma':>7s} {'k(t=0)':>8s} "
           f"{'k(t>0)':>8s} {'±':>6s} {'g*w*v':>8s} {'k/gwv':>7s} "
           f"{'hbar_eff':>9s} {'/Q':>6s}")
    print(hdr)
    for p in prefixes:
        v, g, rows = analyze(p)
        k0 = [r[2] for r in rows if r[1] == 0.0]
        kd = [r[2] for r in rows if r[1] > 0.0]
        k0 = abs(np.mean(k0)) if k0 else float("nan")
        kmean = abs(np.mean(kd))
        kstd = np.std([abs(r[2]) for r in rows if r[1] > 0.0])
        kpred = g * OMEGA * abs(v)
        p_momentum = g * E0 * abs(v)
        hbar_eff = p_momentum / kmean
        name = p.split("/")[-1]
        print(f"{name:8s} {v:8.4f} {g:7.4f} {k0:8.4f} {kmean:8.4f} "
              f"{kstd:6.4f} {kpred:8.4f} {kmean/kpred:7.4f} "
              f"{hbar_eff:9.2f} {hbar_eff/Q0:6.4f}")
        for fi, t, k, ks in rows:
            print(f"    frame {fi} t={t:5.1f}  k={k:+8.4f} ± {ks:.4f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
