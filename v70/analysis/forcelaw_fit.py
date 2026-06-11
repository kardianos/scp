#!/usr/bin/env python3
"""forcelaw_fit.py — measure a_rel(D) per run from sfa_qball_track TSVs and
fit the force-law exponent n in a_C(D) = C / D^n.

Per run: take the two most massive clusters per frame, D(t) = min-image
separation, fit D(t) = D0 + v0*(t-t0) + (a_rel/2)*(t-t0)^2 for t >= t_min
(v0 absorbs the seed settling kick; a_rel is the physics). Residual
bootstrap gives the error bar.

Combining runs: the Coulomb estimator at each D is
    a_C(D) = ( a_rel[same](D) - a_rel[opp](D) ) / 2
which cancels additive systematics (lattice drift, sponge asymmetry) and
leaves Coulomb (sign flips) plus any phase-coherent contact residual
(negligible beyond D~26: F_t ~ e^{-0.564 D}).

Prediction (no free parameters): a_rel = 2 g^2 Q^2 / (4 pi D^2 M)
with g=0.05, Q per-ball Noether charge, M per-ball energy.

Usage:
  forcelaw_fit.py --pairs same1.tsv:opp1.tsv:D1 [same2.tsv:opp2.tsv:D2 ...]
                  [--single run.tsv:D ...] [--tmin 20] [--Q 315.4] [--M 456.3]
                  [--nboot 2000] [--plot out.png]

  --pairs   same-charge and opposite-charge tracker TSVs at nominal D
  --single  same-charge-only runs (reported, not used in the exponent fit)
"""
import sys
import argparse
import numpy as np

G = 0.05


def load_sep(tsv_path):
    d = np.genfromtxt(tsv_path, names=True)
    frames = np.unique(d["frame"]).astype(int)
    t, sep = [], []
    for fi in frames:
        m = d["frame"] == fi
        if m.sum() < 2:
            continue
        order = np.argsort(d["mass"][m])[::-1]
        cx = d["cx"][m][order[:2]]
        cy = d["cy"][m][order[:2]]
        cz = d["cz"][m][order[:2]]
        dx = cx[0] - cx[1]
        dy = cy[0] - cy[1]
        dz = cz[0] - cz[1]
        t.append(d["t"][m][0])
        sep.append(np.sqrt(dx * dx + dy * dy + dz * dz))
    return np.array(t), np.array(sep)


def fit_accel(t, sep, tmin, nboot=2000, rng=None):
    """quadratic fit with v0 term; returns a_rel, sigma_a (bootstrap)."""
    m = t >= tmin
    if m.sum() < 4:
        raise ValueError(f"only {m.sum()} frames with t >= {tmin}")
    tt, ss = t[m] - t[m][0], sep[m]
    A = np.vstack([np.ones_like(tt), tt, 0.5 * tt * tt]).T
    coef, *_ = np.linalg.lstsq(A, ss, rcond=None)
    resid = ss - A @ coef
    if rng is None:
        rng = np.random.default_rng(20260611)
    boots = []
    for _ in range(nboot):
        sb = A @ coef + rng.choice(resid, size=len(resid), replace=True)
        cb, *_ = np.linalg.lstsq(A, sb, rcond=None)
        boots.append(cb[2])
    return coef[2], float(np.std(boots)), resid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", nargs="*", default=[],
                    help="same.tsv:opp.tsv:D")
    ap.add_argument("--single", nargs="*", default=[],
                    help="same.tsv:D (reported only)")
    ap.add_argument("--tmin", type=float, default=20.0)
    ap.add_argument("--Q", type=float, default=315.4,
                    help="per-ball Noether charge (diag Q_phi/2 for pairs)")
    ap.add_argument("--M", type=float, default=456.3,
                    help="per-ball total energy (gscan E_total at omega=1.42)")
    ap.add_argument("--nboot", type=int, default=2000)
    args = ap.parse_args()

    coef_pred = 2 * G * G * args.Q ** 2 / (4 * np.pi * args.M)
    print(f"# Coulomb prediction: a_rel(D) = {coef_pred:.4g} / D^2 "
          f"(g={G}, Q={args.Q}, M={args.M})")

    print(f"\n{'run':>28s} {'D_nom':>6s} {'D_fit':>7s} {'a_rel':>11s} "
          f"{'±':>9s} {'a_pred':>11s} {'ratio':>6s}")

    def report(tsv, Dnom, sign_label):
        t, sep = load_sep(tsv)
        a, sa, resid = fit_accel(t, sep, args.tmin, args.nboot)
        Dfit = float(np.mean(sep[t >= args.tmin]))
        apred = coef_pred / Dfit ** 2
        name = tsv.split("/")[-1].replace("_track.tsv", "")
        print(f"{name:>28s} {Dnom:6.0f} {Dfit:7.2f} {a:+11.3e} "
              f"{sa:9.2e} {apred:11.3e} {a/apred:+6.2f}  rms_resid={np.std(resid):.4f}")
        return Dfit, a, sa

    results = []
    for spec in args.pairs:
        same_tsv, opp_tsv, D = spec.rsplit(":", 2)[0], spec.rsplit(":", 2)[1], float(spec.rsplit(":", 2)[2])
        Ds, a_s, e_s = report(same_tsv, D, "same")
        Do, a_o, e_o = report(opp_tsv, D, "opp")
        aC = 0.5 * (a_s - a_o)
        eC = 0.5 * np.hypot(e_s, e_o)
        Dm = 0.5 * (Ds + Do)
        results.append((Dm, aC, eC))
        print(f"{'-> Coulomb estimator':>28s} {D:6.0f} {Dm:7.2f} {aC:+11.3e} "
              f"{eC:9.2e} {coef_pred/Dm**2:11.3e} {aC/(coef_pred/Dm**2):+6.2f}")

    for spec in args.single:
        tsv, D = spec.rsplit(":", 1)[0], float(spec.rsplit(":", 1)[1])
        report(tsv, D, "single")

    if len(results) >= 3:
        D = np.array([r[0] for r in results])
        a = np.array([r[1] for r in results])
        e = np.array([r[2] for r in results])
        ok = a > 0
        if ok.sum() >= 3:
            w = 1.0 / (e[ok] / a[ok]) ** 2
            # weighted log-log fit a = C / D^n
            X = np.log(D[ok]); Y = np.log(a[ok])
            W = np.diag(w)
            A = np.vstack([np.ones_like(X), -X]).T
            cov = np.linalg.inv(A.T @ W @ A)
            beta = cov @ A.T @ W @ Y
            # parameter errors via bootstrap over points with gaussian noise
            rng = np.random.default_rng(7)
            ns = []
            for _ in range(4000):
                Yb = np.log(np.maximum(a[ok] + rng.normal(0, e[ok]), 1e-12))
                bb = cov @ A.T @ W @ Yb
                ns.append(bb[1])
            n_err = float(np.std(ns))
            Cfit = float(np.exp(beta[0]))
            print(f"\n# exponent fit over {ok.sum()} Coulomb estimator points:")
            print(f"#   n     = {beta[1]:.3f} ± {n_err:.3f}   (Coulomb: 2.000)")
            print(f"#   C     = {Cfit:.4g}            (predicted {coef_pred:.4g}, "
                  f"ratio {Cfit/coef_pred:.3f})")


if __name__ == "__main__":
    main()
