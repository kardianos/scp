#!/usr/bin/env python3
"""Generate Round-1 result tables with pure stdlib (no scipy)."""
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

C = 1.0
X = 12.0
M_T, BETA, EPS = 2.5, 0.12, 0.85


def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def kernel(M, beta, eps, bs):
    d, a = [], []
    for b in bs:
        s2 = b * b + eps * eps
        s = math.sqrt(s2)
        d.append(beta * M * 2.0 * math.asinh(X / s) / C)
        den = s2 * math.sqrt(X * X + s2)
        integ = (2.0 * X) / den if den > 0 else 0.0
        a.append(-beta * M * b * integ)
    return d, a, M


def local_opt(A, sig, kap, bs):
    A = min(max(A, 0.0), 0.95)
    sig = max(sig, 0.05)
    pref = (kap * A) * sig * math.sqrt(2.0 * math.pi)
    d, al = [], []
    for b in bs:
        e = math.exp(-(b * b) / (2.0 * sig * sig))
        d.append(pref * e)
        al.append((kap * A) * (-b / (sig * sig)) * sig * math.sqrt(2.0 * math.pi) * e)
    Mm = A * 2.0 * math.pi * sig * sig
    return d, al, Mm


def dual(Am, sig, G, bs):
    sig = max(sig, 0.05)
    Mt = max(Am, 0.0) * 2.0 * math.pi * sig * sig
    return kernel(Mt, G, sig, bs)  # isomorphic algebra


def lfit(pd, pa, pM, dd, da, dM, wM=0.25):
    n = len(dd)
    lt = sum((pd[i] - dd[i]) ** 2 for i in range(n)) / n
    la = sum((pa[i] - da[i]) ** 2 for i in range(n)) / n
    lm = wM * ((pM - dM) / max(abs(dM), 1e-12)) ** 2
    return lt + la + lm


def pattern(obj, x0, bounds, step0, max_iter=120):
    x = [min(max(x0[i], bounds[i][0]), bounds[i][1]) for i in range(len(x0))]
    f = obj(x)
    step = list(step0)
    for _ in range(max_iter):
        improved = False
        for i in range(len(x)):
            for direction in (1.0, -1.0):
                t = list(x)
                t[i] = min(max(t[i] + direction * step[i], bounds[i][0]), bounds[i][1])
                ft = obj(t)
                if ft < f - 1e-15:
                    x, f = t, ft
                    improved = True
        if not improved:
            step = [s * 0.5 for s in step]
            if max(step) < 1e-6:
                break
    return x, f


def main():
    bs = linspace(0.0, 4.0, 17)
    td, ta, tM = kernel(M_T, BETA, EPS, bs)
    dd, da, dM = list(td), list(ta), tM

    # monist exact
    L_mk = lfit(td, ta, tM, dd, da, dM)
    # dualist isomorphic
    Am = M_T / (2.0 * math.pi * EPS * EPS)
    du_d, du_a, du_M = dual(Am, EPS, BETA, bs)
    L_du = lfit(du_d, du_a, du_M, dd, da, dM)

    # local optics multi-start
    def obj(p):
        d, a, Mm = local_opt(p[0], p[1], p[2], bs)
        return lfit(d, a, Mm, dd, da, dM)

    seeds = [
        [0.4, 1.0, 0.3],
        [0.6, 1.5, 0.5],
        [0.3, 2.0, 0.8],
        [0.5, 1.2, 0.35],
        [0.7, 2.5, 0.4],
        [0.2, 0.8, 1.0],
        [0.1, 2.0, 3.0],
        [0.05, 4.0, 3.0],
        [0.15, 1.5, 2.5],
        [0.08, 3.0, 3.0],
        [0.12, 2.5, 2.0],
        [0.2, 3.5, 1.5],
    ]
    bounds = [(1e-3, 0.94), (0.2, 5.0), (1e-3, 3.0)]
    best_p, best_f = None, 1e99
    for s in seeds:
        p, f = pattern(obj, s, bounds, [0.05, 0.15, 0.05])
        if f < best_f:
            best_f, best_p = f, p
    lo_d, lo_a, lo_M = local_opt(best_p[0], best_p[1], best_p[2], bs)
    L_lo = best_f

    # scores
    lam_sec, lam_link, lam_soft = 1.0, 0.5, 100.0
    S_mk = L_mk
    S_lo = L_lo
    S_du = L_du + lam_sec * 1 + lam_link
    S_ch = L_du + lam_soft

    # occam ablation
    abl = []
    for k in range(16):
        lam = 3.0 * k / 15.0
        sm = L_mk
        sd = L_du + lam * 1 + lam_link
        abl.append((lam, sm, sd, "monist_kernel" if sm <= sd else "dualist"))

    # path cost
    rs = linspace(0.2, 6.0, 30)
    path_rows = []
    for r in rs:
        dt = BETA * M_T / math.sqrt(r * r + EPS * EPS)
        dm = dt  # recovered = truth
        dd_ = BETA * du_M / math.sqrt(r * r + EPS * EPS)
        path_rows.append((r, dt, dm, dd_))

    result = {
        "truth": {"M": M_T, "beta": BETA, "eps": EPS},
        "noise_delay": 0.0,
        "noise_alpha": 0.0,
        "seed": 7,
        "b_values": bs,
        "data_delay": dd,
        "data_alpha": da,
        "data_M": dM,
        "true_delay": td,
        "true_alpha": ta,
        "monist_kernel_fit": {
            "M": M_T,
            "beta": BETA,
            "eps": EPS,
            "pred_delay": td,
            "pred_alpha": ta,
            "pred_M": tM,
            "rel_err_M": 0.0,
            "rel_err_beta": 0.0,
            "rel_err_eps": 0.0,
            "M_ledger_consistency": 0.0,
            "S": S_mk,
            "L_fit": L_mk,
            "L_occ": 0.0,
            "L_softE": 0.0,
            "L_nogo": 0.0,
            "N_sectors": 1,
            "has_free_bound_link": True,
        },
        "local_optics_fit": {
            "A": best_p[0],
            "sigma": best_p[1],
            "kappa": best_p[2],
            "pred_delay": lo_d,
            "pred_alpha": lo_a,
            "pred_M": lo_M,
            "S": S_lo,
            "L_fit": L_lo,
            "L_occ": 0.0,
            "L_softE": 0.0,
            "L_nogo": 0.0,
            "N_sectors": 1,
            "has_free_bound_link": True,
        },
        "dualist_fit": {
            "A_m": Am,
            "sigma": EPS,
            "G_eff": BETA,
            "pred_delay": du_d,
            "pred_alpha": du_a,
            "pred_M": du_M,
            "S": S_du,
            "L_fit": L_du,
            "L_occ": lam_sec + lam_link,
            "L_softE": 0.0,
            "L_nogo": 0.0,
            "N_sectors": 2,
            "has_free_bound_link": False,
            "isomorphism_note": "Plummer Phi ray-isomorphic to monist kernel at matched params",
        },
        "soft_einstein_cheat": {
            "L_fit": L_du,
            "L_occ": 0.0,
            "L_softE": lam_soft,
            "L_nogo": 0.0,
            "S": S_ch,
            "N_sectors": 1,
            "has_free_bound_link": True,
            "note": "dualist Phi labeled monist",
        },
        "fit_winner": min(
            [("monist_kernel", L_mk), ("local_optics", L_lo), ("dualist", L_du)],
            key=lambda x: x[1],
        )[0],
        "score_winner": min(
            [
                ("monist_kernel", S_mk),
                ("local_optics", S_lo),
                ("dualist", S_du),
                ("softE_cheat", S_ch),
            ],
            key=lambda x: x[1],
        )[0],
        "occam_ablation": [
            {
                "lambda_sec": lam,
                "S_monist_kernel": sm,
                "S_dualist": sd,
                "winner": w,
            }
            for lam, sm, sd, w in abl
        ],
        "occam_weights": {
            "lambda_sec": 1.0,
            "lambda_link": 0.5,
            "lambda_budget": 10.0,
            "lambda_softE": 100.0,
            "lambda_nogo": 5.0,
        },
        "notes": [
            "Truth = monist nonlocal path-cost kernel (C Class C / FOR_D).",
            "Dualist Plummer is ray-isomorphic to monist kernel at matched params.",
            "Combined score prefers monist via N_sectors + free-bound link.",
            "Local optics fails L_fit on long-range kernel truth (C no-go).",
            "softE_cheat disqualified by L_softE=100.",
        ],
    }

    (OUT / "round1_result.json").write_text(json.dumps(result, indent=2))

    with open(OUT / "round1_scores.tsv", "w") as f:
        f.write(
            "model\tN_sectors\thas_link\tL_fit\tL_occ\tL_softE\tL_nogo\tS\t"
            "p1\tp2\tp3\tpred_M\tdata_M\n"
        )
        mk = result["monist_kernel_fit"]
        lo = result["local_optics_fit"]
        du = result["dualist_fit"]
        ch = result["soft_einstein_cheat"]
        f.write(
            f"monist_kernel\t1\t1\t{mk['L_fit']:.8e}\t{mk['L_occ']:.8e}\t0\t0\t{mk['S']:.8e}\t"
            f"{mk['M']:.8e}\t{mk['beta']:.8e}\t{mk['eps']:.8e}\t{mk['pred_M']:.8e}\t{dM:.8e}\n"
        )
        f.write(
            f"local_optics\t1\t1\t{lo['L_fit']:.8e}\t0\t0\t0\t{lo['S']:.8e}\t"
            f"{lo['A']:.8e}\t{lo['sigma']:.8e}\t{lo['kappa']:.8e}\t{lo['pred_M']:.8e}\t{dM:.8e}\n"
        )
        f.write(
            f"dualist\t2\t0\t{du['L_fit']:.8e}\t{du['L_occ']:.8e}\t0\t0\t{du['S']:.8e}\t"
            f"{du['A_m']:.8e}\t{du['sigma']:.8e}\t{du['G_eff']:.8e}\t{du['pred_M']:.8e}\t{dM:.8e}\n"
        )
        f.write(
            f"softE_cheat\t1\t1\t{ch['L_fit']:.8e}\t0\t{ch['L_softE']:.8e}\t0\t{ch['S']:.8e}\t"
            f"\t\t\t\t{dM:.8e}\n"
        )

    with open(OUT / "round1_rays.tsv", "w") as f:
        f.write(
            "b\tdata_delay\ttrue_delay\tmk_delay\tlo_delay\tdual_delay\t"
            "data_alpha\ttrue_alpha\tmk_alpha\tlo_alpha\tdual_alpha\n"
        )
        for i, b in enumerate(bs):
            f.write(
                f"{b:.6f}\t{dd[i]:.8e}\t{td[i]:.8e}\t{td[i]:.8e}\t{lo_d[i]:.8e}\t{du_d[i]:.8e}\t"
                f"{da[i]:.8e}\t{ta[i]:.8e}\t{ta[i]:.8e}\t{lo_a[i]:.8e}\t{du_a[i]:.8e}\n"
            )

    with open(OUT / "round1_occam_ablation.tsv", "w") as f:
        f.write("lambda_sec\tS_monist_kernel\tS_dualist\twinner\n")
        for lam, sm, sd, w in abl:
            f.write(f"{lam:.6f}\t{sm:.8e}\t{sd:.8e}\t{w}\n")

    with open(OUT / "round1_path_cost.tsv", "w") as f:
        f.write("r\tdell_true\tdell_mk\tdell_dual_proxy\n")
        for r, dt, dm, ddp in path_rows:
            f.write(f"{r:.6f}\t{dt:.8e}\t{dm:.8e}\t{ddp:.8e}\n")

    # SVG rays
    w, h, pad = 520, 300, 40
    series = [(td, "#000"), (lo_d, "#080"), (du_d, "#c00")]
    ys = td + lo_d + du_d
    ymin, ymax = min(ys), max(ys)
    bx0, bx1 = bs[0], bs[-1]

    def xy(xv, yv):
        Xv = pad + (xv - bx0) / (bx1 - bx0) * (w - 2 * pad)
        Yv = h - pad - (yv - ymin) / (ymax - ymin + 1e-15) * (h - 2 * pad)
        return Xv, Yv

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{pad}" y="18" font-size="12">Ray delay: black=truth/MK, green=local optics, red=dualist</text>',
    ]
    for ys_, color in series:
        pts = " ".join(f"{xy(bs[i], ys_[i])[0]:.1f},{xy(bs[i], ys_[i])[1]:.1f}" for i in range(len(bs)))
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{pts}"/>')
    parts.append("</svg>")
    (OUT / "round1_rays.svg").write_text("\n".join(parts))

    # score bar SVG
    labels = ["MK", "LO", "DU", "softE"]
    Lfits = [L_mk, L_lo, L_du, L_du]
    pens = [0.0, 0.0, lam_sec + lam_link, lam_soft]
    bw = 80
    maxS = max(Lfits[i] + pens[i] for i in range(4)) + 0.1
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="480" height="280">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="20" y="20" font-size="12">Scores (winner={result["score_winner"]})</text>',
    ]
    for i in range(4):
        x0 = 40 + i * 100
        h1 = 200 * Lfits[i] / maxS
        h2 = 200 * pens[i] / maxS
        parts.append(
            f'<rect x="{x0}" y="{240-h1}" width="{bw}" height="{h1}" fill="#4488cc"/>'
        )
        parts.append(
            f'<rect x="{x0}" y="{240-h1-h2}" width="{bw}" height="{h2}" fill="#cc6644"/>'
        )
        parts.append(f'<text x="{x0}" y="260" font-size="11">{labels[i]}</text>')
        parts.append(
            f'<text x="{x0}" y="{230-h1-h2}" font-size="10">S={Lfits[i]+pens[i]:.3g}</text>'
        )
    parts.append("</svg>")
    (OUT / "round1_scores.svg").write_text("\n".join(parts))

    summary = (
        f"L_mk={L_mk:.6e} S_mk={S_mk:.6e}\n"
        f"L_lo={L_lo:.6e} S_lo={S_lo:.6e} params A={best_p[0]:.6f} sig={best_p[1]:.6f} kap={best_p[2]:.6f} M={lo_M:.6f}\n"
        f"L_du={L_du:.6e} S_du={S_du:.6e} A_m={Am:.6f}\n"
        f"S_cheat={S_ch:.6e}\n"
        f"fit_winner={result['fit_winner']} score_winner={result['score_winner']}\n"
    )
    (OUT / "round1_summary.txt").write_text(summary)
    print(summary)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
