#!/usr/bin/env python3
"""
v76 Approach D — reverse-numeric demo (Round 1)

Aligned with C FOR_D / A FOR_D:
  - Loss targets path-cost / ray observables (ℓ-class), not δρ_free as primary.
  - Monist truth = compact lock + nonlocal free-response kernel (Class C).
  - Dualist adversary = matter density + separate long-range gravity potential
    (S6+S9 / S2 scaffolding): two sectors, no free–bound identity.
  - Local-optics monist (n∝ρ_bound, pointwise budget) is a *wrong class* fit
    when truth is long-range kernel — scores poorly on L_fit (C no-go).
  - Combined score penalizes N_sectors, missing free–bound link, soft-Einstein
    cheat, and non-integrable free-energy 1/r under claimed local budget.

Pure stdlib. No scp_sim. No soft-coded Einstein as monist win.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

C_LIGHT = 1.0
RHO0 = 1.0
X_RAY = 12.0  # half-length of Born ray chart path


# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------
@dataclass
class OccamWeights:
    lambda_sec: float = 1.0
    lambda_link: float = 0.5
    lambda_budget: float = 10.0
    lambda_softE: float = 100.0
    lambda_nogo: float = 5.0  # local-budget + claims long-range ρ tail


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def linspace(a: float, b: float, n: int) -> List[float]:
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def trapz(y: Sequence[float], dx: float) -> float:
    if len(y) < 2:
        return 0.0
    s = 0.5 * (y[0] + y[-1])
    for i in range(1, len(y) - 1):
        s += y[i]
    return s * dx


class LCG:
    """Minimal deterministic RNG so offline recomputation matches (no system MT)."""

    def __init__(self, seed: int = 7):
        self.state = seed & 0xFFFFFFFF

    def random(self) -> float:
        # Numerical Recipes LCG
        self.state = (1664525 * self.state + 1013904223) & 0xFFFFFFFF
        return self.state / 4294967296.0

    def gauss(self, sigma: float) -> float:
        u1 = max(self.random(), 1e-15)
        u2 = self.random()
        return sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def mse(a: Sequence[float], b: Sequence[float]) -> float:
    n = len(a)
    return sum((a[i] - b[i]) ** 2 for i in range(n)) / max(n, 1)


# ---------------------------------------------------------------------------
# Monist kernel forward (preferred Class C — one sector)
# Path-cost excess: δℓ(r) = β * M / sqrt(r² + ε²)
# with M = E_star / c². This IS free-response geometry, not T→G.
# ---------------------------------------------------------------------------
def kernel_delay_alpha(
    M: float, beta: float, eps: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    """
    Born delay: ∫ δℓ dx = β M ∫ dx / sqrt(x² + b² + ε²)
              = β M * 2 * asinh(X / s),  s=sqrt(b²+ε²)

    Born alpha ~ ∫ ∂_y δℓ dx:
      ∂_y δℓ = β M * (-1/2)*2 y / (r²+ε²)^{3/2} = - β M y / (r²+ε²)^{3/2}
      ∫ ∂_y δℓ dx at y=b = -β M b ∫ dx / (x²+s²)^{3/2}
                         = -β M b * (2 X) / (s² * sqrt(X²+s²))
    """
    delays: List[float] = []
    alphas: List[float] = []
    X = X_RAY
    for b in b_values:
        s2 = b * b + eps * eps
        s = math.sqrt(s2)
        # delay
        delays.append(beta * M * 2.0 * math.asinh(X / s) / C_LIGHT)
        # alpha
        den = s2 * math.sqrt(X * X + s2)
        integ = (2.0 * X) / den if den > 0 else 0.0
        alphas.append(-beta * M * b * integ)
    # E_star = M c²; report M as bound mass ledger
    return delays, alphas, M


# ---------------------------------------------------------------------------
# Local-optics monist (wrong class for long-range truth; C no-go)
# n-1 = κ A exp(-r²/2σ²)/ρ0 ; compact support in practice
# ---------------------------------------------------------------------------
def local_optics_delay_alpha(
    A: float, sigma: float, kappa: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    A_eff = min(max(A, 0.0), 0.95 * RHO0)
    sigma = max(sigma, 0.05)
    pref = (kappa * A_eff / RHO0) * sigma * math.sqrt(2.0 * math.pi)
    delays, alphas = [], []
    for b in b_values:
        e = math.exp(-(b * b) / (2.0 * sigma * sigma))
        delays.append(pref * e / C_LIGHT)
        alphas.append(
            (kappa * A_eff / RHO0)
            * (-b / (sigma * sigma))
            * sigma
            * math.sqrt(2.0 * math.pi)
            * e
        )
    M = A_eff * 2.0 * math.pi * sigma * sigma / (C_LIGHT ** 2)
    return delays, alphas, M


# ---------------------------------------------------------------------------
# Dualist: matter blob + separate Plummer gravity (2 sectors)
# ---------------------------------------------------------------------------
def dualist_delay_alpha(
    A_m: float,
    sigma: float,
    G_eff: float,
    b_values: Sequence[float],
    n_x: int = 301,
) -> Tuple[List[float], List[float], float]:
    """
    Dualist Plummer gravity sector (independent of free-budget).

    Phi = -G_eff * M_tot / sqrt(r^2 + a^2), a=sigma.
    Closed Born integrals (same algebra as monist kernel with beta→G_eff):
      delay = G M * 2 asinh(X/s) / c
      alpha = -G M b * 2X / (s^2 sqrt(X^2+s^2))

    NOTE: ray map is *isomorphic* to monist kernel when G_eff=beta, a=eps,
    M_tot=M. Discrimination for D3 is then ontology Occam, not pure L_fit.
    (n_x retained for API compat; analytical path used.)
    """
    _ = n_x
    sig = max(sigma, 0.05)
    M_tot = max(A_m, 0.0) * 2.0 * math.pi * sig * sig
    a2 = sig * sig
    X = X_RAY
    delays, alphas = [], []
    for b in b_values:
        s2 = b * b + a2
        s = math.sqrt(s2)
        delays.append(G_eff * M_tot * 2.0 * math.asinh(X / s) / C_LIGHT)
        den = s2 * math.sqrt(X * X + s2)
        integ = (2.0 * X) / den if den > 0 else 0.0
        alphas.append(-G_eff * M_tot * b * integ)
    return delays, alphas, M_tot / (C_LIGHT ** 2)


# ---------------------------------------------------------------------------
# Losses / scores
# ---------------------------------------------------------------------------
def l_fit(
    pred_d: Sequence[float],
    pred_a: Sequence[float],
    pred_M: float,
    data_d: Sequence[float],
    data_a: Sequence[float],
    data_M: float,
    w_t: float = 1.0,
    w_a: float = 1.0,
    w_M: float = 0.25,
) -> float:
    scale_M = max(abs(data_M), 1e-12)
    return (
        w_t * mse(pred_d, data_d)
        + w_a * mse(pred_a, data_a)
        + w_M * ((pred_M - data_M) / scale_M) ** 2
    )


def score_model(
    Lfit: float,
    *,
    n_sectors: int,
    has_free_bound_link: bool,
    claims_soft_einstein: bool,
    claims_local_budget_with_long_range_rho: bool,
    occ: OccamWeights,
) -> dict:
    Locc = occ.lambda_sec * max(n_sectors - 1, 0)
    if not has_free_bound_link:
        Locc += occ.lambda_link
    Lsoft = occ.lambda_softE if claims_soft_einstein else 0.0
    Lnogo = occ.lambda_nogo if claims_local_budget_with_long_range_rho else 0.0
    S = Lfit + Locc + Lsoft + Lnogo
    return {
        "S": S,
        "L_fit": Lfit,
        "L_occ": Locc,
        "L_softE": Lsoft,
        "L_nogo": Lnogo,
        "N_sectors": n_sectors,
        "has_free_bound_link": has_free_bound_link,
    }


# ---------------------------------------------------------------------------
# Pattern search optimizer
# ---------------------------------------------------------------------------
def pattern_search(
    obj: Callable[[List[float]], float],
    x0: List[float],
    bounds: List[Tuple[float, float]],
    step0: List[float],
    max_iter: int = 100,
    shrink: float = 0.5,
    min_step: float = 1e-6,
) -> Tuple[List[float], float]:
    x = [min(max(x0[i], bounds[i][0]), bounds[i][1]) for i in range(len(x0))]
    f = obj(x)
    step = list(step0)
    for _ in range(max_iter):
        improved = False
        for i in range(len(x)):
            for direction in (+1.0, -1.0):
                trial = list(x)
                trial[i] = min(max(trial[i] + direction * step[i], bounds[i][0]), bounds[i][1])
                ft = obj(trial)
                if ft < f - 1e-15:
                    x, f = trial, ft
                    improved = True
        if not improved:
            step = [s * shrink for s in step]
            if max(step) < min_step:
                break
    return x, f


def multi_start(
    obj: Callable[[List[float]], float],
    seeds: List[List[float]],
    bounds: List[Tuple[float, float]],
    step0: List[float],
) -> Tuple[List[float], float]:
    best_p, best_f = None, float("inf")
    for s0 in seeds:
        p, f = pattern_search(obj, s0, bounds, step0)
        if f < best_f:
            best_f, best_p = f, p
    assert best_p is not None
    return best_p, best_f


# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------
@dataclass
class Truth:
    M: float = 2.5
    beta: float = 0.12
    eps: float = 0.85


def run(
    truth: Truth | None = None,
    noise_delay: float = 0.0,
    noise_alpha: float = 0.0,
    seed: int = 7,
    occ: OccamWeights | None = None,
) -> dict:
    truth = truth or Truth()
    occ = occ or OccamWeights()
    rng = LCG(seed)
    b_values = linspace(0.0, 4.0, 17)

    true_d, true_a, true_M = kernel_delay_alpha(truth.M, truth.beta, truth.eps, b_values)
    data_d = [true_d[i] + rng.gauss(noise_delay) for i in range(len(b_values))]
    data_a = [true_a[i] + rng.gauss(noise_alpha) for i in range(len(b_values))]
    data_M = true_M  # unlock/bound ledger known (C: M = E_star/c²)

    # --- fit monist kernel (correct class) ---
    def obj_mk(p):
        M, beta, eps = p
        if M <= 0 or beta <= 0 or eps <= 0.05:
            return 1e9
        d, a, Mm = kernel_delay_alpha(M, beta, eps, b_values)
        return l_fit(d, a, Mm, data_d, data_a, data_M)

    p_mk, L_mk = multi_start(
        obj_mk,
        seeds=[
            [2.5, 0.12, 0.85],
            [2.0, 0.10, 1.0],
            [3.0, 0.15, 0.7],
            [1.5, 0.08, 1.2],
            [2.8, 0.14, 0.9],
            [2.2, 0.11, 0.6],
        ],
        bounds=[(0.1, 10.0), (1e-3, 2.0), (0.1, 3.0)],
        step0=[0.2, 0.02, 0.1],
    )
    d_mk, a_mk, M_mk = kernel_delay_alpha(p_mk[0], p_mk[1], p_mk[2], b_values)
    sc_mk = score_model(
        L_mk,
        n_sectors=1,
        has_free_bound_link=True,
        claims_soft_einstein=False,
        claims_local_budget_with_long_range_rho=False,
        occ=occ,
    )

    # --- fit local-optics monist (wrong class) ---
    def obj_lo(p):
        A, sig, kap = p
        if A <= 0 or sig <= 0.05 or kap <= 0:
            return 1e9
        d, a, Mm = local_optics_delay_alpha(A, sig, kap, b_values)
        return l_fit(d, a, Mm, data_d, data_a, data_M)

    p_lo, L_lo = multi_start(
        obj_lo,
        seeds=[
            [0.4, 1.0, 0.3],
            [0.6, 1.5, 0.5],
            [0.3, 2.0, 0.8],
            [0.5, 1.2, 0.35],
            [0.7, 2.5, 0.4],
            [0.2, 0.8, 1.0],
        ],
        bounds=[(1e-3, 0.94), (0.2, 5.0), (1e-3, 3.0)],
        step0=[0.05, 0.15, 0.05],
    )
    d_lo, a_lo, M_lo = local_optics_delay_alpha(p_lo[0], p_lo[1], p_lo[2], b_values)
    # Local optics + compact budget cannot make long-range 1/r ρ; if someone
    # claims it does, nogo fires. Here the *model class* is local optics for
    # long-range truth → high L_fit; we do NOT auto-apply nogo unless claimed.
    sc_lo = score_model(
        L_lo,
        n_sectors=1,
        has_free_bound_link=True,
        claims_soft_einstein=False,
        claims_local_budget_with_long_range_rho=False,
        occ=occ,
    )

    # --- dualist adversary ---
    def obj_du(p):
        Am, sig, G = p
        if Am <= 0 or sig <= 0.05 or G <= 0:
            return 1e9
        d, a, Mm = dualist_delay_alpha(Am, sig, G, b_values)
        return l_fit(d, a, Mm, data_d, data_a, data_M)

    p_du, L_du = multi_start(
        obj_du,
        seeds=[
            [0.3, 1.0, 0.05],
            [0.4, 0.9, 0.08],
            [0.2, 1.2, 0.04],
            [0.5, 0.85, 0.06],
            [0.35, 1.1, 0.1],
            [0.25, 0.7, 0.12],
        ],
        bounds=[(1e-3, 5.0), (0.2, 4.0), (1e-4, 1.0)],
        step0=[0.05, 0.1, 0.01],
    )
    d_du, a_du, M_du = dualist_delay_alpha(p_du[0], p_du[1], p_du[2], b_values)
    sc_du = score_model(
        L_du,
        n_sectors=2,
        has_free_bound_link=False,
        claims_soft_einstein=False,
        claims_local_budget_with_long_range_rho=False,
        occ=occ,
    )

    # soft-Einstein cheat: dualist map labeled monist
    sc_cheat = score_model(
        L_du,
        n_sectors=1,  # lies about sector count
        has_free_bound_link=True,  # lies about link
        claims_soft_einstein=True,
        claims_local_budget_with_long_range_rho=False,
        occ=occ,
    )
    # fix: cheat still pays softE even if it fakes N_sectors=1
    # Recompute honestly: if it claims monist but uses Phi gravity:
    sc_cheat = {
        **sc_cheat,
        "L_occ": 0.0,  # pretends one sector
        "S": L_du + occ.lambda_softE,
        "note": "dualist Phi labeled monist → only softE penalty applied",
    }

    # ranking
    scores = {
        "monist_kernel": sc_mk["S"],
        "local_optics": sc_lo["S"],
        "dualist": sc_du["S"],
        "softE_cheat": sc_cheat["S"],
    }
    score_winner = min(scores, key=scores.get)
    fit_scores = {
        "monist_kernel": L_mk,
        "local_optics": L_lo,
        "dualist": L_du,
    }
    fit_winner = min(fit_scores, key=fit_scores.get)

    # occam ablation on monist_kernel vs dualist
    abl = []
    for k in range(16):
        lam = 3.0 * k / 15.0
        occ_a = OccamWeights(
            lambda_sec=lam,
            lambda_link=occ.lambda_link,
            lambda_softE=occ.lambda_softE,
            lambda_nogo=occ.lambda_nogo,
        )
        sm = score_model(
            L_mk,
            n_sectors=1,
            has_free_bound_link=True,
            claims_soft_einstein=False,
            claims_local_budget_with_long_range_rho=False,
            occ=occ_a,
        )
        sd = score_model(
            L_du,
            n_sectors=2,
            has_free_bound_link=False,
            claims_soft_einstein=False,
            claims_local_budget_with_long_range_rho=False,
            occ=occ_a,
        )
        abl.append(
            {
                "lambda_sec": lam,
                "S_monist_kernel": sm["S"],
                "S_dualist": sd["S"],
                "winner": "monist_kernel" if sm["S"] <= sd["S"] else "dualist",
            }
        )

    # M consistency: lensing-inferred from kernel amp vs ledger
    # From delay shape, M_hat from monist kernel fit should match data_M
    M_consistency = abs(M_mk - data_M) / max(abs(data_M), 1e-12)

    return {
        "truth": asdict(truth),
        "noise_delay": noise_delay,
        "noise_alpha": noise_alpha,
        "seed": seed,
        "b_values": b_values,
        "data_delay": data_d,
        "data_alpha": data_a,
        "data_M": data_M,
        "true_delay": true_d,
        "true_alpha": true_a,
        "monist_kernel_fit": {
            "M": p_mk[0],
            "beta": p_mk[1],
            "eps": p_mk[2],
            "pred_delay": d_mk,
            "pred_alpha": a_mk,
            "pred_M": M_mk,
            "rel_err_M": abs(p_mk[0] - truth.M) / truth.M,
            "rel_err_beta": abs(p_mk[1] - truth.beta) / truth.beta,
            "rel_err_eps": abs(p_mk[2] - truth.eps) / truth.eps,
            "M_ledger_consistency": M_consistency,
            **sc_mk,
        },
        "local_optics_fit": {
            "A": p_lo[0],
            "sigma": p_lo[1],
            "kappa": p_lo[2],
            "pred_delay": d_lo,
            "pred_alpha": a_lo,
            "pred_M": M_lo,
            **sc_lo,
        },
        "dualist_fit": {
            "A_m": p_du[0],
            "sigma": p_du[1],
            "G_eff": p_du[2],
            "pred_delay": d_du,
            "pred_alpha": a_du,
            "pred_M": M_du,
            **sc_du,
        },
        "soft_einstein_cheat": sc_cheat,
        "fit_winner": fit_winner,
        "score_winner": score_winner,
        "occam_ablation": abl,
        "occam_weights": asdict(occ),
        "notes": [
            "Truth = monist nonlocal path-cost kernel (C Class C / FOR_D).",
            "Dualist = independent matter + Plummer Phi (two sectors).",
            "Local optics = compact n(ρ); expected poor long-range fit (C no-go).",
            "softE_cheat = dualist Phi labeled monist; L_softE disqualifies.",
        ],
    }


def save_outputs(result: dict) -> None:
    with open(OUT / "round1_result.json", "w") as f:
        json.dump(result, f, indent=2)

    mk = result["monist_kernel_fit"]
    lo = result["local_optics_fit"]
    du = result["dualist_fit"]
    ch = result["soft_einstein_cheat"]

    with open(OUT / "round1_scores.tsv", "w") as f:
        f.write(
            "model\tN_sectors\thas_link\tL_fit\tL_occ\tL_softE\tL_nogo\tS\t"
            "p1\tp2\tp3\tpred_M\tdata_M\n"
        )
        f.write(
            f"monist_kernel\t{mk['N_sectors']}\t{int(mk['has_free_bound_link'])}\t"
            f"{mk['L_fit']:.8e}\t{mk['L_occ']:.8e}\t{mk['L_softE']:.8e}\t{mk['L_nogo']:.8e}\t"
            f"{mk['S']:.8e}\t{mk['M']:.8e}\t{mk['beta']:.8e}\t{mk['eps']:.8e}\t"
            f"{mk['pred_M']:.8e}\t{result['data_M']:.8e}\n"
        )
        f.write(
            f"local_optics\t{lo['N_sectors']}\t{int(lo['has_free_bound_link'])}\t"
            f"{lo['L_fit']:.8e}\t{lo['L_occ']:.8e}\t{lo['L_softE']:.8e}\t{lo['L_nogo']:.8e}\t"
            f"{lo['S']:.8e}\t{lo['A']:.8e}\t{lo['sigma']:.8e}\t{lo['kappa']:.8e}\t"
            f"{lo['pred_M']:.8e}\t{result['data_M']:.8e}\n"
        )
        f.write(
            f"dualist\t{du['N_sectors']}\t{int(du['has_free_bound_link'])}\t"
            f"{du['L_fit']:.8e}\t{du['L_occ']:.8e}\t{du['L_softE']:.8e}\t{du['L_nogo']:.8e}\t"
            f"{du['S']:.8e}\t{du['A_m']:.8e}\t{du['sigma']:.8e}\t{du['G_eff']:.8e}\t"
            f"{du['pred_M']:.8e}\t{result['data_M']:.8e}\n"
        )
        f.write(
            f"softE_cheat\t1\t1\t"
            f"{ch['L_fit']:.8e}\t{ch['L_occ']:.8e}\t{ch['L_softE']:.8e}\t{ch.get('L_nogo',0):.8e}\t"
            f"{ch['S']:.8e}\t\t\t\t\t{result['data_M']:.8e}\n"
        )

    b = result["b_values"]
    with open(OUT / "round1_rays.tsv", "w") as f:
        f.write(
            "b\tdata_delay\ttrue_delay\tmk_delay\tlo_delay\tdual_delay\t"
            "data_alpha\ttrue_alpha\tmk_alpha\tlo_alpha\tdual_alpha\n"
        )
        for i in range(len(b)):
            f.write(
                f"{b[i]:.6f}\t"
                f"{result['data_delay'][i]:.8e}\t{result['true_delay'][i]:.8e}\t"
                f"{mk['pred_delay'][i]:.8e}\t{lo['pred_delay'][i]:.8e}\t{du['pred_delay'][i]:.8e}\t"
                f"{result['data_alpha'][i]:.8e}\t{result['true_alpha'][i]:.8e}\t"
                f"{mk['pred_alpha'][i]:.8e}\t{lo['pred_alpha'][i]:.8e}\t{du['pred_alpha'][i]:.8e}\n"
            )

    with open(OUT / "round1_occam_ablation.tsv", "w") as f:
        f.write("lambda_sec\tS_monist_kernel\tS_dualist\twinner\n")
        for row in result["occam_ablation"]:
            f.write(
                f"{row['lambda_sec']:.6f}\t{row['S_monist_kernel']:.8e}\t"
                f"{row['S_dualist']:.8e}\t{row['winner']}\n"
            )

    # path-cost profile sample at discrete r (truth vs recovered kernel)
    rs = linspace(0.2, 6.0, 30)
    with open(OUT / "round1_path_cost.tsv", "w") as f:
        f.write("r\tdell_true\tdell_mk\tdell_dual_proxy\n")
        t = result["truth"]
        for r in rs:
            dell_t = t["beta"] * t["M"] / math.sqrt(r * r + t["eps"] ** 2)
            dell_m = mk["beta"] * mk["M"] / math.sqrt(r * r + mk["eps"] ** 2)
            # dualist proxy: -Phi / c^2 style cost ~ G M / sqrt(r^2+a^2)
            a = du["sigma"]
            dell_d = du["G_eff"] * du["pred_M"] / math.sqrt(r * r + a * a)
            f.write(f"{r:.6f}\t{dell_t:.8e}\t{dell_m:.8e}\t{dell_d:.8e}\n")

    try_plot(result)
    write_svg(result)


def try_plot(result: dict) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"matplotlib skip: {e}")
        return

    b = result["b_values"]
    mk, lo, du = result["monist_kernel_fit"], result["local_optics_fit"], result["dualist_fit"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(b, result["true_delay"], "k-", lw=2, label="truth (kernel)")
    axes[0].plot(b, result["data_delay"], "ko", ms=3, alpha=0.45, label="data")
    axes[0].plot(b, mk["pred_delay"], "b--", label="monist kernel")
    axes[0].plot(b, lo["pred_delay"], "g-.", label="local optics")
    axes[0].plot(b, du["pred_delay"], "r:", label="dualist")
    axes[0].set_xlabel("b")
    axes[0].set_ylabel("delay")
    axes[0].set_title("Ray delay inversion")
    axes[0].legend(fontsize=7)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(b, result["true_alpha"], "k-", lw=2, label="truth")
    axes[1].plot(b, result["data_alpha"], "ko", ms=3, alpha=0.45, label="data")
    axes[1].plot(b, mk["pred_alpha"], "b--", label="monist kernel")
    axes[1].plot(b, lo["pred_alpha"], "g-.", label="local optics")
    axes[1].plot(b, du["pred_alpha"], "r:", label="dualist")
    axes[1].set_xlabel("b")
    axes[1].set_ylabel("alpha")
    axes[1].set_title("Ray deflection (Born)")
    axes[1].legend(fontsize=7)
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "round1_rays.png", dpi=140)
    plt.close(fig)

    labels = ["monist_kernel", "local_optics", "dualist", "softE_cheat"]
    Lfits = [mk["L_fit"], lo["L_fit"], du["L_fit"], result["soft_einstein_cheat"]["L_fit"]]
    penalties = [
        mk["S"] - mk["L_fit"],
        lo["S"] - lo["L_fit"],
        du["S"] - du["L_fit"],
        result["soft_einstein_cheat"]["S"] - result["soft_einstein_cheat"]["L_fit"],
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    x = list(range(4))
    ax.bar(x, Lfits, 0.55, label="L_fit")
    ax.bar(x, penalties, 0.55, bottom=Lfits, label="ontology penalties")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("score")
    ax.set_title(f"Combined score winner={result['score_winner']} | fit={result['fit_winner']}")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "round1_scores.png", dpi=140)
    plt.close(fig)

    abl = result["occam_ablation"]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([r["lambda_sec"] for r in abl], [r["S_monist_kernel"] for r in abl], "b-o", ms=4, label="S monist kernel")
    ax.plot([r["lambda_sec"] for r in abl], [r["S_dualist"] for r in abl], "r-s", ms=4, label="S dualist")
    ax.set_xlabel("λ_sec")
    ax.set_ylabel("S")
    ax.set_title("Occam ablation (kernel vs dualist)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "round1_occam_ablation.png", dpi=140)
    plt.close(fig)

    # path cost radial
    rs, dt, dm, dd = [], [], [], []
    with open(OUT / "round1_path_cost.tsv") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            rs.append(float(parts[0]))
            dt.append(float(parts[1]))
            dm.append(float(parts[2]))
            dd.append(float(parts[3]))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(rs, dt, "k-", lw=2, label="true δℓ(r)")
    ax.plot(rs, dm, "b--", label="recovered monist kernel")
    ax.plot(rs, dd, "r:", label="dualist Φ proxy")
    ax.set_xlabel("r")
    ax.set_ylabel("path-cost excess")
    ax.set_title("ℓ(r) recovery (C FOR_D target)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "round1_path_cost.png", dpi=140)
    plt.close(fig)
    print(f"Wrote PNG figures under {OUT}")


def write_svg(result: dict) -> None:
    b = result["b_values"]
    series = [
        (result["true_delay"], "#000000"),
        (result["monist_kernel_fit"]["pred_delay"], "#0000cc"),
        (result["local_optics_fit"]["pred_delay"], "#008800"),
        (result["dualist_fit"]["pred_delay"], "#cc0000"),
    ]
    w, h, pad = 520, 300, 40
    ys_all = []
    for s, _ in series:
        ys_all.extend(s)
    ymin, ymax = min(ys_all), max(ys_all)
    if abs(ymax - ymin) < 1e-15:
        ymax = ymin + 1.0
    bx0, bx1 = min(b), max(b)

    def xy(xv, yv):
        X = pad + (xv - bx0) / (bx1 - bx0) * (w - 2 * pad)
        Y = h - pad - (yv - ymin) / (ymax - ymin) * (h - 2 * pad)
        return X, Y

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{pad}" y="18" font-size="12">Ray delay (SVG)</text>',
    ]
    for ys, color in series:
        pts = " ".join(f"{xy(b[i], ys[i])[0]:.1f},{xy(b[i], ys[i])[1]:.1f}" for i in range(len(b)))
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{pts}"/>')
    parts.append("</svg>")
    (OUT / "round1_rays.svg").write_text("\n".join(parts))


def main():
    print("=== v76 Approach D invert_demo Round 1 ===")
    result = run()
    save_outputs(result)
    mk = result["monist_kernel_fit"]
    lo = result["local_optics_fit"]
    du = result["dualist_fit"]
    print(f"Truth: M={result['truth']['M']}, beta={result['truth']['beta']}, eps={result['truth']['eps']}")
    print(
        f"Monist kernel: M={mk['M']:.6f} beta={mk['beta']:.6f} eps={mk['eps']:.6f} "
        f"L_fit={mk['L_fit']:.6e} S={mk['S']:.6e}"
    )
    print(
        f"  rel_err M={mk['rel_err_M']:.4e} beta={mk['rel_err_beta']:.4e} eps={mk['rel_err_eps']:.4e}"
    )
    print(f"Local optics: L_fit={lo['L_fit']:.6e} S={lo['S']:.6e} pred_M={lo['pred_M']:.4f}")
    print(
        f"Dualist: L_fit={du['L_fit']:.6e} S={du['S']:.6e} "
        f"A_m={du['A_m']:.4f} sig={du['sigma']:.4f} G={du['G_eff']:.4f}"
    )
    print(f"softE cheat S={result['soft_einstein_cheat']['S']:.6e}")
    print(f"Pure-fit winner: {result['fit_winner']}")
    print(f"Combined-score winner: {result['score_winner']}")
    print(f"Outputs: {OUT}")
    return result


if __name__ == "__main__":
    main()
