#!/usr/bin/env python3
"""
v76 Approach D — Round 2 congruence scorer

Ingest B forward maps (local GRIN + postulated kernel).
Fit monist candidates and NON-isomorphic dualist adversaries.
Score fit + ontology Occam + C necessary-condition checklist.

No scp_sim. Pure stdlib.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)
B_OUT = Path("/home/d/code/scp/v76/work/B/outputs")

C_LIGHT = 1.0
X_RAY = 12.0
RHO0 = 1.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def linspace(a: float, b: float, n: int) -> List[float]:
    if n == 1:
        return [a]
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def trapz(y: Sequence[float], dx: float) -> float:
    if len(y) < 2:
        return 0.0
    s = 0.5 * (y[0] + y[-1])
    for i in range(1, len(y) - 1):
        s += y[i]
    return s * dx


def mse(a: Sequence[float], b: Sequence[float]) -> float:
    n = len(a)
    if n == 0:
        return 0.0
    return sum((a[i] - b[i]) ** 2 for i in range(n)) / n


def pattern_search(
    obj: Callable[[List[float]], float],
    x0: List[float],
    bounds: List[Tuple[float, float]],
    step0: List[float],
    max_iter: int = 100,
) -> Tuple[List[float], float]:
    x = [min(max(x0[i], bounds[i][0]), bounds[i][1]) for i in range(len(x0))]
    f = obj(x)
    step = list(step0)
    for _ in range(max_iter):
        improved = False
        for i in range(len(x)):
            for direction in (+1.0, -1.0):
                trial = list(x)
                trial[i] = min(
                    max(trial[i] + direction * step[i], bounds[i][0]), bounds[i][1]
                )
                ft = obj(trial)
                if ft < f - 1e-15:
                    x, f = trial, ft
                    improved = True
        if not improved:
            step = [s * 0.5 for s in step]
            if max(step) < 1e-6:
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
# Forward maps
# ---------------------------------------------------------------------------
def local_grin_rays(
    A: float, sigma: float, b_values: Sequence[float], n_x: int = 401
) -> Tuple[List[float], List[float], float]:
    """
    B-style optical monism: rho_b = A exp(-r^2/2σ^2) clipped,
    rho_f = rho0 - rho_b, n = rho0 / rho_f.
    Born along x at y=b.
    """
    A = min(max(A, 0.0), 0.94)
    sigma = max(sigma, 0.05)
    xs = linspace(-X_RAY, X_RAY, n_x)
    dx = xs[1] - xs[0]
    delays, alphas = [], []
    for b in b_values:
        n_line = []
        dndy_line = []
        for x in xs:
            r2 = x * x + b * b
            rb = A * math.exp(-r2 / (2.0 * sigma * sigma))
            rb = min(rb, RHO0 - 0.05)
            rf = RHO0 - rb
            n = RHO0 / max(rf, 1e-9)
            # ∂n/∂y = ∂n/∂rb * ∂rb/∂y; n=1/(1-rb/ρ0) wait n=ρ0/(ρ0-rb)
            # dn/drb = ρ0/(ρ0-rb)^2 = n/rf
            drb_dy = rb * (-b / (sigma * sigma))
            dndy = (RHO0 / (rf * rf)) * drb_dy
            n_line.append(n)
            dndy_line.append(dndy)
        delays.append(trapz([n - 1.0 for n in n_line], dx) / C_LIGHT)
        alphas.append(trapz(dndy_line, dx))
    M = A * 2.0 * math.pi * sigma * sigma / (C_LIGHT ** 2)
    return delays, alphas, M


def monist_kernel_rays(
    M: float, beta: float, eps: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    """Path-cost δℓ = β M / sqrt(r²+ε²); Born delay/alpha (same as R1)."""
    delays, alphas = [], []
    X = X_RAY
    for b in b_values:
        s2 = b * b + eps * eps
        s = math.sqrt(s2)
        delays.append(beta * M * 2.0 * math.asinh(X / s) / C_LIGHT)
        den = s2 * math.sqrt(X * X + s2)
        integ = (2.0 * X) / den if den > 0 else 0.0
        alphas.append(-beta * M * b * integ)
    return delays, alphas, M


def dualist_plummer_rays(
    M: float, G: float, a: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    """Isomorphic to monist kernel (R1 degeneracy). Kept for contrast only."""
    return monist_kernel_rays(M, G, a, b_values)


def dualist_log_rays(
    M: float, G: float, a: float, b_values: Sequence[float], n_x: int = 401
) -> Tuple[List[float], List[float], float]:
    """
    NON-isomorphic dualist: soft 2D log potential (attractive)
      Φ = -G M log(sqrt(r²+a²)/a) ≤ 0
    delay = ∫ -Φ dx > 0
    alpha = -G M b * ∫ dx/(x²+b²+a²)   (sign: alpha * b < 0, toward mass)
          = -G M b * (2/s0) arctan(X/s0),  s0=sqrt(b²+a²)
    Large-X: |alpha| → G M π |b|/s0 → G M π (not 1/b) — NOT kernel-isomorphic.
    """
    _ = n_x
    M = max(M, 1e-9)
    G = max(G, 0.0)
    a = max(a, 0.05)
    X = X_RAY
    delays, alphas = [], []
    xs = linspace(-X, X, 401)
    dx = xs[1] - xs[0]
    for b in b_values:
        s0 = math.sqrt(b * b + a * a)
        # delay: ∫ -Φ dx = G M ∫ log(s/a) dx
        vals = [
            G * M * math.log(math.sqrt(x * x + b * b + a * a) / a) for x in xs
        ]
        delays.append(trapz(vals, dx) / C_LIGHT)
        integ = (2.0 / s0) * math.atan(X / s0)
        alphas.append(-G * M * b * integ)
    return delays, alphas, M


def dualist_yukawa_rays(
    M: float, G: float, mu: float, a: float, b_values: Sequence[float], n_x: int = 401
) -> Tuple[List[float], List[float], float]:
    """
    NON-isomorphic dualist: Yukawa Φ = -G M exp(-μ r_soft) / r_soft
    with r_soft = sqrt(r²+a²). Wrong range μ>0 cannot match pure 1/r exterior.
    """
    M = max(M, 1e-9)
    G = max(G, 0.0)
    mu = max(mu, 0.0)
    a = max(a, 0.05)
    xs = linspace(-X_RAY, X_RAY, n_x)
    dx = xs[1] - xs[0]
    delays, alphas = [], []
    for b in b_values:
        phi_line = []
        dphidy_line = []
        for x in xs:
            r2 = x * x + b * b
            rs = math.sqrt(r2 + a * a)
            e = math.exp(-mu * rs)
            phi = -G * M * e / rs
            # ∂y Φ: d/dy [e/rs] = e * (-μ ∂rs/∂y * rs - ∂rs/∂y)/rs² wait
            # d(e/rs)/dy = e*(-μ)*(y/rs)/rs + e*(-1/rs²)*(y/rs)
            #            = -e y / rs² * (μ + 1/rs)
            dphidy = G * M * e * b / (rs * rs) * (mu + 1.0 / rs)
            # Φ = -G M e/rs ⇒ ∂y Φ = -G M ∂y(e/rs) = -G M (-e y/rs² (μ+1/rs))
            #                 = G M e y /rs² (μ+1/rs)
            # yes dphidy as above with y=b
            phi_line.append(phi)
            dphidy_line.append(dphidy)
        delays.append(trapz([-p for p in phi_line], dx) / C_LIGHT)
        # Φ = -GMe/rs, ∂yΦ > 0 for b>0; want alpha < 0 toward mass
        alphas.append(-trapz(dphidy_line, dx))
    return delays, alphas, M


def dualist_point_mass_1overb(
    amp: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    """
    Dualist weak formula defl = -amp/b only (no free-bound).
    Delay proxy ~ amp * 2 * asinh(X/|b|) with soft floor.
    Free G*M product; M not identified with ledger unless forced.
    """
    delays, alphas = [], []
    for b in b_values:
        bb = b if abs(b) > 0.15 else (0.15 if b >= 0 else -0.15)
        alphas.append(-amp / bb)
        s = abs(bb)
        delays.append(amp * 2.0 * math.asinh(X_RAY / s))
    return delays, alphas, float("nan")  # M free / undefined


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
@dataclass
class Occam:
    lambda_sec: float = 1.0
    lambda_link: float = 0.5
    lambda_softE: float = 100.0
    lambda_M: float = 0.5  # |M_pred - M_ledger| / M_ledger
    lambda_nogo: float = 2.0  # local budget claims long-range
    lambda_postulate: float = 0.3  # postulated kernel not dynamically derived


@dataclass
class ModelScore:
    name: str
    L_fit: float
    L_occ: float
    L_M: float
    L_softE: float
    L_extra: float
    S: float
    N_sectors: int
    has_free_bound_link: bool
    params: Dict
    pred_M: float
    claims_soft_einstein: bool = False
    is_postulated_kernel: bool = False
    notes: str = ""


def l_fit_rays(
    pred_d: Sequence[float],
    pred_a: Sequence[float],
    data_d: Sequence[float],
    data_a: Sequence[float],
    w_t: float = 1.0,
    w_a: float = 1.0,
) -> float:
    # normalize by data scale so delay/alpha comparable
    scale_d = max(math.sqrt(mse(data_d, [0.0] * len(data_d))), 1e-6)
    scale_a = max(math.sqrt(mse(data_a, [0.0] * len(data_a))), 1e-6)
    nd = [(pred_d[i] - data_d[i]) / scale_d for i in range(len(data_d))]
    na = [(pred_a[i] - data_a[i]) / scale_a for i in range(len(data_a))]
    return w_t * mse(nd, [0.0] * len(nd)) + w_a * mse(na, [0.0] * len(na))


def score_model(
    name: str,
    Lfit: float,
    pred_M: float,
    M_ledger: float,
    *,
    n_sectors: int,
    has_link: bool,
    softE: bool,
    postulated: bool,
    long_range_local_claim: bool,
    params: Dict,
    occ: Occam,
    notes: str = "",
) -> ModelScore:
    Locc = occ.lambda_sec * max(n_sectors - 1, 0)
    if not has_link:
        Locc += occ.lambda_link
    LM = 0.0
    if pred_M == pred_M and M_ledger > 0:  # not NaN
        LM = occ.lambda_M * abs(pred_M - M_ledger) / M_ledger
    Lsoft = occ.lambda_softE if softE else 0.0
    Lextra = 0.0
    if postulated:
        Lextra += occ.lambda_postulate
    if long_range_local_claim:
        Lextra += occ.lambda_nogo
    S = Lfit + Locc + LM + Lsoft + Lextra
    return ModelScore(
        name=name,
        L_fit=Lfit,
        L_occ=Locc,
        L_M=LM,
        L_softE=Lsoft,
        L_extra=Lextra,
        S=S,
        N_sectors=n_sectors,
        has_free_bound_link=has_link,
        params=params,
        pred_M=pred_M if pred_M == pred_M else -1.0,
        claims_soft_einstein=softE,
        is_postulated_kernel=postulated,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Ingest B
# ---------------------------------------------------------------------------
@dataclass
class RayMap:
    source: str
    tag: str
    b: List[float]
    delay: List[float]
    alpha: List[float]  # deflection_rad
    M_ledger: float
    free_deficit_core: float
    gravity_solver: Optional[str]
    design: str
    meta: Dict = field(default_factory=dict)


def load_b_local(path: Path = B_OUT / "results.json") -> List[RayMap]:
    with open(path) as f:
        data = json.load(f)
    maps = []
    for case in data.get("cases", []):
        if case.get("tag") == "vacuum":
            continue
        rays = case["rays"]
        # use positive-b half for cleaner fit (odd/even symmetry)
        b, dly, alf = [], [], []
        for r in rays:
            if r["b"] > 0:
                b.append(float(r["b"]))
                dly.append(float(r["delay"]))
                alf.append(float(r["deflection_rad"]))
        maps.append(
            RayMap(
                source=str(path),
                tag=case["tag"],
                b=b,
                delay=dly,
                alpha=alf,
                M_ledger=float(case["stats"]["m_bound"]),
                free_deficit_core=float(case["stats"]["free_deficit_core"]),
                gravity_solver=data.get("gravity_solver"),
                design=data.get("design", "local"),
                meta={
                    "A": case.get("A"),
                    "sigma": case.get("sigma"),
                    "channel": "local_GRIN",
                    "sector_tag": "1-sector-candidate",
                    "dynamical_free_response": False,
                },
            )
        )
    return maps


def load_b_kernel(path: Path = B_OUT / "results_kernel.json") -> List[RayMap]:
    with open(path) as f:
        data = json.load(f)
    rays = data["weak_field_Born_point_mass"]["rays"]
    b, dly, alf = [], [], []
    alpha_M = float(data["weak_field_Born_point_mass"]["alpha_M"])
    for r in rays:
        if r["b"] > 0:
            bb = float(r["b"])
            b.append(bb)
            # B kernel file may omit delay; synthesize Shapiro-like ∫Φ for Φ=αM/r
            # delay = αM * 2 * asinh(X/|b|) soft
            s = max(abs(bb), 0.15)
            dly.append(alpha_M * 2.0 * math.asinh(X_RAY / s))
            alf.append(float(r["deflection_rad"]))
    return [
        RayMap(
            source=str(path),
            tag="kernel_A0.35_s1.2",
            b=b,
            delay=dly,
            alpha=alf,
            M_ledger=float(data["M_bound_continuous"]),
            free_deficit_core=float(data["free_deficit_core"]),
            gravity_solver=data.get("gravity_solver"),
            design=data.get("design", "kernel"),
            meta={
                "alpha": data.get("alpha"),
                "A": data.get("A"),
                "sigma": data.get("sigma"),
                "channel": "postulated_kernel",
                "sector_tag": "1-sector-candidate-postulated",
                "dynamical_free_response": False,
                "exterior_Phi": data.get("exterior_Phi_vs_monopole"),
                "contrast_r6": data.get("contrast_local_vs_kernel_at_r6"),
            },
        )
    ]


def load_b_round2_tsv(path: Path = B_OUT / "round2_rays.tsv") -> List[RayMap]:
    """Ingest B Round-2 tagged ray exports (M1/M2/M3/R1)."""
    if not path.exists():
        return []
    # group by mechanism
    from collections import defaultdict

    groups: Dict[str, list] = defaultdict(list)
    sector_of: Dict[str, str] = {}
    m_of: Dict[str, float] = {}
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            row = dict(zip(header, parts))
            mech = row["mechanism"]
            groups[mech].append(row)
            sector_of[mech] = row["sector_tag"]
            m_of[mech] = float(row["m_ledger"])

    # free deficit from companion file if present
    free_def = 0.1513163
    fd_path = B_OUT / "round2_free_deficit.tsv"
    if fd_path.exists():
        with open(fd_path) as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and "core" in parts[0].lower():
                    try:
                        free_def = float(parts[1])
                    except ValueError:
                        pass
                    break

    out: List[RayMap] = []
    for mech, rows in groups.items():
        b, dly, alf = [], [], []
        for r in rows:
            bb = float(r["b"])
            if bb > 0:
                b.append(bb)
                dly.append(float(r["delay"]))
                alf.append(float(r["deflection_rad"]))
        if not b:
            continue
        sector = sector_of[mech]
        is_dyn = "free_laplace" in mech or "M2" in mech
        is_post = "postulated" in mech or "R1" in mech
        channel = (
            "dynamical_free_response"
            if is_dyn
            else (
                "postulated_kernel"
                if is_post
                else (
                    "dualist_poisson"
                    if "dualist" in sector or "poisson" in mech.lower()
                    else "local_GRIN"
                )
            )
        )
        out.append(
            RayMap(
                source=str(path),
                tag=mech,
                b=b,
                delay=dly,
                alpha=alf,
                M_ledger=m_of[mech],
                free_deficit_core=free_def,
                gravity_solver=(
                    "poisson_2d"
                    if "poisson" in mech.lower()
                    else (
                        "analytic_postulated"
                        if is_post
                        else None
                    )
                ),
                design=f"B_R2_{mech}",
                meta={
                    "channel": channel,
                    "sector_tag": sector,
                    "dynamical_free_response": is_dyn,
                    "ground_truth_sector": sector,
                    "mechanism": mech,
                },
            )
        )
    return out


def load_b_dynamical() -> List[RayMap]:
    """Round-2 TSV exports preferred; legacy JSON fallbacks."""
    out = load_b_round2_tsv()
    if out:
        return out
    candidates = [
        B_OUT / "results_dynamical.json",
        B_OUT / "results_free_response.json",
        B_OUT / "results_r2.json",
        B_OUT / "round2_result.json",
    ]
    for p in candidates:
        if not p.exists():
            continue
        # round2_result.json is meta-only; skip if no ray list
        with open(p) as f:
            data = json.load(f)
        rays = data.get("rays") or data.get("weak_field_Born_point_mass", {}).get("rays")
        if not rays:
            continue
        b, dly, alf = [], [], []
        for r in rays:
            if float(r.get("b", 0)) > 0:
                b.append(float(r["b"]))
                dly.append(float(r.get("delay", r.get("born_delay", 0.0))))
                alf.append(float(r.get("deflection_rad", r.get("born_defl_rad", 0.0))))
        out.append(
            RayMap(
                source=str(p),
                tag=data.get("tag", p.stem),
                b=b,
                delay=dly,
                alpha=alf,
                M_ledger=float(
                    data.get("M_bound_continuous", data.get("m_bound", data.get("M", 1.0)))
                ),
                free_deficit_core=float(data.get("free_deficit_core", 0.0)),
                gravity_solver=data.get("gravity_solver"),
                design=data.get("design", "dynamical"),
                meta={
                    "channel": "dynamical_free_response",
                    "sector_tag": data.get("sector_tag", "1-sector-candidate"),
                    "dynamical_free_response": True,
                },
            )
        )
    return out


# ---------------------------------------------------------------------------
# Fit each model class to a RayMap
# ---------------------------------------------------------------------------
def fit_all_models(rm: RayMap, occ: Occam) -> List[ModelScore]:
    b, data_d, data_a = rm.b, rm.delay, rm.alpha
    Mled = rm.M_ledger
    scores: List[ModelScore] = []

    # --- monist local GRIN ---
    def obj_lo(p):
        d, a, _ = local_grin_rays(p[0], p[1], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_lo, L_lo = multi_start(
        obj_lo,
        seeds=[
            [0.35, 1.2],
            [0.15, 1.2],
            [0.5, 1.0],
            [0.2, 1.5],
            [0.4, 0.9],
            [0.7, 1.2],
        ],
        bounds=[(1e-3, 0.94), (0.2, 4.0)],
        step0=[0.05, 0.1],
    )
    d_lo, a_lo, M_lo = local_grin_rays(p_lo[0], p_lo[1], b)
    scores.append(
        score_model(
            "monist_local_GRIN",
            L_lo,
            M_lo,
            Mled,
            n_sectors=1,
            has_link=True,
            softE=False,
            postulated=False,
            long_range_local_claim=False,
            params={"A": p_lo[0], "sigma": p_lo[1]},
            occ=occ,
            notes="B-style n=ρ0/ρ_free; compact GRIN",
        )
    )

    # --- monist postulated kernel ---
    def obj_mk(p):
        d, a, _ = monist_kernel_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_mk, L_mk = multi_start(
        obj_mk,
        seeds=[
            [Mled, 0.08, 0.85],
            [Mled, 0.12, 1.0],
            [Mled, 0.05, 0.5],
            [Mled * 0.8, 0.1, 1.2],
            [Mled * 1.2, 0.15, 0.7],
            [2.5, 0.12, 0.85],
        ],
        bounds=[(0.1, 20.0), (1e-3, 2.0), (0.05, 4.0)],
        step0=[0.2, 0.02, 0.1],
    )
    d_mk, a_mk, M_mk = monist_kernel_rays(p_mk[0], p_mk[1], p_mk[2], b)
    # force ledger comparison: also score with M fixed to ledger
    d_mk_fix, a_mk_fix, _ = monist_kernel_rays(Mled, p_mk[1], p_mk[2], b)
    L_mk_led = l_fit_rays(d_mk_fix, a_mk_fix, data_d, data_a)
    scores.append(
        score_model(
            "monist_kernel_postulated",
            min(L_mk, L_mk_led),
            p_mk[0],
            Mled,
            n_sectors=1,
            has_link=True,
            softE=False,
            postulated=True,
            long_range_local_claim=False,
            params={"M": p_mk[0], "beta": p_mk[1], "eps": p_mk[2], "L_fit_free_M": L_mk, "L_fit_M_ledger": L_mk_led},
            occ=occ,
            notes="Class-C path-cost kernel; β postulated not dynamical (λ_postulate)",
        )
    )

    # --- monist dynamical free-response (if this map is tagged dynamical) ---
    if rm.meta.get("dynamical_free_response"):
        # same kernel form but NO postulate penalty
        scores.append(
            score_model(
                "monist_dynamical_free_response",
                min(L_mk, L_mk_led),
                p_mk[0],
                Mled,
                n_sectors=1,
                has_link=True,
                softE=False,
                postulated=False,
                long_range_local_claim=False,
                params={"M": p_mk[0], "beta": p_mk[1], "eps": p_mk[2]},
                occ=occ,
                notes="Dynamical free-response export from B (no postulate penalty)",
            )
        )

    # --- dualist Plummer (isomorphic — should tie on kernel data) ---
    def obj_pl(p):
        d, a, _ = dualist_plummer_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_pl, L_pl = multi_start(
        obj_pl,
        seeds=[
            [Mled, 0.08, 0.85],
            [Mled, 0.12, 1.0],
            [Mled, 0.05, 0.5],
            [3.0, 0.1, 1.0],
        ],
        bounds=[(0.1, 20.0), (1e-3, 2.0), (0.05, 4.0)],
        step0=[0.2, 0.02, 0.1],
    )
    scores.append(
        score_model(
            "dualist_plummer_isomorphic",
            L_pl,
            p_pl[0],
            Mled,
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            long_range_local_claim=False,
            params={"M": p_pl[0], "G": p_pl[1], "a": p_pl[2]},
            occ=occ,
            notes="R1 degeneracy control: may tie L_fit on kernel data",
        )
    )

    # --- dualist 2D log (NON-isomorphic) ---
    def obj_log(p):
        d, a, _ = dualist_log_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_log, L_log = multi_start(
        obj_log,
        seeds=[
            [Mled, 0.05, 0.5],
            [Mled, 0.1, 1.0],
            [Mled, 0.2, 0.8],
            [2.0, 0.08, 1.2],
            [4.0, 0.15, 0.6],
            [Mled, 0.3, 2.0],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.05, 5.0)],
        step0=[0.2, 0.03, 0.15],
    )
    scores.append(
        score_model(
            "dualist_log_2D",
            L_log,
            p_log[0],
            Mled,
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            long_range_local_claim=False,
            params={"M": p_log[0], "G": p_log[1], "a": p_log[2]},
            occ=occ,
            notes="NON-isomorphic dualist adversary (log potential)",
        )
    )

    # --- dualist Yukawa (NON-isomorphic, wrong range) ---
    def obj_yk(p):
        d, a, _ = dualist_yukawa_rays(p[0], p[1], p[2], p[3], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_yk, L_yk = multi_start(
        obj_yk,
        seeds=[
            [Mled, 0.1, 0.5, 0.5],
            [Mled, 0.2, 1.0, 0.8],
            [Mled, 0.08, 0.2, 1.0],
            [Mled, 0.15, 0.8, 0.3],
            [3.0, 0.12, 0.3, 0.85],
            [Mled, 0.05, 0.05, 0.5],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.0, 3.0), (0.05, 4.0)],
        step0=[0.2, 0.03, 0.1, 0.1],
    )
    scores.append(
        score_model(
            "dualist_yukawa",
            L_yk,
            p_yk[0],
            Mled,
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            long_range_local_claim=False,
            params={"M": p_yk[0], "G": p_yk[1], "mu": p_yk[2], "a": p_yk[3]},
            occ=occ,
            notes="NON-isomorphic dualist (Yukawa wrong range)",
        )
    )

    # --- dualist free amp 1/b (G free, M not ledger-tied) ---
    def obj_1b(p):
        d, a, _ = dualist_point_mass_1overb(p[0], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_1b, L_1b = multi_start(
        obj_1b,
        seeds=[[0.25], [0.5], [0.1], [1.0], [0.05], [0.8]],
        bounds=[(1e-4, 5.0)],
        step0=[0.05],
    )
    # M undefined → use NaN path: assign pred_M far from ledger unless amp implies M
    scores.append(
        score_model(
            "dualist_free_GM_1overb",
            L_1b,
            float("nan"),
            Mled,
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            long_range_local_claim=False,
            params={"amp_GM": p_1b[0]},
            occ=occ,
            notes="Free G*M product; no free-bound ledger (M_consistency N/A → no L_M but no link)",
        )
    )

    # --- softE cheat: dualist log labeled monist ---
    scores.append(
        score_model(
            "softE_dualist_log_as_monist",
            L_log,
            p_log[0],
            Mled,
            n_sectors=1,  # lies
            has_link=True,  # lies
            softE=True,
            postulated=False,
            long_range_local_claim=False,
            params={"M": p_log[0], "G": p_log[1], "a": p_log[2]},
            occ=occ,
            notes="Dualist log Phi labeled monist → λ_softE",
        )
    )
    # Fix softE score: only softE penalty on L_fit (honest cheat accounting)
    # Recompute S as L_fit + softE only (as if it fakes sectors)
    cheat = scores[-1]
    cheat.L_occ = 0.0
    cheat.L_M = 0.0
    cheat.L_extra = 0.0
    cheat.S = cheat.L_fit + occ.lambda_softE

    return scores


# ---------------------------------------------------------------------------
# C necessary-condition checklist (numeric filters)
# ---------------------------------------------------------------------------
def nc_checklist(rm: RayMap, winner: ModelScore, all_scores: List[ModelScore]) -> Dict:
    """Heuristic pass/fail against C NC package for this map + winning model."""
    channel = rm.meta.get("channel", "")
    # exterior ray falloff: |alpha|(b=4) / |alpha|(b=1)
    try:
        i1 = rm.b.index(1.0) if 1.0 in rm.b else 0
        i4 = rm.b.index(4.0) if 4.0 in rm.b else -1
        a1 = abs(rm.alpha[i1])
        a4 = abs(rm.alpha[i4])
        ratio_4_1 = a4 / max(a1, 1e-12)
    except Exception:
        ratio_4_1 = float("nan")

    # monopole class: pure 1/b has ratio (1/4)/(1/1)=0.25; compact GRIN << 0.25
    long_range = ratio_4_1 > 0.15  # generous

    monist_kernel = next(
        (s for s in all_scores if s.name == "monist_kernel_postulated"), None
    )
    monist_local = next(
        (s for s in all_scores if s.name == "monist_local_GRIN"), None
    )
    dual_log = next((s for s in all_scores if s.name == "dualist_log_2D"), None)

    checks = {
        "NC_free_deficit_positive": rm.free_deficit_core > 0.01,
        "NC_M_ledger_positive": rm.M_ledger > 0,
        "NC_rays_nonzero": max(abs(a) for a in rm.alpha) > 1e-4,
        "NC_long_range_path_cost": long_range,
        "NC_alpha_ratio_b4_over_b1": ratio_4_1,
        "NC_winner_one_sector": winner.N_sectors == 1,
        "NC_winner_has_link": winner.has_free_bound_link,
        "NC_winner_not_softE": not winner.claims_soft_einstein,
        "NC_kernel_beats_local_on_long_range": (
            monist_kernel is not None
            and monist_local is not None
            and long_range
            and monist_kernel.L_fit < monist_local.L_fit
        ),
        "NC_local_beats_kernel_on_compact": (
            monist_kernel is not None
            and monist_local is not None
            and (not long_range)
            and monist_local.L_fit < monist_kernel.L_fit
        ),
        "NC_noniso_dualist_worse_fit_than_monist_kernel": (
            monist_kernel is not None
            and dual_log is not None
            and monist_kernel.L_fit + 1e-6 < dual_log.L_fit
        ),
        "NC_gravity_solver_null": rm.gravity_solver is None,
        "NC_dynamical_free_response": bool(rm.meta.get("dynamical_free_response")),
        "NC_postulated_kernel_only": channel == "postulated_kernel",
        "channel": channel,
        "winner": winner.name,
        "winner_S": winner.S,
        "winner_L_fit": winner.L_fit,
    }
    # overall "C weak-field monist package" pass for THIS map
    checks["C_package_local_only"] = (
        checks["NC_free_deficit_positive"]
        and checks["NC_rays_nonzero"]
        and checks["NC_winner_one_sector"]
        and not checks["NC_long_range_path_cost"]
    )
    checks["C_package_long_range_candidate"] = (
        checks["NC_free_deficit_positive"]
        and checks["NC_long_range_path_cost"]
        and checks["NC_winner_one_sector"]
        and checks["NC_winner_has_link"]
        and checks["NC_winner_not_softE"]
        and checks["NC_gravity_solver_null"]
    )
    # Honest monist theory+numerics: needs dynamical free response
    checks["C_package_monist_complete"] = (
        checks["C_package_long_range_candidate"]
        and checks["NC_dynamical_free_response"]
        and not checks["NC_postulated_kernel_only"]
    )
    return checks


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run() -> dict:
    occ = Occam()
    maps: List[RayMap] = []
    maps.extend(load_b_local())
    maps.extend(load_b_kernel())
    maps.extend(load_b_dynamical())

    # Prefer primary lock for local (A=0.35)
    primary_local = [m for m in maps if m.tag.startswith("lock_A0.35")]
    primary_kernel = [m for m in maps if "kernel" in m.tag]
    dynamical = [m for m in maps if m.meta.get("dynamical_free_response")]
    eval_maps = primary_local[:1] + primary_kernel[:1] + dynamical

    # also score weak local for scaling check
    weak_local = [m for m in maps if "A0.15" in m.tag]
    if weak_local:
        eval_maps.append(weak_local[0])

    report = {
        "round": 2,
        "occam": asdict(occ),
        "maps_found": [
            {
                "tag": m.tag,
                "channel": m.meta.get("channel"),
                "M_ledger": m.M_ledger,
                "free_deficit_core": m.free_deficit_core,
                "dynamical": m.meta.get("dynamical_free_response"),
                "source": m.source,
            }
            for m in maps
        ],
        "evaluations": [],
        "dynamical_free_response_present": len(dynamical) > 0,
    }

    all_rows = []
    for rm in eval_maps:
        scores = fit_all_models(rm, occ)
        scores_sorted = sorted(scores, key=lambda s: s.S)
        winner = scores_sorted[0]
        fit_winner = min(scores, key=lambda s: s.L_fit)
        nc = nc_checklist(rm, winner, scores)
        ev = {
            "map_tag": rm.tag,
            "channel": rm.meta.get("channel"),
            "M_ledger": rm.M_ledger,
            "free_deficit_core": rm.free_deficit_core,
            "score_winner": winner.name,
            "fit_winner": fit_winner.name,
            "scores": [asdict(s) for s in scores_sorted],
            "nc_checklist": nc,
        }
        report["evaluations"].append(ev)
        for s in scores_sorted:
            all_rows.append(
                {
                    "map": rm.tag,
                    "channel": rm.meta.get("channel"),
                    "model": s.name,
                    "L_fit": s.L_fit,
                    "L_occ": s.L_occ,
                    "L_M": s.L_M,
                    "L_softE": s.L_softE,
                    "L_extra": s.L_extra,
                    "S": s.S,
                    "N_sectors": s.N_sectors,
                    "has_link": int(s.has_free_bound_link),
                    "pred_M": s.pred_M,
                    "M_ledger": rm.M_ledger,
                }
            )

    # global verdict fields
    report["summary"] = make_summary(report)
    return report


def make_summary(report: dict) -> dict:
    dyn = report["dynamical_free_response_present"]
    evals = report["evaluations"]
    kernel_ev = next(
        (e for e in evals if e.get("channel") == "postulated_kernel"), None
    )
    local_ev = next(
        (e for e in evals if e.get("channel") == "local_GRIN"), None
    )
    dyn_ev = next(
        (e for e in evals if e.get("channel") == "dynamical_free_response"), None
    )

    def best_monist(ev):
        if not ev:
            return None
        mon = [
            s
            for s in ev["scores"]
            if s["name"].startswith("monist") and not s.get("claims_soft_einstein")
        ]
        return min(mon, key=lambda s: s["S"]) if mon else None

    def best_noniso_dualist(ev):
        if not ev:
            return None
        du = [
            s
            for s in ev["scores"]
            if s["name"] in ("dualist_log_2D", "dualist_yukawa", "dualist_free_GM_1overb")
        ]
        return min(du, key=lambda s: s["S"]) if du else None

    out = {
        "dynamical_present": dyn,
        "local_score_winner": local_ev["score_winner"] if local_ev else None,
        "kernel_score_winner": kernel_ev["score_winner"] if kernel_ev else None,
        "dynamical_score_winner": dyn_ev["score_winner"] if dyn_ev else None,
        "local_best_monist": best_monist(local_ev),
        "kernel_best_monist": best_monist(kernel_ev),
        "kernel_best_noniso_dualist": best_noniso_dualist(kernel_ev),
        "local_nc": local_ev["nc_checklist"] if local_ev else None,
        "kernel_nc": kernel_ev["nc_checklist"] if kernel_ev else None,
        "goal2_closer": False,
        "goal2_note": "",
    }

    # Goal condition (2): theory+numerics congruent workable
    # Closer if: non-iso dualist loses on kernel map AND monist kernel wins S
    # Complete only if dynamical present and wins
    if kernel_ev:
        km = best_monist(kernel_ev)
        kd = best_noniso_dualist(kernel_ev)
        if km and kd and km["S"] < kd["S"] and km["L_fit"] <= kd["L_fit"] + 0.05:
            out["goal2_closer"] = True
            out["goal2_note"] = (
                "Monist kernel beats non-isomorphic dualist on combined S "
                "(and competitive L_fit) for B kernel map."
            )
        if km and kd and km["L_fit"] + 1e-9 < kd["L_fit"]:
            out["goal2_note"] += " Pure L_fit also prefers monist kernel over log/Yukawa dualist."
    if dyn and dyn_ev:
        dm = best_monist(dyn_ev)
        if dm and dm["S"] < 1.0 and dyn_ev["nc_checklist"].get("C_package_monist_complete"):
            out["goal2_closer"] = True
            out["goal2_note"] += " Dynamical free-response map present and C_package_monist_complete."
        else:
            out["goal2_note"] += " Dynamical map present but incomplete NC package."
    else:
        out["goal2_note"] += (
            " No dynamical free-response export from B yet — "
            "postulated kernel still carries λ_postulate; goal (2) NOT met."
        )
    return out


def save(report: dict) -> None:
    with open(OUT / "round2_result.json", "w") as f:
        json.dump(report, f, indent=2)

    # scores TSV
    with open(OUT / "round2_scores.tsv", "w") as f:
        f.write(
            "map\tchannel\tmodel\tL_fit\tL_occ\tL_M\tL_softE\tL_extra\tS\t"
            "N_sectors\thas_link\tpred_M\tM_ledger\n"
        )
        for ev in report["evaluations"]:
            for s in ev["scores"]:
                f.write(
                    f"{ev['map_tag']}\t{ev['channel']}\t{s['name']}\t"
                    f"{s['L_fit']:.8e}\t{s['L_occ']:.8e}\t{s['L_M']:.8e}\t"
                    f"{s['L_softE']:.8e}\t{s['L_extra']:.8e}\t{s['S']:.8e}\t"
                    f"{s['N_sectors']}\t{int(s['has_free_bound_link'])}\t"
                    f"{s['pred_M']:.8e}\t{ev['M_ledger']:.8e}\n"
                )

    # winners TSV
    with open(OUT / "round2_winners.tsv", "w") as f:
        f.write(
            "map\tchannel\tscore_winner\tfit_winner\twinner_S\twinner_L_fit\t"
            "long_range\tC_long_range_candidate\tC_monist_complete\n"
        )
        for ev in report["evaluations"]:
            nc = ev["nc_checklist"]
            win = next(s for s in ev["scores"] if s["name"] == ev["score_winner"])
            f.write(
                f"{ev['map_tag']}\t{ev['channel']}\t{ev['score_winner']}\t"
                f"{ev['fit_winner']}\t{win['S']:.6e}\t{win['L_fit']:.6e}\t"
                f"{nc.get('NC_long_range_path_cost')}\t"
                f"{nc.get('C_package_long_range_candidate')}\t"
                f"{nc.get('C_package_monist_complete')}\n"
            )

    # NC checklist TSV
    with open(OUT / "round2_nc_checklist.tsv", "w") as f:
        f.write("map\tcheck\tvalue\n")
        for ev in report["evaluations"]:
            for k, v in ev["nc_checklist"].items():
                f.write(f"{ev['map_tag']}\t{k}\t{v}\n")

    # rays comparison for primary maps
    with open(OUT / "round2_rays_ingest.tsv", "w") as f:
        f.write("map\tb\tdelay\talpha\tM_ledger\n")
        for ev in report["evaluations"]:
            # re-load from maps_found is insufficient; re-parse from scores only
            pass
    # rewrite rays from fresh load
    maps = load_b_local() + load_b_kernel() + load_b_dynamical()
    with open(OUT / "round2_rays_ingest.tsv", "w") as f:
        f.write("map\tchannel\tb\tdelay\talpha\tM_ledger\tfree_deficit\n")
        for m in maps:
            if m.tag.startswith("lock_A0.7"):
                continue  # skip strong Born-caveat case for clean table
            for i, b in enumerate(m.b):
                f.write(
                    f"{m.tag}\t{m.meta.get('channel')}\t{b:.4f}\t"
                    f"{m.delay[i]:.8e}\t{m.alpha[i]:.8e}\t"
                    f"{m.M_ledger:.6f}\t{m.free_deficit_core:.6f}\n"
                )

    sm = report["summary"]
    lines = [
        "v76 Approach D Round 2 — congruence score summary",
        f"dynamical_free_response_present = {sm['dynamical_present']}",
        f"local_score_winner = {sm['local_score_winner']}",
        f"kernel_score_winner = {sm['kernel_score_winner']}",
        f"dynamical_score_winner = {sm['dynamical_score_winner']}",
        f"goal2_closer = {sm['goal2_closer']}",
        f"goal2_note = {sm['goal2_note']}",
        "",
    ]
    if sm.get("kernel_best_monist"):
        km = sm["kernel_best_monist"]
        lines.append(
            f"kernel best monist: {km['name']} L_fit={km['L_fit']:.4e} S={km['S']:.4e}"
        )
    if sm.get("kernel_best_noniso_dualist"):
        kd = sm["kernel_best_noniso_dualist"]
        lines.append(
            f"kernel best noniso dualist: {kd['name']} L_fit={kd['L_fit']:.4e} S={kd['S']:.4e}"
        )
    if sm.get("local_best_monist"):
        lm = sm["local_best_monist"]
        lines.append(
            f"local best monist: {lm['name']} L_fit={lm['L_fit']:.4e} S={lm['S']:.4e}"
        )
    lines.append("")
    lines.append("Per-map winners:")
    for ev in report["evaluations"]:
        lines.append(
            f"  {ev['map_tag']} ({ev['channel']}): score={ev['score_winner']} "
            f"fit={ev['fit_winner']} long_range={ev['nc_checklist'].get('NC_long_range_path_cost')}"
        )
    text = "\n".join(lines) + "\n"
    (OUT / "round2_summary.txt").write_text(text)
    print(text)


def main():
    print("=== v76 Approach D congruence_score_r2 ===")
    report = run()
    save(report)
    print(f"Wrote under {OUT}")
    return report


if __name__ == "__main__":
    main()
