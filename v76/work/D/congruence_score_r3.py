#!/usr/bin/env python3
"""
v76 Approach D — Round 3 congruence scorer (3D free-response)

Targets (O-003):
  - Monist 3D free Laplace / free-capacity Green ~ 1/r
  - Dualist 3D Poisson adversary (sector_tag dualist_2sector)
  - C weak-field multipole target: path-cost / deflection ~ M/b class

Ingest:
  - work/B/outputs/round3_*.tsv / round3_result.json when present
  - else synthetic 3D free-Green truth + dualist Poisson twin (documented)

No scp_sim. Pure stdlib.
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)
B_OUT = Path("/home/d/code/scp/v76/work/B/outputs")

C_LIGHT = 1.0
X_RAY = 20.0  # longer chart path for 3D weak-field tails


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
    max_iter: int = 80,
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
# 3D path-cost / potential forward maps
# ---------------------------------------------------------------------------
def soft_r(x: float, y: float, z: float, eps: float) -> float:
    return math.sqrt(x * x + y * y + z * z + eps * eps)


def monist_3d_green_rays(
    M: float,
    gamma: float,
    eps: float,
    b_values: Sequence[float],
    n_x: int = 401,
) -> Tuple[List[float], List[float], List[float], float]:
    """
    Monist 3D free-response (F1/M2-3D class):
      ψ(r) = gamma * M / r_soft     (3D Laplace Green exterior)
      path-cost excess δℓ = ψ  (γ absorbs coupling)
    Born along x at impact b (y=b, z=0):
      delay = ∫ ψ dx
      alpha = ∫ ∂y ψ dx   with sign convention alpha*b < 0 toward mass
    Analytical:
      ∫ dx / R = 2 asinh(X / s), s=sqrt(b²+ε²)
      ∂y ψ = -gamma M y / R^3
      ∫ ∂y ψ dx = -gamma M b * 2 X / (s² sqrt(X²+s²))
    """
    M = max(M, 1e-12)
    gamma = max(gamma, 0.0)
    eps = max(eps, 1e-4)
    X = X_RAY
    delays, alphas, path_at_b = [], [], []
    for b in b_values:
        s2 = b * b + eps * eps
        s = math.sqrt(s2)
        delays.append(gamma * M * 2.0 * math.asinh(X / s) / C_LIGHT)
        den = s2 * math.sqrt(X * X + s2)
        integ = (2.0 * X) / den if den > 0 else 0.0
        # ∂y ψ < 0 for b>0 if ψ>0? ψ = gM/R, ∂yψ = -gM b/R^3 < 0 for b>0
        # monist R2 used alpha < 0 for b>0: alpha = ∫ ∂y ψ
        alphas.append(gamma * M * (-b) * integ)
        # exterior path cost sample at r=|b| in transverse plane
        path_at_b.append(gamma * M / math.sqrt(b * b + eps * eps))
    return delays, alphas, path_at_b, M


def dualist_3d_poisson_rays(
    M: float,
    G: float,
    eps: float,
    b_values: Sequence[float],
) -> Tuple[List[float], List[float], List[float], float]:
    """
    Dualist 3D Newtonian / Poisson potential (2 sectors):
      Φ = -G M / r_soft
    Observables from −Φ as path-cost proxy (Shapiro-like):
      same algebra as monist_3d_green with gamma ↔ G
    ISOMORPHIC to monist 3D free Green on rays — Occam + tags separate.
    """
    # Map G→gamma with positive path cost: δℓ = -Φ = G M / r
    d, a, p, Mm = monist_3d_green_rays(M, G, eps, b_values)
    return d, a, p, Mm


def dualist_2d_log_rays(
    M: float, G: float, a: float, b_values: Sequence[float]
) -> Tuple[List[float], List[float], float]:
    """NON-isomorphic in 3D context: 2D log multipole (wrong dimension)."""
    M = max(M, 1e-12)
    G = max(G, 0.0)
    a = max(a, 0.05)
    X = X_RAY
    xs = linspace(-X, X, 401)
    dx = xs[1] - xs[0]
    delays, alphas = [], []
    for b in b_values:
        s0 = math.sqrt(b * b + a * a)
        vals = [
            G * M * math.log(math.sqrt(x * x + b * b + a * a) / a) for x in xs
        ]
        delays.append(trapz(vals, dx) / C_LIGHT)
        integ = (2.0 / s0) * math.atan(X / s0)
        alphas.append(-G * M * b * integ)
    return delays, alphas, M


def dualist_yukawa_3d(
    M: float, G: float, mu: float, eps: float, b_values: Sequence[float], n_x: int = 401
) -> Tuple[List[float], List[float], float]:
    """NON-iso if μ>0: Yukawa screened 3D potential."""
    M = max(M, 1e-12)
    G = max(G, 0.0)
    mu = max(mu, 0.0)
    eps = max(eps, 1e-4)
    xs = linspace(-X_RAY, X_RAY, n_x)
    dx = xs[1] - xs[0]
    delays, alphas = [], []
    for b in b_values:
        phi_line, dphidy_line = [], []
        for x in xs:
            rs = math.sqrt(x * x + b * b + eps * eps)
            e = math.exp(-mu * rs)
            phi = -G * M * e / rs
            # ∂y Φ = G M e (μ + 1/rs) * y / rs²
            dphidy = G * M * e * (mu + 1.0 / rs) * b / (rs * rs)
            phi_line.append(phi)
            dphidy_line.append(dphidy)
        delays.append(trapz([-p for p in phi_line], dx) / C_LIGHT)
        alphas.append(-trapz(dphidy_line, dx))
    return delays, alphas, M


def local_grin_compact(
    A: float, sigma: float, b_values: Sequence[float], n_x: int = 301
) -> Tuple[List[float], List[float], float]:
    """Compact wrong-class control (2D-style GRIN residual)."""
    RHO0 = 1.0
    A = min(max(A, 0.0), 0.94)
    sigma = max(sigma, 0.05)
    xs = linspace(-X_RAY, X_RAY, n_x)
    dx = xs[1] - xs[0]
    delays, alphas = [], []
    for b in b_values:
        n_line, dndy_line = [], []
        for x in xs:
            r2 = x * x + b * b
            rb = min(A * math.exp(-r2 / (2.0 * sigma * sigma)), RHO0 - 0.05)
            rf = RHO0 - rb
            n = RHO0 / max(rf, 1e-9)
            drb_dy = rb * (-b / (sigma * sigma))
            dndy = (RHO0 / (rf * rf)) * drb_dy
            n_line.append(n)
            dndy_line.append(dndy)
        delays.append(trapz([n - 1.0 for n in n_line], dx))
        alphas.append(trapz(dndy_line, dx))
    M = A * 2.0 * math.pi * sigma * sigma
    return delays, alphas, M


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
@dataclass
class Occam:
    lambda_sec: float = 1.0
    lambda_link: float = 0.5
    lambda_softE: float = 100.0
    lambda_M: float = 0.5
    lambda_postulate: float = 0.3
    lambda_wrong_multipole: float = 0.5  # claim 1/r but ratio far from 0.25
    lambda_gravity_solver: float = 1.0  # monist claim + gravity_solver set


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
    multipole_ratio: float
    notes: str = ""


def l_fit_rays(
    pred_d: Sequence[float],
    pred_a: Sequence[float],
    data_d: Sequence[float],
    data_a: Sequence[float],
) -> float:
    scale_d = max(math.sqrt(mse(data_d, [0.0] * len(data_d))), 1e-6)
    scale_a = max(math.sqrt(mse(data_a, [0.0] * len(data_a))), 1e-6)
    nd = [(pred_d[i] - data_d[i]) / scale_d for i in range(len(data_d))]
    na = [(pred_a[i] - data_a[i]) / scale_a for i in range(len(data_a))]
    return mse(nd, [0.0] * len(nd)) + mse(na, [0.0] * len(na))


def multipole_ratio(b: Sequence[float], alpha: Sequence[float]) -> float:
    """|α|(b≈4) / |α|(b≈1); ~0.25 for 1/r, ~0.7 for log, <<0.1 compact."""
    def nearest(target):
        i = min(range(len(b)), key=lambda j: abs(b[j] - target))
        return abs(alpha[i])

    a1 = nearest(1.0)
    a4 = nearest(4.0)
    return a4 / max(a1, 1e-12)


def score_model(
    name: str,
    Lfit: float,
    pred_M: float,
    M_ledger: float,
    multipole_r: float,
    *,
    n_sectors: int,
    has_link: bool,
    softE: bool,
    postulated: bool,
    claims_1r: bool,
    gravity_solver_present: bool,
    params: Dict,
    occ: Occam,
    notes: str = "",
) -> ModelScore:
    Locc = occ.lambda_sec * max(n_sectors - 1, 0)
    if not has_link:
        Locc += occ.lambda_link
    LM = 0.0
    if pred_M == pred_M and M_ledger > 0:
        LM = occ.lambda_M * abs(pred_M - M_ledger) / M_ledger
    Lsoft = occ.lambda_softE if softE else 0.0
    Lextra = 0.0
    if postulated:
        Lextra += occ.lambda_postulate
    if claims_1r and abs(multipole_r - 0.25) > 0.12:
        Lextra += occ.lambda_wrong_multipole
    # monist 1-sector claim while data/model carries gravity_solver
    if n_sectors == 1 and gravity_solver_present and "dualist" not in name:
        # only penalize if this is a monist-labeled model of dualist machinery
        if softE or "poisson" in name.lower():
            Lextra += occ.lambda_gravity_solver
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
        multipole_ratio=multipole_r,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Ray maps
# ---------------------------------------------------------------------------
@dataclass
class RayMap:
    source: str
    tag: str
    b: List[float]
    delay: List[float]
    alpha: List[float]
    M_ledger: float
    free_deficit_core: float
    gravity_solver: Optional[str]
    design: str
    meta: Dict = field(default_factory=dict)


def poll_b_round3() -> List[str]:
    """Return list of round3-related files currently in B outputs."""
    if not B_OUT.exists():
        return []
    found = []
    for name in sorted(os.listdir(B_OUT)):
        if name.startswith("round3") or "3d" in name.lower() or "3D" in name:
            found.append(str(B_OUT / name))
    return found


def load_b_round3_tsv(path: Path = B_OUT / "round3_rays.tsv") -> List[RayMap]:
    if not path.exists():
        return []
    groups: Dict[str, list] = defaultdict(list)
    sector_of: Dict[str, str] = {}
    origin_of: Dict[str, str] = {}
    m_of: Dict[str, float] = {}
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            row = dict(zip(header, parts))
            # B R3 schema: sector_tag, phi_origin, b, deflection_rad, delay, m_ledger
            # older: mechanism column
            sector = row.get("sector_tag", "unknown")
            origin = row.get("phi_origin", row.get("mechanism", "unknown"))
            mech = row.get("mechanism") or f"{sector}__{origin}"
            groups[mech].append(row)
            sector_of[mech] = sector
            origin_of[mech] = origin
            m_of[mech] = float(row.get("m_ledger", row.get("M_ledger", 1.0)))
    free_def = 0.15
    fd_path = B_OUT / "round3_free_deficit.tsv"
    if fd_path.exists():
        with open(fd_path) as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and "free_deficit" in parts[0]:
                    try:
                        free_def = float(parts[1])
                    except ValueError:
                        pass
    # also read round3_result.json for free deficit / gates
    meta_json = load_b_round3_json()
    if meta_json.get("free_deficit_core") is not None:
        free_def = float(meta_json["free_deficit_core"])

    out = []
    for mech, rows in groups.items():
        b, dly, alf = [], [], []
        for r in rows:
            bb = float(r.get("b", 0))
            if bb > 0:
                b.append(bb)
                dly.append(float(r.get("delay", 0)))
                alf.append(float(r.get("deflection_rad", r.get("alpha", 0))))
        if not b:
            continue
        sector = sector_of[mech]
        origin = origin_of[mech]
        is_dual = "dualist" in sector
        is_monist = ("monist" in sector) and not is_dual
        is_post = "postulat" in origin.lower() or "hand" in origin.lower()
        is_free_origin = (
            "free" in origin.lower()
            or "capacity" in origin.lower()
            or "green" in origin.lower()
            or "relaxation" in origin.lower()
        )
        gravity = None
        if is_dual or "poisson" in origin.lower():
            gravity = "poisson_3d_tagged_dualist"
        if is_monist and is_free_origin:
            gravity = None
        channel = (
            "monist_3d_free"
            if is_monist and not is_post
            else (
                "dualist_3d_poisson"
                if is_dual
                else ("postulated_kernel" if is_post else "other")
            )
        )
        tag = f"B_R3_{origin}" if origin != "unknown" else mech
        out.append(
            RayMap(
                source=str(path),
                tag=tag,
                b=b,
                delay=dly,
                alpha=alf,
                M_ledger=m_of[mech],
                free_deficit_core=free_def,
                gravity_solver=gravity,
                design=f"B_R3_{origin}",
                meta={
                    "channel": channel,
                    "sector_tag": sector,
                    "phi_origin": origin,
                    "dynamical_free_response": is_monist and is_free_origin,
                    "dimension": 3,
                    "ground_truth_sector": sector,
                    "synthetic": False,
                    "b_gates": meta_json.get("gates"),
                    "b_verdict": meta_json.get("verdict"),
                },
            )
        )
    return out


def load_b_round3_json(path: Path = B_OUT / "round3_result.json") -> Dict:
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def synthetic_3d_maps(
    M: float = 3.1673,
    gamma: float = 0.08,
    eps: float = 0.25,
    free_deficit: float = 0.1513,
) -> List[RayMap]:
    """
    Synthetic 3D free-Green monist truth + dualist Poisson twin + controls.
    Used when B has not exported round3 yet (O-003 allows).
    """
    b = linspace(0.5, 4.0, 8)
    # Monist 3D free Green truth
    d_m, a_m, p_m, _ = monist_3d_green_rays(M, gamma, eps, b)
    monist = RayMap(
        source="synthetic_3d_free_green",
        tag="SYN_M2_3D_free_Laplace",
        b=list(b),
        delay=list(d_m),
        alpha=list(a_m),
        M_ledger=M,
        free_deficit_core=free_deficit,
        gravity_solver=None,
        design="synthetic_monist_3d_free_green",
        meta={
            "channel": "monist_3d_free",
            "sector_tag": "monist_1sector",
            "dynamical_free_response": True,
            "phi_origin": "free_relaxation_3d_green",
            "dimension": 3,
            "ground_truth_sector": "monist_1sector",
            "synthetic": True,
            "gamma_true": gamma,
            "eps_true": eps,
            "path_cost_samples": p_m,
        },
    )
    # Dualist 3D Poisson twin (same rays)
    d_d, a_d, p_d, _ = dualist_3d_poisson_rays(M, gamma, eps, b)
    dualist = RayMap(
        source="synthetic_3d_poisson",
        tag="SYN_M3_3D_poisson",
        b=list(b),
        delay=list(d_d),
        alpha=list(a_d),
        M_ledger=M,
        free_deficit_core=free_deficit,
        gravity_solver="poisson_3d",
        design="synthetic_dualist_3d_poisson",
        meta={
            "channel": "dualist_3d_poisson",
            "sector_tag": "dualist_2sector",
            "dynamical_free_response": False,
            "phi_origin": "poisson_solve",
            "dimension": 3,
            "ground_truth_sector": "dualist_2sector",
            "synthetic": True,
            "G_true": gamma,
            "eps_true": eps,
        },
    )
    # 2D log control (wrong dimension for 3D target)
    d_l, a_l, _ = dualist_2d_log_rays(M, 0.12, 0.8, b)
    log2d = RayMap(
        source="synthetic_2d_log_control",
        tag="SYN_2D_log_control",
        b=list(b),
        delay=list(d_l),
        alpha=list(a_l),
        M_ledger=M,
        free_deficit_core=free_deficit,
        gravity_solver="poisson_2d",
        design="synthetic_2d_log",
        meta={
            "channel": "dualist_2d_log",
            "sector_tag": "dualist_2sector",
            "dynamical_free_response": False,
            "dimension": 2,
            "ground_truth_sector": "dualist_2sector",
            "synthetic": True,
        },
    )
    return [monist, dualist, log2d]


# ---------------------------------------------------------------------------
# Fit suite on one map
# ---------------------------------------------------------------------------
def fit_all(rm: RayMap, occ: Occam) -> List[ModelScore]:
    b, data_d, data_a = rm.b, rm.delay, rm.alpha
    Mled = rm.M_ledger
    mp_data = multipole_ratio(b, data_a)
    scores: List[ModelScore] = []

    def add(sc: ModelScore):
        scores.append(sc)

    # --- monist 3D free Green ---
    def obj_m3(p):
        d, a, _, _ = monist_3d_green_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_m3, L_m3 = multi_start(
        obj_m3,
        seeds=[
            [Mled, 0.08, 0.25],
            [Mled, 0.1, 0.5],
            [Mled, 0.05, 0.15],
            [Mled, 0.12, 0.4],
            [Mled * 0.9, 0.08, 0.3],
            [Mled, 0.15, 0.2],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.05, 3.0)],
        step0=[0.15, 0.015, 0.08],
    )
    d_m3, a_m3, _, M_m3 = monist_3d_green_rays(p_m3[0], p_m3[1], p_m3[2], b)
    mp_m3 = multipole_ratio(b, a_m3)
    is_dyn = bool(rm.meta.get("dynamical_free_response"))
    add(
        score_model(
            "monist_3d_free_green",
            L_m3,
            M_m3,
            Mled,
            mp_m3,
            n_sectors=1,
            has_link=True,
            softE=False,
            postulated=not is_dyn and rm.meta.get("channel") == "postulated_kernel",
            claims_1r=True,
            gravity_solver_present=False,
            params={"M": p_m3[0], "gamma": p_m3[1], "eps": p_m3[2]},
            occ=occ,
            notes="3D free Green / F1 path-cost; 1 sector",
        )
    )
    if is_dyn:
        # alias without postulate for dynamical maps
        sc = scores[-1]
        sc.name = "monist_3d_dynamical_free"
        sc.notes = "dynamical 3D free-response (no postulate penalty)"
        sc.L_extra = max(0.0, sc.L_extra - occ.lambda_postulate)
        sc.S = sc.L_fit + sc.L_occ + sc.L_M + sc.L_softE + sc.L_extra

    # --- dualist 3D Poisson (iso to monist 3D green) ---
    def obj_p3(p):
        d, a, _, _ = dualist_3d_poisson_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_p3, L_p3 = multi_start(
        obj_p3,
        seeds=[
            [Mled, 0.08, 0.25],
            [Mled, 0.1, 0.5],
            [Mled, 0.05, 0.15],
            [Mled, 0.12, 0.4],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.05, 3.0)],
        step0=[0.15, 0.015, 0.08],
    )
    d_p3, a_p3, _, M_p3 = dualist_3d_poisson_rays(p_p3[0], p_p3[1], p_p3[2], b)
    add(
        score_model(
            "dualist_3d_poisson",
            L_p3,
            M_p3,
            Mled,
            multipole_ratio(b, a_p3),
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            claims_1r=True,
            gravity_solver_present=True,
            params={"M": p_p3[0], "G": p_p3[1], "eps": p_p3[2]},
            occ=occ,
            notes="3D Poisson adversary; ray-isomorphic to monist 3D green",
        )
    )

    # --- dualist 2D log (non-iso for 3D 1/r truth) ---
    def obj_log(p):
        d, a, _ = dualist_2d_log_rays(p[0], p[1], p[2], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_log, L_log = multi_start(
        obj_log,
        seeds=[
            [Mled, 0.05, 0.5],
            [Mled, 0.1, 1.0],
            [Mled, 0.2, 0.8],
            [Mled, 0.08, 0.3],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.05, 5.0)],
        step0=[0.15, 0.02, 0.12],
    )
    d_log, a_log, _ = dualist_2d_log_rays(p_log[0], p_log[1], p_log[2], b)
    add(
        score_model(
            "dualist_2d_log",
            L_log,
            p_log[0],
            Mled,
            multipole_ratio(b, a_log),
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            claims_1r=False,
            gravity_solver_present=True,
            params={"M": p_log[0], "G": p_log[1], "a": p_log[2]},
            occ=occ,
            notes="NON-iso vs 3D 1/r (wrong dimension multipole)",
        )
    )

    # --- Yukawa ---
    def obj_yk(p):
        d, a, _ = dualist_yukawa_3d(p[0], p[1], p[2], p[3], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_yk, L_yk = multi_start(
        obj_yk,
        seeds=[
            [Mled, 0.08, 0.0, 0.25],
            [Mled, 0.1, 0.3, 0.5],
            [Mled, 0.12, 0.8, 0.3],
            [Mled, 0.08, 0.1, 0.2],
        ],
        bounds=[(0.1, 20.0), (1e-4, 2.0), (0.0, 3.0), (0.05, 3.0)],
        step0=[0.15, 0.02, 0.08, 0.08],
    )
    d_yk, a_yk, _ = dualist_yukawa_3d(p_yk[0], p_yk[1], p_yk[2], p_yk[3], b)
    add(
        score_model(
            "dualist_yukawa_3d",
            L_yk,
            p_yk[0],
            Mled,
            multipole_ratio(b, a_yk),
            n_sectors=2,
            has_link=False,
            softE=False,
            postulated=False,
            claims_1r=False,
            gravity_solver_present=True,
            params={"M": p_yk[0], "G": p_yk[1], "mu": p_yk[2], "eps": p_yk[3]},
            occ=occ,
            notes="screened; μ→0 near-iso to 3D 1/r",
        )
    )

    # --- compact GRIN control ---
    def obj_lo(p):
        d, a, _ = local_grin_compact(p[0], p[1], b)
        return l_fit_rays(d, a, data_d, data_a)

    p_lo, L_lo = multi_start(
        obj_lo,
        seeds=[[0.35, 1.2], [0.2, 1.5], [0.5, 1.0], [0.15, 2.0]],
        bounds=[(1e-3, 0.94), (0.2, 5.0)],
        step0=[0.05, 0.1],
    )
    d_lo, a_lo, M_lo = local_grin_compact(p_lo[0], p_lo[1], b)
    add(
        score_model(
            "monist_local_GRIN_compact",
            L_lo,
            M_lo,
            Mled,
            multipole_ratio(b, a_lo),
            n_sectors=1,
            has_link=True,
            softE=False,
            postulated=False,
            claims_1r=False,
            gravity_solver_present=False,
            params={"A": p_lo[0], "sigma": p_lo[1]},
            occ=occ,
            notes="compact control; should lose on 1/r truth",
        )
    )

    # --- softE: dualist 3D Poisson labeled monist ---
    add(
        score_model(
            "softE_3d_poisson_as_monist",
            L_p3,
            M_p3,
            Mled,
            multipole_ratio(b, a_p3),
            n_sectors=1,
            has_link=True,
            softE=True,
            postulated=False,
            claims_1r=True,
            gravity_solver_present=True,
            params={"M": p_p3[0], "G": p_p3[1], "eps": p_p3[2]},
            occ=occ,
            notes="Poisson labeled monist → λ_softE",
        )
    )
    # fix softE S to L_fit + softE only (honest cheat)
    cheat = scores[-1]
    cheat.L_occ = 0.0
    cheat.L_M = 0.0
    cheat.L_extra = 0.0
    cheat.S = L_p3 + occ.lambda_softE

    return scores


def nc_checklist_3d(rm: RayMap, winner: ModelScore, scores: List[ModelScore]) -> Dict:
    mp = multipole_ratio(rm.b, rm.alpha)
    long_range_1r = abs(mp - 0.25) < 0.08
    long_range_any = mp > 0.12
    monist_3d = next(
        (s for s in scores if "monist_3d" in s.name or s.name == "monist_3d_dynamical_free"),
        None,
    )
    dual_3d = next((s for s in scores if s.name == "dualist_3d_poisson"), None)
    dual_log = next((s for s in scores if s.name == "dualist_2d_log"), None)
    compact = next((s for s in scores if "local_GRIN" in s.name), None)

    checks = {
        "data_multipole_ratio_b4_b1": mp,
        "NC_free_deficit_positive": rm.free_deficit_core > 0.01,
        "NC_M_ledger_positive": rm.M_ledger > 0,
        "NC_rays_nonzero": max(abs(a) for a in rm.alpha) > 1e-4,
        "NC_long_range_any": long_range_any,
        "NC_long_range_1r_class": long_range_1r,
        "NC_winner_one_sector": winner.N_sectors == 1,
        "NC_winner_has_link": winner.has_free_bound_link,
        "NC_winner_not_softE": winner.L_softE < 1.0,
        "NC_gravity_solver_null_on_monist_data": (
            rm.gravity_solver is None
            if "monist" in rm.meta.get("ground_truth_sector", "")
            else True
        ),
        "NC_dynamical_free_response": bool(rm.meta.get("dynamical_free_response")),
        "NC_dimension_3": rm.meta.get("dimension") == 3,
        "NC_monist_3d_beats_2d_log_fit": (
            monist_3d is not None
            and dual_log is not None
            and monist_3d.L_fit + 1e-9 < dual_log.L_fit
        ),
        "NC_monist_3d_beats_compact_fit": (
            monist_3d is not None
            and compact is not None
            and monist_3d.L_fit + 1e-9 < compact.L_fit
        ),
        "NC_iso_poisson_L_fit_tie": (
            monist_3d is not None
            and dual_3d is not None
            and abs(monist_3d.L_fit - dual_3d.L_fit) < 0.05
        ),
        "NC_monist_beats_poisson_on_S": (
            monist_3d is not None
            and dual_3d is not None
            and monist_3d.S < dual_3d.S
        ),
        "ground_truth_sector": rm.meta.get("ground_truth_sector"),
        "winner": winner.name,
        "winner_S": winner.S,
        "winner_L_fit": winner.L_fit,
        "synthetic": bool(rm.meta.get("synthetic")),
    }
    # C package for 3D monist 1/r
    checks["C_package_3d_1r_candidate"] = (
        checks["NC_free_deficit_positive"]
        and checks["NC_long_range_1r_class"]
        and checks["NC_winner_one_sector"]
        and checks["NC_winner_has_link"]
        and checks["NC_winner_not_softE"]
        and checks["NC_dimension_3"]
        and checks["NC_monist_beats_poisson_on_S"]
    )
    checks["C_package_monist_complete_partial"] = (
        checks["C_package_3d_1r_candidate"]
        and checks["NC_dynamical_free_response"]
        and checks["NC_monist_3d_beats_2d_log_fit"]
    )
    # Full complete still needs inertia triad etc.
    checks["C_package_monist_complete_full"] = False  # inertia not in R3 maps
    checks["false_positive_risk"] = (
        rm.meta.get("ground_truth_sector") == "dualist_2sector"
        and winner.N_sectors == 1
        and "monist" in winner.name
    )
    return checks


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run() -> dict:
    occ = Occam()
    b_files = poll_b_round3()
    maps = load_b_round3_tsv()
    b_json = load_b_round3_json()
    used_synthetic = False
    if not maps:
        maps = synthetic_3d_maps()
        used_synthetic = True

    report = {
        "round": 3,
        "occam": asdict(occ),
        "b_round3_files_polled": b_files,
        "used_synthetic": used_synthetic,
        "b_round3_json_keys": list(b_json.keys()) if b_json else [],
        "maps_found": [
            {
                "tag": m.tag,
                "channel": m.meta.get("channel"),
                "sector_tag": m.meta.get("sector_tag"),
                "M_ledger": m.M_ledger,
                "synthetic": m.meta.get("synthetic"),
                "dimension": m.meta.get("dimension"),
                "dynamical": m.meta.get("dynamical_free_response"),
                "multipole_ratio": multipole_ratio(m.b, m.alpha),
            }
            for m in maps
        ],
        "evaluations": [],
    }

    for rm in maps:
        scores = fit_all(rm, occ)
        scores_sorted = sorted(scores, key=lambda s: s.S)
        winner = scores_sorted[0]
        fit_winner = min(scores, key=lambda s: s.L_fit)
        nc = nc_checklist_3d(rm, winner, scores)
        report["evaluations"].append(
            {
                "map_tag": rm.tag,
                "channel": rm.meta.get("channel"),
                "ground_truth_sector": rm.meta.get("ground_truth_sector"),
                "M_ledger": rm.M_ledger,
                "free_deficit_core": rm.free_deficit_core,
                "score_winner": winner.name,
                "fit_winner": fit_winner.name,
                "scores": [asdict(s) for s in scores_sorted],
                "nc_checklist": nc,
            }
        )

    report["summary"] = make_summary(report, used_synthetic)
    return report


def make_summary(report: dict, used_synthetic: bool) -> dict:
    evals = report["evaluations"]
    monist_ev = next(
        (
            e
            for e in evals
            if e.get("ground_truth_sector") == "monist_1sector"
            or e.get("channel") == "monist_3d_free"
        ),
        None,
    )
    dual_ev = next(
        (
            e
            for e in evals
            if e.get("ground_truth_sector") == "dualist_2sector"
            and "poisson" in e.get("map_tag", "").lower()
        ),
        None,
    )
    log_ev = next(
        (e for e in evals if "log" in e.get("map_tag", "").lower()),
        None,
    )

    def get_score(ev, name_substr):
        if not ev:
            return None
        for s in ev["scores"]:
            if name_substr in s["name"]:
                return s
        return None

    out = {
        "used_synthetic": used_synthetic,
        "monist_map": monist_ev["map_tag"] if monist_ev else None,
        "monist_score_winner": monist_ev["score_winner"] if monist_ev else None,
        "monist_fit_winner": monist_ev["fit_winner"] if monist_ev else None,
        "monist_nc": monist_ev["nc_checklist"] if monist_ev else None,
        "dualist_map": dual_ev["map_tag"] if dual_ev else None,
        "dualist_score_winner": dual_ev["score_winner"] if dual_ev else None,
        "dualist_false_positive": (
            dual_ev["nc_checklist"].get("false_positive_risk") if dual_ev else None
        ),
        "partial_goal2": False,
        "goal2_met": False,
        "goal2_note": "",
    }

    if monist_ev:
        nc = monist_ev["nc_checklist"]
        m3 = get_score(monist_ev, "monist_3d")
        d3 = get_score(monist_ev, "dualist_3d_poisson")
        dlog = get_score(monist_ev, "dualist_2d_log")
        out["monist_3d_L_fit"] = m3["L_fit"] if m3 else None
        out["monist_3d_S"] = m3["S"] if m3 else None
        out["dualist_3d_L_fit"] = d3["L_fit"] if d3 else None
        out["dualist_3d_S"] = d3["S"] if d3 else None
        out["dualist_log_L_fit"] = dlog["L_fit"] if dlog else None
        # Partial goal 2: 3D monist 1/r candidate passes C_package_monist_complete_partial
        if nc.get("C_package_monist_complete_partial"):
            out["partial_goal2"] = True
            out["goal2_note"] = (
                "Monist 3D free-response map passes fit+Occam+C 1/r multipole "
                "candidate package (dynamical, 1 sector, beats 2D log on L_fit, "
                "beats dualist Poisson on S). FULL goal (2) still needs inertia "
                "triad + non-synthetic B export + theory congruence J4 closed."
            )
        elif nc.get("C_package_3d_1r_candidate"):
            out["partial_goal2"] = True
            out["goal2_note"] = (
                "3D 1/r monist candidate package PASS on fit+Occam multipole; "
                "dynamical flag or B export still weak — partial only."
            )
        else:
            out["goal2_note"] = (
                "Monist 3D map did not clear C_package_3d_1r_candidate."
            )
        if used_synthetic:
            out["goal2_note"] += (
                " DATA ARE SYNTHETIC (B round3 not present at D run) — "
                "partial goal (2) is conditional on B confirming same multipole "
                "with monist_1sector free-relaxation tags."
            )
            # synthetic alone cannot fully claim partial goal(2) for project
            out["partial_goal2_project_claim"] = False
            out["partial_goal2_synthetic_ready"] = out["partial_goal2"]
        else:
            out["partial_goal2_project_claim"] = out["partial_goal2"]
            out["partial_goal2_synthetic_ready"] = False

    if dual_ev and dual_ev["nc_checklist"].get("false_positive_risk"):
        out["goal2_note"] += (
            " WARNING: dualist 3D Poisson data still false-awards monist 1-sector "
            "on rays+Occam alone (iso Green) — sector_tag/gravity_solver mandatory."
        )

    out["goal2_met"] = False  # full goal never without triad + real B + theory
    return out


def save(report: dict) -> None:
    with open(OUT / "round3_result.json", "w") as f:
        json.dump(report, f, indent=2)

    with open(OUT / "round3_scores.tsv", "w") as f:
        f.write(
            "map\tchannel\tground_truth\tmodel\tL_fit\tL_occ\tL_M\tL_softE\t"
            "L_extra\tS\tN_sectors\thas_link\tpred_M\tM_ledger\tmultipole_ratio\n"
        )
        for ev in report["evaluations"]:
            for s in ev["scores"]:
                f.write(
                    f"{ev['map_tag']}\t{ev['channel']}\t{ev['ground_truth_sector']}\t"
                    f"{s['name']}\t{s['L_fit']:.8e}\t{s['L_occ']:.8e}\t{s['L_M']:.8e}\t"
                    f"{s['L_softE']:.8e}\t{s['L_extra']:.8e}\t{s['S']:.8e}\t"
                    f"{s['N_sectors']}\t{int(s['has_free_bound_link'])}\t"
                    f"{s['pred_M']:.8e}\t{ev['M_ledger']:.8e}\t"
                    f"{s['multipole_ratio']:.6f}\n"
                )

    with open(OUT / "round3_winners.tsv", "w") as f:
        f.write(
            "map\tchannel\tground_truth\tscore_winner\tfit_winner\twinner_S\t"
            "winner_L_fit\tmultipole_data\tC_3d_1r\tC_partial_complete\t"
            "false_positive_risk\tsynthetic\n"
        )
        for ev in report["evaluations"]:
            nc = ev["nc_checklist"]
            win = next(s for s in ev["scores"] if s["name"] == ev["score_winner"])
            f.write(
                f"{ev['map_tag']}\t{ev['channel']}\t{ev['ground_truth_sector']}\t"
                f"{ev['score_winner']}\t{ev['fit_winner']}\t{win['S']:.6e}\t"
                f"{win['L_fit']:.6e}\t{nc.get('data_multipole_ratio_b4_b1')}\t"
                f"{nc.get('C_package_3d_1r_candidate')}\t"
                f"{nc.get('C_package_monist_complete_partial')}\t"
                f"{nc.get('false_positive_risk')}\t{nc.get('synthetic')}\n"
            )

    with open(OUT / "round3_nc_checklist.tsv", "w") as f:
        f.write("map\tcheck\tvalue\n")
        for ev in report["evaluations"]:
            for k, v in ev["nc_checklist"].items():
                f.write(f"{ev['map_tag']}\t{k}\t{v}\n")

    # rays export
    with open(OUT / "round3_rays_ingest.tsv", "w") as f:
        f.write("map\tchannel\tb\tdelay\talpha\tM_ledger\tsynthetic\n")
        # reconstruct from evaluations is incomplete; re-poll maps
        maps = load_b_round3_tsv()
        if not maps:
            maps = synthetic_3d_maps()
        for m in maps:
            for i, b in enumerate(m.b):
                f.write(
                    f"{m.tag}\t{m.meta.get('channel')}\t{b:.4f}\t"
                    f"{m.delay[i]:.8e}\t{m.alpha[i]:.8e}\t{m.M_ledger:.6f}\t"
                    f"{m.meta.get('synthetic')}\n"
                )

    sm = report["summary"]
    lines = [
        "v76 Approach D Round 3 — 3D free-response congruence",
        f"used_synthetic = {sm.get('used_synthetic')}",
        f"b_round3_files = {report.get('b_round3_files_polled')}",
        f"monist_map = {sm.get('monist_map')}",
        f"monist_score_winner = {sm.get('monist_score_winner')}",
        f"monist_3d L_fit={sm.get('monist_3d_L_fit')} S={sm.get('monist_3d_S')}",
        f"dualist_3d L_fit={sm.get('dualist_3d_L_fit')} S={sm.get('dualist_3d_S')}",
        f"dualist_log L_fit={sm.get('dualist_log_L_fit')}",
        f"dualist_false_positive = {sm.get('dualist_false_positive')}",
        f"partial_goal2 = {sm.get('partial_goal2')}",
        f"partial_goal2_project_claim = {sm.get('partial_goal2_project_claim')}",
        f"goal2_met = {sm.get('goal2_met')}",
        f"goal2_note = {sm.get('goal2_note')}",
        "",
        "Per-map winners:",
    ]
    for ev in report["evaluations"]:
        lines.append(
            f"  {ev['map_tag']}: score={ev['score_winner']} fit={ev['fit_winner']} "
            f"1r={ev['nc_checklist'].get('NC_long_range_1r_class')} "
            f"fp={ev['nc_checklist'].get('false_positive_risk')}"
        )
    text = "\n".join(lines) + "\n"
    (OUT / "round3_summary.txt").write_text(text)
    print(text)


def main():
    print("=== v76 Approach D congruence_score_r3 (3D) ===")
    print("B round3 poll:", poll_b_round3())
    report = run()
    save(report)
    print(f"Wrote under {OUT}")
    return report


if __name__ == "__main__":
    main()
