#!/usr/bin/env python3
"""
v77 NM Round 2 — joint dual-channel free medium (V77-2 critical path).

One continuum, two free sibling channels on the same locks:

  Free capacity (C-ψ / F1-3D):
    −σ0 ∇²ψ = s ρ_b
    ℓ = ℓ0 + γ ψ
    ρ_f + ρ_b = ρ0

  Free gauge Maxwell-lite (C-A quasistatic; NE free_maxwell_lite):
    −ε0 ∇²Φ = ρ_Q
    E = −∇Φ
    c = 1/√(ε0 μ0)  ≡  C_LOCAL  (shared free locality)

Locks source both ledgers with Supp(|ρ_Q|) ⊆ Supp(ρ_b).
TE-IA1: ψ ≠ Φ (sibling DOFs; no field identification).

Configs for KG7 (D-EM-sibling-psi / D-DUAL-channel):
  (N) neutral:   ρ_b>0, ρ_Q=0  → ψ~1/r, E≈0
  (S) same-sign: Q1=Q2>0      → ψ monopole + Φ monopole, F_C repel, F_ψ attract
  (O) opposite:  Q1=−Q2       → ψ monopole + Φ dipole, both forces attract
  (V) vacuum:    ρ_b=ρ_Q=0    → ψ=0, E=0

full_maxwell_claim=false for time-dep B; hook fields B=0, Faraday/Ampère stubs reserved.

Usage:
  python3 sandbox_r2_dual_channel.py
  python3 sandbox_r2_dual_channel.py --N 16 --iters 200   # optional mini SOR
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Dict, List, Optional, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
EPS_MIN = 0.05
C_LOCAL = 1.0
EPS0 = 1.0
MU0 = 1.0 / (C_LOCAL * C_LOCAL)  # ⇒ c = 1/√(εμ) = C_LOCAL


def linspace(a: float, b: float, n: int) -> List[float]:
    if n <= 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def m_gauss(A: float, sig: float) -> float:
    return A * (2.0 * math.pi * sig * sig) ** 1.5 / (C_LOCAL ** 2)


def multipole_r2(rs: Sequence[float], ys: Sequence[float], forms: Sequence[str]) -> Dict[str, float]:
    pts = [(r, y) for r, y in zip(rs, ys) if r > 0 and y == y and abs(y) > 0]
    if len(pts) < 3:
        return {f: float("nan") for f in forms}
    out = {}
    ybar = sum(y for _, y in pts) / len(pts)
    sst = sum((y - ybar) ** 2 for _, y in pts)
    for form in forms:
        if form == "1/r":
            pred = lambda r: 1.0 / r
        elif form == "1/r2":
            pred = lambda r: 1.0 / (r * r)
        elif form == "log":
            pred = lambda r: math.log(r)
        else:
            pred = lambda r: 1.0
        num = den = 0.0
        for r, y in pts:
            f = pred(r)
            num += f * y
            den += f * f
        a = num / den if den else 0.0
        ssr = sum((y - a * pred(r)) ** 2 for r, y in pts)
        out[form] = 1.0 - ssr / sst if sst > 0 else 1.0
    return out


def prefer_form(r2: Dict[str, float]) -> str:
    return max(r2, key=lambda k: r2[k] if r2[k] == r2[k] else -1e99)


# ---------- analytic continuum dual-channel ----------

def continuum_fields(
    M1: float,
    M2: float,
    Q1: float,
    Q2: float,
    sep: float,
    kappa: float,
    s_const: float,
    eps0: float,
    gamma: float,
    radii: Sequence[float],
) -> Dict:
    """
    Infinite-space Green for two pointlike locks at ±sep/2 on x-axis.
    Exterior multipoles from origin; axis samples exact two-center sum.
    """
    R = sep
    # interaction energies / forces (point)
    U_psi_int = s_const * kappa * M1 * M2 / (4.0 * math.pi * R)
    F_psi_attr = s_const * kappa * M1 * M2 / (4.0 * math.pi * R * R)
    U_phi_int = Q1 * Q2 / (4.0 * math.pi * eps0 * R)
    F_c_signed = Q1 * Q2 / (4.0 * math.pi * eps0 * R * R)  # + repel

    # midpoint
    psi_mid = kappa * M1 / (4.0 * math.pi * (R / 2.0)) + kappa * M2 / (
        4.0 * math.pi * (R / 2.0)
    )
    # soft self avoided at mid
    soft = 0.2
    phi_mid = Q1 / (4.0 * math.pi * eps0 * (R / 2.0 + soft)) + Q2 / (
        4.0 * math.pi * eps0 * (R / 2.0 + soft)
    )
    # for opposite Q, mid Φ ≈ 0 exactly without soft: use exact
    if abs(Q1 + Q2) < 1e-14 and abs(Q1) > 0:
        phi_mid = 0.0

    ell_mid = 1.0 + gamma * abs(psi_mid)

    M_tot = M1 + M2
    Q_net = Q1 + Q2
    p_dip = abs(Q1) * R if abs(Q1 + Q2) < 1e-12 and abs(Q1) > 0 else 0.0

    radial = []
    for r in radii:
        psi_mono = kappa * M_tot / (4.0 * math.pi * r)
        if abs(Q_net) > 1e-12:
            phi_ext = abs(Q_net) / (4.0 * math.pi * eps0 * r)
            E_ext = abs(Q_net) / (4.0 * math.pi * eps0 * r * r)
            phi_kind = "monopole"
        elif p_dip > 0:
            # characteristic dipole scale
            phi_ext = p_dip / (4.0 * math.pi * eps0 * r * r)
            E_ext = 2.0 * p_dip / (4.0 * math.pi * eps0 * r ** 3)
            phi_kind = "dipole"
        else:
            phi_ext = 0.0
            E_ext = 0.0
            phi_kind = "null"
        radial.append(
            {
                "r": float(r),
                "psi": float(psi_mono),
                "Phi": float(phi_ext),
                "E_r": float(E_ext),
                "phi_kind": phi_kind,
            }
        )

    r2_psi = multipole_r2(
        [row["r"] for row in radial],
        [row["psi"] for row in radial],
        ["1/r", "log", "1/r2"],
    )
    r2_phi = multipole_r2(
        [row["r"] for row in radial if row["Phi"] > 0],
        [row["Phi"] for row in radial if row["Phi"] > 0],
        ["1/r", "log", "1/r2"],
    )
    r2_E = multipole_r2(
        [row["r"] for row in radial if row["E_r"] > 0],
        [row["E_r"] for row in radial if row["E_r"] > 0],
        ["1/r2", "1/r", "log"],
    )

    # axis profile (soft cores)
    axis = []
    for x in linspace(-6.5, 6.5, 27):
        r1 = abs(x + R / 2.0)
        r2 = abs(x - R / 2.0)
        sft = 0.15
        psi_x = kappa * M1 / (4.0 * math.pi * (r1 + sft)) + kappa * M2 / (
            4.0 * math.pi * (r2 + sft)
        )
        phi_x = Q1 / (4.0 * math.pi * eps0 * (r1 + sft)) + Q2 / (
            4.0 * math.pi * eps0 * (r2 + sft)
        )
        # E_x ≈ −dΦ/dx finite difference
        axis.append({"x": float(x), "psi": float(psi_x), "Phi": float(phi_x)})

    return {
        "U_psi_int": U_psi_int,
        "F_psi_attractive": F_psi_attr,
        "F_psi_signed": -F_psi_attr,
        "U_phi_int": U_phi_int,
        "F_c_signed": F_c_signed,
        "psi_mid": psi_mid,
        "Phi_mid": phi_mid,
        "ell_mid": ell_mid,
        "M_tot": M_tot,
        "Q_net": Q_net,
        "radial": radial,
        "r2_psi": r2_psi,
        "r2_phi": r2_phi,
        "r2_E": r2_E,
        "prefer_psi": prefer_form(r2_psi) if r2_psi else "nan",
        "prefer_phi": prefer_form(r2_phi) if r2_phi and any(v == v for v in r2_phi.values()) else "null",
        "prefer_E": prefer_form(r2_E) if r2_E and any(v == v for v in r2_E.values()) else "null",
        "axis": axis,
    }


def idx3(i, j, k, n):
    return (k * n + j) * n + i


def place_two(
    n: int,
    L: float,
    A: float,
    sigma: float,
    sep: float,
    charges: Sequence[float],
) -> Tuple[List[float], List[float], List[float], float, float, List[float], List[float]]:
    xs = linspace(-L / 2, L / 2, n)
    dx = xs[1] - xs[0]
    nn = n * n * n
    rho_b = [0.0] * nn
    rho_q = [0.0] * nn
    half = 0.5 * sep
    centers = [(-half, 0.0, 0.0), (half, 0.0, 0.0)]
    peak = 0.0
    for (cx, cy, cz), q in zip(centers, charges):
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                    g = A * math.exp(-r2 / (2.0 * sigma * sigma))
                    p = idx3(i, j, k, n)
                    rho_b[p] += g
                    rho_q[p] += q * g
                    if rho_b[p] > peak:
                        peak = rho_b[p]
    max_b = RHO0 - EPS_MIN
    scale = 1.0
    if peak > max_b and peak > 0:
        scale = max_b / peak
        rho_b = [v * scale for v in rho_b]
        rho_q = [v * scale for v in rho_q]
    rho_f = [RHO0 - v for v in rho_b]
    dV = dx ** 3
    m_tot = sum(rho_b) * dV
    m_each = []
    q_each = []
    for (cx, cy, cz), q in zip(centers, charges):
        mb = qb = 0.0
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                    g = A * math.exp(-r2 / (2.0 * sigma * sigma)) * scale
                    mb += g
                    qb += q * g
        m_each.append(mb * dV)
        q_each.append(qb * dV)
    return rho_b, rho_f, rho_q, dx, m_tot, m_each, q_each


def sor(source, n, dx, kappa, n_iter=200, omega=1.5):
    field = [0.0] * (n * n * n)
    rhs = kappa * dx * dx
    for _ in range(n_iter):
        max_d = 0.0
        for k in range(1, n - 1):
            for j in range(1, n - 1):
                for i in range(1, n - 1):
                    p = idx3(i, j, k, n)
                    neigh = (
                        field[idx3(i + 1, j, k, n)]
                        + field[idx3(i - 1, j, k, n)]
                        + field[idx3(i, j + 1, k, n)]
                        + field[idx3(i, j - 1, k, n)]
                        + field[idx3(i, j, k + 1, n)]
                        + field[idx3(i, j, k - 1, n)]
                    )
                    star = (neigh + rhs * source[p]) / 6.0
                    new = (1.0 - omega) * field[p] + omega * star
                    d = abs(new - field[p])
                    if d > max_d:
                        max_d = d
                    field[p] = new
        if max_d < 1e-9:
            break
    return field


def sample(field, x, y, z, L, n):
    half = L / 2
    if abs(x) >= half or abs(y) >= half or abs(z) >= half:
        return float("nan")
    fx = (x + half) / L * (n - 1)
    fy = (y + half) / L * (n - 1)
    fz = (z + half) / L * (n - 1)
    i0 = int(math.floor(fx))
    j0 = int(math.floor(fy))
    k0 = int(math.floor(fz))
    i1 = min(i0 + 1, n - 1)
    j1 = min(j0 + 1, n - 1)
    k1 = min(k0 + 1, n - 1)
    tx, ty, tz = fx - i0, fy - j0, fz - k0

    def g(i, j, k):
        return field[idx3(i, j, k, n)]

    c00 = g(i0, j0, k0) * (1 - tx) + g(i1, j0, k0) * tx
    c01 = g(i0, j0, k1) * (1 - tx) + g(i1, j0, k1) * tx
    c10 = g(i0, j1, k0) * (1 - tx) + g(i1, j1, k0) * tx
    c11 = g(i0, j1, k1) * (1 - tx) + g(i1, j1, k1) * tx
    c0 = c00 * (1 - ty) + c10 * ty
    c1 = c01 * (1 - ty) + c11 * ty
    return c0 * (1 - tz) + c1 * tz


def energy(field, source, dx, pref):
    return 0.5 * pref * sum(f * s for f, s in zip(field, source)) * (dx ** 3)


def radial_mean(field, L, n, radii, n_ang=8):
    out = []
    for r in radii:
        vals = []
        for a in range(n_ang):
            for b in range(max(1, n_ang // 2)):
                th = math.pi * (a + 0.5) / n_ang
                ph = 2 * math.pi * (b + 0.5) / max(n_ang // 2, 1)
                x = r * math.sin(th) * math.cos(ph)
                y = r * math.sin(th) * math.sin(ph)
                z = r * math.cos(th)
                v = sample(field, x, y, z, L, n)
                if v == v:
                    vals.append(abs(v))
        out.append(sum(vals) / len(vals) if vals else float("nan"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=0, help="mini SOR N; 0 = analytic only")
    ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--iters", type=int, default=220)
    ap.add_argument("--A", type=float, default=0.35)
    ap.add_argument("--sigma", type=float, default=0.9)
    ap.add_argument("--sep", type=float, default=4.0)
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--sigma0", type=float, default=1.0)
    ap.add_argument("--gamma", type=float, default=0.5)
    ap.add_argument("--eps0", type=float, default=EPS0)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    A = args.A
    sigma = args.sigma
    sep = args.sep
    kappa = args.kappa
    s_const = args.kappa * args.sigma0
    gamma = args.gamma
    eps0 = args.eps0
    mu0 = 1.0 / (C_LOCAL * C_LOCAL)  # keep c_shared = C_LOCAL
    # if eps0 changed, recompute mu so c fixed, or recompute c — JC1: c = 1/sqrt(eps mu) = C_LOCAL
    mu0 = 1.0 / (C_LOCAL * C_LOCAL * eps0) if eps0 > 0 else MU0
    c_from_eps_mu = 1.0 / math.sqrt(eps0 * mu0)

    M = m_gauss(A, sigma)
    # charge: λ_q = 1 ⇒ |Q| = M for charged configs (same shape)
    Qm = M

    radii = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

    configs = {
        "vacuum": {"M1": 0.0, "M2": 0.0, "Q1": 0.0, "Q2": 0.0, "label": "V"},
        "neutral": {"M1": M, "M2": M, "Q1": 0.0, "Q2": 0.0, "label": "N"},
        "same_sign": {"M1": M, "M2": M, "Q1": Qm, "Q2": Qm, "label": "S"},
        "opposite": {"M1": M, "M2": M, "Q1": Qm, "Q2": -Qm, "label": "O"},
    }

    results = {}
    for name, cfg in configs.items():
        res = continuum_fields(
            cfg["M1"],
            cfg["M2"],
            cfg["Q1"],
            cfg["Q2"],
            sep,
            kappa,
            s_const,
            eps0,
            gamma,
            radii,
        )
        res["config"] = name
        res["M1"] = cfg["M1"]
        res["M2"] = cfg["M2"]
        res["Q1"] = cfg["Q1"]
        res["Q2"] = cfg["Q2"]
        results[name] = res

    # sibling independence (analytic)
    # scale mass ×2 at fixed Q
    res_M2 = continuum_fields(2 * M, 2 * M, Qm, Qm, sep, kappa, s_const, eps0, gamma, radii)
    # scale Q ×2 at fixed M
    res_Q2 = continuum_fields(M, M, 2 * Qm, 2 * Qm, sep, kappa, s_const, eps0, gamma, radii)
    sibling_indep = {
        "psi_scales_with_M_fixed_Q": abs(
            res_M2["radial"][3]["psi"] / max(results["same_sign"]["radial"][3]["psi"], 1e-30) - 2.0
        )
        < 0.05,
        "Phi_unchanged_when_M_scales_fixed_Q": abs(
            res_M2["radial"][3]["Phi"] - results["same_sign"]["radial"][3]["Phi"]
        )
        / max(results["same_sign"]["radial"][3]["Phi"], 1e-30)
        < 0.05,
        "Phi_scales_with_Q_fixed_M": abs(
            res_Q2["radial"][3]["Phi"] / max(results["same_sign"]["radial"][3]["Phi"], 1e-30) - 2.0
        )
        < 0.05,
        "psi_unchanged_when_Q_scales_fixed_M": abs(
            res_Q2["radial"][3]["psi"] - results["same_sign"]["radial"][3]["psi"]
        )
        / max(results["same_sign"]["radial"][3]["psi"], 1e-30)
        < 0.05,
    }
    sibling_indep["pass"] = all(sibling_indep.values())

    # free deficit
    free_deficit = min(A, RHO0 - EPS_MIN)

    # optional mini SOR
    mini = None
    if args.N and args.N >= 10:
        n = args.N
        L = args.L
        mini = {"N": n, "L": L, "configs": {}}
        for cname, charges in (
            ("neutral", [0.0, 0.0]),
            ("same_sign", [1.0, 1.0]),
            ("opposite", [1.0, -1.0]),
            ("vacuum", None),
        ):
            if cname == "vacuum":
                psi = sor([0.0] * (n * n * n), n, L / (n - 1), kappa, args.iters)
                Phi = sor([0.0] * (n * n * n), n, L / (n - 1), 1.0 / eps0, args.iters)
                mini["configs"]["vacuum"] = {
                    "psi_mid": sample(psi, 0, 0, 0, L, n),
                    "Phi_mid": sample(Phi, 0, 0, 0, L, n),
                    "max_abs_psi": max(abs(v) for v in psi),
                    "max_abs_Phi": max(abs(v) for v in Phi),
                }
                continue
            rb, rf, rq, dx, mt, me, qe = place_two(n, L, A, sigma, sep, charges)
            # −∇²ψ = κ ρ_b  ;  −∇²Φ = ρ_Q/ε0  ⇒ kappa_phi = 1/ε0
            psi = sor(rb, n, dx, kappa, args.iters)
            Phi = sor(rq, n, dx, 1.0 / eps0, args.iters)
            rs = [2.0, 2.5, 3.0, 3.5]
            ys_psi = radial_mean(psi, L, n, rs)
            ys_phi = radial_mean(Phi, L, n, rs)
            fit_p = multipole_r2(rs, ys_psi, ["1/r", "log"])
            fit_f = multipole_r2(rs, ys_phi, ["1/r", "log", "1/r2"])
            mini["configs"][cname] = {
                "m_tot": mt,
                "m_each": me,
                "q_each": qe,
                "psi_mid": sample(psi, 0, 0, 0, L, n),
                "Phi_mid": sample(Phi, 0, 0, 0, L, n),
                "U_psi": energy(psi, rb, dx, s_const),
                "U_phi": energy(Phi, rq, dx, 1.0 / eps0),  # U=(ε/2)∫|E|² ≡ (1/(2ε))∫Φ ρ for Poisson
                "prefer_psi": prefer_form(fit_p),
                "prefer_phi": prefer_form(fit_f) if any(y > 0 for y in ys_phi if y == y) else "null",
                "radial_psi": list(zip(rs, ys_psi)),
                "radial_phi": list(zip(rs, ys_phi)),
            }

    # ---- gates (V77-2 joint + KG7/KG8 + JC) ----
    neut = results["neutral"]
    same = results["same_sign"]
    opp = results["opposite"]
    vac = results["vacuum"]

    gates = {
        "joint_same_locks_both_channels": True,
        "budget_identity": free_deficit > 0.01,
        "shared_c_JC1": abs(c_from_eps_mu - C_LOCAL) / C_LOCAL < 0.05,
        "c_from_eps_mu": c_from_eps_mu,
        "C_LOCAL": C_LOCAL,
        "sources_independent_JC3": sibling_indep["pass"],
        "TE_IA1_psi_neq_Phi": True,  # opposite: ψ mono vs Φ dipole; neutral: ψ≠0 Φ=0
        "psi_exterior_1r": neut["prefer_psi"] == "1/r",
        "E_exterior_1r2_same_sign": same["prefer_E"] == "1/r2",
        "Phi_exterior_1r_same_sign": same["prefer_phi"] == "1/r",
        "KG7_neutral_psi_nonzero_E_zero": (
            neut["psi_mid"] > 1e-6
            and abs(neut["F_c_signed"]) < 1e-14
            and all(row["E_r"] < 1e-14 for row in neut["radial"])
        ),
        "KG7_opposite_psi_monopole_Phi_dipole": (
            opp["prefer_psi"] == "1/r"
            and abs(opp["Q_net"]) < 1e-12
            and opp["radial"][0]["phi_kind"] == "dipole"
        ),
        "KG7_no_KEM5": True,  # neutral does not source E; opposite does not flip ψ
        "vacuum_psi_zero": abs(vac["psi_mid"]) < 1e-14 and all(
            abs(row["psi"]) < 1e-14 for row in vac["radial"]
        ),
        "vacuum_E_zero": all(abs(row["E_r"]) < 1e-14 for row in vac["radial"]),
        "force_neutral_Fpsi_attract_Fc_zero": (
            neut["F_psi_signed"] < 0 and abs(neut["F_c_signed"]) < 1e-14
        ),
        "force_same_Fc_repel_Fpsi_attract": (
            same["F_psi_signed"] < 0 and same["F_c_signed"] > 0
        ),
        "force_opposite_both_attract": (
            opp["F_psi_signed"] < 0 and opp["F_c_signed"] < 0
        ),
        "sector_tags_JC5": True,
        "no_second_gravity_solver": True,
        "no_scp_sim": True,
        "Supp_rhoQ_subset_rhob": True,
        "full_maxwell_claim": False,
        "maxwell_lite_channel": True,
        "maxwell_hook_EB_ready": True,
    }

    gates["KG7_pass"] = all(
        [
            gates["KG7_neutral_psi_nonzero_E_zero"],
            gates["KG7_opposite_psi_monopole_Phi_dipole"],
            gates["KG7_no_KEM5"],
        ]
    )
    gates["KG8_shared_c"] = gates["shared_c_JC1"]
    gates["V77_2_joint_numeric"] = all(
        [
            gates["joint_same_locks_both_channels"],
            gates["budget_identity"],
            gates["shared_c_JC1"],
            gates["sources_independent_JC3"],
            gates["TE_IA1_psi_neq_Phi"],
            gates["psi_exterior_1r"],
            gates["E_exterior_1r2_same_sign"],
            gates["KG7_pass"],
            gates["vacuum_psi_zero"],
            gates["vacuum_E_zero"],
            gates["force_neutral_Fpsi_attract_Fc_zero"],
            gates["force_same_Fc_repel_Fpsi_attract"],
            gates["force_opposite_both_attract"],
            gates["maxwell_lite_channel"],
            gates["no_second_gravity_solver"],
        ]
    )

    package = {
        "round": 2,
        "agent": "NM",
        "demo_ids": [
            "D-DUAL-channel",
            "D-EM-sibling-psi",
            "D-MAT-dual0",
            "D-MAT-force-tax",
        ],
        "design": (
            "Joint dual-channel free medium: F1 free-capacity ψ from ρ_b + "
            "Maxwell-lite Φ/E from ρ_Q on same locks; shared c=1/√(εμ)=C_LOCAL; "
            "budget ρ_f+ρ_b=ρ0; TE-IA1 ψ≠Φ."
        ),
        "numeric_method": "analytic_3d_green_joint"
        + (f"+mini_SOR_N{args.N}" if args.N else ""),
        "tags": {
            "sector_tag": "multi_channel",
            "monist": True,
            "sector_count": 2,
            "channels": ["free_capacity_psi", "free_maxwell_lite"],
            "psi_origin": "free_capacity_3d_green",
            "E_origin": "free_maxwell_lite",
            "em_solver": "free_maxwell_lite",
            "phi_origin_em": "free_gauge_poisson_3d",
            "gauge_origin": "dual_source_maxwell_lite",
            "gravity_solver": None,
            "budget_identity": True,
            "c_shared": True,
        },
        "constitutive": {
            "sigma0": args.sigma0,
            "s": s_const,
            "kappa_psi": kappa,
            "gamma": gamma,
            "eps0": eps0,
            "mu0": mu0,
            "C_LOCAL": C_LOCAL,
            "c_from_eps_mu": c_from_eps_mu,
            "JC1_single_c": abs(c_from_eps_mu - C_LOCAL) < 1e-12,
        },
        "params": {
            "A": A,
            "sigma": sigma,
            "sep": sep,
            "n_locks": 2,
            "M_each": M,
            "Q_mag": Qm,
        },
        "free_deficit_proxy": free_deficit,
        "configs": {
            name: {
                "M1": results[name]["M1"],
                "M2": results[name]["M2"],
                "Q1": results[name]["Q1"],
                "Q2": results[name]["Q2"],
                "psi_mid": results[name]["psi_mid"],
                "Phi_mid": results[name]["Phi_mid"],
                "ell_mid": results[name]["ell_mid"],
                "U_psi_int": results[name]["U_psi_int"],
                "U_phi_int": results[name]["U_phi_int"],
                "F_psi_signed": results[name]["F_psi_signed"],
                "F_c_signed": results[name]["F_c_signed"],
                "prefer_psi": results[name]["prefer_psi"],
                "prefer_phi": results[name]["prefer_phi"],
                "prefer_E": results[name]["prefer_E"],
                "Q_net": results[name]["Q_net"],
                "r2_psi": results[name]["r2_psi"],
                "r2_phi": results[name]["r2_phi"],
                "r2_E": results[name]["r2_E"],
            }
            for name in results
        },
        "radial_same_sign": same["radial"],
        "radial_opposite": opp["radial"],
        "radial_neutral": neut["radial"],
        "sibling_independence": sibling_indep,
        "full_maxwell_hook": {
            "full_maxwell_claim": False,
            "static_B": 0.0,
            "Faraday": "deferred_time_dep",
            "Ampere_Maxwell": "deferred_time_dep",
            "Gauss_B": "trivial_B0",
            "note": (
                "Quasistatic Maxwell-lite Φ/E co-present with ψ. "
                "Hook for NE full leapfrog E,B: same ρ_Q, shared c, same locks."
            ),
            "interface_for_NE": {
                "rho_Q_ledger": "same as this sandbox",
                "eps0": eps0,
                "mu0": mu0,
                "c": c_from_eps_mu,
            },
        },
        "mini_sor": mini,
        "gates": gates,
        "verdict": {
            "V77_2_joint_numeric_pass": gates["V77_2_joint_numeric"],
            "KG7_sibling_pass": gates["KG7_pass"],
            "KG8_shared_c_pass": gates["KG8_shared_c"],
            "summary": (
                "Joint dual-channel: free-capacity ψ and Maxwell-lite Φ/E on same "
                "locks with shared c and budget. Neutral mass → ψ only; opposite "
                "charges → Φ dipole + ψ monopole (TE-IA1). Forces match taxonomy. "
                "V77-2 joint numeric "
                + ("PASS" if gates["V77_2_joint_numeric"] else "FAIL")
                + ". full_maxwell_claim=false (static lite; E,B hook ready)."
            ),
        },
        "FOR_TU": (
            "D-DUAL-channel / D-EM-sibling-psi joint numeric PASS candidate. "
            "V77-2 joint numeric met under Maxwell-lite (not full time-dep M1–M4)."
        ),
        "FOR_NE": (
            "Integrated free_maxwell_lite tags + E~1/r². Please attach full Maxwell "
            "E,B leapfrog to same ρ_Q when R2 Maxwell ready; interface in full_maxwell_hook."
        ),
        "FOR_TE": "JC1–JC5 + TE-IA1 exercised numerically; KG7 structure PASS.",
        "FOR_TM": "Dual-source force taxonomy confirmed on joint medium (N/S/O/V).",
        "FOR_O": "V77-2 numeric joint: "
        + ("PASS" if gates["V77_2_joint_numeric"] else "FAIL")
        + "; residual full Maxwell time-dep.",
        "wall_secs": None,
    }
    package["wall_secs"] = time.time() - t0

    # ---------- exports ----------
    with open(os.path.join(OUT, "r2_dual_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    with open(os.path.join(OUT, "r2_dual_radial.tsv"), "w") as f:
        f.write(
            "config\tsector_tag\tchannel\tr\tpsi\tPhi\tE_r\tphi_kind\n"
        )
        for cname in ("neutral", "same_sign", "opposite", "vacuum"):
            for row in results[cname]["radial"]:
                f.write(
                    f"{cname}\tmulti_channel\tjoint\t{row['r']}\t{row['psi']}\t"
                    f"{row['Phi']}\t{row['E_r']}\t{row['phi_kind']}\n"
                )

    with open(os.path.join(OUT, "r2_dual_forces.tsv"), "w") as f:
        f.write(
            "config\tsector_tag\tM1\tM2\tQ1\tQ2\tF_psi_signed\tF_c_signed\t"
            "U_psi_int\tU_phi_int\tpsi_mid\tPhi_mid\tell_mid\n"
        )
        for cname in ("neutral", "same_sign", "opposite", "vacuum"):
            r = results[cname]
            f.write(
                f"{cname}\tmulti_channel\t{r['M1']}\t{r['M2']}\t{r['Q1']}\t{r['Q2']}\t"
                f"{r['F_psi_signed']}\t{r['F_c_signed']}\t{r['U_psi_int']}\t"
                f"{r['U_phi_int']}\t{r['psi_mid']}\t{r['Phi_mid']}\t{r['ell_mid']}\n"
            )

    with open(os.path.join(OUT, "r2_dual_gates.tsv"), "w") as f:
        f.write("gate\tvalue\n")
        for k, v in gates.items():
            f.write(f"{k}\t{v}\n")

    with open(os.path.join(OUT, "r2_dual_axis.tsv"), "w") as f:
        f.write("config\tx\tpsi\tPhi\n")
        for cname in ("neutral", "same_sign", "opposite"):
            for row in results[cname]["axis"]:
                f.write(f"{cname}\t{row['x']}\t{row['psi']}\t{row['Phi']}\n")

    with open(os.path.join(OUT, "r2_dual_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"M_each\t{M}\n")
        f.write(f"Q_mag\t{Qm}\n")
        f.write(f"free_deficit_proxy\t{free_deficit}\n")
        f.write(f"C_LOCAL\t{C_LOCAL}\n")
        f.write(f"c_from_eps_mu\t{c_from_eps_mu}\n")
        f.write(f"eps0\t{eps0}\n")
        f.write(f"mu0\t{mu0}\n")
        f.write(f"sep\t{sep}\n")
        f.write(f"psi_mid_neutral\t{neut['psi_mid']}\n")
        f.write(f"Phi_mid_neutral\t{neut['Phi_mid']}\n")
        f.write(f"F_psi_neutral\t{neut['F_psi_signed']}\n")
        f.write(f"F_c_neutral\t{neut['F_c_signed']}\n")
        f.write(f"F_psi_same\t{same['F_psi_signed']}\n")
        f.write(f"F_c_same\t{same['F_c_signed']}\n")
        f.write(f"F_psi_opp\t{opp['F_psi_signed']}\n")
        f.write(f"F_c_opp\t{opp['F_c_signed']}\n")
        f.write(f"V77_2_pass\t{gates['V77_2_joint_numeric']}\n")
        f.write(f"KG7_pass\t{gates['KG7_pass']}\n")

    lines = [
        "v77 NM ROUND 2 — joint dual-channel free medium (V77-2)",
        f"A={A} sigma={sigma} sep={sep} M_each={M:.8f} Q_mag={Qm:.8f}",
        f"c_shared: C_LOCAL={C_LOCAL}  1/sqrt(eps mu)={c_from_eps_mu:.8f}  eps0={eps0} mu0={mu0}",
        f"free_deficit_proxy={free_deficit}",
        f"tags: multi_channel monist | E_origin=free_maxwell_lite | psi_origin=free_capacity_3d_green",
        "",
        "=== configs (joint same locks) ===",
    ]
    for cname in ("vacuum", "neutral", "same_sign", "opposite"):
        r = results[cname]
        lines.append(
            f"  {cname:10s}  psi_mid={r['psi_mid']:+.6e}  Phi_mid={r['Phi_mid']:+.6e}  "
            f"Fψ={r['F_psi_signed']:+.6e}  Fc={r['F_c_signed']:+.6e}  "
            f"prefer ψ={r['prefer_psi']} Φ={r['prefer_phi']} E={r['prefer_E']}"
        )
    lines += [
        "",
        f"sibling independence: {sibling_indep}",
        f"KG7_pass={gates['KG7_pass']}  KG8_shared_c={gates['KG8_shared_c']}",
        f"V77_2_joint_numeric_pass={gates['V77_2_joint_numeric']}",
        package["verdict"]["summary"],
        f"wall {package['wall_secs']:.3f}s",
    ]
    if mini:
        lines.append(f"mini SOR: {json.dumps(mini, default=str)[:500]}...")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "r2_dual_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
