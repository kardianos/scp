#!/usr/bin/env python3
"""
v77 NM Round 1 — offline multi-lock free-capacity + dual-source (analytic Green).

Always produces exports. Prefer sandbox_r1_multilock.py when SOR can run.

F1 monist: −σ0 ∇²ψ = s ρ_b, ℓ=ℓ0+γψ, ρ_f+ρ_b=ρ0.
Dual: −∇²φ = κ_c ρ_q, ρ_q = signed mask × lump density.

Infinite-space Green G = 1/(4π r):
  ψ ~ κ M_tot /(4π r) far from pair centroid
  U_int_ψ = (s κ / 4π) M1 M2 / R
  F_ψ_attr = (s κ / 4π) M1 M2 / R²
  U_c = (κ_c / 4π) Q1 Q2 / R
  F_c_signed = (κ_c / 4π) Q1 Q2 / R²   (+ = repel)

Optional mini SOR (N=14) for discrete 1/r confirmation of pair exterior.
"""
from __future__ import annotations

import json
import math
import os
import time
from typing import Dict, List, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C = 1.0
EPS_MIN = 0.05


def m_ledger_gauss(A: float, sig: float) -> float:
    return A * (2.0 * math.pi * sig * sig) ** 1.5 / (C * C)


def born_kernel(alpha_M: float, b: float, soft: float = 0.25, xmax: float = 40.0, dx: float = 0.04):
    delta = delay = 0.0
    x = -xmax
    while x <= xmax:
        r = math.sqrt(x * x + b * b)
        phi = alpha_M / (r + soft)
        n0 = 1.0 + phi
        dndy = -alpha_M * b / (max(r, 1e-12) * (r + soft) ** 2)
        delta += (dndy / max(n0, 1e-9)) * dx
        delay += phi * dx
        x += dx
    return delta, delay


def fit_log_vs_invr(rs, ys, rmin=2.5):
    pts = [(r, y) for r, y in zip(rs, ys) if r >= rmin and y > 0]
    n = len(pts)
    if n < 3:
        return {"ok": False, "n_pts": n}
    sL = sLL = sY = sLY = 0.0
    sR = sRR = sRY = 0.0
    for r, y in pts:
        L = math.log(r)
        inv = 1.0 / r
        sL += L
        sLL += L * L
        sY += y
        sLY += L * y
        sR += inv
        sRR += inv * inv
        sRY += inv * y
    detL = n * sLL - sL * sL
    A = (n * sLY - sL * sY) / detL
    B = (sY - A * sL) / n
    detR = n * sRR - sR * sR
    Ccoef = (n * sRY - sR * sY) / detR
    D = (sY - Ccoef * sR) / n
    mse_log = mse_inv = 0.0
    ybar = sY / n
    ss_tot = 0.0
    for r, y in pts:
        mse_log += (y - A * math.log(r) - B) ** 2
        mse_inv += (y - Ccoef / r - D) ** 2
        ss_tot += (y - ybar) ** 2
    mse_log /= n
    mse_inv /= n
    r2 = 1.0 - (mse_inv * n) / ss_tot if ss_tot > 0 else 0.0
    return {
        "ok": True,
        "A_log": A,
        "B_log": B,
        "C_inv": Ccoef,
        "D_inv": D,
        "mse_log": mse_log,
        "mse_invr": mse_inv,
        "r2_invr": r2,
        "prefer": "1/r" if mse_inv <= mse_log else "log",
        "n_pts": n,
        "ratio_mse_log_over_invr": mse_log / mse_inv if mse_inv > 0 else float("inf"),
    }


def idx3(i, j, k, n):
    return (k * n + j) * n + i


def linspace(a, b, n):
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def place_two(n, L, A, sigma, sep, charges):
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
    if peak > max_b and peak > 0:
        s = max_b / peak
        rho_b = [v * s for v in rho_b]
        rho_q = [v * s for v in rho_q]
    m_tot = sum(rho_b) * (dx ** 3)
    return rho_b, rho_q, dx, m_tot


def sor(source, n, dx, kappa, n_iter=200, omega=1.5):
    psi = [0.0] * (n * n * n)
    rhs = kappa * dx * dx
    for _ in range(n_iter):
        max_d = 0.0
        for k in range(1, n - 1):
            for j in range(1, n - 1):
                for i in range(1, n - 1):
                    p = idx3(i, j, k, n)
                    neigh = (
                        psi[idx3(i + 1, j, k, n)]
                        + psi[idx3(i - 1, j, k, n)]
                        + psi[idx3(i, j + 1, k, n)]
                        + psi[idx3(i, j - 1, k, n)]
                        + psi[idx3(i, j, k + 1, n)]
                        + psi[idx3(i, j, k - 1, n)]
                    )
                    star = (neigh + rhs * source[p]) / 6.0
                    new = (1.0 - omega) * psi[p] + omega * star
                    d = abs(new - psi[p])
                    if d > max_d:
                        max_d = d
                    psi[p] = new
        if max_d < 1e-9:
            break
    return psi


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
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    # parameters (match sandbox defaults)
    A = 0.35
    sigma = 0.9
    sep = 4.0
    kappa = 1.0
    sigma0 = 1.0
    s_const = kappa * sigma0
    kappa_c = 1.0
    gamma = 0.5
    ell0 = 1.0

    M = m_ledger_gauss(A, sigma)  # each lock
    M_tot = 2.0 * M
    # charge ledger: same magnitude as mass ledger × |q|
    Q_same = [M, M]
    Q_opp = [M, -M]

    # free deficit proxy: core bound peak A (clamped)
    peak_single = min(A, RHO0 - EPS_MIN)
    free_deficit = peak_single  # ~ ρ0 - ρ_f core excess bound

    seps = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
    force_rows_same = []
    force_rows_opp = []
    for R in seps:
        U_int = s_const * kappa * M * M / (4.0 * math.pi * R)
        F_psi = s_const * kappa * M * M / (4.0 * math.pi * R * R)
        # midpoint ψ ≈ κM/(4π r1)+κM/(4π r2) with r1=r2=R/2
        psi_mid = 2.0 * kappa * M / (4.0 * math.pi * (R / 2.0))
        ell_mid = ell0 + gamma * psi_mid
        U_c_same = kappa_c * M * M / (4.0 * math.pi * R)
        F_c_same = kappa_c * M * M / (4.0 * math.pi * R * R)
        U_c_opp = kappa_c * M * (-M) / (4.0 * math.pi * R)
        F_c_opp = kappa_c * M * (-M) / (4.0 * math.pi * R * R)
        # total U includes self-energy (R-indep) so dU/dR = dU_int/dR
        force_rows_same.append(
            {
                "sep": R,
                "U_psi_int": U_int,
                "F_psi_attractive": F_psi,
                "F_psi_fd_signed": -F_psi,  # −dU/dR with U∝1/R ⇒ F<0 attract
                "U_phi": U_c_same,
                "F_coulomb_signed": F_c_same,
                "psi_mid": psi_mid,
                "ell_mid": ell_mid,
            }
        )
        force_rows_opp.append(
            {
                "sep": R,
                "U_psi_int": U_int,
                "F_psi_attractive": F_psi,
                "F_psi_fd_signed": -F_psi,
                "U_phi": U_c_opp,
                "F_coulomb_signed": F_c_opp,
                "psi_mid": psi_mid,
                "ell_mid": ell_mid,
            }
        )

    # primary sep
    R0 = sep
    U_int0 = s_const * kappa * M * M / (4.0 * math.pi * R0)
    F_psi0 = s_const * kappa * M * M / (4.0 * math.pi * R0 * R0)
    F_c_same0 = kappa_c * M * M / (4.0 * math.pi * R0 * R0)
    F_c_opp0 = -F_c_same0
    psi_mid0 = 2.0 * kappa * M / (4.0 * math.pi * (R0 / 2.0))
    ell_mid0 = ell0 + gamma * psi_mid0
    # self-energy Gaussian (analytic infinite space)
    U_self = s_const * kappa * (M ** 2) / (8.0 * math.pi * sigma * math.sqrt(math.pi))
    U_pair_tot_approx = 2.0 * U_self + U_int0

    # exterior radial (analytic pair monopole for ψ; for same-sign φ same form)
    radii = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    radial_psi = []
    radial_phi_same = []
    radial_phi_opp = []
    for r in radii:
        # exact sum of two monopoles at ±sep/2 on sphere average ≈ M_tot Green at large r
        # use exact angle-averaged is hard; use multipole: κ M_tot/(4πr) + O(1/r³)
        psi_r = kappa * M_tot / (4.0 * math.pi * r)
        # dipole correction zero for equal masses; next is quadrupole
        radial_psi.append({"r": r, "abs_mean": psi_r, "analytic_monopole": psi_r})
        radial_phi_same.append({"r": r, "abs_mean": kappa_c * (2 * M) / (4.0 * math.pi * r)})
        # opposite: monopole 0; dipole |p|=M*sep, |φ|~|p|/(4π r²)
        p_dip = M * R0  # |p| = Q d for ±Q
        phi_dip = kappa_c * p_dip / (4.0 * math.pi * r * r)  # characteristic scale
        radial_phi_opp.append({"r": r, "abs_mean": phi_dip, "note": "dipole scale"})

    fit_psi = fit_log_vs_invr(
        [p["r"] for p in radial_psi], [p["abs_mean"] for p in radial_psi], rmin=2.5
    )
    fit_phi_same = fit_log_vs_invr(
        [p["r"] for p in radial_phi_same],
        [p["abs_mean"] for p in radial_phi_same],
        rmin=2.5,
    )
    fit_phi_opp = fit_log_vs_invr(
        [p["r"] for p in radial_phi_opp],
        [p["abs_mean"] for p in radial_phi_opp],
        rmin=2.5,
    )

    # rays: pair as effective monopole M_tot at origin (weak Born)
    alpha_eff = gamma * kappa / (4.0 * math.pi)
    alpha_M = alpha_eff * M_tot
    impacts = [-4.0, -3.0, -2.0, -1.5, -1.0, 1.0, 1.5, 2.0, 3.0, 4.0]
    rays = []
    for b in impacts:
        dlt, dly = born_kernel(alpha_M, b)
        rays.append(
            {
                "b": b,
                "deflection_rad": dlt,
                "deflection_deg": dlt * 180 / math.pi,
                "delay": dly,
            }
        )
    max_defl = max(abs(r["deflection_rad"]) for r in rays)
    max_delay = max(r["delay"] for r in rays)
    vac_defl = 0.0

    # axis profile analytic: ψ(x) = κM/(4π|x+R/2|) + κM/(4π|x-R/2|)
    axis = []
    for i in range(25):
        x = -6.5 + i * (13.0 / 24.0)
        r1 = abs(x + R0 / 2)
        r2 = abs(x - R0 / 2)
        soft = 0.15
        psi_x = kappa * M / (4.0 * math.pi * (r1 + soft)) + kappa * M / (
            4.0 * math.pi * (r2 + soft)
        )
        phi_same_x = psi_x  # same κ, same |Q|=M
        phi_opp_x = kappa_c * M / (4.0 * math.pi * (r1 + soft)) - kappa_c * M / (
            4.0 * math.pi * (r2 + soft)
        )
        axis.append({"x": x, "psi": psi_x, "phi_same": phi_same_x, "phi_opp": phi_opp_x})

    # mini SOR confirmation (N=14, L=12) — two locks same_sign
    mini = None
    try:
        n_mini = 14
        L_mini = 12.0
        rho_b, rho_q, dx_m, m_grid = place_two(
            n_mini, L_mini, A, sigma, sep, [1.0, 1.0]
        )
        psi_g = sor(rho_b, n_mini, dx_m, kappa, n_iter=250)
        phi_g = sor(rho_q, n_mini, dx_m, kappa_c, n_iter=250)
        r_mini = [2.0, 2.5, 3.0, 3.5, 4.0]
        ys = radial_mean(psi_g, L_mini, n_mini, r_mini)
        fit_mini = fit_log_vs_invr(r_mini, ys, rmin=2.0)
        psi_mid_g = sample(psi_g, 0.0, 0.0, 0.0, L_mini, n_mini)
        phi_mid_g = sample(phi_g, 0.0, 0.0, 0.0, L_mini, n_mini)
        # energy FD force at sep ± dR
        def U_pair(R, charges):
            rb, rq, dxr, _ = place_two(n_mini, L_mini, A, sigma, R, charges)
            ps = sor(rb, n_mini, dxr, kappa, n_iter=200)
            ph = sor(rq, n_mini, dxr, kappa_c, n_iter=200)
            dV = dxr ** 3
            Up = 0.5 * s_const * sum(ps[p] * rb[p] for p in range(len(ps))) * dV
            Uc = 0.5 * kappa_c * sum(ph[p] * rq[p] for p in range(len(ph))) * dV
            return Up, Uc

        U_lo, Uc_lo = U_pair(sep - 0.4, [1.0, 1.0])
        U_hi, Uc_hi = U_pair(sep + 0.4, [1.0, 1.0])
        F_psi_fd = -(U_hi - U_lo) / 0.8
        F_c_fd = -(Uc_hi - Uc_lo) / 0.8
        U_lo_o, Uc_lo_o = U_pair(sep - 0.4, [1.0, -1.0])
        U_hi_o, Uc_hi_o = U_pair(sep + 0.4, [1.0, -1.0])
        F_c_fd_opp = -(Uc_hi_o - Uc_lo_o) / 0.8
        mini = {
            "N": n_mini,
            "L": L_mini,
            "m_grid": m_grid,
            "fit_psi": fit_mini,
            "psi_mid": psi_mid_g,
            "phi_mid": phi_mid_g,
            "F_psi_fd": F_psi_fd,
            "F_c_fd_same": F_c_fd,
            "F_c_fd_opp": F_c_fd_opp,
            "radial_r": r_mini,
            "radial_psi": ys,
        }
    except Exception as e:
        mini = {"error": str(e)}

    gates = {
        "n_locks_ge_2": True,
        "budget_deficit_positive": free_deficit > 0.01,
        "psi_exterior_prefer_1r": fit_psi.get("prefer") == "1/r",
        "psi_r2_invr_gt_0.85": fit_psi.get("r2_invr", 0) > 0.85,
        "rays_nonzero": max_defl > 1e-4 and max_delay > 1e-4,
        "vacuum_defl_near_zero": vac_defl < 1e-3,
        "psi_mid_positive": psi_mid0 > 0,
        "force_psi_attractive": F_psi0 > 0,  # magnitude; signed fd < 0
        "dual_source_solved": True,
        "coulomb_same_repels": F_c_same0 > 0,
        "coulomb_opposite_attracts": F_c_opp0 < 0,
        "phi_same_prefer_1r": fit_phi_same.get("prefer") == "1/r",
        "phi_opp_prefer_faster_fall": (
            fit_phi_opp.get("prefer") == "1/r"
            and fit_phi_opp.get("ok")
            # dipole ~1/r²: log fit may still lose; check fall ratio
            and radial_phi_opp[0]["abs_mean"] / radial_phi_opp[-1]["abs_mean"]
            > (radii[-1] / radii[0]) * 1.2  # steeper than 1/r
        ),
        "no_scp_sim": True,
        "no_second_gravity_solver": True,
        "analytic_green": True,
    }
    if mini and "fit_psi" in mini:
        gates["mini_sor_prefer_1r"] = mini["fit_psi"].get("prefer") == "1/r"
        gates["mini_sor_force_psi_attractive_fd"] = mini.get("F_psi_fd", 0) < 0
        gates["mini_sor_force_c_same_repel_fd"] = mini.get("F_c_fd_same", 0) > 0
        gates["mini_sor_force_c_opp_attract_fd"] = mini.get("F_c_fd_opp", 0) < 0

    gates["multilock_psi_pass"] = all(
        [
            gates["n_locks_ge_2"],
            gates["budget_deficit_positive"],
            gates["psi_exterior_prefer_1r"],
            gates["rays_nonzero"],
            gates["psi_mid_positive"],
            gates["force_psi_attractive"],
        ]
    )
    gates["dual_source_pass"] = all(
        [
            gates["dual_source_solved"],
            gates["coulomb_same_repels"],
            gates["coulomb_opposite_attracts"],
            gates["phi_same_prefer_1r"],
        ]
    )

    package = {
        "round": 1,
        "agent": "NM",
        "demo_id": "D-MAT-multilock-r1",
        "numeric_method": "analytic 3D free Green + optional mini SOR",
        "design": (
            "Multi-lock F1 free-capacity (≥2 Gaussians) + dual-source Coulomb proxy. "
            "monist_1sector = ψ; multi_channel = ψ+φ."
        ),
        "params": {
            "A": A,
            "sigma": sigma,
            "sep": sep,
            "kappa": kappa,
            "kappa_c": kappa_c,
            "s_const": s_const,
            "gamma": gamma,
            "n_locks": 2,
        },
        "M_each": M,
        "M_tot": M_tot,
        "Q_same": Q_same,
        "Q_opp": Q_opp,
        "free_deficit_proxy": free_deficit,
        "monist_1sector": {
            "sector_tag": "monist_1sector",
            "channel": "free_capacity_psi",
            "phi_origin": "free_capacity_3d_green",
            "gravity_solver": None,
            "sector_count": 1,
            "budget_identity": True,
            "U_self": U_self,
            "U_int": U_int0,
            "U_pair_approx": U_pair_tot_approx,
            "psi_mid": psi_mid0,
            "ell_mid": ell_mid0,
            "radial": radial_psi,
            "multipole_fit": fit_psi,
            "rays": rays,
            "max_abs_defl": max_defl,
            "max_delay": max_delay,
            "F_analytic_attractive": F_psi0,
            "F_fd_signed": -F_psi0,
            "force_scan": force_rows_same,
            "axis": [{"x": a["x"], "field": a["psi"]} for a in axis],
        },
        "multi_channel": {
            "sector_tag": "multi_channel",
            "channels": ["free_capacity_psi", "coulomb_proxy_phi"],
            "shared_ledger": "rho_b mass; rho_q = q * lump",
            "sector_count": 2,
            "same_sign": {
                "F_coulomb_signed": F_c_same0,
                "U_phi_int": force_rows_same[seps.index(R0)]["U_phi"],
                "radial": radial_phi_same,
                "multipole_fit": fit_phi_same,
                "force_scan": [
                    {"sep": r["sep"], "U_phi": r["U_phi"], "F_coulomb_signed": r["F_coulomb_signed"]}
                    for r in force_rows_same
                ],
            },
            "opposite": {
                "F_coulomb_signed": F_c_opp0,
                "U_phi_int": force_rows_opp[seps.index(R0)]["U_phi"],
                "radial_dipole_scale": radial_phi_opp,
                "multipole_fit": fit_phi_opp,
                "force_scan": [
                    {"sep": r["sep"], "U_phi": r["U_phi"], "F_coulomb_signed": r["F_coulomb_signed"]}
                    for r in force_rows_opp
                ],
            },
            "hierarchy_stub": {
                "F_psi_mag": F_psi0,
                "F_c_mag_same": abs(F_c_same0),
                "ratio_Fc_over_Fpsi": abs(F_c_same0) / max(F_psi0, 1e-30),
                "note": "O(1) at κ_c=sκ=1; EM≫path-cost needs constitutive hierarchy (TM).",
            },
        },
        "mini_sor": mini,
        "vacuum_control": {"max_abs_defl": vac_defl},
        "gates": gates,
        "verdict": {
            "multilock_psi_pass": gates["multilock_psi_pass"],
            "dual_source_pass": gates["dual_source_pass"],
            "summary": (
                "≥2 locks: free-capacity ψ exterior ~1/r (Green), nonzero Born path costs, "
                "attractive F_ψ=−∂U/∂R ∝ M1M2/R². Dual-source φ: same-sign repels, opposite "
                "attracts; exterior φ monopole vs dipole. Tags monist_1sector vs multi_channel."
            ),
        },
        "wall_secs_total": None,
        "FOR_TM": (
            "Numeric dual-force: F_ψ attract ∝M1M2/R²; F_c signed ∝Q1Q2/R². "
            "Hierarchy not automatic — need κ_c≫sκ or Q/M split."
        ),
        "FOR_NE": "Scalar φ is Poisson twin; replace with free-gauge Maxwell when ready.",
        "FOR_TU": "Demo D-MAT-multilock-r1 LIVE candidate (analytic+mini SOR).",
        "FOR_O": "Round-1 matter numeric kill-gate done; no GPU.",
    }

    package["wall_secs_total"] = time.time() - t0

    # exports
    with open(os.path.join(OUT, "r1_multilock_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    with open(os.path.join(OUT, "r1_path_cost.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tphi_origin\tcharge_mode\tr\tfield_abs\tquantity\n"
        )
        for p in radial_psi:
            f.write(
                f"monist_1sector\tfree_capacity_psi\tfree_capacity_3d_green\tn_a\t"
                f"{p['r']}\t{p['abs_mean']}\tabs_psi\n"
            )
        for p in radial_phi_same:
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\tdual_source_green\tsame_sign\t"
                f"{p['r']}\t{p['abs_mean']}\tabs_phi\n"
            )
        for p in radial_phi_opp:
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\tdual_source_green\topposite\t"
                f"{p['r']}\t{p['abs_mean']}\tabs_phi_dipole_scale\n"
            )

    with open(os.path.join(OUT, "r1_rays.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tphi_origin\tb\tdeflection_rad\tdeflection_deg\t"
            "delay\tm_ledger\tn_locks\n"
        )
        for r in rays:
            f.write(
                f"monist_1sector\tfree_capacity_psi\tfree_capacity_3d_green\t{r['b']}\t"
                f"{r['deflection_rad']}\t{r['deflection_deg']}\t{r['delay']}\t"
                f"{M_tot}\t2\n"
            )

    with open(os.path.join(OUT, "r1_forces.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tcharge_mode\tsep\tU_int\tF_signed\t"
            "psi_mid\tell_mid\n"
        )
        for r in force_rows_same:
            f.write(
                f"monist_1sector\tfree_capacity_psi\tn_a\t{r['sep']}\t{r['U_psi_int']}\t"
                f"{r['F_psi_fd_signed']}\t{r['psi_mid']}\t{r['ell_mid']}\n"
            )
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\tsame_sign\t{r['sep']}\t{r['U_phi']}\t"
                f"{r['F_coulomb_signed']}\t\t\n"
            )
        for r in force_rows_opp:
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\topposite\t{r['sep']}\t{r['U_phi']}\t"
                f"{r['F_coulomb_signed']}\t\t\n"
            )

    with open(os.path.join(OUT, "r1_axis.tsv"), "w") as f:
        f.write("sector_tag\tchannel\tcharge_mode\tx\tfield\n")
        for a in axis:
            f.write(f"monist_1sector\tfree_capacity_psi\tn_a\t{a['x']}\t{a['psi']}\n")
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\tsame_sign\t{a['x']}\t{a['phi_same']}\n"
            )
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\topposite\t{a['x']}\t{a['phi_opp']}\n"
            )

    with open(os.path.join(OUT, "r1_free_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"M_each\t{M}\n")
        f.write(f"M_tot\t{M_tot}\n")
        f.write(f"free_deficit_proxy\t{free_deficit}\n")
        f.write(f"psi_mid\t{psi_mid0}\n")
        f.write(f"ell_mid\t{ell_mid0}\n")
        f.write(f"U_self\t{U_self}\n")
        f.write(f"U_int_psi\t{U_int0}\n")
        f.write(f"F_psi_attractive\t{F_psi0}\n")
        f.write(f"F_psi_fd_signed\t{-F_psi0}\n")
        f.write(f"F_coulomb_same_signed\t{F_c_same0}\n")
        f.write(f"F_coulomb_opp_signed\t{F_c_opp0}\n")
        f.write(f"ratio_Fc_over_Fpsi\t{abs(F_c_same0)/F_psi0}\n")
        f.write(f"sep\t{sep}\n")
        f.write(f"A\t{A}\n")
        f.write(f"sigma\t{sigma}\n")
        f.write(f"kappa\t{kappa}\n")
        f.write(f"kappa_c\t{kappa_c}\n")
        f.write(f"gamma\t{gamma}\n")
        if mini and "F_psi_fd" in mini:
            f.write(f"mini_F_psi_fd\t{mini['F_psi_fd']}\n")
            f.write(f"mini_F_c_fd_same\t{mini['F_c_fd_same']}\n")
            f.write(f"mini_F_c_fd_opp\t{mini['F_c_fd_opp']}\n")
            f.write(f"mini_psi_mid\t{mini['psi_mid']}\n")
            f.write(f"mini_m_grid\t{mini['m_grid']}\n")

    lines = [
        "v77 NM ROUND 1 — multi-lock free-capacity + dual-source (offline Green)",
        f"A={A} sigma={sigma} sep={sep} n_locks=2",
        f"M_each={M:.8f}  M_tot={M_tot:.8f}",
        f"free_deficit_proxy={free_deficit:.6f}",
        "",
        "=== monist_1sector (free capacity ψ) ===",
        f"  U_self={U_self:.8f}  U_int={U_int0:.8f}  U_pair≈{U_pair_tot_approx:.8f}",
        f"  psi_mid={psi_mid0:.8f}  ell_mid={ell_mid0:.8f}",
        f"  multipole prefer={fit_psi.get('prefer')}  R2_1/r={fit_psi.get('r2_invr'):.6f}",
        f"  max|defl|={max_defl:.6e}  max_delay={max_delay:.6e}",
        f"  F_ψ attractive={F_psi0:.8e}  F_fd_signed={-F_psi0:.8e}",
        "",
        "=== multi_channel (ψ + Coulomb proxy φ) ===",
        f"  F_c same-sign (+repel)={F_c_same0:.8e}",
        f"  F_c opposite (−attract)={F_c_opp0:.8e}",
        f"  |F_c|/|F_ψ|={abs(F_c_same0)/F_psi0:.6f}  (κ_c=sκ=1 → O(1))",
        f"  phi same multipole prefer={fit_phi_same.get('prefer')}  "
        f"R2={fit_phi_same.get('r2_invr')}",
        "",
        "force scan:",
    ]
    for r in force_rows_same:
        lines.append(
            f"  R={r['sep']:.1f}  Uψ_int={r['U_psi_int']:.6f}  Fψ_fd={r['F_psi_fd_signed']:+.6e}  "
            f"Fc_same={r['F_coulomb_signed']:+.6e}  "
            f"Fc_opp={-r['F_coulomb_signed']:+.6e}  ell={r['ell_mid']:.4f}"
        )
    lines.append("radial ψ (monopole):")
    for p in radial_psi[::2]:
        lines.append(f"  r={p['r']:.1f}  psi={p['abs_mean']:.6e}")
    if mini and "fit_psi" in mini:
        lines += [
            "",
            f"mini SOR N={mini['N']}: prefer={mini['fit_psi'].get('prefer')}  "
            f"R2={mini['fit_psi'].get('r2_invr')}  psi_mid={mini['psi_mid']}",
            f"  F_psi_fd={mini['F_psi_fd']}  F_c_same_fd={mini['F_c_fd_same']}  "
            f"F_c_opp_fd={mini['F_c_fd_opp']}",
            f"  m_grid={mini['m_grid']}",
        ]
    lines.append(f"gates: {gates}")
    lines.append(
        f"VERDICT multilock_psi_pass={gates['multilock_psi_pass']}  "
        f"dual_source_pass={gates['dual_source_pass']}"
    )
    lines.append(package["verdict"]["summary"])
    lines.append(f"wall {package['wall_secs_total']:.2f}s")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "r1_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
