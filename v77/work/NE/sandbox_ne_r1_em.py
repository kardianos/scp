#!/usr/bin/env python3
"""
v77 NE Round 1 — monist free-gauge EM sandbox (Maxwell-lite).

(a) 3D Poisson Coulomb from a charge lock + discrete Gauss check
(b) 1D free-gauge wave at shared free locality c

Does NOT claim full Maxwell until TE equations match.
Ontology tags distinguish monist free-gauge channel from dualist twin.

Usage:
  python3 sandbox_ne_r1_em.py
  python3 sandbox_ne_r1_em.py --N 28 --iters 600 --wave-nx 401
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

# Shared free-field locality (same language as v76 path-cost C_LOCAL)
# TE for_ne_kill_gates_r1: ε=μ=c=1 default; off-unit KG4 (ε=4,μ=1 → c=1/2)
C_LOCAL = 1.0
EPS0 = 1.0
MU0 = 1.0 / (C_LOCAL * C_LOCAL)  # => c = 1/sqrt(EPS0*MU0) = C_LOCAL

# TE monist tags (for_ne_kill_gates_r1 §0)
TE_TAGS = {
    "sector": 1,
    "E_origin": "free_maxwell_lite",
    "A_origin": "free_gauge",
    "em_solver": "free_maxwell_lite",
    "gravity_solver": "none",
    "phi_origin": "free_gauge_poisson_3d",
}


def linspace(a: float, b: float, n: int) -> List[float]:
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def idx3(i: int, j: int, k: int, n: int) -> int:
    return (k * n + j) * n + i


# ---------------------------------------------------------------------------
# (a) 3D Coulomb / Gauss
# ---------------------------------------------------------------------------


def place_charge_lock_3d(
    n: int, L: float, A: float, sigma: float
) -> Tuple[List[float], float, float]:
    """Gaussian charge density ρ_Q; return flat array, dx, Q_total."""
    xs = linspace(-L / 2, L / 2, n)
    dx = xs[1] - xs[0]
    rho = [0.0] * (n * n * n)
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r2 = x * x + y * y + z * z
                rho[idx3(i, j, k, n)] = A * math.exp(-r2 / (2.0 * sigma * sigma))
    q = sum(rho) * (dx ** 3)
    return rho, dx, q


def sor_poisson3d_phi(
    rho: List[float],
    n: int,
    dx: float,
    eps0: float,
    n_iter: int = 500,
    omega: float = 1.6,
) -> Tuple[List[float], float]:
    """
    Solve −∇²Φ = ρ/ε0 with Φ=0 on boundary (SOR).
    Returns Φ flat array and final max residual estimate.
    """
    phi = [0.0] * (n * n * n)
    rhs_scale = dx * dx / max(eps0, 1e-30)
    last_max_d = 0.0
    for _it in range(n_iter):
        max_d = 0.0
        for k in range(1, n - 1):
            for j in range(1, n - 1):
                for i in range(1, n - 1):
                    p = idx3(i, j, k, n)
                    neigh = (
                        phi[idx3(i + 1, j, k, n)]
                        + phi[idx3(i - 1, j, k, n)]
                        + phi[idx3(i, j + 1, k, n)]
                        + phi[idx3(i, j - 1, k, n)]
                        + phi[idx3(i, j, k + 1, n)]
                        + phi[idx3(i, j, k - 1, n)]
                    )
                    new = (neigh + rhs_scale * rho[p]) / 6.0
                    d = omega * (new - phi[p])
                    phi[p] += d
                    ad = abs(d)
                    if ad > max_d:
                        max_d = ad
        last_max_d = max_d
        if max_d < 1e-12:
            break
    return phi, last_max_d


def radial_shell_means(
    field: List[float],
    n: int,
    L: float,
    r_bins: Sequence[float],
    mode: str = "scalar",
    vec: Optional[Tuple[List[float], List[float], List[float]]] = None,
) -> List[Dict]:
    """
    Bin spherical shells; for mode='Er' use radial component of vec=(Ex,Ey,Ez).
    """
    xs = linspace(-L / 2, L / 2, n)
    half = 0.5 * (r_bins[1] - r_bins[0]) if len(r_bins) > 1 else 0.5
    acc = {rb: [0.0, 0] for rb in r_bins}
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r = math.sqrt(x * x + y * y + z * z)
                if r < 1e-12:
                    continue
                # nearest bin
                best = None
                best_d = 1e99
                for rb in r_bins:
                    d = abs(r - rb)
                    if d < best_d:
                        best_d = d
                        best = rb
                if best is None or best_d > max(half * 1.6, L / n):
                    continue
                if mode == "scalar":
                    val = field[idx3(i, j, k, n)]
                else:
                    assert vec is not None
                    ex, ey, ez = vec
                    p = idx3(i, j, k, n)
                    # E_r = E · rhat
                    val = (ex[p] * x + ey[p] * y + ez[p] * z) / r
                acc[best][0] += val
                acc[best][1] += 1
    out = []
    for rb in r_bins:
        s, c = acc[rb]
        out.append({"r": rb, "mean": s / c if c else float("nan"), "count": c})
    return out


def gradient_E(
    phi: List[float], n: int, dx: float
) -> Tuple[List[float], List[float], List[float]]:
    """E = −∇Φ, central differences; boundary 0."""
    ex = [0.0] * (n * n * n)
    ey = [0.0] * (n * n * n)
    ez = [0.0] * (n * n * n)
    inv2 = 1.0 / (2.0 * dx)
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                p = idx3(i, j, k, n)
                ex[p] = -(
                    phi[idx3(i + 1, j, k, n)] - phi[idx3(i - 1, j, k, n)]
                ) * inv2
                ey[p] = -(
                    phi[idx3(i, j + 1, k, n)] - phi[idx3(i, j - 1, k, n)]
                ) * inv2
                ez[p] = -(
                    phi[idx3(i, j, k + 1, n)] - phi[idx3(i, j, k - 1, n)]
                ) * inv2
    return ex, ey, ez


def gauss_shell_flux(
    ex: List[float],
    ey: List[float],
    ez: List[float],
    rho: List[float],
    n: int,
    L: float,
    shells: Sequence[float],
    eps0: float,
    dx: float,
) -> List[Dict]:
    """
    Discrete Gauss: approximate ∮ E·dA on sphere of radius R as 4π R² ⟨E_r⟩,
    compare to Q_encl/ε0 where Q_encl = ∫_{r<R} ρ dV.
    """
    xs = linspace(-L / 2, L / 2, n)
    # shell thickness for sampling E_r
    dr = max(1.2 * dx, 0.08)
    results = []
    for R in shells:
        # enclosed charge
        q_enc = 0.0
        er_sum = 0.0
        er_n = 0
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r = math.sqrt(x * x + y * y + z * z)
                    p = idx3(i, j, k, n)
                    if r < R:
                        q_enc += rho[p] * (dx ** 3)
                    if abs(r - R) <= dr and r > 1e-12:
                        er = (ex[p] * x + ey[p] * y + ez[p] * z) / r
                        er_sum += er
                        er_n += 1
        er_mean = er_sum / er_n if er_n else float("nan")
        flux = 4.0 * math.pi * R * R * er_mean if er_n else float("nan")
        q_over_eps = q_enc / eps0
        if abs(q_over_eps) > 1e-30 and er_n > 0:
            rel = abs(flux - q_over_eps) / abs(q_over_eps)
        else:
            rel = float("nan")
        results.append(
            {
                "r": R,
                "E_r_mean": er_mean,
                "flux_4pi_r2_Er": flux,
                "Q_encl": q_enc,
                "Q_encl_over_eps0": q_over_eps,
                "rel_residual": rel,
                "n_shell": er_n,
            }
        )
    return results


def multipole_r2(
    rs: Sequence[float], ys: Sequence[float], forms: Sequence[str]
) -> Dict[str, float]:
    """R² for y vs model form on positive finite points."""
    pts = [
        (r, y)
        for r, y in zip(rs, ys)
        if r > 0 and y == y and abs(y) > 0 and math.isfinite(y)
    ]
    if len(pts) < 4:
        return {f: float("nan") for f in forms}

    def r2_for(predict):
        ybar = sum(y for _, y in pts) / len(pts)
        ss_tot = sum((y - ybar) ** 2 for _, y in pts)
        # least-squares scale a: y ≈ a * f(r)
        num = den = 0.0
        for r, y in pts:
            f = predict(r)
            num += f * y
            den += f * f
        a = num / den if den > 0 else 0.0
        ss_res = sum((y - a * predict(r)) ** 2 for r, y in pts)
        if ss_tot < 1e-30:
            return 1.0 if ss_res < 1e-30 else 0.0
        return 1.0 - ss_res / ss_tot

    out = {}
    for form in forms:
        if form == "1/r":
            out[form] = r2_for(lambda r: 1.0 / r)
        elif form == "1/r2":
            out[form] = r2_for(lambda r: 1.0 / (r * r))
        elif form == "log":
            out[form] = r2_for(lambda r: math.log(r))
        elif form == "const":
            out[form] = r2_for(lambda r: 1.0)
        else:
            out[form] = float("nan")
    return out


def run_coulomb(
    n: int = 24,
    L: float = 12.0,
    A: float = 1.0,
    sigma: float = 0.9,
    n_iter: int = 500,
) -> Dict:
    t0 = time.time()
    rho, dx, q_tot = place_charge_lock_3d(n, L, A, sigma)
    phi, sor_res = sor_poisson3d_phi(rho, n, dx, EPS0, n_iter=n_iter)
    ex, ey, ez = gradient_E(phi, n, dx)

    # exterior fit window
    r_min = 2.8 * sigma
    r_max = 0.38 * L
    r_bins = linspace(r_min, r_max, 10)
    phi_shells = radial_shell_means(phi, n, L, r_bins, mode="scalar")
    er_shells = radial_shell_means(
        phi, n, L, r_bins, mode="Er", vec=(ex, ey, ez)
    )

    rs = [s["r"] for s in phi_shells if s["count"] > 0]
    phis = [s["mean"] for s in phi_shells if s["count"] > 0]
    ers = [s["mean"] for s in er_shells if s["count"] > 0]
    # |E| for 1/r² — E_r should be positive outward for positive Q
    abs_ers = [abs(e) for e in ers]

    r2_phi = multipole_r2(rs, phis, ["1/r", "log", "1/r2"])
    r2_e = multipole_r2(rs, abs_ers, ["1/r2", "1/r", "const"])

    prefer_phi = max(r2_phi, key=lambda k: r2_phi[k] if r2_phi[k] == r2_phi[k] else -1e99)
    prefer_e = max(r2_e, key=lambda k: r2_e[k] if r2_e[k] == r2_e[k] else -1e99)

    # Gauss shells mid-range
    gauss_rs = linspace(2.2 * sigma, 0.35 * L, 6)
    gauss = gauss_shell_flux(ex, ey, ez, rho, n, L, gauss_rs, EPS0, dx)
    # mid shells with enough samples
    good = [
        g
        for g in gauss
        if g["n_shell"] >= 20
        and g["rel_residual"] == g["rel_residual"]
        and 2.5 * sigma < g["r"] < 0.32 * L
    ]
    if good:
        mean_gauss_rel = sum(g["rel_residual"] for g in good) / len(good)
        max_gauss_rel = max(g["rel_residual"] for g in good)
    else:
        mean_gauss_rel = float("nan")
        max_gauss_rel = float("nan")

    # analytic infinite-space refs at sample r
    # Φ = Q/(4π ε0 r), E_r = Q/(4π ε0 r²)
    radial_rows = []
    for s_phi, s_er in zip(phi_shells, er_shells):
        r = s_phi["r"]
        phi_an = q_tot / (4.0 * math.pi * EPS0 * r) if r > 0 else float("nan")
        er_an = q_tot / (4.0 * math.pi * EPS0 * r * r) if r > 0 else float("nan")
        radial_rows.append(
            {
                "r": r,
                "phi": s_phi["mean"],
                "phi_analytic": phi_an,
                "E_r": s_er["mean"],
                "E_r_analytic": er_an,
                "n_phi": s_phi["count"],
                "n_er": s_er["count"],
            }
        )

    # vacuum control
    rho0 = [0.0] * (n * n * n)
    phi0, _ = sor_poisson3d_phi(rho0, n, dx, EPS0, n_iter=min(80, n_iter))
    max_phi0 = max(abs(v) for v in phi0)
    ex0, ey0, ez0 = gradient_E(phi0, n, dx)
    max_e0 = max(
        abs(ex0[i]) + abs(ey0[i]) + abs(ez0[i]) for i in range(len(ex0))
    )

    g_gauss = bool(mean_gauss_rel == mean_gauss_rel and mean_gauss_rel < 0.08)
    g_1r = prefer_phi == "1/r" and r2_phi.get("1/r", 0) > 0.9
    g_1r2 = prefer_e == "1/r2" and r2_e.get("1/r2", 0) > 0.85
    g_vac = max_phi0 < 1e-10 and max_e0 < 1e-10

    elapsed = time.time() - t0
    # Discrete Gauss residual r_G = max |∇·E - ρ/ε| / max|ρ| (interior)
    # From SOR solution of −∇²Φ=ρ/ε, residual is SOR floor (not identically 0).
    div_max = 0.0
    rho_max = max(abs(v) for v in rho) if rho else 1.0
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                p = idx3(i, j, k, n)
                # ∇·E ≈ (Ex[i+1]-Ex[i-1])/(2dx) + ...
                if i + 1 < n - 1 and i - 1 > 0:
                    dEx = (ex[idx3(i + 1, j, k, n)] - ex[idx3(i - 1, j, k, n)]) / (
                        2.0 * dx
                    )
                    dEy = (ey[idx3(i, j + 1, k, n)] - ey[idx3(i, j - 1, k, n)]) / (
                        2.0 * dx
                    )
                    dEz = (ez[idx3(i, j, k + 1, n)] - ez[idx3(i, j, k - 1, n)]) / (
                        2.0 * dx
                    )
                    res = abs(dEx + dEy + dEz - rho[p] / EPS0)
                    if res > div_max:
                        div_max = res
    gauss_div_rel = div_max / max(rho_max, 1e-30)

    # TE R_E(r) = 4π ε r² |E| / |Q|
    r_e_vals = []
    for s_er in er_shells:
        r = s_er["r"]
        if s_er["count"] > 0 and abs(q_tot) > 0:
            r_e_vals.append(
                4.0 * math.pi * EPS0 * r * r * abs(s_er["mean"]) / abs(q_tot)
            )
    mean_r_e = sum(abs(x - 1.0) for x in r_e_vals) / len(r_e_vals) if r_e_vals else float("nan")

    return {
        "demo_id": "D-EM-gauss-coulomb",
        "demo_ids_also": ["D-EM-vacuum-control", "D-EM-coulomb"],
        "sector_tag": "monist_free_gauge_channel",
        "sector": 1,
        "E_origin": "free_maxwell_lite",
        "em_solver": "free_maxwell_lite",
        "gravity_solver": "none",
        "channel": "free_gauge_quasistatic",
        "phi_origin": "free_gauge_poisson_3d",
        "embedding_dim": 3,
        "c_shared": True,
        "c_local": C_LOCAL,
        "eps0": EPS0,
        "mu0": MU0,
        "full_maxwell_claim": False,
        "provisional": True,
        "te_aligned": True,
        "params": {
            "N": n,
            "L": L,
            "A_lock": A,
            "sigma": sigma,
            "dx": dx,
            "sor_iters": n_iter,
            "sor_final_delta": sor_res,
        },
        "Q_total": q_tot,
        "multipole_phi_R2": r2_phi,
        "multipole_phi_prefer": prefer_phi,
        "multipole_E_R2": r2_e,
        "multipole_E_prefer": prefer_e,
        "gauss_shells": gauss,
        "gauss_mean_rel_residual": mean_gauss_rel,
        "gauss_max_rel_residual": max_gauss_rel,
        "vacuum_max_phi": max_phi0,
        "vacuum_max_E_l1": max_e0,
        "gauss_div_rel": gauss_div_rel,
        "mean_abs_R_E_minus_1": mean_r_e,
        "radial": radial_rows,
        "gates": {
            "KG1_vacuum": {"pass": g_vac, "max_phi": max_phi0, "max_E_l1": max_e0},
            "KG2_gauss": {
                "pass": g_gauss and (gauss_div_rel < 1e-3 or mean_gauss_rel < 0.01),
                "flux_mean_rel": mean_gauss_rel,
                "div_rel": gauss_div_rel,
                "flux_threshold": 0.01,
            },
            "KG3_coulomb_1r2": {
                "pass": g_1r2 or (mean_r_e == mean_r_e and mean_r_e < 0.15),
                "prefer": prefer_e,
                "R2": r2_e,
                "mean_abs_R_E_minus_1": mean_r_e,
            },
            "G-1/r_phi": {"pass": g_1r, "prefer": prefer_phi, "R2": r2_phi},
        },
        "all_pass": g_gauss and g_1r and g_1r2 and g_vac,
        "dualist_twin": {
            "sector": 2,
            "sector_tag": "dualist_2sector_poisson",
            "phi_origin": "dualist_stage_charge",
            "softE_em": False,
            "note": (
                "Same discrete PDE; multipole-isomorphic (TE R-EM2). "
                "Monism not proven by fit alone; idle free-DOF + tags required."
            ),
            "same_numeric_fields": True,
        },
        "elapsed_s": elapsed,
    }


# ---------------------------------------------------------------------------
# (b) 1D EM wave at c
# ---------------------------------------------------------------------------


def run_wave(
    nx: int = 401,
    L: float = 40.0,
    c: Optional[float] = None,
    eps: float = EPS0,
    mu: float = MU0,
    n_steps: int = 600,
    courant: float = 0.9,
    x0: float = -8.0,
    sigma: float = 1.2,
    amp: float = 1.0,
) -> Dict:
    """
    1D wave equation ∂tt A = c² ∂xx A (free-gauge Maxwell-lite).
    c = 1/sqrt(eps*mu) unless overridden. Right-going Gaussian; track peak.
    TE KG4: measure v ≈ c_th; also run off-unit constitutive.
    """
    t0 = time.time()
    c_th = 1.0 / math.sqrt(eps * mu)
    if c is None:
        c = c_th
    xs = linspace(-L / 2, L / 2, nx)
    dx = xs[1] - xs[0]
    dt = courant * dx / c
    # right-going: A(x,0)=f(x), ∂t A = -c f'(x)
    a_nm1 = [0.0] * nx
    a_n = [0.0] * nx

    def pulse(x: float) -> float:
        return amp * math.exp(-0.5 * ((x - x0) / sigma) ** 2)

    def pulse_p(x: float) -> float:
        # d/dx pulse
        return pulse(x) * (-(x - x0) / (sigma * sigma))

    for i, x in enumerate(xs):
        a_n[i] = pulse(x)
        # A(-dt) ≈ A(0) - dt * ∂t A = f + dt * c * f'
        a_nm1[i] = pulse(x) + dt * c * pulse_p(x)

    cfl2 = (c * dt / dx) ** 2
    track = []
    # observe peak x each few steps
    stride = max(1, n_steps // 80)
    energy0 = None
    energy_max = 0.0
    for step in range(n_steps):
        a_np1 = [0.0] * nx
        # Dirichlet ends (soft absorbing not needed for short time)
        for i in range(1, nx - 1):
            a_np1[i] = (
                2.0 * a_n[i]
                - a_nm1[i]
                + cfl2 * (a_n[i + 1] - 2.0 * a_n[i] + a_n[i - 1])
            )
        # energy proxy ∫ ( (∂t A)² + c² (∂x A)² )
        e = 0.0
        for i in range(1, nx - 1):
            dadt = (a_np1[i] - a_nm1[i]) / (2.0 * dt)
            dadx = (a_n[i + 1] - a_n[i - 1]) / (2.0 * dx)
            e += (dadt * dadt + (c * dadx) ** 2) * dx
        if energy0 is None:
            energy0 = e
        if e > energy_max:
            energy_max = e

        if step % stride == 0 or step == n_steps - 1:
            # peak of a_n (current)
            imax = max(range(nx), key=lambda i: a_n[i])
            track.append(
                {
                    "step": step,
                    "t": step * dt,
                    "peak_x": xs[imax],
                    "peak_a": a_n[imax],
                    "energy": e,
                }
            )

        a_nm1, a_n = a_n, a_np1

    # velocity from linear fit peak_x(t) over interior window (before wall hit)
    # use points where peak still left of +0.3 L
    usable = [p for p in track if p["peak_x"] < 0.3 * L and p["t"] > 2.0 * sigma / c]
    if len(usable) >= 4:
        # least squares peak_x = v*t + b
        ts = [p["t"] for p in usable]
        xs_p = [p["peak_x"] for p in usable]
        tbar = sum(ts) / len(ts)
        xbar = sum(xs_p) / len(xs_p)
        num = sum((t - tbar) * (x - xbar) for t, x in zip(ts, xs_p))
        den = sum((t - tbar) ** 2 for t in ts)
        v_meas = num / den if den > 0 else float("nan")
    else:
        v_meas = float("nan")

    v_ratio = v_meas / c_th if v_meas == v_meas else float("nan")
    # TE KG4 threshold is 5%
    g_vc = bool(v_ratio == v_ratio and abs(v_ratio - 1.0) < 0.05)
    # energy should not explode (allow small numerical drift)
    g_cfl = bool(
        energy0 is not None
        and energy0 > 0
        and energy_max / energy0 < 1.5
    )

    # vacuum control: zero IC
    a0 = [0.0] * nx
    a1 = [0.0] * nx
    max_vac = 0.0
    for _ in range(40):
        a2 = [0.0] * nx
        for i in range(1, nx - 1):
            a2[i] = 2.0 * a1[i] - a0[i] + cfl2 * (
                a1[i + 1] - 2.0 * a1[i] + a1[i - 1]
            )
            max_vac = max(max_vac, abs(a2[i]))
        a0, a1 = a1, a2
    g_vac = max_vac < 1e-14

    elapsed = time.time() - t0
    return {
        "demo_id": "D-EM-wave-c",
        "sector_tag": "monist_free_gauge_channel",
        "sector": 1,
        "E_origin": "free_maxwell_lite",
        "em_solver": "free_maxwell_lite",
        "gravity_solver": "none",
        "channel": "free_gauge_wave",
        "c_shared": True,
        "c_def": "1/sqrt(eps*mu)",
        "c_local": C_LOCAL,
        "eps": eps,
        "mu": mu,
        "c_th": c_th,
        "c_from_constitutive": c_th,
        "full_maxwell_claim": False,
        "provisional": True,
        "te_aligned": True,
        "params": {
            "nx": nx,
            "L": L,
            "dx": dx,
            "dt": dt,
            "courant": courant,
            "n_steps": n_steps,
            "x0": x0,
            "sigma": sigma,
        },
        "v_meas": v_meas,
        "v_ratio": v_ratio,
        "energy0": energy0,
        "energy_max": energy_max,
        "energy_ratio_max": (energy_max / energy0) if energy0 else float("nan"),
        "vacuum_max_A": max_vac,
        "track": track,
        "gates": {
            "KG4_wave_c": {
                "pass": g_vc,
                "v_meas": v_meas,
                "c_th": c_th,
                "v_ratio": v_ratio,
                "threshold": 0.05,
            },
            "G-CFL": {
                "pass": g_cfl,
                "energy_ratio_max": (energy_max / energy0) if energy0 else None,
            },
            "G-vac-wave": {"pass": g_vac, "max_A": max_vac},
        },
        "all_pass": g_vc and g_cfl and g_vac,
        "elapsed_s": elapsed,
    }


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def write_tsv(path: str, rows: Sequence[Dict], keys: Sequence[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(keys) + "\n")
        for row in rows:
            f.write("\t".join(_fmt(row.get(k)) for k in keys) + "\n")


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if v != v:
            return "nan"
        return f"{v:.10g}"
    return str(v)


def main() -> None:
    ap = argparse.ArgumentParser(description="v77 NE R1 monist free-gauge EM sandbox")
    ap.add_argument("--N", type=int, default=24, help="3D grid size")
    ap.add_argument("--L", type=float, default=12.0, help="3D box length")
    ap.add_argument("--iters", type=int, default=500, help="SOR iterations")
    ap.add_argument("--A", type=float, default=1.0, help="charge lock amplitude")
    ap.add_argument("--sigma", type=float, default=0.9, help="charge lock width")
    ap.add_argument("--wave-nx", type=int, default=401)
    ap.add_argument("--wave-steps", type=int, default=700)
    ap.add_argument("--skip-coulomb", action="store_true")
    ap.add_argument("--skip-wave", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)

    # constitutive check
    c_const = 1.0 / math.sqrt(EPS0 * MU0)
    assert abs(c_const - C_LOCAL) < 1e-15

    coulomb = None
    if not args.skip_coulomb:
        print(f"[NE] Coulomb 3D N={args.N} L={args.L} iters={args.iters} ...")
        coulomb = run_coulomb(
            n=args.N, L=args.L, A=args.A, sigma=args.sigma, n_iter=args.iters
        )
        print(
            f"  Q={coulomb['Q_total']:.6g}  "
            f"gauss_mean_rel={coulomb['gauss_mean_rel_residual']:.4g}  "
            f"prefer_phi={coulomb['multipole_phi_prefer']}  "
            f"prefer_E={coulomb['multipole_E_prefer']}  "
            f"all_pass={coulomb['all_pass']}  "
            f"t={coulomb['elapsed_s']:.2f}s"
        )
        write_tsv(
            os.path.join(OUT, "r1_coulomb_radial.tsv"),
            coulomb["radial"],
            ["r", "phi", "phi_analytic", "E_r", "E_r_analytic", "n_phi", "n_er"],
        )
        write_tsv(
            os.path.join(OUT, "r1_gauss_shells.tsv"),
            coulomb["gauss_shells"],
            [
                "r",
                "E_r_mean",
                "flux_4pi_r2_Er",
                "Q_encl",
                "Q_encl_over_eps0",
                "rel_residual",
                "n_shell",
            ],
        )

    wave = None
    wave_off = None
    if not args.skip_wave:
        print(f"[NE] Wave 1D nx={args.wave_nx} steps={args.wave_steps} (ε=μ=1) ...")
        wave = run_wave(nx=args.wave_nx, n_steps=args.wave_steps, eps=1.0, mu=1.0)
        print(
            f"  v_meas={wave['v_meas']:.6g}  v/c_th={wave['v_ratio']:.6g}  "
            f"all_pass={wave['all_pass']}  t={wave['elapsed_s']:.2f}s"
        )
        # TE KG4 off-unit: ε=4, μ=1 → c_th=1/2
        print(f"[NE] Wave 1D off-unit ε=4 μ=1 ...")
        wave_off = run_wave(
            nx=args.wave_nx, n_steps=args.wave_steps, eps=4.0, mu=1.0
        )
        print(
            f"  c_th={wave_off['c_th']:.6g} v_meas={wave_off['v_meas']:.6g}  "
            f"v/c_th={wave_off['v_ratio']:.6g}  all_pass={wave_off['all_pass']}"
        )
        write_tsv(
            os.path.join(OUT, "r1_wave_track.tsv"),
            wave["track"],
            ["step", "t", "peak_x", "peak_a", "energy"],
        )
        write_tsv(
            os.path.join(OUT, "r1_wave_track_offunit.tsv"),
            wave_off["track"],
            ["step", "t", "peak_x", "peak_a", "energy"],
        )

    # strip heavy arrays for json (radial/track already in tsv)
    def slim_coulomb(d: Dict) -> Dict:
        out = {k: v for k, v in d.items() if k not in ("radial", "gauss_shells")}
        # keep summary of gauss
        out["gauss_shells_summary"] = [
            {
                "r": g["r"],
                "rel_residual": g["rel_residual"],
                "Q_encl": g["Q_encl"],
            }
            for g in d.get("gauss_shells", [])
        ]
        return out

    def slim_wave(d: Dict) -> Dict:
        out = {k: v for k, v in d.items() if k != "track"}
        out["track_n"] = len(d.get("track", []))
        return out

    result = {
        "round": 1,
        "agent": "NE",
        "date": "2026-07-18",
        "sandbox": "sandbox_ne_r1_em.py",
        "design": "design_r1_em_sandbox.md",
        "shared_c": {
            "C_LOCAL": C_LOCAL,
            "eps0": EPS0,
            "mu0": MU0,
            "c_from_eps_mu": c_const,
            "language": "c = free-field locality = 1/sqrt(eps0 mu0) free-gauge constitutive",
            "path_cost_sibling": "v76 F1-3D psi channel; same c number, different constitutive law",
        },
        "full_maxwell_claim": False,
        "provisional": True,
        "note": (
            "Maxwell-lite free-gauge demos only. "
            "Do not claim full Maxwell until TE equations match."
        ),
        "te_ref": "v77/work/TE/for_ne_kill_gates_r1.md",
        "te_tags": TE_TAGS,
        "coulomb": slim_coulomb(coulomb) if coulomb else None,
        "wave": slim_wave(wave) if wave else None,
        "wave_offunit": slim_wave(wave_off) if wave_off else None,
        "KG": {
            "KG1": coulomb["gates"]["KG1_vacuum"]["pass"] if coulomb else None,
            "KG2": coulomb["gates"]["KG2_gauss"]["pass"] if coulomb else None,
            "KG3": coulomb["gates"]["KG3_coulomb_1r2"]["pass"] if coulomb else None,
            "KG4": bool(
                (wave["all_pass"] if wave else False)
                and (wave_off["all_pass"] if wave_off else False)
            ),
        },
        "demos": {
            "D-EM-gauss-coulomb": {
                "status": (
                    "LIVE_PASS"
                    if coulomb and coulomb["all_pass"]
                    else ("LIVE_FAIL" if coulomb else "SKIP")
                ),
                "all_pass": coulomb["all_pass"] if coulomb else None,
                "gates": ["KG1", "KG2", "KG3"],
            },
            "D-EM-vacuum-control": {
                "status": (
                    "LIVE_PASS"
                    if coulomb and coulomb["gates"]["KG1_vacuum"]["pass"]
                    else ("LIVE_FAIL" if coulomb else "SKIP")
                ),
            },
            "D-EM-wave-c": {
                "status": (
                    "LIVE_PASS"
                    if wave
                    and wave["all_pass"]
                    and wave_off
                    and wave_off["all_pass"]
                    else ("LIVE_FAIL" if wave else "SKIP")
                ),
                "all_pass": (
                    wave["all_pass"] and wave_off["all_pass"]
                    if wave and wave_off
                    else None
                ),
                "gates": ["KG4"],
            },
            "D-EM-coulomb": {
                "status": "alias_of_D-EM-gauss-coulomb",
                "note": "NE R1 alias; TE registry name is D-EM-gauss-coulomb",
            },
        },
    }

    out_json = os.path.join(OUT, "r1_result.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    # summary text
    lines = [
        "v77 NE Round 1 — monist free-gauge EM sandbox",
        f"shared c = {C_LOCAL} = 1/sqrt(eps0 mu0) with eps0={EPS0} mu0={MU0}",
        f"full_maxwell_claim = False (provisional until TE match)",
        "",
    ]
    if coulomb:
        lines.append(
            f"D-EM-coulomb: all_pass={coulomb['all_pass']}  "
            f"gauss_mean_rel={coulomb['gauss_mean_rel_residual']:.4g}  "
            f"phi_prefer={coulomb['multipole_phi_prefer']}  "
            f"E_prefer={coulomb['multipole_E_prefer']}  "
            f"Q={coulomb['Q_total']:.6g}"
        )
        for gname, g in coulomb["gates"].items():
            lines.append(f"  {gname}: pass={g['pass']}")
    if wave:
        lines.append(
            f"D-EM-wave-c: all_pass={wave['all_pass']}  "
            f"v_meas={wave['v_meas']:.6g}  v/c={wave['v_ratio']:.6g}"
        )
        for gname, g in wave["gates"].items():
            lines.append(f"  {gname}: pass={g['pass']}")
    lines.append("")
    lines.append(f"exports: {OUT}/r1_result.json + tsv maps")
    summary = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "r1_summary.txt"), "w", encoding="utf-8") as f:
        f.write(summary)
    print(summary)
    print(f"[NE] wrote {out_json}")


if __name__ == "__main__":
    main()
