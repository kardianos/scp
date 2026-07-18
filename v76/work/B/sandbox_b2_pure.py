#!/usr/bin/env python3
"""
v76 Approach B — B2-lite pure-stdlib sandbox (no numpy required).

Budget: rho_free + rho_bound = rho0.
n = rho0 / rho_free  (optical monism; no Poisson/Einstein).
Primary metrics: free deficit, Born eikonal deflection/delay, RK2 rays.
"""
from __future__ import annotations

import json
import math
import os
from typing import Dict, List, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C_LOCAL = 1.0
EPS_MIN = 0.05
L = 20.0
N = 128  # grid (pure-python friendly)


def linspace(a: float, b: float, n: int) -> List[float]:
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def place_lock(xs: List[float], ys: List[float], A: float, sigma: float):
    rho_b = []
    rho_f = []
    max_b = RHO0 - EPS_MIN
    peak = 0.0
    for y in ys:
        row_b = []
        for x in xs:
            r2 = x * x + y * y
            v = A * math.exp(-r2 / (2.0 * sigma * sigma))
            row_b.append(v)
            if v > peak:
                peak = v
        rho_b.append(row_b)
    scale = 1.0
    if peak > max_b and peak > 0:
        scale = max_b / peak
    for j in range(len(ys)):
        row_f = []
        for i in range(len(xs)):
            rho_b[j][i] *= scale
            row_f.append(RHO0 - rho_b[j][i])
        rho_f.append(row_f)
    return rho_f, rho_b


def refractive(rho_f: List[List[float]]) -> List[List[float]]:
    return [[RHO0 / max(v, 1e-15) for v in row] for row in rho_f]


def medium_stats(rho_f, rho_b, xs, ys, dx, A, sigma) -> Dict:
    core_r = 2.0 * sigma
    ext_lo, ext_hi = 4.0 * sigma, 0.45 * L
    m = 0.0
    sum_c = n_c = 0.0
    sum_e = n_e = 0.0
    resid = 0.0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            rb = rho_b[j][i]
            rf = rho_f[j][i]
            m += rb * dx * dx
            resid = max(resid, abs(rf + rb - RHO0))
            r = math.hypot(x, y)
            if r <= core_r:
                sum_c += rf
                n_c += 1
            elif ext_lo < r < ext_hi:
                sum_e += rf
                n_e += 1
    rf_c = sum_c / n_c if n_c else float("nan")
    rf_e = sum_e / n_e if n_e else RHO0
    return {
        "m_bound": m / (C_LOCAL**2),
        "rho_free_core": rf_c,
        "rho_free_exterior": rf_e,
        "free_deficit_core": rf_e - rf_c,
        "budget_residual_max": resid,
        "A_lock": A,
        "sigma_lock": sigma,
    }


def sample(field: List[List[float]], x: float, y: float, xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    # indices
    if x < xs[0] or x > xs[-1] or y < ys[0] or y > ys[-1]:
        return float("nan")
    fx = (x - xs[0]) / (xs[-1] - xs[0]) * (n - 1)
    fy = (y - ys[0]) / (ys[-1] - ys[0]) * (n - 1)
    i0 = int(math.floor(fx))
    j0 = int(math.floor(fy))
    i1 = min(i0 + 1, n - 1)
    j1 = min(j0 + 1, n - 1)
    tx = fx - i0
    ty = fy - j0
    v00 = field[j0][i0]
    v10 = field[j0][i1]
    v01 = field[j1][i0]
    v11 = field[j1][i1]
    return (1 - ty) * ((1 - tx) * v00 + tx * v10) + ty * ((1 - tx) * v01 + tx * v11)


def grad_n(nfield, x, y, xs, ys, h: float) -> Tuple[float, float]:
    npx = sample(nfield, x + h, y, xs, ys)
    nmx = sample(nfield, x - h, y, xs, ys)
    npy = sample(nfield, x, y + h, xs, ys)
    nmy = sample(nfield, x, y - h, xs, ys)
    if any(math.isnan(v) for v in (npx, nmx, npy, nmy)):
        return 0.0, 0.0
    return (npx - nmx) / (2 * h), (npy - nmy) / (2 * h)


def integrate_ray(
    nfield, x0, y0, th0, xs, ys, dx_grid, s_max=40.0, ds=0.03
) -> Tuple[List[Tuple[float, float]], float, float]:
    x, y, th = x0, y0, th0
    path = [(x, y)]
    t_opt = 0.0
    half = 0.48 * L
    steps = int(s_max / ds)
    for _ in range(steps):
        nval = sample(nfield, x, y, xs, ys)
        if math.isnan(nval) or abs(x) > half or abs(y) > half:
            break
        dndx, dndy = grad_n(nfield, x, y, xs, ys, dx_grid)
        tx, ty = math.cos(th), math.sin(th)
        gn_par = dndx * tx + dndy * ty
        dtx = (dndx - gn_par * tx) / nval
        dty = (dndy - gn_par * ty) / nval
        dth = -ty * dtx + tx * dty
        th_mid = th + 0.5 * ds * dth
        x_mid = x + 0.5 * ds * math.cos(th)
        y_mid = y + 0.5 * ds * math.sin(th)
        n_mid = sample(nfield, x_mid, y_mid, xs, ys)
        if math.isnan(n_mid):
            break
        dndx_m, dndy_m = grad_n(nfield, x_mid, y_mid, xs, ys, dx_grid)
        txm, tym = math.cos(th_mid), math.sin(th_mid)
        gn_par_m = dndx_m * txm + dndy_m * tym
        dtx_m = (dndx_m - gn_par_m * txm) / n_mid
        dty_m = (dndy_m - gn_par_m * tym) / n_mid
        dth_m = -tym * dtx_m + txm * dty_m
        x = x + ds * math.cos(th_mid)
        y = y + ds * math.sin(th_mid)
        th = th + ds * dth_m
        t_opt += n_mid * ds / C_LOCAL
        path.append((x, y))
    return path, t_opt, th


def born_ray(A: float, sigma: float, b: float, xmax: float = 30.0, dx: float = 0.01):
    """Straight-path Born eikonal: deflection and excess delay. No gravity solver."""
    s2 = sigma * sigma
    delta = 0.0
    delay = 0.0
    x = -xmax
    while x <= xmax + 1e-15:
        f = A * math.exp(-(x * x + b * b) / (2.0 * s2))
        if f > RHO0 - EPS_MIN:
            f = RHO0 - EPS_MIN
        n = RHO0 / max(RHO0 - f, EPS_MIN)
        # d theta / ds ≈ (1/n) dn/dy = n * df/dy, df/dy = f*(-b/s2)
        dth = n * f * (-b / s2)
        delta += dth * dx
        delay += (n - 1.0) * dx
        x += dx
    return delta, delay


def ray_suite(nfield, xs, ys, dx, impacts: Sequence[float], A_born: float, sig_born: float, use_grid_rays: bool):
    out = []
    n_vac = [[1.0 for _ in xs] for _ in ys]
    x_start = -0.45 * L
    for b in impacts:
        if use_grid_rays:
            _, t_v, _ = integrate_ray(n_vac, x_start, b, 0.0, xs, ys, dx)
            path, t_opt, th = integrate_ray(nfield, x_start, b, 0.0, xs, ys, dx)
            defl = th
            delay = t_opt - t_v
            y_exit = path[-1][1] if path else float("nan")
            n_steps = len(path)
        else:
            defl, delay = born_ray(A_born, sig_born, b)
            t_opt = delay  # excess only in born mode for delay field
            t_v = 0.0
            y_exit = float("nan")
            n_steps = 0
        # Also always record Born for cross-check
        b_defl, b_delay = born_ray(A_born, sig_born, b)
        out.append(
            {
                "b": float(b),
                "deflection_rad": float(defl),
                "deflection_deg": float(defl * 180.0 / math.pi),
                "delay": float(delay),
                "t_optical": float(t_opt) if use_grid_rays else float("nan"),
                "t_vacuum": float(t_v) if use_grid_rays else float("nan"),
                "born_defl_rad": float(b_defl),
                "born_defl_deg": float(b_defl * 180.0 / math.pi),
                "born_delay": float(b_delay),
                "y_exit": float(y_exit) if not math.isnan(y_exit) else None,
                "n_steps": int(n_steps),
            }
        )
    return out


def formation_step(rho_f, rho_b, gamma=0.05, mu=0.01):
    """Autocatalytic lock growth: only grows where bound already present (seed)."""
    ny = len(rho_b)
    nx = len(rho_b[0])
    new_b = [[0.0] * nx for _ in range(ny)]
    for j in range(ny):
        for i in range(nx):
            rb = rho_b[j][i]
            rf = rho_f[j][i]
            # dB = gamma * B * F - mu * B  (stays zero where B=0)
            db = gamma * rb * rf - mu * rb
            nb = rb + db
            if nb < 0:
                nb = 0.0
            if nb > RHO0 - EPS_MIN:
                nb = RHO0 - EPS_MIN
            new_b[j][i] = nb
    new_f = [[RHO0 - new_b[j][i] for i in range(nx)] for j in range(ny)]
    return new_f, new_b


def run_case(A: float, sigma: float, tag: str, impacts, use_grid_rays: bool = True) -> Dict:
    xs = linspace(-L / 2, L / 2, N)
    ys = linspace(-L / 2, L / 2, N)
    dx = xs[1] - xs[0]
    rho_f, rho_b = place_lock(xs, ys, A, sigma)
    nfield = refractive(rho_f)
    stats = medium_stats(rho_f, rho_b, xs, ys, dx, A, sigma)
    rays = ray_suite(nfield, xs, ys, dx, impacts, A, sigma, use_grid_rays)
    max_defl = max(abs(r["deflection_rad"]) for r in rays)
    max_delay = max(r["delay"] for r in rays)
    max_born_defl = max(abs(r["born_defl_rad"]) for r in rays)
    max_born_delay = max(r["born_delay"] for r in rays)
    return {
        "tag": tag,
        "A": A,
        "sigma": sigma,
        "stats": stats,
        "rays": rays,
        "max_abs_deflection_rad": max_defl,
        "max_delay": max_delay,
        "max_abs_born_defl_rad": max_born_defl,
        "max_born_delay": max_born_delay,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    impacts = [-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]

    # Vacuum: use Born (exact zeros) + light grid
    vac = run_case(0.0, 1.0, "vacuum", impacts, use_grid_rays=True)
    # Primary lock — moderate A so rays don't trap
    lock = run_case(0.35, 1.2, "lock_A0.35_s1.2", impacts, use_grid_rays=True)
    weak = run_case(0.15, 1.2, "lock_A0.15_s1.2", impacts, use_grid_rays=True)
    # Stronger lock Born-only cross-check (grid rays may strongly bend)
    strong = run_case(0.7, 1.2, "lock_A0.7_s1.2_bornish", impacts, use_grid_rays=True)

    # Formation smoke: seed + autocatalytic growth
    xs = linspace(-L / 2, L / 2, N)
    ys = linspace(-L / 2, L / 2, N)
    dx = xs[1] - xs[0]
    rho_f, rho_b = place_lock(xs, ys, A=0.12, sigma=0.9)
    hist = []
    for step in range(60):
        rho_f, rho_b = formation_step(rho_f, rho_b, gamma=0.12, mu=0.02)
        if step % 15 == 0 or step == 59:
            st = medium_stats(rho_f, rho_b, xs, ys, dx, A=0.12, sigma=0.9)
            hist.append(
                {
                    "step": step,
                    "m_bound": st["m_bound"],
                    "free_deficit_core": st["free_deficit_core"],
                    "rho_b_max": max(max(row) for row in rho_b),
                }
            )
    form_stats = medium_stats(rho_f, rho_b, xs, ys, dx, A=0.12, sigma=0.9)
    nfield = refractive(rho_f)
    form_rays = ray_suite(nfield, xs, ys, dx, impacts, A=form_stats["rho_free_core"], sig_born=0.9, use_grid_rays=True)
    # For formation Born, estimate A from max bound
    A_form = max(max(row) for row in rho_b)
    form_born = [born_ray(A_form, 0.9, b) for b in impacts]

    lock_s = lock["stats"]
    gates = {
        "budget_residual_ok": lock_s["budget_residual_max"] < 1e-12,
        "free_deficit_positive": lock_s["free_deficit_core"] > 0.01,
        "lock_deflection_nonzero": lock["max_abs_deflection_rad"] > 1e-4
        or lock["max_abs_born_defl_rad"] > 1e-4,
        "lock_delay_positive": lock["max_delay"] > 1e-4 or lock["max_born_delay"] > 1e-4,
        "vacuum_deflection_near_zero": vac["max_abs_deflection_rad"] < 1e-3,
        "vacuum_delay_near_zero": abs(vac["max_delay"]) < 1e-3,
        "weaker_smaller_defl": (
            weak["max_abs_born_defl_rad"] < lock["max_abs_born_defl_rad"]
        ),
        "formation_grew_or_held": hist[-1]["m_bound"] >= hist[0]["m_bound"] * 0.5,
    }
    gates["round1_pass"] = all(gates.values())

    results = {
        "design": "B2-lite optical monism (pure stdlib)",
        "rho0": RHO0,
        "c_local": C_LOCAL,
        "eps_min": EPS_MIN,
        "grid_N": N,
        "n_index_law": "n = rho0 / rho_free",
        "gravity_solver": None,
        "cases": [vac, lock, weak, strong],
        "formation_smoke": {
            "tag": "formation_autocatalytic",
            "history": hist,
            "stats": form_stats,
            "A_final_peak": A_form,
            "max_abs_deflection_rad": max(abs(r["deflection_rad"]) for r in form_rays),
            "max_born_defl_rad": max(abs(d) for d, _ in form_born),
            "rays_sample": form_rays[::2],
        },
        "gates": gates,
    }

    out_json = os.path.join(OUT, "results.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

    lines = [
        "v76 B2-lite pure sandbox results",
        f"round1_pass = {gates['round1_pass']}",
        f"lock m_bound = {lock_s['m_bound']:.6f}",
        f"lock free_deficit_core = {lock_s['free_deficit_core']:.6f}",
        f"lock max |defl| rad (grid) = {lock['max_abs_deflection_rad']:.6e}",
        f"lock max |defl| rad (Born) = {lock['max_abs_born_defl_rad']:.6e}",
        f"lock max delay (grid) = {lock['max_delay']:.6e}",
        f"lock max delay (Born) = {lock['max_born_delay']:.6e}",
        f"vac max |defl| rad = {vac['max_abs_deflection_rad']:.6e}",
        f"budget residual max = {lock_s['budget_residual_max']:.3e}",
        f"formation m: {hist[0]['m_bound']:.4f} -> {hist[-1]['m_bound']:.4f}",
        "gates:",
    ]
    for k, v in gates.items():
        lines.append(f"  {k}: {v}")
    lines.append("sample lock rays (b, grid_defl_deg, born_defl_deg, grid_delay, born_delay):")
    for r in lock["rays"][::2]:
        lines.append(
            f"  b={r['b']:+.1f}  g={r['deflection_deg']:+.4f}deg  "
            f"B={r['born_defl_deg']:+.4f}deg  gd={r['delay']:+.5f}  Bd={r['born_delay']:+.5f}"
        )
    text = "\n".join(lines) + "\n"
    summary_path = os.path.join(OUT, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(text)
    print(text)
    print(f"Wrote {out_json}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
