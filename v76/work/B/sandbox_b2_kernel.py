#!/usr/bin/env python3
"""
v76 Approach B — B2-response kernel sandbox (Round 1 extension).

Addresses C no-go: local n(ρ_free)+pointwise budget ⇒ compact GRIN only,
not long-range path-cost monopole.

Here: compact ρ_bound lock + free-medium path-cost response
  Φ(x) = α ∫ ρ_bound(x') / (|x-x'| + ε) dA'
  n(x) = 1 + Φ(x)     (chart optical factor; local c still seed)
Rays from n only — no Poisson/Einstein *solver as second sector*; Φ is the
explicit free-response kernel (C Class preferred shape). Dualist residue:
kernel form is postulated (Round 2: derive from free dynamics).

Also retains budget ρ_free + ρ_bound = ρ0 for ledger/deficit bookkeeping;
rays use Φ[ρ_bound], not local n(ρ_free), for the long-range channel.
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
L = 24.0
N = 96  # coarser: O(N^4) naive convolution; use N=96 carefully
# Soft core for 1/R
SOFT = 0.15


def linspace(a: float, b: float, n: int) -> List[float]:
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def place_lock(xs, ys, A, sigma):
    rho_b = []
    peak = 0.0
    max_b = RHO0 - EPS_MIN
    for y in ys:
        row = []
        for x in xs:
            v = A * math.exp(-(x * x + y * y) / (2.0 * sigma * sigma))
            row.append(v)
            peak = max(peak, v)
        rho_b.append(row)
    if peak > max_b and peak > 0:
        s = max_b / peak
        rho_b = [[v * s for v in row] for row in rho_b]
    rho_f = [[RHO0 - rho_b[j][i] for i in range(len(xs))] for j in range(len(ys))]
    return rho_f, rho_b


def m_bound(rho_b, dx):
    return sum(sum(row) for row in rho_b) * dx * dx


def free_response_phi(rho_b, xs, ys, dx, alpha: float) -> List[List[float]]:
    """Φ = α ∫ ρ_b / (R+soft) dA  — free path-cost kernel (C preferred class)."""
    ny, nx = len(ys), len(xs)
    phi = [[0.0] * nx for _ in range(ny)]
    # downsample source for speed: use every p-th cell as multipole lumps
    p = 2
    src = []
    for j in range(0, ny, p):
        for i in range(0, nx, p):
            mass = 0.0
            for jj in range(j, min(j + p, ny)):
                for ii in range(i, min(i + p, nx)):
                    mass += rho_b[jj][ii]
            mass *= dx * dx
            if mass > 1e-12:
                # center of block
                cx = xs[min(i + p // 2, nx - 1)]
                cy = ys[min(j + p // 2, ny - 1)]
                src.append((cx, cy, mass))
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            s = 0.0
            for cx, cy, mass in src:
                R = math.hypot(x - cx, y - cy)
                s += mass / (R + SOFT)
            phi[j][i] = alpha * s
    return phi


def sample(field, x, y, xs, ys):
    n = len(xs)
    if x < xs[0] or x > xs[-1] or y < ys[0] or y > ys[-1]:
        return float("nan")
    fx = (x - xs[0]) / (xs[-1] - xs[0]) * (n - 1)
    fy = (y - ys[0]) / (ys[-1] - ys[0]) * (n - 1)
    i0 = int(math.floor(fx))
    j0 = int(math.floor(fy))
    i1 = min(i0 + 1, n - 1)
    j1 = min(j0 + 1, n - 1)
    tx, ty = fx - i0, fy - j0
    return (
        (1 - ty) * ((1 - tx) * field[j0][i0] + tx * field[j0][i1])
        + ty * ((1 - tx) * field[j1][i0] + tx * field[j1][i1])
    )


def grad(field, x, y, xs, ys, h):
    return (
        (sample(field, x + h, y, xs, ys) - sample(field, x - h, y, xs, ys)) / (2 * h),
        (sample(field, x, y + h, xs, ys) - sample(field, x, y - h, xs, ys)) / (2 * h),
    )


def integrate_ray(nfield, x0, y0, th0, xs, ys, h, s_max=50.0, ds=0.04):
    x, y, th = x0, y0, th0
    t_opt = 0.0
    half = 0.48 * L
    for _ in range(int(s_max / ds)):
        nval = sample(nfield, x, y, xs, ys)
        if math.isnan(nval) or abs(x) > half or abs(y) > half:
            break
        dndx, dndy = grad(nfield, x, y, xs, ys, h)
        if math.isnan(dndx):
            break
        tx, ty = math.cos(th), math.sin(th)
        gn = dndx * tx + dndy * ty
        dtx = (dndx - gn * tx) / max(nval, 1e-9)
        dty = (dndy - gn * ty) / max(nval, 1e-9)
        dth = -ty * dtx + tx * dty
        th_m = th + 0.5 * ds * dth
        xm = x + 0.5 * ds * math.cos(th)
        ym = y + 0.5 * ds * math.sin(th)
        nm = sample(nfield, xm, ym, xs, ys)
        if math.isnan(nm):
            break
        dndx_m, dndy_m = grad(nfield, xm, ym, xs, ys, h)
        txm, tym = math.cos(th_m), math.sin(th_m)
        gn_m = dndx_m * txm + dndy_m * tym
        dtx_m = (dndx_m - gn_m * txm) / max(nm, 1e-9)
        dty_m = (dndy_m - gn_m * tym) / max(nm, 1e-9)
        dth_m = -tym * dtx_m + txm * dty_m
        x += ds * math.cos(th_m)
        y += ds * math.sin(th_m)
        th += ds * dth_m
        t_opt += nm * ds
    return t_opt, th


def born_from_phi_along_line(phi_func, b, xmax=40.0, dx=0.05):
    """Straight-path Born on n=1+Φ: dθ/ds=(1/n)∂n/∂y, delay=∫(n-1)dx."""
    delta = delay = 0.0
    x = -xmax
    while x <= xmax:
        # finite-diff ∂Φ/∂y
        n0 = 1.0 + phi_func(x, b)
        npy = 1.0 + phi_func(x, b + dx)
        nmy = 1.0 + phi_func(x, b - dx)
        dndy = (npy - nmy) / (2 * dx)
        dth = dndy / max(n0, 1e-9)
        delta += dth * dx
        delay += (n0 - 1.0) * dx
        x += dx
    return delta, delay


def analytic_phi_point(x, y, M, alpha, soft=SOFT):
    """Point-mass kernel Φ = α M / (R+soft) — monopole exterior."""
    return alpha * M / (math.hypot(x, y) + soft)


def main():
    os.makedirs(OUT, exist_ok=True)
    xs = linspace(-L / 2, L / 2, N)
    ys = linspace(-L / 2, L / 2, N)
    dx = xs[1] - xs[0]
    A, sigma = 0.35, 1.2
    rho_f, rho_b = place_lock(xs, ys, A, sigma)
    M = m_bound(rho_b, dx)
    # Choose alpha so weak-field scale is visible
    alpha = 0.08
    phi = free_response_phi(rho_b, xs, ys, dx, alpha)
    nfield = [[1.0 + phi[j][i] for i in range(N)] for j in range(N)]

    # Exterior Φ vs 1/r check at sample radii
    exterior = []
    for r in [3.0, 4.0, 5.0, 6.0, 8.0]:
        ph = sample(phi, r, 0.0, xs, ys)
        # monopole expectation α M / r
        mono = alpha * M / r
        exterior.append(
            {
                "r": r,
                "Phi": ph,
                "alpha_M_over_r": mono,
                "ratio_Phi_to_monopole": ph / mono if mono else float("nan"),
            }
        )

    impacts = [-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    # Point-mass Born (clean monopole) for comparison
    rays_point = []
    for b in impacts:
        d, delay = born_from_phi_along_line(
            lambda x, y, M=M, a=alpha: analytic_phi_point(x, y, M, a), b
        )
        rays_point.append(
            {
                "b": b,
                "deflection_rad": d,
                "deflection_deg": d * 180 / math.pi,
                "delay": delay,
            }
        )

    # Grid rays
    rays_grid = []
    h = dx
    x_start = -0.45 * L
    for b in impacts:
        t_opt, th = integrate_ray(nfield, x_start, b, 0.0, xs, ys, h)
        # vacuum
        nvac = [[1.0] * N for _ in range(N)]
        t_v, _ = integrate_ray(nvac, x_start, b, 0.0, xs, ys, h)
        rays_grid.append(
            {
                "b": b,
                "deflection_rad": th,
                "deflection_deg": th * 180 / math.pi,
                "delay": t_opt - t_v,
            }
        )

    # Local-optics GRIN for contrast (same lock, n=rho0/rho_free)
    n_local = [[RHO0 / max(rho_f[j][i], EPS_MIN) for i in range(N)] for j in range(N)]
    # sample Phi_local-equivalent: n-1 at exterior r=6
    n_ext = sample(n_local, 6.0, 0.0, xs, ys)
    phi_ext = sample(phi, 6.0, 0.0, xs, ys)

    # free deficit
    core_r = 2.4
    sum_c = n_c = sum_e = n_e = 0.0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            r = math.hypot(x, y)
            if r <= core_r:
                sum_c += rho_f[j][i]
                n_c += 1
            elif 4.0 < r < 0.45 * L:
                sum_e += rho_f[j][i]
                n_e += 1
    deficit = (sum_e / n_e) - (sum_c / n_c)

    max_defl_pt = max(abs(r["deflection_rad"]) for r in rays_point)
    max_defl_g = max(abs(r["deflection_rad"]) for r in rays_grid)

    results = {
        "design": "B2-response free kernel Φ=α∫ρ_b/(R+ε)",
        "gravity_solver": None,
        "alpha": alpha,
        "M_bound": M,
        "A": A,
        "sigma": sigma,
        "free_deficit_core": deficit,
        "exterior_Phi_vs_monopole": exterior,
        "rays_point_mass_Born": rays_point,
        "rays_grid": rays_grid,
        "max_abs_defl_point_Born": max_defl_pt,
        "max_abs_defl_grid": max_defl_g,
        "contrast_at_r6": {
            "n_local_optics": n_ext,
            "n_local_minus_1": n_ext - 1.0,
            "Phi_kernel": phi_ext,
            "note": "local optics n-1 ~0 at r=6 (compact); kernel Phi remains O(alpha M/r)",
        },
        "gates": {
            "free_deficit_positive": deficit > 0.01,
            "kernel_ray_defl_nonzero": max_defl_pt > 1e-4 or max_defl_g > 1e-4,
            "exterior_Phi_near_monopole": all(
                0.5 < e["ratio_Phi_to_monopole"] < 1.5 for e in exterior if e["r"] >= 4
            ),
            "local_optics_exterior_small": (n_ext - 1.0) < 0.05,
        },
    }
    results["gates"]["kernel_round1_interesting"] = all(
        results["gates"][k]
        for k in (
            "free_deficit_positive",
            "kernel_ray_defl_nonzero",
            "exterior_Phi_near_monopole",
            "local_optics_exterior_small",
        )
    )

    path = os.path.join(OUT, "results_kernel.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)

    # Also write human summary fragment
    lines = [
        "B2-response kernel results",
        f"M_bound = {M:.6f}",
        f"free_deficit_core = {deficit:.6f}",
        f"max |defl| point-Born = {max_defl_pt:.6e}",
        f"max |defl| grid = {max_defl_g:.6e}",
        f"n_local(r=6)-1 = {n_ext-1:.6e}  (should be ~0)",
        f"Phi_kernel(r=6) = {phi_ext:.6e}  (should be ~ alpha M/r = {alpha*M/6:.6e})",
        "exterior Phi/monopole ratios:",
    ]
    for e in exterior:
        lines.append(
            f"  r={e['r']}: Phi={e['Phi']:.5f}  mono={e['alpha_M_over_r']:.5f}  "
            f"ratio={e['ratio_Phi_to_monopole']:.3f}"
        )
    lines.append("sample point-Born rays:")
    for r in rays_point[::2]:
        lines.append(
            f"  b={r['b']:+.1f}  defl={r['deflection_deg']:+.4f}deg  delay={r['delay']:+.5f}"
        )
    lines.append(f"gates: {results['gates']}")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "summary_kernel.txt"), "w") as f:
        f.write(text)
    print(text)
    print("Wrote", path)


if __name__ == "__main__":
    main()
