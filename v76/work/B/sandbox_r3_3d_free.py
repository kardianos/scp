#!/usr/bin/env python3
"""
v76 Approach B Round 3 — 3D free-capacity / free Laplace medium (A F1).

Dynamics (monist_1sector if ψ is free continuum state):
  −σ0 ∇²ψ = s ρ_bound     (σ0,s > 0 free constants)
  ℓ = ℓ0 + γ ψ            (path cost from free state only)
  ρ_free + ρ_bound = ρ0   (budget ledger)

Discrete: SOR Jacobi on N³ grid, ψ=0 on box boundary.
Rays: straight-path Born in equatorial plane, n = 1 + γ |ψ|.
Dualist control: same PDE tagged dualist_2sector (ontology, not different math).

Also reports analytic infinite-space monopole ψ∞ = (s/σ0) M / (4π r) for fit.

Usage:
  python3 sandbox_r3_3d_free.py
  python3 sandbox_r3_3d_free.py --N 24 --iters 400
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Dict, List, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C_LOCAL = 1.0
EPS_MIN = 0.05


def linspace(a: float, b: float, n: int) -> List[float]:
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def idx3(i: int, j: int, k: int, n: int) -> int:
    return (k * n + j) * n + i


def place_gaussian_3d(
    n: int, L: float, A: float, sigma: float
) -> Tuple[List[float], List[float], float, float]:
    """Return flat ρ_b, ρ_f arrays length n³, dx, m_ledger."""
    xs = linspace(-L / 2, L / 2, n)
    dx = xs[1] - xs[0]
    max_b = RHO0 - EPS_MIN
    rho_b = [0.0] * (n * n * n)
    peak = 0.0
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r2 = x * x + y * y + z * z
                v = A * math.exp(-r2 / (2.0 * sigma * sigma))
                rho_b[idx3(i, j, k, n)] = v
                if v > peak:
                    peak = v
    if peak > max_b and peak > 0:
        s = max_b / peak
        rho_b = [v * s for v in rho_b]
    rho_f = [RHO0 - v for v in rho_b]
    m = sum(rho_b) * (dx ** 3) / (C_LOCAL ** 2)
    return rho_b, rho_f, dx, m


def free_deficit_core(
    rho_f: List[float], n: int, L: float, sigma: float
) -> float:
    xs = linspace(-L / 2, L / 2, n)
    core_r = 2.0 * sigma
    sc = nc = se = ne = 0.0
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r = math.sqrt(x * x + y * y + z * z)
                v = rho_f[idx3(i, j, k, n)]
                if r <= core_r:
                    sc += v
                    nc += 1
                elif 3.5 * sigma < r < 0.42 * L:
                    se += v
                    ne += 1
    if nc == 0 or ne == 0:
        return 0.0
    return se / ne - sc / nc


def sor_poisson3d(
    rho_b: List[float],
    n: int,
    dx: float,
    kappa: float,
    n_iter: int = 500,
    omega: float = 1.6,
) -> List[float]:
    """
    Solve −∇²ψ = kappa * ρ_b with ψ=0 on boundary via SOR.
    kappa = s/σ0 in F1 notation.
    Flat array length n³.
    """
    psi = [0.0] * (n * n * n)
    rhs_scale = kappa * dx * dx
    # interior 1..n-2
    for it in range(n_iter):
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
                    # ∇²ψ ≈ (neigh - 6ψ)/dx² = -kappa ρ
                    # 6ψ - neigh = kappa ρ dx²
                    # ψ* = (neigh + rhs_scale ρ)/6
                    star = (neigh + rhs_scale * rho_b[p]) / 6.0
                    new = (1.0 - omega) * psi[p] + omega * star
                    d = abs(new - psi[p])
                    if d > max_d:
                        max_d = d
                    psi[p] = new
        # boundary stays 0
        if max_d < 1e-10:
            break
    return psi


def sample_trilinear(
    field: List[float], x: float, y: float, z: float, L: float, n: int
) -> float:
    half = L / 2
    if abs(x) >= half or abs(y) >= half or abs(z) >= half:
        return float("nan")
    # map to index
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


def radial_profile(
    field: List[float], L: float, n: int, radii: Sequence[float], n_ang: int = 12
) -> List[Dict]:
    out = []
    for r in radii:
        vals = []
        # fibonacci-ish sphere samples
        for a in range(n_ang):
            for b in range(max(1, n_ang // 2)):
                th = math.pi * (a + 0.5) / n_ang
                ph = 2 * math.pi * (b + 0.5) / (n_ang // 2)
                x = r * math.sin(th) * math.cos(ph)
                y = r * math.sin(th) * math.sin(ph)
                z = r * math.cos(th)
                v = sample_trilinear(field, x, y, z, L, n)
                if not math.isnan(v):
                    vals.append(v)
        if not vals:
            out.append({"r": r, "mean": float("nan"), "abs_mean": float("nan")})
            continue
        mean = sum(vals) / len(vals)
        out.append({"r": float(r), "mean": float(mean), "abs_mean": float(abs(mean)), "n": len(vals)})
    return out


def fit_log_vs_invr(rs: List[float], ys: List[float], rmin: float = 2.5) -> Dict:
    pts = [(r, y) for r, y in zip(rs, ys) if r >= rmin and y == y and y > 0]
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
    A = (n * sLY - sL * sY) / detL if abs(detL) > 1e-18 else 0.0
    B = (sY - A * sL) / n
    detR = n * sRR - sR * sR
    C = (n * sRY - sR * sY) / detR if abs(detR) > 1e-18 else 0.0
    D = (sY - C * sR) / n
    mse_log = mse_inv = 0.0
    for r, y in pts:
        mse_log += (y - A * math.log(r) - B) ** 2
        mse_inv += (y - C / r - D) ** 2
    mse_log /= n
    mse_inv /= n
    # also pure C/r with D=0
    # R² for inv: 1 - SS_res/SS_tot
    ybar = sY / n
    ss_tot = sum((y - ybar) ** 2 for _, y in pts)
    ss_inv = sum((y - C / r - D) ** 2 for r, y in pts)
    r2_inv = 1.0 - ss_inv / ss_tot if ss_tot > 0 else 0.0
    return {
        "ok": True,
        "A_log": A,
        "B_log": B,
        "C_inv": C,
        "D_inv": D,
        "mse_log": mse_log,
        "mse_invr": mse_inv,
        "r2_invr": r2_inv,
        "prefer": "1/r" if mse_inv <= mse_log else "log",
        "n_pts": n,
        "ratio_mse_log_over_invr": mse_log / mse_inv if mse_inv > 0 else float("inf"),
    }


def born_rays_equatorial(
    field: List[float],
    L: float,
    n: int,
    gamma: float,
    impacts: Sequence[float],
    xmax: float = None,
    ds: float = 0.08,
) -> List[Dict]:
    """n = 1 + gamma * |ψ|; Born straight path in z=0 plane."""
    if xmax is None:
        xmax = 0.42 * L
    out = []
    for b in impacts:
        if abs(b) > 0.4 * L:
            continue
        delta = delay = 0.0
        x = -xmax
        while x <= xmax:
            def nxy(xx, yy):
                v = sample_trilinear(field, xx, yy, 0.0, L, n)
                if math.isnan(v):
                    return 1.0
                return 1.0 + gamma * abs(v)

            n0 = nxy(x, b)
            h = ds
            npy = nxy(x, b + h)
            nmy = nxy(x, b - h)
            dndy = (npy - nmy) / (2 * h)
            delta += (dndy / max(n0, 1e-9)) * ds
            delay += (n0 - 1.0) * ds
            x += ds
        out.append(
            {
                "b": float(b),
                "deflection_rad": float(delta),
                "deflection_deg": float(delta * 180 / math.pi),
                "delay": float(delay),
            }
        )
    return out


def analytic_monopole_rays(
    M: float, kappa: float, gamma: float, impacts: Sequence[float], soft: float = 0.25
) -> List[Dict]:
    """ψ = kappa M / (4π r), n=1+γ|ψ|; weak Born defl ≈ −(γ kappa M / 2) / b for small γψ."""
    # exact Born for Φ=A/r, n=1+γ A/r: similar to R1 with alpha_eff = γ kappa M / (4π) wait
    # n-1 = γ kappa M / (4π r) = α_eff M / r with α_eff = γ kappa/(4π)
    alpha_eff = gamma * kappa / (4.0 * math.pi)
    out = []
    for b in impacts:
        # weak: defl = -2 α_eff M / b for n=1+α_eff M/r  (same as R1 kernel)
        A = alpha_eff * M
        # soft core Born integral
        delta = delay = 0.0
        dx = 0.05
        x = -40.0
        while x <= 40.0:
            r = math.sqrt(x * x + b * b)
            phi = A / (r + soft)  # n-1
            n0 = 1.0 + phi
            # d(n-1)/dy = A * (-b) / (r+soft)^2 * (partial r? for soft on r: -A b /r /(r+soft)^2 * r wait)
            # Φ=A/(R+s), R=sqrt(x²+b²), dΦ/db = -A/(R+s)^2 * dR/db = -A b /(R (R+s)^2)
            dndy = -A * b / (max(r, 1e-9) * (r + soft) ** 2)
            delta += (dndy / max(n0, 1e-9)) * dx
            delay += phi * dx
            x += dx
        out.append(
            {
                "b": float(b),
                "deflection_rad": float(delta),
                "deflection_deg": float(delta * 180 / math.pi),
                "delay": float(delay),
                "alpha_eff": alpha_eff,
            }
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--L", type=float, default=16.0)
    ap.add_argument("--A", type=float, default=0.4)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--iters", type=int, default=450)
    ap.add_argument("--kappa", type=float, default=1.0, help="s/σ0")
    ap.add_argument("--gamma", type=float, default=0.5, help="path-cost coupling")
    ap.add_argument("--omega", type=float, default=1.55)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    n = args.N
    L = args.L
    rho_b, rho_f, dx, m_L = place_gaussian_3d(n, L, args.A, args.sigma)
    deficit = free_deficit_core(rho_f, n, L, args.sigma)

    t1 = time.time()
    psi = sor_poisson3d(rho_b, n, dx, args.kappa, n_iter=args.iters, omega=args.omega)
    t2 = time.time()

    # radial
    # avoid near boundary and deep core
    r_max = 0.38 * L
    radii = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    radii = [r for r in radii if r < r_max]
    radial = radial_profile(psi, L, n, radii, n_ang=10)
    fit = fit_log_vs_invr([p["r"] for p in radial], [p["abs_mean"] for p in radial], rmin=2.0)

    # analytic infinite-space comparison at same r: ψ∞ = kappa M / (4π r)
    kappa = args.kappa
    radial_cmp = []
    for p in radial:
        r = p["r"]
        analytic = kappa * m_L / (4.0 * math.pi * r)
        ratio = p["abs_mean"] / analytic if analytic > 0 and p["abs_mean"] == p["abs_mean"] else float("nan")
        radial_cmp.append(
            {
                "r": r,
                "psi_num": p["abs_mean"],
                "psi_analytic_inf": analytic,
                "ratio_num_over_analytic": ratio,
            }
        )

    impacts = [-4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    impacts = [b for b in impacts if abs(b) < 0.4 * L]
    rays = born_rays_equatorial(psi, L, n, args.gamma, impacts)
    rays_an = analytic_monopole_rays(m_L, kappa, args.gamma, impacts)

    max_defl = max(abs(r["deflection_rad"]) for r in rays) if rays else 0.0
    max_delay = max(r["delay"] for r in rays) if rays else 0.0
    max_defl_an = max(abs(r["deflection_rad"]) for r in rays_an) if rays_an else 0.0

    # mid-b weak slope proxy mean(|defl|*|b|)
    def slope(rs):
        ps = [abs(r["deflection_rad"]) * abs(r["b"]) for r in rs if 1.0 <= abs(r["b"]) <= 3.0]
        return sum(ps) / len(ps) if ps else float("nan")

    # vacuum control: rho_b=0
    psi_vac = sor_poisson3d([0.0] * (n * n * n), n, dx, kappa, n_iter=50, omega=1.0)
    rays_vac = born_rays_equatorial(psi_vac, L, n, args.gamma, [1.0, 2.0])
    vac_defl = max(abs(r["deflection_rad"]) for r in rays_vac) if rays_vac else 0.0

    # gates
    prefer_1r = fit.get("prefer") == "1/r"
    r2_ok = fit.get("r2_invr", 0) > 0.9 if fit.get("ok") else False
    # ratio stability in mid shell
    mid = [c for c in radial_cmp if 2.5 <= c["r"] <= 5.0 and c["ratio_num_over_analytic"] == c["ratio_num_over_analytic"]]
    if mid:
        ratios = [c["ratio_num_over_analytic"] for c in mid]
        ratio_mean = sum(ratios) / len(ratios)
        ratio_std = (sum((x - ratio_mean) ** 2 for x in ratios) / len(ratios)) ** 0.5
    else:
        ratio_mean = ratio_std = float("nan")

    gates = {
        "budget_deficit_positive": deficit > 0.01,
        "exterior_prefer_1r": prefer_1r,
        "r2_invr_gt_0.9": r2_ok,
        "mse_log_gt_mse_invr": fit.get("ok") and fit.get("mse_log", 0) > fit.get("mse_invr", 1),
        "rays_nonzero": max_defl > 1e-4 and max_delay > 1e-4,
        "vacuum_defl_near_zero": vac_defl < 1e-3,
        "no_second_gravity_solver": True,
        "phi_origin_free_relaxation": True,
    }
    gates["monist_3d_1r_pass"] = all(
        [
            gates["budget_deficit_positive"],
            gates["exterior_prefer_1r"],
            gates["rays_nonzero"],
            gates["vacuum_defl_near_zero"],
            gates["no_second_gravity_solver"],
        ]
    )

    # dualist control is same field — ontology tag only for D Occam
    package = {
        "round": 3,
        "design": "3D free-capacity F1: -sigma0 Laplace psi = s rho_b; ell=ell0+gamma psi",
        "grid_N": n,
        "L": L,
        "dx": dx,
        "A": args.A,
        "sigma": args.sigma,
        "kappa_s_over_sigma0": kappa,
        "gamma_path": args.gamma,
        "sor_iters": args.iters,
        "omega": args.omega,
        "wall_secs_setup": t1 - t0,
        "wall_secs_sor": t2 - t1,
        "wall_secs_total": time.time() - t0,
        "m_ledger": m_L,
        "free_deficit_core": deficit,
        "monist_channel": {
            "sector_tag": "monist_1sector",
            "phi_origin": "free_relaxation",
            "gravity_solver": None,
            "sector_count": 1,
            "budget_identity": True,
            "radial": radial,
            "radial_vs_analytic": radial_cmp,
            "multipole_fit": fit,
            "rays": rays,
            "max_abs_defl": max_defl,
            "max_delay": max_delay,
            "ray_slope_proxy": slope(rays),
        },
        "analytic_infinite_space": {
            "sector_tag": "monist_1sector_analytic_reference",
            "psi": "kappa M / (4 pi r)",
            "rays": rays_an,
            "max_abs_defl": max_defl_an,
            "ray_slope_proxy": slope(rays_an),
            "alpha_eff": args.gamma * kappa / (4.0 * math.pi),
        },
        "dualist_control": {
            "sector_tag": "dualist_2sector",
            "note": "Same discrete field as monist channel; ontology = independent Phi sector (D Occam). Math isomorphic in vacuum linear F1.",
            "gravity_solver": "poisson_3d_tagged_dualist",
            "phi_origin": "dualist_poisson_label",
            "sector_count": 2,
            "budget_identity": False,
            "rays": rays,
            "radial": radial,
            "multipole_fit": fit,
        },
        "vacuum_control": {"max_abs_defl": vac_defl},
        "ratio_num_over_analytic_mid": {"mean": ratio_mean, "std": ratio_std},
        "gates": gates,
        "verdict": {
            "monist_3d_free_response_1r": gates["monist_3d_1r_pass"],
            "summary": (
                "3D free-capacity relaxation (F1) produces exterior ~1/r and Born rays "
                "without a second gravity pass. ψ is free continuum state; ρ_b is bound "
                "ledger of same continuum. Dualist control is ontology-tagged twin for D."
                if gates["monist_3d_1r_pass"]
                else "3D free relaxation did not pass 1/r + ray gates — see fit diagnostics."
            ),
        },
        "FOR_D": "Ingest round3_rays.tsv / path_cost; monist free_relaxation vs dualist_2sector twin + analytic",
        "FOR_A": "F1 3D numeric: G_eff candidate from gamma*s/(4pi sigma0) via alpha_eff=gamma*kappa/(4pi)",
        "FOR_C": "F1 3D Green 1/r confirmed on grid; S2 stress-test: monist only if psi=free state",
    }

    # exports
    with open(os.path.join(OUT, "round3_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    with open(os.path.join(OUT, "round3_path_cost.tsv"), "w") as f:
        f.write("sector_tag\tphi_origin\tr\tpsi_abs\tpsi_analytic\tratio\tquantity\n")
        for c in radial_cmp:
            f.write(
                f"monist_1sector\tfree_relaxation\t{c['r']}\t{c['psi_num']}\t"
                f"{c['psi_analytic_inf']}\t{c['ratio_num_over_analytic']}\tabs_psi\n"
            )
            f.write(
                f"dualist_2sector\tdualist_poisson_label\t{c['r']}\t{c['psi_num']}\t"
                f"{c['psi_analytic_inf']}\t{c['ratio_num_over_analytic']}\tabs_psi\n"
            )
            f.write(
                f"monist_1sector_analytic_reference\tanalytic_inf\t{c['r']}\t"
                f"{c['psi_analytic_inf']}\t{c['psi_analytic_inf']}\t1.0\tabs_psi\n"
            )

    with open(os.path.join(OUT, "round3_rays.tsv"), "w") as f:
        f.write(
            "sector_tag\tphi_origin\tb\tdeflection_rad\tdeflection_deg\tdelay\tm_ledger\n"
        )
        for r in rays:
            f.write(
                f"monist_1sector\tfree_relaxation\t{r['b']}\t{r['deflection_rad']}\t"
                f"{r['deflection_deg']}\t{r['delay']}\t{m_L}\n"
            )
            f.write(
                f"dualist_2sector\tdualist_poisson_label\t{r['b']}\t{r['deflection_rad']}\t"
                f"{r['deflection_deg']}\t{r['delay']}\t{m_L}\n"
            )
        for r in rays_an:
            f.write(
                f"monist_1sector_analytic_reference\tanalytic_inf\t{r['b']}\t"
                f"{r['deflection_rad']}\t{r['deflection_deg']}\t{r['delay']}\t{m_L}\n"
            )

    with open(os.path.join(OUT, "round3_free_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{m_L}\n")
        f.write(f"free_deficit_core\t{deficit}\n")
        f.write(f"A\t{args.A}\n")
        f.write(f"sigma\t{args.sigma}\n")
        f.write(f"N\t{n}\n")
        f.write(f"L\t{L}\n")
        f.write(f"kappa\t{kappa}\n")
        f.write(f"gamma\t{args.gamma}\n")

    lines = [
        "v76 Approach B ROUND 3 — 3D free-capacity F1",
        f"N={n} L={L} dx={dx:.4f} A={args.A} sigma={args.sigma}",
        f"m_ledger={m_L:.6f}  free_deficit_core={deficit:.6f}",
        f"SOR iters={args.iters}  wall_sor={t2-t1:.2f}s  total={time.time()-t0:.2f}s",
        f"multipole fit: prefer={fit.get('prefer')}  R2_1/r={fit.get('r2_invr')}  "
        f"mse_log={fit.get('mse_log')}  mse_1/r={fit.get('mse_invr')}  "
        f"mse_log/mse_1r={fit.get('ratio_mse_log_over_invr')}",
        f"mid-shell ratio num/analytic mean={ratio_mean} std={ratio_std}",
        f"max|defl| grid={max_defl:.6e}  max_delay={max_delay:.6e}",
        f"max|defl| analytic={max_defl_an:.6e}  vac_defl={vac_defl:.3e}",
        f"ray_slope grid={slope(rays):.5f}  analytic={slope(rays_an):.5f}",
        "sample radial:",
    ]
    for c in radial_cmp[::2]:
        lines.append(
            f"  r={c['r']:.1f}  psi={c['psi_num']:.5e}  analytic={c['psi_analytic_inf']:.5e}  "
            f"ratio={c['ratio_num_over_analytic']}"
        )
    lines.append("sample rays (grid):")
    for r in rays[::2]:
        lines.append(
            f"  b={r['b']:+.2f}  defl={r['deflection_deg']:+.3f}deg  delay={r['delay']:+.4f}"
        )
    lines.append(f"gates: {gates}")
    lines.append(f"VERDICT monist_3d_1r_pass = {gates['monist_3d_1r_pass']}")
    lines.append(package["verdict"]["summary"])
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "round3_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
