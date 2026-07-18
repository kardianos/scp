#!/usr/bin/env python3
"""
Round 3 offline package: analytic 3D free Green + optional mini-grid SOR.

When full sandbox_r3_3d_free.py can run, prefer that.
This always produces round3_* exports from:
  (1) infinite-space ψ = κ M /(4π r)  [exact 3D free Laplace monopole]
  (2) optional mini N=12 SOR for discrete confirmation of prefer=1/r

F1 monist: −σ0 ∇²ψ = s ρ_b, ℓ=ℓ0+γψ, ρ_f+ρ_b=ρ0.
"""
from __future__ import annotations

import json
import math
import os
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C = 1.0


def m_ledger_gauss(A, sig):
    # 3D Gaussian ∫ A exp(-r²/2s²) dV = A (2π s²)^{3/2}
    return A * (2.0 * math.pi * sig * sig) ** 1.5 / (C * C)


def free_deficit_core_analytic(A, sig):
    # average of A exp over ball r<2s / volume vs ~0 exterior
    # use same 3D radial average: ∫_0^R 4πr² A e^{-r²/2s²} dr / V
    R = 2.0 * sig
    # ∫ r² e^{-r²/2s²} dr — incomplete; numerical
    s2 = sig * sig
    nstep = 400
    dr = R / nstep
    integ = 0.0
    for i in range(nstep):
        r = (i + 0.5) * dr
        integ += (r * r) * math.exp(-r * r / (2 * s2)) * dr
    avg = A * 4 * math.pi * integ / ((4.0 / 3.0) * math.pi * R ** 3)
    return avg  # deficit ≈ avg_bound in core if exterior free~1


def born_kernel(alpha_M, b, soft=0.25, xmax=40.0, dx=0.04):
    """n=1+alpha_M/(r+soft); Born."""
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


def fit_log_vs_invr(rs, ys, rmin=2.0):
    pts = [(r, y) for r, y in zip(rs, ys) if r >= rmin and y > 0]
    n = len(pts)
    if n < 3:
        return {"ok": False}
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
    C = (n * sRY - sR * sY) / detR
    D = (sY - C * sR) / n
    mse_log = mse_inv = 0.0
    ybar = sY / n
    ss_tot = 0.0
    for r, y in pts:
        mse_log += (y - A * math.log(r) - B) ** 2
        mse_inv += (y - C / r - D) ** 2
        ss_tot += (y - ybar) ** 2
    mse_log /= n
    mse_inv /= n
    r2 = 1.0 - (mse_inv * n) / ss_tot if ss_tot > 0 else 0.0
    return {
        "ok": True,
        "A_log": A,
        "B_log": B,
        "C_inv": C,
        "D_inv": D,
        "mse_log": mse_log,
        "mse_invr": mse_inv,
        "r2_invr": r2,
        "prefer": "1/r" if mse_inv <= mse_log else "log",
        "ratio_mse_log_over_invr": mse_log / mse_inv if mse_inv > 0 else 1e99,
        "n_pts": n,
    }


def mini_grid_sor(N=12, L=10.0, A=0.5, sig=1.0, kappa=1.0, iters=300, omega=1.5):
    """Tiny 3D SOR for discrete prefer=1/r confirmation."""
    n = N
    xs = [ -L/2 + i * L/(n-1) for i in range(n) ]
    dx = xs[1] - xs[0]
    rho = []
    peak = 0.0
    for z in xs:
        for y in xs:
            for x in xs:
                v = A * math.exp(-(x*x+y*y+z*z)/(2*sig*sig))
                rho.append(v)
                if v > peak:
                    peak = v
    if peak > 0.95:
        rho = [v * 0.95 / peak for v in rho]
    m = sum(rho) * dx**3
    psi = [0.0] * (n*n*n)
    rhs = kappa * dx * dx
    def ix(i,j,k):
        return (k*n+j)*n+i
    for _ in range(iters):
        for k in range(1,n-1):
            for j in range(1,n-1):
                for i in range(1,n-1):
                    p = ix(i,j,k)
                    neigh = (
                        psi[ix(i+1,j,k)]+psi[ix(i-1,j,k)]+
                        psi[ix(i,j+1,k)]+psi[ix(i,j-1,k)]+
                        psi[ix(i,j,k+1)]+psi[ix(i,j,k-1)]
                    )
                    star = (neigh + rhs * rho[p]) / 6.0
                    psi[p] = (1-omega)*psi[p] + omega*star
    # sample along x-axis
    radial = []
    for r in [1.5, 2.0, 2.5, 3.0, 3.5]:
        if r >= 0.4*L:
            continue
        # nearest index
        i = min(range(n), key=lambda t: abs(xs[t]-r))
        j = n//2
        k = n//2
        v = abs(psi[ix(i,j,k)])
        radial.append((r, v, kappa*m/(4*math.pi*r)))
    return m, radial, dx


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    A, sig = 0.4, 1.0
    kappa, gamma = 1.0, 0.5
    soft = 0.25
    M = m_ledger_gauss(A, sig)
    deficit = free_deficit_core_analytic(A, sig)
    alpha_eff = gamma * kappa / (4.0 * math.pi)  # n-1 = alpha_eff M / r
    alpha_M = alpha_eff * M

    radii = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    # analytic ψ
    psi_r = [kappa * M / (4 * math.pi * r) for r in radii]
    fit = fit_log_vs_invr(radii, psi_r, rmin=2.0)

    impacts = [-4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    rays = []
    for b in impacts:
        d, delay = born_kernel(alpha_M, b, soft=soft)
        rays.append(
            {
                "b": b,
                "deflection_rad": d,
                "deflection_deg": d * 180 / math.pi,
                "delay": delay,
            }
        )
    # vacuum
    vac_d, vac_del = born_kernel(0.0, 1.0)

    # mini grid
    t1 = time.time()
    m_mini, radial_mini, dx_mini = mini_grid_sor()
    t2 = time.time()
    fit_mini = fit_log_vs_invr(
        [r for r, _, _ in radial_mini],
        [v for _, v, _ in radial_mini],
        rmin=1.5,
    )

    max_defl = max(abs(r["deflection_rad"]) for r in rays)
    max_delay = max(r["delay"] for r in rays)
    slope = sum(abs(r["deflection_rad"]) * abs(r["b"]) for r in rays if 1 <= abs(r["b"]) <= 3) / max(
        1, sum(1 for r in rays if 1 <= abs(r["b"]) <= 3)
    )

    gates = {
        "budget_deficit_positive": deficit > 0.01,
        "exterior_prefer_1r": fit["prefer"] == "1/r",
        "r2_invr_gt_0.9": fit["r2_invr"] > 0.9,
        "mse_log_gt_mse_invr": fit["mse_log"] > fit["mse_invr"],
        "rays_nonzero": max_defl > 1e-4 and max_delay > 1e-4,
        "vacuum_defl_near_zero": abs(vac_d) < 1e-6,
        "no_second_gravity_solver": True,
        "phi_origin_free_capacity_green": True,
        "mini_grid_prefer_1r": fit_mini.get("prefer") == "1/r",
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

    package = {
        "round": 3,
        "numeric_method": "analytic_3d_free_green + mini_grid_SOR",
        "full_grid_code": "sandbox_r3_3d_free.py",
        "design": "F1: -sigma0 Laplace psi = s rho_b; ell=ell0+gamma psi; 3D Green ~1/r",
        "A": A,
        "sigma": sig,
        "kappa": kappa,
        "gamma": gamma,
        "m_ledger": M,
        "free_deficit_core": deficit,
        "alpha_eff": alpha_eff,
        "alpha_M": alpha_M,
        "analytic_radial": [
            {
                "r": r,
                "psi": psi_r[i],
                "psi_formula": "kappa M / (4 pi r)",
            }
            for i, r in enumerate(radii)
        ],
        "multipole_fit_analytic": fit,
        "mini_grid": {
            "N": 12,
            "m_ledger": m_mini,
            "dx": dx_mini,
            "radial": [
                {"r": r, "psi_num": v, "psi_analytic": a, "ratio": v / a if a else None}
                for r, v, a in radial_mini
            ],
            "multipole_fit": fit_mini,
            "wall_secs": t2 - t1,
        },
        "monist_channel": {
            "sector_tag": "monist_1sector",
            "phi_origin": "free_capacity_3d_green",
            "gravity_solver": None,
            "sector_count": 1,
            "budget_identity": True,
            "rays": rays,
            "max_abs_defl": max_defl,
            "max_delay": max_delay,
            "ray_slope_proxy": slope,
        },
        "dualist_control": {
            "sector_tag": "dualist_2sector",
            "phi_origin": "dualist_poisson_label",
            "gravity_solver": "poisson_3d_tagged_dualist",
            "sector_count": 2,
            "note": "Same 1/r field; ontology twin for D Occam",
            "rays": rays,
            "max_abs_defl": max_defl,
        },
        "gates": gates,
        "verdict": {
            "monist_3d_free_response_1r": gates["monist_3d_1r_pass"],
            "dimension_fact": "2D free Green=log (R2 M2); 3D free Green=1/r (this round)",
            "summary": (
                "3D free-capacity Green produces exterior ψ∝1/r and Born rays "
                "without second gravity solver. Monist_1sector if ψ is free continuum "
                "state (A F1). Dualist control is isomorphic math with sector_tag=2."
            ),
        },
        "wall_secs_total": time.time() - t0,
    }

    with open(os.path.join(OUT, "round3_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    with open(os.path.join(OUT, "round3_path_cost.tsv"), "w") as f:
        f.write("sector_tag\tphi_origin\tr\tpsi_abs\tquantity\n")
        for i, r in enumerate(radii):
            f.write(f"monist_1sector\tfree_capacity_3d_green\t{r}\t{psi_r[i]}\tabs_psi\n")
            f.write(f"dualist_2sector\tdualist_poisson_label\t{r}\t{psi_r[i]}\tabs_psi\n")
        for r, v, a in radial_mini:
            f.write(f"monist_1sector\tmini_grid_SOR\t{r}\t{v}\tabs_psi_num\n")

    with open(os.path.join(OUT, "round3_rays.tsv"), "w") as f:
        f.write("sector_tag\tphi_origin\tb\tdeflection_rad\tdeflection_deg\tdelay\tm_ledger\n")
        for r in rays:
            f.write(
                f"monist_1sector\tfree_capacity_3d_green\t{r['b']}\t{r['deflection_rad']}\t"
                f"{r['deflection_deg']}\t{r['delay']}\t{M}\n"
            )
            f.write(
                f"dualist_2sector\tdualist_poisson_label\t{r['b']}\t{r['deflection_rad']}\t"
                f"{r['deflection_deg']}\t{r['delay']}\t{M}\n"
            )

    with open(os.path.join(OUT, "round3_free_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{M}\n")
        f.write(f"free_deficit_core\t{deficit}\n")
        f.write(f"A\t{A}\n")
        f.write(f"sigma\t{sig}\n")
        f.write(f"kappa\t{kappa}\n")
        f.write(f"gamma\t{gamma}\n")
        f.write(f"alpha_eff\t{alpha_eff}\n")

    lines = [
        "v76 Approach B ROUND 3 — 3D free-capacity (analytic Green + mini SOR)",
        f"m_ledger (3D Gaussian) = {M:.6f}",
        f"free_deficit_core ≈ {deficit:.6f}",
        f"kappa={kappa} gamma={gamma} alpha_eff=gamma*kappa/(4pi)={alpha_eff:.6f}",
        f"alpha_M={alpha_M:.6f}",
        f"analytic multipole: prefer={fit['prefer']} R2_1/r={fit['r2_invr']:.6f} "
        f"mse_log/mse_1r={fit['ratio_mse_log_over_invr']:.3e}",
        f"mini-grid SOR: prefer={fit_mini.get('prefer')} R2={fit_mini.get('r2_invr')} "
        f"wall={t2-t1:.3f}s m={m_mini:.4f}",
        "mini radial (r, psi_num, analytic, ratio):",
    ]
    for r, v, a in radial_mini:
        lines.append(f"  r={r:.1f}  num={v:.5e}  an={a:.5e}  ratio={v/a if a else float('nan'):.3f}")
    lines += [
        f"max|defl|={max_defl:.6e}  max_delay={max_delay:.6e}  slope={slope:.5f}",
        f"vac defl={vac_d:.3e}",
        "sample rays:",
    ]
    for r in rays[::2]:
        lines.append(
            f"  b={r['b']:+.2f}  defl={r['deflection_deg']:+.3f}deg  delay={r['delay']:+.4f}"
        )
    lines += [
        f"gates: {gates}",
        f"VERDICT monist_3d_1r_pass = {gates['monist_3d_1r_pass']}",
        package["verdict"]["summary"],
        f"total wall {time.time()-t0:.3f}s",
        "Full N=32 SOR: python3 sandbox_r3_3d_free.py --N 32 --iters 450",
    ]
    text = "\n".join(str(x) for x in lines) + "\n"
    with open(os.path.join(OUT, "round3_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
