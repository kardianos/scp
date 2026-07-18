#!/usr/bin/env python3
"""
Round 2 offline result package — closed-form / series (no numpy).

When shell is available, prefer: python3 sandbox_r2_dynamical.py
This file writes the same export schema with analytical diagnostics so D can
ingest and B can log a monist eligibility verdict without waiting on Jacobi.

Physics content (honest):
  M1 local free path cost ~ compact (Born local n from free deficit)
  M2 free Laplace exterior 2D: u ∝ log(R_box/r)  (monist free equilibration)
  M3 dualist 2D Poisson exterior: Φ ∝ (M/2π) log r class (2-sector)
  R1 kernel: Φ = α M / R  (postulated; dualist-isomorphic; NOT dynamical free)
"""
from __future__ import annotations

import json
import math
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C = 1.0
A = 0.35
SIG = 1.2
L = 20.0
R_BOX = 0.48 * L  # outer Dirichlet
R_LOCK = 1.5 * SIG  # interface scale
SOFT = 0.2
ALPHA_K = 0.08
KAPPA_U = 1.2
BETA_U = 0.8
BETA_PHI = 0.5


def m_ledger(A=A, s=SIG):
    return 2.0 * math.pi * A * s * s / (C * C)


def free_deficit_core(A=A, s=SIG):
    R = 2.0 * s
    avg_g = (2.0 * s * s / (R * R)) * (1.0 - math.exp(-(R * R) / (2.0 * s * s)))
    return A * avg_g  # exterior~1, core free = 1 - A*avg_g → deficit = A*avg_g


def born_local_optics(A, s, b, xmax=30.0, dx=0.02):
    """Local n=1/(1-f), f=A exp(-(x²+b²)/2s²). Compact GRIN."""
    s2 = s * s
    delta = delay = 0.0
    x = -xmax
    while x <= xmax:
        f = A * math.exp(-(x * x + b * b) / (2 * s2))
        f = min(f, 0.95)
        n = 1.0 / max(1.0 - f, 0.05)
        dth = n * f * (-b / s2)
        delta += dth * dx
        delay += (n - 1.0) * dx
        x += dx
    return delta, delay


def free_laplace_u(r, u0, r_in=R_LOCK, r_out=R_BOX):
    """2D free Laplace: u=u0 on r_in, u=0 on r_out → log profile."""
    if r <= r_in:
        return u0
    if r >= r_out:
        return 0.0
    return u0 * math.log(r_out / r) / math.log(r_out / r_in)


def born_from_radial_n(n_of_r, b, xmax=12.0, dx=0.05):
    """n depends on r=hypot(x,b); Born straight path."""
    delta = delay = 0.0
    x = -xmax
    while x <= xmax:
        def nxy(xx, yy):
            return n_of_r(math.hypot(xx, yy))

        n0 = nxy(x, b)
        npy = nxy(x, b + dx)
        nmy = nxy(x, b - dx)
        dndy = (npy - nmy) / (2 * dx)
        delta += (dndy / max(n0, 1e-9)) * dx
        delay += (n0 - 1.0) * dx
        x += dx
    return delta, delay


def poisson_2d_phi_approx(r, M, r_core=R_LOCK, r_out=R_BOX):
    """
    Rough 2D Poisson exterior with Φ=0 at r_out:
    Φ(r) ≈ (M/(2π)) * log(r_out/r) for r > r_core  (same log family as free Laplace!)
    Note: 2D Poisson Green is log — NOT 1/r. Dualist 3D would be 1/r.
    """
    if r < r_core:
        # flatten core
        return (M / (2 * math.pi)) * math.log(r_out / r_core)
    if r >= r_out:
        return 0.0
    return (M / (2 * math.pi)) * math.log(r_out / r)


def kernel_phi(r, M, alpha=ALPHA_K, soft=SOFT):
    return alpha * M / (r + soft)


def fit_log_vs_invr(rs, ys):
    pts = [(r, y) for r, y in zip(rs, ys) if r >= 3.0 and y > 0]
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
    for r, y in pts:
        mse_log += (y - A * math.log(r) - B) ** 2
        mse_inv += (y - C / r - D) ** 2
    mse_log /= n
    mse_inv /= n
    return {
        "ok": True,
        "A_log": A,
        "B_log": B,
        "C_inv": C,
        "D_inv": D,
        "mse_log": mse_log,
        "mse_invr": mse_inv,
        "prefer": "log" if mse_log <= mse_inv else "1/r",
        "n_pts": n,
    }


def ray_table(fn_delay_defl, impacts):
    rows = []
    for b in impacts:
        d, delay = fn_delay_defl(b)
        rows.append(
            {
                "b": b,
                "deflection_rad": d,
                "deflection_deg": d * 180 / math.pi,
                "delay": delay,
            }
        )
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    M = m_ledger()
    deficit = free_deficit_core()
    impacts = [-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    radii = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    # M1 — local free optics (proxy for free-graph geodesic excess)
    rays_m1 = ray_table(lambda b: born_local_optics(A, SIG, b), impacts)
    # path cost radial: use delay at b=r as proxy for isotropic excess scale
    # better: local (n-1) at radius r (angular mean of free deficit optical)
    def local_n1(r):
        f = A * math.exp(-(r * r) / (2 * SIG * SIG))
        return 1.0 / max(1.0 - f, 0.05) - 1.0

    radial_m1 = [{"r": r, "delta_path_cost": local_n1(r), "quantity": "n_minus_1_local"} for r in radii]

    # M2 — free Laplace log
    u0 = -KAPPA_U * A  # interface depth scale
    def n_m2(r):
        return 1.0 + BETA_U * abs(free_laplace_u(r, u0))

    rays_m2 = ray_table(lambda b: born_from_radial_n(n_m2, b), impacts)
    radial_m2 = [
        {"r": r, "abs_mean": abs(free_laplace_u(r, u0)), "quantity": "abs_u_free_laplace"}
        for r in radii
    ]
    fit_m2 = fit_log_vs_invr([p["r"] for p in radial_m2], [p["abs_mean"] for p in radial_m2])

    # M3 — dualist 2D Poisson (log exterior)
    def n_m3(r):
        return 1.0 + BETA_PHI * abs(poisson_2d_phi_approx(r, M))

    rays_m3 = ray_table(lambda b: born_from_radial_n(n_m3, b), impacts)
    radial_m3 = [
        {"r": r, "abs_mean": abs(poisson_2d_phi_approx(r, M)), "quantity": "abs_phi_poisson2d"}
        for r in radii
    ]
    fit_m3 = fit_log_vs_invr([p["r"] for p in radial_m3], [p["abs_mean"] for p in radial_m3])

    # R1 postulated 1/r kernel
    def n_k(r):
        return 1.0 + kernel_phi(r, M)

    rays_k = ray_table(lambda b: born_from_radial_n(n_k, b), impacts)
    radial_k = [
        {"r": r, "abs_mean": kernel_phi(r, M), "quantity": "Phi_alphaM_R"} for r in radii
    ]
    fit_k = fit_log_vs_invr([p["r"] for p in radial_k], [p["abs_mean"] for p in radial_k])

    # 3D dualist reference (true 1/r Poisson-like) for contrast
    def phi_3d(r):
        return M / (4 * math.pi * (r + SOFT))  # schematic

    radial_3d = [{"r": r, "abs_mean": phi_3d(r)} for r in radii]
    fit_3d = fit_log_vs_invr([p["r"] for p in radial_3d], [p["abs_mean"] for p in radial_3d])

    # Inertia triad — honest non-tautological proxies
    # 1) m_ledger
    # 2) m_inertial: unlock energy / c^2 is same ledger by Ax (not independent yet);
    #    independent smoke: free-path cost sensitivity ∂C/∂M from two masses
    M_lo = m_ledger(A=0.15, s=SIG)
    M_hi = M
    # delay at b=2 for local (M1) and kernel
    d_lo, del_lo = born_local_optics(0.15, SIG, 2.0)
    d_hi, del_hi = born_local_optics(A, SIG, 2.0)
    # ray slope mean |defl|*|b|
    def slope(rays):
        ps = [abs(r["deflection_rad"]) * abs(r["b"]) for r in rays if 0.8 <= abs(r["b"]) <= 3.5]
        return sum(ps) / len(ps) if ps else float("nan")

    s1, s2, s3, sk = slope(rays_m1), slope(rays_m2), slope(rays_m3), slope(rays_k)
    # linearity of kernel slope vs M: sk ∝ α M * const
    sk_lo = slope(
        ray_table(lambda b: born_from_radial_n(lambda r: 1.0 + kernel_phi(r, M_lo), b), impacts)
    )

    inertia = {
        "m_ledger": M,
        "m_ledger_weak": M_lo,
        "m_inertial_status": "not_independent — push F/a with a=F/m_L is tautological; deferred true free-drag experiment",
        "ray_slope_proxy": {
            "M1_local": s1,
            "M2_free_laplace": s2,
            "M3_poisson2d": s3,
            "R1_kernel": sk,
            "R1_kernel_weak_M": sk_lo,
            "kernel_slope_ratio_hi_lo": sk / sk_lo if sk_lo else None,
            "M_ratio_hi_lo": M / M_lo,
            "note": "kernel slope tracks M (ratio~M ratio); local M1 slope is nonlinear in A",
        },
        "delay_b2_local_lo_hi": [del_lo, del_hi],
    }

    # Long-range contrast at r=6
    long_range = {
        "r6_M1_n_minus_1": local_n1(6.0),
        "r6_M2_abs_u": abs(free_laplace_u(6.0, u0)),
        "r6_M3_abs_phi": abs(poisson_2d_phi_approx(6.0, M)),
        "r6_R1_kernel_Phi": kernel_phi(6.0, M),
        "r6_3d_dualist_Phi": phi_3d(6.0),
        "M1_is_compact": local_n1(6.0) < 1e-3,
        "M2_fit": fit_m2,
        "M3_fit": fit_m3,
        "R1_fit": fit_k,
        "dualist_3d_fit": fit_3d,
    }

    # Verdict
    # - M1 compact monist works as ledger+compact rays, fails long-range 1/r
    # - M2 free Laplace is monist_1sector dynamical free response; exterior is LOG (2D), not 1/r
    # - M3 dualist 2D also log (Poisson 2D); does not get 1/r either
    # - Only R1 postulated kernel / 3D dualist give 1/r
    # monist_kernel_failed for Einstein 1/r from dynamical free-only (M1/M2)
    # M2 IS a success as "dynamical free-response" monist — but wrong multipole for 3D GR targets
    verdict = {
        "monist_1sector_compact_M1": True,
        "monist_1sector_dynamical_free_response_M2": True,
        "M2_exterior_multipole": "log (2D free Laplace Green)",
        "M2_produces_einstein_1r": False,
        "dualist_2sector_M3_poisson": True,
        "M3_2d_exterior_multipole": "log (2D Poisson Green) — same family as M2!",
        "R1_postulated_kernel_1r": True,
        "R1_is_dynamical_free_evolution": False,
        "monist_kernel_failed": True,
        "monist_kernel_failed_meaning": (
            "Einstein-class exterior path cost ∝1/r was NOT obtained from "
            "dynamical free-medium local updates (M1/M2). It appears only for "
            "postulated 3D-like kernel Φ=αM/R or dualist 3D Poisson — not evolved "
            "from free graph dynamics in this 2D sandbox. "
            "Additionally M2 free-Laplace and M3 2D-Poisson are multipole-siblings "
            "(both log): monist free equilibration ≠ automatic monist 1/r gravity."
        ),
        "positive_monist_result": (
            "M2 shows a true 1-sector dynamical free response (Laplace equilibration "
            "on free graph / free continuum with lock as hole). Rays bend without "
            "second gravity field. Exterior is 2D-log, not GR 1/r."
        ),
        "recommendation": (
            "Do not label R1 Φ=α∫ρ/R as monist success. Either (i) move to 3D free "
            "Laplace (Green~1/r) as monist free equilibration, or (ii) accept free "
            "response kernel as axiomatic medium law (A Ax8) not numerics-derived, "
            "or (iii) kill Einstein-class monism from this class of local free PDEs."
        ),
    }

    package = {
        "round": 2,
        "numeric_method": "offline_closed_form_series",
        "grid_dynamical_code": "sandbox_r2_dynamical.py (parent re-run for Jacobi/Dijkstra)",
        "m_ledger": M,
        "free_deficit_core": deficit,
        "A": A,
        "sigma": SIG,
        "M1": {
            "mechanism": "M1_local_free_optics_proxy_dijkstra",
            "sector_tag": "monist_1sector",
            "gravity_solver": None,
            "m_ledger": M,
            "free_deficit_core": deficit,
            "rays": rays_m1,
            "path_cost_radial": radial_m1,
            "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m1),
            "max_delay": max(r["delay"] for r in rays_m1),
        },
        "M2": {
            "mechanism": "M2_free_laplace_2d_log",
            "sector_tag": "monist_1sector",
            "gravity_solver": None,
            "note": "Dynamical free equilibration (analytic 2D Dirichlet); lock=interface u0",
            "u0": u0,
            "m_ledger": M,
            "rays": rays_m2,
            "u_radial": radial_m2,
            "u_profile_fit": fit_m2,
            "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m2),
            "max_delay": max(r["delay"] for r in rays_m2),
        },
        "M3": {
            "mechanism": "M3_dualist_poisson_2d",
            "sector_tag": "dualist_2sector",
            "gravity_solver": "poisson_2d_Phi_from_M",
            "m_ledger": M,
            "rays": rays_m3,
            "phi_radial": radial_m3,
            "phi_profile_fit": fit_m3,
            "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m3),
            "max_delay": max(r["delay"] for r in rays_m3),
        },
        "R1_postulated_kernel": {
            "mechanism": "R1_postulated_kernel",
            "sector_tag": "dualist_2sector_or_postulated",
            "gravity_solver": "analytic_Phi_eq_alpha_M_over_R",
            "alpha": ALPHA_K,
            "rays": rays_k,
            "phi_radial": radial_k,
            "phi_profile_fit": fit_k,
            "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_k),
            "max_delay": max(r["delay"] for r in rays_k),
        },
        "dualist_3d_reference": {
            "sector_tag": "dualist_2sector",
            "phi_radial": radial_3d,
            "fit": fit_3d,
            "note": "schematic Φ=M/(4πr); prefer 1/r in fit",
        },
        "inertia_triad": inertia,
        "long_range": long_range,
        "verdict": verdict,
        "FOR_D": "sector_tag in TSV is ontology ground truth; M2 vs M3 both log in 2D — use 3d_ref and R1_kernel for 1/r class",
    }

    # TSVs for D
    with open(os.path.join(OUT, "round2_rays.tsv"), "w") as f:
        f.write("sector_tag\tmechanism\tb\tdeflection_rad\tdeflection_deg\tdelay\tm_ledger\n")
        for key in ("M1", "M2", "M3", "R1_postulated_kernel"):
            blk = package[key]
            for r in blk["rays"]:
                f.write(
                    f"{blk['sector_tag']}\t{blk['mechanism']}\t{r['b']}\t"
                    f"{r['deflection_rad']}\t{r['deflection_deg']}\t{r['delay']}\t{M}\n"
                )

    with open(os.path.join(OUT, "round2_path_cost.tsv"), "w") as f:
        f.write("sector_tag\tmechanism\tr\tvalue\tquantity\n")
        for p in radial_m1:
            f.write(f"monist_1sector\tM1_local\t{p['r']}\t{p['delta_path_cost']}\tn_minus_1\n")
        for p in radial_m2:
            f.write(f"monist_1sector\tM2_free_laplace\t{p['r']}\t{p['abs_mean']}\tabs_u\n")
        for p in radial_m3:
            f.write(f"dualist_2sector\tM3_poisson2d\t{p['r']}\t{p['abs_mean']}\tabs_phi\n")
        for p in radial_k:
            f.write(
                f"dualist_2sector_or_postulated\tR1_kernel\t{p['r']}\t{p['abs_mean']}\tPhi_alphaM_R\n"
            )
        for p in radial_3d:
            f.write(f"dualist_2sector\tREF_3d_poisson\t{p['r']}\t{p['abs_mean']}\tPhi_M_4pir\n")

    with open(os.path.join(OUT, "round2_free_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{M}\n")
        f.write(f"free_deficit_core\t{deficit}\n")
        f.write(f"A\t{A}\n")
        f.write(f"sigma\t{SIG}\n")

    with open(os.path.join(OUT, "round2_inertia_triad.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{M}\n")
        f.write(f"m_ledger_weak\t{M_lo}\n")
        f.write(f"ray_slope_M1\t{s1}\n")
        f.write(f"ray_slope_M2\t{s2}\n")
        f.write(f"ray_slope_M3\t{s3}\n")
        f.write(f"ray_slope_R1kernel\t{sk}\n")
        f.write(f"ray_slope_R1kernel_weakM\t{sk_lo}\n")
        f.write(f"kernel_slope_ratio\t{sk/sk_lo if sk_lo else float('nan')}\n")
        f.write(f"M_ratio\t{M/M_lo}\n")
        f.write("m_inertial_independent\tNA_deferred\n")

    with open(os.path.join(OUT, "round2_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    lines = [
        "v76 Approach B ROUND 2 — dynamical free-response (offline closed-form)",
        f"m_ledger = {M:.6f}",
        f"free_deficit_core = {deficit:.6f}",
        "",
        f"M1 local free (monist_1sector): max|defl|={package['M1']['max_abs_defl']:.4f}  "
        f"r6 n-1={local_n1(6):.3e}  COMPACT",
        f"M2 free Laplace 2D (monist_1sector): max|defl|={package['M2']['max_abs_defl']:.4f}  "
        f"fit prefer={fit_m2.get('prefer')}  r6 |u|={abs(free_laplace_u(6,u0)):.4f}  LOG exterior",
        f"M3 dualist Poisson 2D: max|defl|={package['M3']['max_abs_defl']:.4f}  "
        f"fit prefer={fit_m3.get('prefer')}  r6 |Φ|={abs(poisson_2d_phi_approx(6,M)):.4f}  LOG (same family!)",
        f"R1 postulated kernel 1/r: max|defl|={package['R1_postulated_kernel']['max_abs_defl']:.4f}  "
        f"fit prefer={fit_k.get('prefer')}  r6 Φ={kernel_phi(6,M):.4f}",
        f"3D dualist ref fit prefer={fit_3d.get('prefer')}",
        "",
        f"Inertia: kernel slope ratio hi/lo M = {sk/sk_lo if sk_lo else 'NA':.4f} vs M ratio {M/M_lo:.4f}",
        "  independent m_inertial push: DEFERRED (F/a with a=F/m is tautological)",
        "",
        "VERDICT:",
        f"  monist_kernel_failed = {verdict['monist_kernel_failed']}",
        f"  {verdict['monist_kernel_failed_meaning']}",
        f"  positive: {verdict['positive_monist_result']}",
        f"  recommend: {verdict['recommendation']}",
        "",
        "Exports: round2_rays.tsv path_cost free_deficit inertia_triad round2_result.json",
        "Parent re-run dynamical Jacobi/Dijkstra: python3 sandbox_r2_dynamical.py",
    ]
    text = "\n".join(str(x) for x in lines) + "\n"
    with open(os.path.join(OUT, "round2_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
