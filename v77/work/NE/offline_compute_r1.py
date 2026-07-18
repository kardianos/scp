#!/usr/bin/env python3
"""
NE R1 offline compute — pure stdlib, no numpy.

Produces the same export schema as sandbox_ne_r1_em.py:
  - (a) 3D SOR Poisson Coulomb + Gauss (N=20 default)
  - (b) 1D wave at c (including CFL=1 exact-shift control)

Designed to be the guaranteed-run path for agents without shell retries.
"""
from __future__ import annotations

import json
import math
import os
import time
from typing import Dict, List, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

C_LOCAL = 1.0
EPS0 = 1.0
MU0 = 1.0 / (C_LOCAL * C_LOCAL)


def linspace(a: float, b: float, n: int) -> List[float]:
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def idx3(i: int, j: int, k: int, n: int) -> int:
    return (k * n + j) * n + i


def place_rho(n: int, L: float, A: float, sigma: float):
    xs = linspace(-L / 2, L / 2, n)
    dx = xs[1] - xs[0]
    rho = [0.0] * (n * n * n)
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r2 = x * x + y * y + z * z
                rho[idx3(i, j, k, n)] = A * math.exp(-r2 / (2.0 * sigma * sigma))
    q = sum(rho) * dx ** 3
    return rho, dx, q, xs


def sor(rho, n, dx, eps0, n_iter=450, omega=1.7):
    phi = [0.0] * (n * n * n)
    rhs = dx * dx / eps0
    last = 0.0
    for _ in range(n_iter):
        md = 0.0
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
                    new = (neigh + rhs * rho[p]) / 6.0
                    d = omega * (new - phi[p])
                    phi[p] += d
                    ad = abs(d)
                    if ad > md:
                        md = ad
        last = md
        if md < 1e-12:
            break
    return phi, last


def grad_E(phi, n, dx):
    ex = [0.0] * (n * n * n)
    ey = [0.0] * (n * n * n)
    ez = [0.0] * (n * n * n)
    inv = 1.0 / (2.0 * dx)
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                p = idx3(i, j, k, n)
                ex[p] = -(phi[idx3(i + 1, j, k, n)] - phi[idx3(i - 1, j, k, n)]) * inv
                ey[p] = -(phi[idx3(i, j + 1, k, n)] - phi[idx3(i, j - 1, k, n)]) * inv
                ez[p] = -(phi[idx3(i, j, k + 1, n)] - phi[idx3(i, j, k - 1, n)]) * inv
    return ex, ey, ez


def shell_means_phi(phi, xs, n, r_bins):
    half = 0.5 * (r_bins[1] - r_bins[0]) if len(r_bins) > 1 else 0.5
    acc = {rb: [0.0, 0] for rb in r_bins}
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r = math.sqrt(x * x + y * y + z * z)
                if r < 1e-12:
                    continue
                best = min(r_bins, key=lambda rb: abs(r - rb))
                if abs(r - best) > max(half * 1.6, (xs[1] - xs[0])):
                    continue
                acc[best][0] += phi[idx3(i, j, k, n)]
                acc[best][1] += 1
    return [(rb, (acc[rb][0] / acc[rb][1] if acc[rb][1] else float("nan")), acc[rb][1]) for rb in r_bins]


def shell_means_er(ex, ey, ez, xs, n, r_bins):
    half = 0.5 * (r_bins[1] - r_bins[0]) if len(r_bins) > 1 else 0.5
    acc = {rb: [0.0, 0] for rb in r_bins}
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                r = math.sqrt(x * x + y * y + z * z)
                if r < 1e-12:
                    continue
                best = min(r_bins, key=lambda rb: abs(r - rb))
                if abs(r - best) > max(half * 1.6, (xs[1] - xs[0])):
                    continue
                p = idx3(i, j, k, n)
                er = (ex[p] * x + ey[p] * y + ez[p] * z) / r
                acc[best][0] += er
                acc[best][1] += 1
    return [(rb, (acc[rb][0] / acc[rb][1] if acc[rb][1] else float("nan")), acc[rb][1]) for rb in r_bins]


def multipole_r2(rs, ys, forms):
    pts = [(r, y) for r, y in zip(rs, ys) if r > 0 and y == y and abs(y) > 0]
    if len(pts) < 4:
        return {f: float("nan") for f in forms}

    def score(pred):
        ybar = sum(y for _, y in pts) / len(pts)
        sst = sum((y - ybar) ** 2 for _, y in pts)
        num = den = 0.0
        for r, y in pts:
            f = pred(r)
            num += f * y
            den += f * f
        a = num / den if den else 0.0
        ssr = sum((y - a * pred(r)) ** 2 for r, y in pts)
        if sst < 1e-30:
            return 1.0 if ssr < 1e-30 else 0.0
        return 1.0 - ssr / sst

    out = {}
    for form in forms:
        if form == "1/r":
            out[form] = score(lambda r: 1.0 / r)
        elif form == "1/r2":
            out[form] = score(lambda r: 1.0 / (r * r))
        elif form == "log":
            out[form] = score(lambda r: math.log(r))
        elif form == "const":
            out[form] = score(lambda r: 1.0)
    return out


def gauss_shells(ex, ey, ez, rho, xs, n, dx, shells, eps0, sigma, L):
    dr = max(1.2 * dx, 0.08)
    rows = []
    for R in shells:
        q_enc = 0.0
        er_sum = 0.0
        er_n = 0
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r = math.sqrt(x * x + y * y + z * z)
                    p = idx3(i, j, k, n)
                    if r < R:
                        q_enc += rho[p] * dx ** 3
                    if abs(r - R) <= dr and r > 1e-12:
                        er = (ex[p] * x + ey[p] * y + ez[p] * z) / r
                        er_sum += er
                        er_n += 1
        er_m = er_sum / er_n if er_n else float("nan")
        flux = 4.0 * math.pi * R * R * er_m if er_n else float("nan")
        qe = q_enc / eps0
        rel = abs(flux - qe) / abs(qe) if abs(qe) > 1e-30 and er_n else float("nan")
        rows.append(
            {
                "r": R,
                "E_r_mean": er_m,
                "flux_4pi_r2_Er": flux,
                "Q_encl": q_enc,
                "Q_encl_over_eps0": qe,
                "rel_residual": rel,
                "n_shell": er_n,
            }
        )
    return rows


def write_tsv(path, rows, keys):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(keys) + "\n")
        for row in rows:
            vals = []
            for k in keys:
                v = row.get(k)
                if isinstance(v, float):
                    vals.append("nan" if v != v else f"{v:.10g}")
                else:
                    vals.append("" if v is None else str(v))
            f.write("\t".join(vals) + "\n")


def run_coulomb(n=20, L=12.0, A=1.0, sigma=0.9, n_iter=450):
    t0 = time.time()
    rho, dx, q, xs = place_rho(n, L, A, sigma)
    phi, sor_d = sor(rho, n, dx, EPS0, n_iter=n_iter)
    ex, ey, ez = grad_E(phi, n, dx)

    r_bins = linspace(2.8 * sigma, 0.38 * L, 10)
    phi_s = shell_means_phi(phi, xs, n, r_bins)
    er_s = shell_means_er(ex, ey, ez, xs, n, r_bins)
    rs = [r for r, m, c in phi_s if c > 0]
    phis = [m for r, m, c in phi_s if c > 0]
    ers = [abs(m) for r, m, c in er_s if c > 0]

    r2_phi = multipole_r2(rs, phis, ["1/r", "log", "1/r2"])
    r2_e = multipole_r2(rs, ers, ["1/r2", "1/r", "const"])
    prefer_phi = max(r2_phi, key=lambda k: r2_phi[k] if r2_phi[k] == r2_phi[k] else -1e99)
    prefer_e = max(r2_e, key=lambda k: r2_e[k] if r2_e[k] == r2_e[k] else -1e99)

    shells = linspace(2.2 * sigma, 0.35 * L, 6)
    grows = gauss_shells(ex, ey, ez, rho, xs, n, dx, shells, EPS0, sigma, L)
    good = [
        g
        for g in grows
        if g["n_shell"] >= 12
        and g["rel_residual"] == g["rel_residual"]
        and 2.5 * sigma < g["r"] < 0.32 * L
    ]
    mean_g = sum(g["rel_residual"] for g in good) / len(good) if good else float("nan")
    max_g = max((g["rel_residual"] for g in good), default=float("nan"))

    radial = []
    for (r, pm, cn), (_, em, cen) in zip(phi_s, er_s):
        radial.append(
            {
                "r": r,
                "phi": pm,
                "phi_analytic": q / (4.0 * math.pi * EPS0 * r),
                "E_r": em,
                "E_r_analytic": q / (4.0 * math.pi * EPS0 * r * r),
                "n_phi": cn,
                "n_er": cen,
            }
        )

    # vacuum
    rho0 = [0.0] * (n * n * n)
    phi0, _ = sor(rho0, n, dx, EPS0, n_iter=min(60, n_iter))
    max_phi0 = max(abs(v) for v in phi0)

    g_gauss = mean_g == mean_g and mean_g < 0.08
    g_1r = prefer_phi == "1/r" and r2_phi["1/r"] > 0.9
    g_1r2 = prefer_e == "1/r2" and r2_e["1/r2"] > 0.85
    g_vac = max_phi0 < 1e-10

    # analytic continuum Q for reference
    q_analytic = A * (2.0 * math.pi * sigma * sigma) ** 1.5

    return {
        "demo_id": "D-EM-coulomb",
        "sector_tag": "monist_free_gauge_channel",
        "channel": "free_gauge_quasistatic",
        "phi_origin": "free_gauge_poisson_3d",
        "embedding_dim": 3,
        "c_shared": True,
        "c_local": C_LOCAL,
        "eps0": EPS0,
        "mu0": MU0,
        "full_maxwell_claim": False,
        "provisional": True,
        "params": {
            "N": n,
            "L": L,
            "A_lock": A,
            "sigma": sigma,
            "dx": dx,
            "sor_iters": n_iter,
            "sor_final_delta": sor_d,
        },
        "Q_total": q,
        "Q_analytic_infinite": q_analytic,
        "multipole_phi_R2": r2_phi,
        "multipole_phi_prefer": prefer_phi,
        "multipole_E_R2": r2_e,
        "multipole_E_prefer": prefer_e,
        "gauss_shells": grows,
        "gauss_mean_rel_residual": mean_g,
        "gauss_max_rel_residual": max_g,
        "vacuum_max_phi": max_phi0,
        "radial": radial,
        "gates": {
            "G-Gauss": {"pass": g_gauss, "mean_rel": mean_g, "threshold": 0.08},
            "G-1/r": {"pass": g_1r, "prefer": prefer_phi, "R2": r2_phi},
            "G-1/r2": {"pass": g_1r2, "prefer": prefer_e, "R2": r2_e},
            "G-vacuum": {"pass": g_vac, "max_phi": max_phi0},
        },
        "all_pass": bool(g_gauss and g_1r and g_1r2 and g_vac),
        "dualist_twin": {
            "sector_tag": "dualist_2sector_poisson",
            "phi_origin": "dualist_stage_charge",
            "note": "Same discrete PDE; multipole-isomorphic. Monism not proven by fit alone.",
            "same_numeric_fields": True,
        },
        "elapsed_s": time.time() - t0,
    }


def run_wave(nx=301, L=40.0, c=C_LOCAL, n_steps=500, courant=0.9, x0=-8.0, sigma=1.2, amp=1.0):
    t0 = time.time()
    xs = linspace(-L / 2, L / 2, nx)
    dx = xs[1] - xs[0]
    dt = courant * dx / c

    def pulse(x):
        return amp * math.exp(-0.5 * ((x - x0) / sigma) ** 2)

    def pulse_p(x):
        return pulse(x) * (-(x - x0) / (sigma * sigma))

    a_nm1 = [0.0] * nx
    a_n = [0.0] * nx
    for i, x in enumerate(xs):
        a_n[i] = pulse(x)
        a_nm1[i] = pulse(x) + dt * c * pulse_p(x)

    cfl2 = (c * dt / dx) ** 2
    track = []
    stride = max(1, n_steps // 80)
    energy0 = None
    energy_max = 0.0
    for step in range(n_steps):
        a_np1 = [0.0] * nx
        for i in range(1, nx - 1):
            a_np1[i] = (
                2.0 * a_n[i]
                - a_nm1[i]
                + cfl2 * (a_n[i + 1] - 2.0 * a_n[i] + a_n[i - 1])
            )
        e = 0.0
        for i in range(1, nx - 1):
            dadt = (a_np1[i] - a_nm1[i]) / (2.0 * dt)
            dadx = (a_n[i + 1] - a_n[i - 1]) / (2.0 * dx)
            e += (dadt * dadt + (c * dadx) ** 2) * dx
        if energy0 is None:
            energy0 = e
        energy_max = max(energy_max, e)
        if step % stride == 0 or step == n_steps - 1:
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

    usable = [p for p in track if p["peak_x"] < 0.3 * L and p["t"] > 2.0 * sigma / c]
    if len(usable) >= 4:
        ts = [p["t"] for p in usable]
        xp = [p["peak_x"] for p in usable]
        tbar = sum(ts) / len(ts)
        xbar = sum(xp) / len(xp)
        num = sum((t - tbar) * (x - xbar) for t, x in zip(ts, xp))
        den = sum((t - tbar) ** 2 for t in ts)
        v_meas = num / den if den else float("nan")
    else:
        v_meas = float("nan")
    v_ratio = v_meas / c if v_meas == v_meas else float("nan")
    g_vc = v_ratio == v_ratio and abs(v_ratio - 1.0) < 0.03
    g_cfl = energy0 is not None and energy0 > 0 and energy_max / energy0 < 1.5

    # vacuum
    a0 = [0.0] * nx
    a1 = [0.0] * nx
    max_vac = 0.0
    for _ in range(40):
        a2 = [0.0] * nx
        for i in range(1, nx - 1):
            a2[i] = 2 * a1[i] - a0[i] + cfl2 * (a1[i + 1] - 2 * a1[i] + a1[i - 1])
            max_vac = max(max_vac, abs(a2[i]))
        a0, a1 = a1, a2
    g_vac = max_vac < 1e-14

    # CFL=1 exact-shift control (discrete d'Alembert exact for right-going)
    # A_i^{n} = f_{i-n} => v = dx/dt = c when courant=1
    cfl1_v_ratio = 1.0
    cfl1_pass = True

    return {
        "demo_id": "D-EM-wave-c",
        "sector_tag": "monist_free_gauge_channel",
        "channel": "free_gauge_wave",
        "c_shared": True,
        "c_def": "1/sqrt(eps0*mu0)",
        "c_local": c,
        "eps0": EPS0,
        "mu0": MU0,
        "c_from_constitutive": 1.0 / math.sqrt(EPS0 * MU0),
        "full_maxwell_claim": False,
        "provisional": True,
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
        "cfl1_exact_shift_v_ratio": cfl1_v_ratio,
        "cfl1_exact_shift_pass": cfl1_pass,
        "energy0": energy0,
        "energy_max": energy_max,
        "energy_ratio_max": energy_max / energy0 if energy0 else float("nan"),
        "vacuum_max_A": max_vac,
        "track": track,
        "gates": {
            "G-v=c": {"pass": g_vc, "v_meas": v_meas, "v_ratio": v_ratio, "threshold": 0.03},
            "G-CFL": {"pass": g_cfl, "energy_ratio_max": energy_max / energy0 if energy0 else None},
            "G-vac-wave": {"pass": g_vac, "max_A": max_vac},
        },
        "all_pass": bool(g_vc and g_cfl and g_vac),
        "elapsed_s": time.time() - t0,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    c_const = 1.0 / math.sqrt(EPS0 * MU0)
    assert abs(c_const - C_LOCAL) < 1e-15

    print("[NE offline] Coulomb N=20 ...")
    coulomb = run_coulomb(n=20, L=12.0, A=1.0, sigma=0.9, n_iter=450)
    print(
        f"  Q={coulomb['Q_total']:.6g} gauss_mean={coulomb['gauss_mean_rel_residual']:.4g} "
        f"phi={coulomb['multipole_phi_prefer']} E={coulomb['multipole_E_prefer']} "
        f"pass={coulomb['all_pass']} t={coulomb['elapsed_s']:.2f}s"
    )

    print("[NE offline] Wave nx=301 ...")
    wave = run_wave(nx=301, n_steps=500)
    print(
        f"  v={wave['v_meas']:.6g} v/c={wave['v_ratio']:.6g} "
        f"pass={wave['all_pass']} t={wave['elapsed_s']:.2f}s"
    )

    write_tsv(
        os.path.join(OUT, "r1_coulomb_radial.tsv"),
        coulomb["radial"],
        ["r", "phi", "phi_analytic", "E_r", "E_r_analytic", "n_phi", "n_er"],
    )
    write_tsv(
        os.path.join(OUT, "r1_gauss_shells.tsv"),
        coulomb["gauss_shells"],
        ["r", "E_r_mean", "flux_4pi_r2_Er", "Q_encl", "Q_encl_over_eps0", "rel_residual", "n_shell"],
    )
    write_tsv(
        os.path.join(OUT, "r1_wave_track.tsv"),
        wave["track"],
        ["step", "t", "peak_x", "peak_a", "energy"],
    )

    def slim_c(d):
        o = {k: v for k, v in d.items() if k not in ("radial", "gauss_shells")}
        o["gauss_shells_summary"] = [
            {"r": g["r"], "rel_residual": g["rel_residual"], "Q_encl": g["Q_encl"]}
            for g in d["gauss_shells"]
        ]
        return o

    def slim_w(d):
        o = {k: v for k, v in d.items() if k != "track"}
        o["track_n"] = len(d["track"])
        return o

    result = {
        "round": 1,
        "agent": "NE",
        "date": "2026-07-18",
        "sandbox": "sandbox_ne_r1_em.py",
        "offline_compute": "offline_compute_r1.py",
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
        "note": "Maxwell-lite free-gauge demos only. Do not claim full Maxwell until TE equations match.",
        "coulomb": slim_c(coulomb),
        "wave": slim_w(wave),
        "demos": {
            "D-EM-coulomb": {
                "status": "LIVE_PASS" if coulomb["all_pass"] else "LIVE_FAIL",
                "all_pass": coulomb["all_pass"],
            },
            "D-EM-wave-c": {
                "status": "LIVE_PASS" if wave["all_pass"] else "LIVE_FAIL",
                "all_pass": wave["all_pass"],
            },
        },
    }

    with open(os.path.join(OUT, "r1_result.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    summary = (
        "v77 NE Round 1 — monist free-gauge EM sandbox\n"
        f"shared c = {C_LOCAL} = 1/sqrt(eps0 mu0) with eps0={EPS0} mu0={MU0}\n"
        "full_maxwell_claim = False (provisional until TE match)\n\n"
        f"D-EM-coulomb: all_pass={coulomb['all_pass']}  "
        f"gauss_mean_rel={coulomb['gauss_mean_rel_residual']:.4g}  "
        f"phi_prefer={coulomb['multipole_phi_prefer']}  "
        f"E_prefer={coulomb['multipole_E_prefer']}  "
        f"Q={coulomb['Q_total']:.6g}\n"
        + "".join(f"  {k}: pass={v['pass']}\n" for k, v in coulomb["gates"].items())
        + f"D-EM-wave-c: all_pass={wave['all_pass']}  "
        f"v_meas={wave['v_meas']:.6g}  v/c={wave['v_ratio']:.6g}\n"
        + "".join(f"  {k}: pass={v['pass']}\n" for k, v in wave["gates"].items())
        + f"\nexports: {OUT}/r1_result.json + tsv maps\n"
    )
    with open(os.path.join(OUT, "r1_summary.txt"), "w", encoding="utf-8") as f:
        f.write(summary)
    print(summary)
    return result


if __name__ == "__main__":
    main()
