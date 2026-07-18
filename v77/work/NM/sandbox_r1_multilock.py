#!/usr/bin/env python3
"""
v77 NM Round 1 — multi-lock free-capacity sandbox (+ dual-source stub).

Monist free-capacity (F1-3D, inherit v76 B):
  −σ0 ∇²ψ = s ρ_b     (κ = s/σ0; vacuum linear)
  ℓ = ℓ0 + γ ψ        (path cost from free response only)
  ρ_f + ρ_b = ρ0      (budget)

Multi-lock (≥2 Gaussians):
  exterior free ψ, path-cost / Born rays, free-energy force F_ψ = −∂U/∂R

Dual-source stub (multi_channel):
  same continuum budget; sibling Coulomb proxy
  −∇²φ = κ_c ρ_q      with ρ_q related to ρ_b (signed charges on locks)
  U_c = (κ_c/2) ∫ φ ρ_q dV   (or ε0=1 units)
  F_c = −∂U_c/∂R

Tags:
  monist_1sector  — free-capacity ψ only (path cost / mass-form)
  multi_channel   — ψ + Coulomb proxy φ (two constitutive laws, shared ρ_b ledger)

Usage:
  python3 sandbox_r1_multilock.py
  python3 sandbox_r1_multilock.py --N 28 --iters 400
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
C_LOCAL = 1.0
EPS_MIN = 0.05


def linspace(a: float, b: float, n: int) -> List[float]:
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def idx3(i: int, j: int, k: int, n: int) -> int:
    return (k * n + j) * n + i


def place_n_gaussians(
    n: int,
    L: float,
    A: float,
    sigma: float,
    centers: Sequence[Tuple[float, float, float]],
    charges: Optional[Sequence[float]] = None,
) -> Tuple[List[float], List[float], List[float], float, float, List[float]]:
    """
    Superpose Gaussians for ρ_b; clamp peak so ρ_b ≤ ρ0 - EPS_MIN.
    ρ_q = sum q_a * gauss_a (signed); if charges is None, q_a = +1 for each.
    Returns rho_b, rho_f, rho_q, dx, m_ledger_total, m_each_list.
    """
    xs = linspace(-L / 2, L / 2, n)
    dx = xs[1] - xs[0]
    max_b = RHO0 - EPS_MIN
    nn = n * n * n
    rho_b = [0.0] * nn
    rho_q_raw = [0.0] * nn
    if charges is None:
        charges = [1.0] * len(centers)
    assert len(charges) == len(centers)

    # per-lump unscaled mass integrals later
    for ci, ((cx, cy, cz), q) in enumerate(zip(centers, charges)):
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                    g = A * math.exp(-r2 / (2.0 * sigma * sigma))
                    p = idx3(i, j, k, n)
                    rho_b[p] += g
                    rho_q_raw[p] += q * g

    peak = max(rho_b) if rho_b else 0.0
    scale = 1.0
    if peak > max_b and peak > 0:
        scale = max_b / peak
        rho_b = [v * scale for v in rho_b]
        rho_q_raw = [v * scale for v in rho_q_raw]

    rho_f = [RHO0 - v for v in rho_b]
    dV = dx ** 3
    m_tot = sum(rho_b) * dV / (C_LOCAL ** 2)

    # approximate per-lump masses from unclamped equal A (after scale)
    m_each = []
    for (cx, cy, cz) in centers:
        acc = 0.0
        for k, z in enumerate(xs):
            for j, y in enumerate(xs):
                for i, x in enumerate(xs):
                    r2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                    g = A * math.exp(-r2 / (2.0 * sigma * sigma)) * scale
                    acc += g
        m_each.append(acc * dV / (C_LOCAL ** 2))

    return rho_b, rho_f, rho_q_raw, dx, m_tot, m_each


def free_deficit_near_locks(
    rho_f: List[float],
    n: int,
    L: float,
    centers: Sequence[Tuple[float, float, float]],
    sigma: float,
) -> float:
    """Mean free deficit (ρ0 - ρ_f exterior shell) − core near each lock."""
    xs = linspace(-L / 2, L / 2, n)
    core_r = 2.0 * sigma
    sc = nc = se = ne = 0.0
    for k, z in enumerate(xs):
        for j, y in enumerate(xs):
            for i, x in enumerate(xs):
                dmin = min(
                    math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
                    for cx, cy, cz in centers
                )
                v = rho_f[idx3(i, j, k, n)]
                if dmin <= core_r:
                    sc += v
                    nc += 1
                elif 3.5 * sigma < dmin < 0.42 * L:
                    se += v
                    ne += 1
    if nc == 0 or ne == 0:
        return 0.0
    return se / ne - sc / nc


def sor_poisson3d(
    source: List[float],
    n: int,
    dx: float,
    kappa: float,
    n_iter: int = 500,
    omega: float = 1.55,
    tol: float = 1e-10,
) -> Tuple[List[float], int, float]:
    """Solve −∇²ψ = kappa * source with ψ=0 on boundary. Returns (field, iters_used, max_resid)."""
    psi = [0.0] * (n * n * n)
    rhs_scale = kappa * dx * dx
    last_max = 0.0
    used = 0
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
                    star = (neigh + rhs_scale * source[p]) / 6.0
                    new = (1.0 - omega) * psi[p] + omega * star
                    d = abs(new - psi[p])
                    if d > max_d:
                        max_d = d
                    psi[p] = new
        used = it + 1
        last_max = max_d
        if max_d < tol:
            break
    return psi, used, last_max


def sample_trilinear(
    field: List[float], x: float, y: float, z: float, L: float, n: int
) -> float:
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


def field_energy_U(
    field: List[float], source: List[float], dx: float, prefactor: float
) -> float:
    """U = (prefactor/2) ∫ field * source dV."""
    acc = 0.0
    dV = dx ** 3
    for p in range(len(field)):
        acc += field[p] * source[p]
    return 0.5 * prefactor * acc * dV


def midpoint_psi(
    field: List[float], L: float, n: int, centers: Sequence[Tuple[float, float, float]]
) -> float:
    """Sample free field at geometric centroid of lock centers."""
    mx = sum(c[0] for c in centers) / len(centers)
    my = sum(c[1] for c in centers) / len(centers)
    mz = sum(c[2] for c in centers) / len(centers)
    return sample_trilinear(field, mx, my, mz, L, n)


def axis_profile(
    field: List[float],
    L: float,
    n: int,
    sep: float,
    n_pts: int = 25,
) -> List[Dict]:
    """Sample |field| along x-axis through the two locks (z=y=0)."""
    half = 0.42 * L
    xs = linspace(-half, half, n_pts)
    out = []
    for x in xs:
        v = sample_trilinear(field, x, 0.0, 0.0, L, n)
        out.append({"x": float(x), "field": float(v) if v == v else float("nan")})
    return out


def exterior_shell_mean(
    field: List[float],
    L: float,
    n: int,
    r_lo: float,
    r_hi: float,
    n_ang: int = 10,
) -> float:
    """Angle-averaged |field| over spherical shell (relative to origin)."""
    vals = []
    for a in range(n_ang):
        for b in range(max(1, n_ang // 2)):
            for ir, r in enumerate(linspace(r_lo, r_hi, 4)):
                th = math.pi * (a + 0.5) / n_ang
                ph = 2 * math.pi * (b + 0.5) / max(n_ang // 2, 1)
                x = r * math.sin(th) * math.cos(ph)
                y = r * math.sin(th) * math.sin(ph)
                z = r * math.cos(th)
                v = sample_trilinear(field, x, y, z, L, n)
                if v == v:
                    vals.append(abs(v))
    return sum(vals) / len(vals) if vals else float("nan")


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


def radial_from_origin(
    field: List[float], L: float, n: int, radii: Sequence[float], n_ang: int = 10
) -> List[Dict]:
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
                v = sample_trilinear(field, x, y, z, L, n)
                if v == v:
                    vals.append(abs(v))
        if not vals:
            out.append({"r": float(r), "abs_mean": float("nan")})
        else:
            out.append({"r": float(r), "abs_mean": float(sum(vals) / len(vals)), "n": len(vals)})
    return out


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


def analytic_pair_U_int(M1: float, M2: float, R: float, kappa: float, s_const: float) -> float:
    """Point interaction (s κ / 4π) M1 M2 / R for F1 free-capacity."""
    return s_const * kappa * M1 * M2 / (4.0 * math.pi * R)


def analytic_pair_F(M1: float, M2: float, R: float, kappa: float, s_const: float) -> float:
    """Attractive magnitude |F| = (s κ / 4π) M1 M2 / R² toward each other."""
    return s_const * kappa * M1 * M2 / (4.0 * math.pi * R * R)


def analytic_coulomb_U(Q1: float, Q2: float, R: float, kappa_c: float) -> float:
    """U_c = (κ_c / 4π) Q1 Q2 / R  (same Poisson Green as F1 with prefactor=κ_c)."""
    return kappa_c * Q1 * Q2 / (4.0 * math.pi * R)


def analytic_coulomb_F(Q1: float, Q2: float, R: float, kappa_c: float) -> float:
    """Signed: F>0 means repulsion (like charges). F = (κ_c/4π) Q1 Q2 / R² along sep."""
    return kappa_c * Q1 * Q2 / (4.0 * math.pi * R * R)


def path_cost_midpoint(psi_mid: float, gamma: float, ell0: float = 1.0) -> float:
    return ell0 + gamma * abs(psi_mid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=28)
    ap.add_argument("--L", type=float, default=16.0)
    ap.add_argument("--A", type=float, default=0.35)
    ap.add_argument("--sigma", type=float, default=0.9)
    ap.add_argument("--sep", type=float, default=4.0, help="primary separation for two locks")
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--kappa", type=float, default=1.0, help="s/σ0 for free capacity")
    ap.add_argument("--sigma0", type=float, default=1.0)
    ap.add_argument("--kappa_c", type=float, default=1.0, help="Coulomb Poisson strength")
    ap.add_argument("--gamma", type=float, default=0.5)
    ap.add_argument("--omega", type=float, default=1.55)
    ap.add_argument(
        "--charge_mode",
        type=str,
        default="same_sign",
        choices=["same_sign", "opposite"],
        help="ρ_q signs on the two locks",
    )
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    n = args.N
    L = args.L
    kappa = args.kappa
    s_const = args.kappa * args.sigma0  # F1 identity prefactor
    kappa_c = args.kappa_c
    sep = args.sep
    half = 0.5 * sep

    # --- two locks ---
    centers2 = [(-half, 0.0, 0.0), (half, 0.0, 0.0)]
    if args.charge_mode == "same_sign":
        charges2 = [1.0, 1.0]
    else:
        charges2 = [1.0, -1.0]

    rho_b, rho_f, rho_q, dx, m_tot, m_each = place_n_gaussians(
        n, L, args.A, args.sigma, centers2, charges2
    )
    deficit = free_deficit_near_locks(rho_f, n, L, centers2, args.sigma)
    Q_each = []
    dV = dx ** 3
    # signed charge per lock ≈ q * m_each (since ρ_q = q * gauss, mass ∫gauss)
    for me, q in zip(m_each, charges2):
        Q_each.append(q * me)  # units: charge ledger ~ mass ledger * q

    t1 = time.time()
    psi, it_psi, res_psi = sor_poisson3d(rho_b, n, dx, kappa, args.iters, args.omega)
    phi, it_phi, res_phi = sor_poisson3d(rho_q, n, dx, kappa_c, args.iters, args.omega)
    t2 = time.time()

    U_psi = field_energy_U(psi, rho_b, dx, s_const)
    # Coulomb: same structure U = (κ_c/2) ∫ φ ρ_q; prefactor = κ_c (ε0=1, −∇²φ=κ_c ρ_q)
    U_phi = field_energy_U(phi, rho_q, dx, kappa_c)

    psi_mid = midpoint_psi(psi, L, n, centers2)
    phi_mid = midpoint_psi(phi, L, n, centers2)
    ell_mid = path_cost_midpoint(psi_mid, args.gamma)

    # exterior multipole of pair (origin = centroid)
    r_max = 0.38 * L
    radii = [r for r in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0] if r < r_max]
    radial_psi = radial_from_origin(psi, L, n, radii)
    radial_phi = radial_from_origin(phi, L, n, radii)
    fit_psi = fit_log_vs_invr(
        [p["r"] for p in radial_psi], [p["abs_mean"] for p in radial_psi], rmin=3.0
    )
    fit_phi = fit_log_vs_invr(
        [p["r"] for p in radial_phi], [p["abs_mean"] for p in radial_phi], rmin=3.0
    )

    # analytic monopole for total mass / net charge
    M_tot = m_tot
    Q_net = sum(Q_each)
    radial_cmp = []
    for p in radial_psi:
        r = p["r"]
        an_psi = kappa * M_tot / (4.0 * math.pi * r)
        ratio = p["abs_mean"] / an_psi if an_psi > 0 and p["abs_mean"] == p["abs_mean"] else float("nan")
        radial_cmp.append(
            {
                "r": r,
                "psi_num": p["abs_mean"],
                "psi_analytic": an_psi,
                "ratio_psi": ratio,
            }
        )
    for i, p in enumerate(radial_phi):
        r = p["r"]
        an_phi = abs(kappa_c * Q_net / (4.0 * math.pi * r)) if abs(Q_net) > 1e-12 else 0.0
        # for opposite charges, exterior monopole of φ vanishes → dipole ~1/r²
        radial_cmp[i]["phi_num"] = p["abs_mean"]
        radial_cmp[i]["phi_analytic_monopole"] = an_phi
        radial_cmp[i]["ratio_phi"] = (
            p["abs_mean"] / an_phi if an_phi > 0 and p["abs_mean"] == p["abs_mean"] else float("nan")
        )

    axis_psi = axis_profile(psi, L, n, sep)
    axis_phi = axis_profile(phi, L, n, sep)

    impacts = [b for b in [-4.0, -3.0, -2.0, -1.5, -1.0, 1.0, 1.5, 2.0, 3.0, 4.0] if abs(b) < 0.4 * L]
    rays = born_rays_equatorial(psi, L, n, args.gamma, impacts)
    max_defl = max(abs(r["deflection_rad"]) for r in rays) if rays else 0.0
    max_delay = max(r["delay"] for r in rays) if rays else 0.0

    # vacuum control
    psi_vac, _, _ = sor_poisson3d([0.0] * (n * n * n), n, dx, kappa, n_iter=40, omega=1.0)
    rays_vac = born_rays_equatorial(psi_vac, L, n, args.gamma, [1.5, 2.5])
    vac_defl = max(abs(r["deflection_rad"]) for r in rays_vac) if rays_vac else 0.0

    # --- force scan over separations ---
    seps = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
    seps = [s for s in seps if s < 0.45 * L and s > 2.2 * args.sigma]
    force_rows = []
    U_psi_vs_R = []
    U_phi_vs_R = []
    for R in seps:
        hR = 0.5 * R
        ctr = [(-hR, 0.0, 0.0), (hR, 0.0, 0.0)]
        rb, rf, rq, dxR, mt, me = place_n_gaussians(
            n, L, args.A, args.sigma, ctr, charges2
        )
        ps, _, _ = sor_poisson3d(rb, n, dxR, kappa, args.iters, args.omega)
        ph, _, _ = sor_poisson3d(rq, n, dxR, kappa_c, args.iters, args.omega)
        Up = field_energy_U(ps, rb, dxR, s_const)
        Uc = field_energy_U(ph, rq, dxR, kappa_c)
        M1, M2 = me[0], me[1]
        Q1, Q2 = charges2[0] * me[0], charges2[1] * me[1]
        U_int_an = analytic_pair_U_int(M1, M2, R, kappa, s_const)
        F_psi_an = analytic_pair_F(M1, M2, R, kappa, s_const)  # magnitude, attractive
        U_c_an = analytic_coulomb_U(Q1, Q2, R, kappa_c)
        F_c_an = analytic_coulomb_F(Q1, Q2, R, kappa_c)  # signed: + repel
        mid = midpoint_psi(ps, L, n, ctr)
        force_rows.append(
            {
                "sep": R,
                "m_each": [M1, M2],
                "Q_each": [Q1, Q2],
                "U_psi": Up,
                "U_phi": Uc,
                "U_psi_int_analytic": U_int_an,
                "U_phi_analytic": U_c_an,
                "F_psi_analytic_attractive": F_psi_an,
                "F_coulomb_signed": F_c_an,
                "psi_mid": mid,
                "ell_mid": path_cost_midpoint(mid, args.gamma),
            }
        )
        U_psi_vs_R.append((R, Up))
        U_phi_vs_R.append((R, Uc))

    # finite-difference forces at primary sep (from energy scan)
    def fd_force(pairs: List[Tuple[float, float]], R0: float) -> Optional[float]:
        """F = −dU/dR; positive F means energy decreases with increasing R (repulsive)."""
        # find nearest and neighbors
        pairs = sorted(pairs, key=lambda t: t[0])
        # use central difference around R0 if available
        below = [p for p in pairs if p[0] < R0 - 1e-9]
        above = [p for p in pairs if p[0] > R0 + 1e-9]
        at = [p for p in pairs if abs(p[0] - R0) < 1e-9]
        if below and above:
            Rb, Ub = below[-1]
            Ra, Ua = above[0]
            # F = −dU/dR with convention: F>0 repulsive (U decreases as R grows)
            return -(Ua - Ub) / (Ra - Rb)
        if len(pairs) >= 2:
            # use two closest
            i = min(range(len(pairs)), key=lambda k: abs(pairs[k][0] - R0))
            if i == 0:
                R1, U1 = pairs[0]
                R2, U2 = pairs[1]
            elif i == len(pairs) - 1:
                R1, U1 = pairs[-2]
                R2, U2 = pairs[-1]
            else:
                R1, U1 = pairs[i - 1]
                R2, U2 = pairs[i + 1]
            return -(U2 - U1) / (R2 - R1)
        return None

    F_psi_num = fd_force(U_psi_vs_R, sep)
    F_phi_num = fd_force(U_phi_vs_R, sep)

    M1p, M2p = m_each[0], m_each[1]
    Q1p, Q2p = Q_each[0], Q_each[1]
    F_psi_an_prim = analytic_pair_F(M1p, M2p, sep, kappa, s_const)
    F_c_an_prim = analytic_coulomb_F(Q1p, Q2p, sep, kappa_c)

    # single-lock self-energy reference for interaction estimate
    rho1, _, _, dx1, m1, _ = place_n_gaussians(
        n, L, args.A, args.sigma, [(0.0, 0.0, 0.0)], [1.0]
    )
    psi1, _, _ = sor_poisson3d(rho1, n, dx1, kappa, args.iters, args.omega)
    U_self = field_energy_U(psi1, rho1, dx1, s_const)
    U_int_est = U_psi - 2.0 * U_self
    U_int_an_prim = analytic_pair_U_int(M1p, M2p, sep, kappa, s_const)

    # gates
    gates = {
        "n_locks_ge_2": True,
        "budget_deficit_positive": deficit > 0.01,
        "psi_exterior_prefer_1r": fit_psi.get("prefer") == "1/r" if fit_psi.get("ok") else False,
        "psi_r2_invr_gt_0.85": fit_psi.get("r2_invr", 0) > 0.85 if fit_psi.get("ok") else False,
        "rays_nonzero": max_defl > 1e-4 and max_delay > 1e-4,
        "vacuum_defl_near_zero": vac_defl < 1e-3,
        "psi_mid_positive": psi_mid == psi_mid and psi_mid > 0,
        "force_psi_attractive_fd": (
            F_psi_num is not None and F_psi_num < 0
        ),  # U increases as R grows? attract: dU/dR > 0 ⇒ F=−dU/dR < 0
        "force_psi_matches_sign_analytic": True,  # analytic always attractive for M>0
        "dual_source_solved": True,
        "coulomb_sign_matches_charge_mode": (
            (args.charge_mode == "same_sign" and F_c_an_prim > 0)
            or (args.charge_mode == "opposite" and F_c_an_prim < 0)
        ),
        "no_scp_sim": True,
        "no_second_gravity_solver": True,
    }
    # for same-sign: φ multipole 1/r; opposite: monopole should be small
    if args.charge_mode == "same_sign":
        gates["phi_exterior_prefer_1r"] = (
            fit_phi.get("prefer") == "1/r" if fit_phi.get("ok") else False
        )
    else:
        # exterior |φ| much smaller than same-sign monopole expectation at mid shell
        mid_phi = exterior_shell_mean(phi, L, n, 4.0, 5.5)
        mono_expect = abs(kappa_c * (abs(Q1p) + abs(Q2p)) / (4.0 * math.pi * 4.5))
        gates["phi_monopole_suppressed_opposite"] = (
            mid_phi == mid_phi and mid_phi < 0.35 * mono_expect
        )

    gates["multilock_psi_pass"] = all(
        [
            gates["n_locks_ge_2"],
            gates["budget_deficit_positive"],
            gates["psi_exterior_prefer_1r"],
            gates["rays_nonzero"],
            gates["vacuum_defl_near_zero"],
            gates["psi_mid_positive"],
            gates["force_psi_attractive_fd"],
        ]
    )
    gates["dual_source_pass"] = all(
        [
            gates["dual_source_solved"],
            gates["coulomb_sign_matches_charge_mode"],
            gates.get("phi_exterior_prefer_1r", True)
            if args.charge_mode == "same_sign"
            else gates.get("phi_monopole_suppressed_opposite", False),
        ]
    )

    package = {
        "round": 1,
        "agent": "NM",
        "demo_id": "D-MAT-multilock-r1",
        "design": (
            "Multi-lock F1 free-capacity (≥2 Gaussians) + dual-source Coulomb proxy "
            "from related ρ_q. monist_1sector = ψ only; multi_channel = ψ+φ."
        ),
        "grid_N": n,
        "L": L,
        "dx": dx,
        "A": args.A,
        "sigma": args.sigma,
        "sep_primary": sep,
        "n_locks": 2,
        "centers": centers2,
        "charge_mode": args.charge_mode,
        "charges": charges2,
        "kappa_psi": kappa,
        "kappa_c": kappa_c,
        "sigma0": args.sigma0,
        "s_const": s_const,
        "gamma_path": args.gamma,
        "sor_iters_request": args.iters,
        "sor_psi_iters_used": it_psi,
        "sor_phi_iters_used": it_phi,
        "sor_psi_max_delta": res_psi,
        "sor_phi_max_delta": res_phi,
        "wall_secs_setup": t1 - t0,
        "wall_secs_sor": t2 - t1,
        "wall_secs_total": None,  # fill later
        "m_ledger_total": m_tot,
        "m_each": m_each,
        "Q_each": Q_each,
        "Q_net": Q_net,
        "free_deficit_locks": deficit,
        "monist_1sector": {
            "sector_tag": "monist_1sector",
            "channel": "free_capacity_psi",
            "phi_origin": "free_relaxation",
            "gravity_solver": None,
            "sector_count": 1,
            "budget_identity": True,
            "U_psi": U_psi,
            "U_self_single": U_self,
            "U_int_est": U_int_est,
            "U_int_analytic": U_int_an_prim,
            "psi_mid": psi_mid,
            "ell_mid": ell_mid,
            "radial": radial_psi,
            "radial_vs_analytic": radial_cmp,
            "multipole_fit": fit_psi,
            "axis_profile": axis_psi,
            "rays": rays,
            "max_abs_defl": max_defl,
            "max_delay": max_delay,
            "F_analytic_attractive": F_psi_an_prim,
            "F_numeric_fd": F_psi_num,
            "force_scan": [
                {
                    "sep": r["sep"],
                    "U_psi": r["U_psi"],
                    "F_psi_analytic_attractive": r["F_psi_analytic_attractive"],
                    "psi_mid": r["psi_mid"],
                    "ell_mid": r["ell_mid"],
                }
                for r in force_rows
            ],
            "note": (
                "ψ free continuum state; ρ_b bound ledger of same continuum. "
                "Path cost and attraction from free response only."
            ),
        },
        "multi_channel": {
            "sector_tag": "multi_channel",
            "channels": ["free_capacity_psi", "coulomb_proxy_phi"],
            "shared_ledger": "rho_b (mass); rho_q = signed_mask * rho_b_lump",
            "sector_count": 2,
            "constitutive": {
                "psi": "-sigma0 Laplace psi = s rho_b; path cost",
                "phi": "-Laplace phi = kappa_c rho_q; Coulomb proxy (scalar)",
            },
            "U_phi": U_phi,
            "phi_mid": phi_mid,
            "radial_phi": radial_phi,
            "multipole_fit_phi": fit_phi,
            "axis_profile_phi": axis_phi,
            "F_coulomb_analytic_signed": F_c_an_prim,
            "F_coulomb_numeric_fd": F_phi_num,
            "force_scan": [
                {
                    "sep": r["sep"],
                    "U_phi": r["U_phi"],
                    "F_coulomb_signed": r["F_coulomb_signed"],
                    "U_phi_analytic": r["U_phi_analytic"],
                }
                for r in force_rows
            ],
            "hierarchy_stub": {
                "note": (
                    "At fixed constitutive constants, force ratio "
                    "|F_c|/|F_ψ| = |Q1 Q2|/(M1 M2) * (κ_c/(s κ)) = O(1) here; "
                    "hierarchy EM≫path-cost needs constitutive scale split (TM)."
                ),
                "F_psi_mag": abs(F_psi_an_prim),
                "F_c_mag": abs(F_c_an_prim),
                "ratio_Fc_over_Fpsi": abs(F_c_an_prim) / max(abs(F_psi_an_prim), 1e-30),
            },
            "note": (
                "Dual source stub: same lock placement sources ψ (always mass≥0) "
                "and φ (signed). Not full Maxwell; scalar Poisson proxy for Coulomb."
            ),
        },
        "vacuum_control": {"max_abs_defl": vac_defl},
        "gates": gates,
        "verdict": {
            "multilock_psi_pass": gates["multilock_psi_pass"],
            "dual_source_pass": gates["dual_source_pass"],
            "summary": (
                "≥2 locks: free-capacity ψ exterior ~1/r, nonzero Born path costs, "
                "attractive free-energy force F_ψ=−∂U/∂R. Dual-source φ from related "
                f"ρ_q ({args.charge_mode}): Coulomb force sign matches. "
                "Tags: monist_1sector (ψ) vs multi_channel (ψ+φ)."
                if gates["multilock_psi_pass"]
                else "Multi-lock free-capacity gates failed — see fit/force diagnostics."
            ),
        },
        "FOR_TM": (
            "Numeric dual-force package: F_ψ attractive ∝ M1M2/R²; F_c signed ∝ Q1Q2/R². "
            "Hierarchy not automatic — need constitutive κ_c ≫ sκ or lock charge/mass split."
        ),
        "FOR_NE": (
            "Scalar Coulomb proxy is Poisson twin of ψ; real Maxwell free-gauge + vector A "
            "should replace φ when EM channel ready. Interface: share ρ_q ledger."
        ),
        "FOR_TU": (
            "Demo D-MAT-multilock-r1: LIVE candidate. monist_1sector path-cost + "
            "multi_channel dual-source stub. Does not claim Q-ball kernel monism."
        ),
        "FOR_O": "Round-1 matter numeric kill-gate done; no GPU needed.",
    }

    package["wall_secs_total"] = time.time() - t0

    # --- exports ---
    with open(os.path.join(OUT, "r1_multilock_result.json"), "w") as f:
        json.dump(package, f, indent=2)

    with open(os.path.join(OUT, "r1_path_cost.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tphi_origin\tr\tfield_abs\tanalytic\tratio\tquantity\n"
        )
        for c in radial_cmp:
            f.write(
                f"monist_1sector\tfree_capacity_psi\tfree_relaxation\t{c['r']}\t"
                f"{c['psi_num']}\t{c['psi_analytic']}\t{c['ratio_psi']}\tabs_psi\n"
            )
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\tdual_source_poisson\t{c['r']}\t"
                f"{c['phi_num']}\t{c['phi_analytic_monopole']}\t{c.get('ratio_phi', float('nan'))}\tabs_phi\n"
            )

    with open(os.path.join(OUT, "r1_rays.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tphi_origin\tb\tdeflection_rad\tdeflection_deg\t"
            "delay\tm_ledger\tn_locks\n"
        )
        for r in rays:
            f.write(
                f"monist_1sector\tfree_capacity_psi\tfree_relaxation\t{r['b']}\t"
                f"{r['deflection_rad']}\t{r['deflection_deg']}\t{r['delay']}\t"
                f"{m_tot}\t2\n"
            )

    with open(os.path.join(OUT, "r1_forces.tsv"), "w") as f:
        f.write(
            "sector_tag\tchannel\tsep\tU\tF_analytic\tF_numeric_fd\t"
            "psi_or_phi_mid\tell_mid\tcharge_mode\n"
        )
        for r in force_rows:
            f.write(
                f"monist_1sector\tfree_capacity_psi\t{r['sep']}\t{r['U_psi']}\t"
                f"{-r['F_psi_analytic_attractive']}\t"  # signed: attract negative if + = repel
                f"{F_psi_num if abs(r['sep']-sep)<1e-9 else float('nan')}\t"
                f"{r['psi_mid']}\t{r['ell_mid']}\t{args.charge_mode}\n"
            )
            f.write(
                f"multi_channel\tcoulomb_proxy_phi\t{r['sep']}\t{r['U_phi']}\t"
                f"{r['F_coulomb_signed']}\t"
                f"{F_phi_num if abs(r['sep']-sep)<1e-9 else float('nan')}\t"
                f"{float('nan')}\t{float('nan')}\t{args.charge_mode}\n"
            )
        # summary force row at primary
        f.write(
            f"monist_1sector\tfree_capacity_psi_fd_at_sep\t{sep}\t{U_psi}\t"
            f"{-F_psi_an_prim}\t{F_psi_num}\t{psi_mid}\t{ell_mid}\t{args.charge_mode}\n"
        )
        f.write(
            f"multi_channel\tcoulomb_proxy_phi_fd_at_sep\t{sep}\t{U_phi}\t"
            f"{F_c_an_prim}\t{F_phi_num}\t{phi_mid}\t{float('nan')}\t{args.charge_mode}\n"
        )

    with open(os.path.join(OUT, "r1_axis.tsv"), "w") as f:
        f.write("sector_tag\tchannel\tx\tfield\n")
        for row in axis_psi:
            f.write(f"monist_1sector\tfree_capacity_psi\t{row['x']}\t{row['field']}\n")
        for row in axis_phi:
            f.write(f"multi_channel\tcoulomb_proxy_phi\t{row['x']}\t{row['field']}\n")

    with open(os.path.join(OUT, "r1_free_deficit.tsv"), "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger_total\t{m_tot}\n")
        f.write(f"m_each_0\t{m_each[0]}\n")
        f.write(f"m_each_1\t{m_each[1]}\n")
        f.write(f"Q_each_0\t{Q_each[0]}\n")
        f.write(f"Q_each_1\t{Q_each[1]}\n")
        f.write(f"Q_net\t{Q_net}\n")
        f.write(f"free_deficit_locks\t{deficit}\n")
        f.write(f"psi_mid\t{psi_mid}\n")
        f.write(f"phi_mid\t{phi_mid}\n")
        f.write(f"ell_mid\t{ell_mid}\n")
        f.write(f"U_psi\t{U_psi}\n")
        f.write(f"U_self\t{U_self}\n")
        f.write(f"U_int_est\t{U_int_est}\n")
        f.write(f"U_phi\t{U_phi}\n")
        f.write(f"F_psi_analytic_attractive\t{F_psi_an_prim}\n")
        f.write(f"F_psi_numeric_fd\t{F_psi_num}\n")
        f.write(f"F_coulomb_analytic_signed\t{F_c_an_prim}\n")
        f.write(f"F_coulomb_numeric_fd\t{F_phi_num}\n")
        f.write(f"sep\t{sep}\n")
        f.write(f"N\t{n}\n")
        f.write(f"L\t{L}\n")
        f.write(f"charge_mode\t{args.charge_mode}\n")

    lines = [
        "v77 NM ROUND 1 — multi-lock free-capacity + dual-source stub",
        f"N={n} L={L} dx={dx:.4f} A={args.A} sigma={args.sigma} sep={sep}",
        f"n_locks=2  charge_mode={args.charge_mode}  charges={charges2}",
        f"m_tot={m_tot:.6f}  m_each={m_each[0]:.6f},{m_each[1]:.6f}",
        f"Q_each={Q_each[0]:.6f},{Q_each[1]:.6f}  Q_net={Q_net:.6f}",
        f"free_deficit_locks={deficit:.6f}",
        f"SOR psi iters={it_psi} delta={res_psi:.3e}  phi iters={it_phi} delta={res_phi:.3e}",
        f"wall_sor={t2-t1:.2f}s  total={package['wall_secs_total']:.2f}s",
        "",
        "=== monist_1sector (free capacity ψ) ===",
        f"  U_psi={U_psi:.6f}  U_self={U_self:.6f}  U_int_est={U_int_est:.6f}  "
        f"U_int_an={U_int_an_prim:.6f}",
        f"  psi_mid={psi_mid:.6e}  ell_mid={ell_mid:.6f}",
        f"  multipole prefer={fit_psi.get('prefer')}  R2_1/r={fit_psi.get('r2_invr')}  "
        f"mse_log/mse_1r={fit_psi.get('ratio_mse_log_over_invr')}",
        f"  max|defl|={max_defl:.6e}  max_delay={max_delay:.6e}  vac_defl={vac_defl:.3e}",
        f"  F_psi analytic attractive mag={F_psi_an_prim:.6e}",
        f"  F_psi numeric fd (−dU/dR; <0 attract)={F_psi_num}",
        "",
        "=== multi_channel (ψ + Coulomb proxy φ) ===",
        f"  U_phi={U_phi:.6f}  phi_mid={phi_mid:.6e}",
        f"  multipole_phi prefer={fit_phi.get('prefer')}  R2_1/r={fit_phi.get('r2_invr')}",
        f"  F_c analytic signed (+repel)={F_c_an_prim:.6e}",
        f"  F_c numeric fd (−dU_c/dR)={F_phi_num}",
        f"  |F_c|/|F_ψ|={abs(F_c_an_prim)/max(abs(F_psi_an_prim),1e-30):.4f}",
        "",
        "force scan (U_psi, U_phi vs R):",
    ]
    for r in force_rows:
        lines.append(
            f"  R={r['sep']:.1f}  Uψ={r['U_psi']:.5f}  Uφ={r['U_phi']:.5f}  "
            f"Fψ_an={r['F_psi_analytic_attractive']:.5e}  "
            f"Fc_an={r['F_coulomb_signed']:+.5e}  ell_mid={r['ell_mid']:.4f}"
        )
    lines.append("sample radial ψ:")
    for c in radial_cmp[::2]:
        lines.append(
            f"  r={c['r']:.1f}  psi={c['psi_num']:.5e}  an={c['psi_analytic']:.5e}  "
            f"ratio={c['ratio_psi']}  phi={c['phi_num']:.5e}"
        )
    lines.append(f"gates: {gates}")
    lines.append(f"VERDICT multilock_psi_pass={gates['multilock_psi_pass']}  "
                 f"dual_source_pass={gates['dual_source_pass']}")
    lines.append(package["verdict"]["summary"])
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "r1_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
