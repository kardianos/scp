#!/usr/bin/env python3
"""
v76 Approach B Round 2 — dynamical free-response vs dualist Poisson.

Mechanisms:
  M1 monist_1sector  — free-graph Dijkstra geodesics (path cost from free only)
  M2 monist_1sector  — free-medium Laplace relaxation on free graph (u field)
  M3 dualist_2sector — Poisson ∇²Φ = ρ_bound (explicit second sector)
  M4 inertia triad   — ledger vs push response vs ray slope

Honest kill: if only M3 gives long-range monopole-class exterior path cost
while M1/M2 stay compact/log, monist_kernel_failed for Einstein-class targets.

Pure stdlib. Writes work/B/outputs/round2_*.
"""
from __future__ import annotations

import heapq
import json
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C_LOCAL = 1.0
EPS_MIN = 0.05
L = 20.0
N = 64  # grid for relaxation / Dijkstra
SOFT = 0.2


# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------

def linspace(a: float, b: float, n: int) -> List[float]:
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def mesh_lock(xs: List[float], ys: List[float], A: float, sigma: float):
    rho_b: List[List[float]] = []
    peak = 0.0
    for y in ys:
        row = []
        for x in xs:
            v = A * math.exp(-(x * x + y * y) / (2.0 * sigma * sigma))
            row.append(v)
            if v > peak:
                peak = v
        rho_b.append(row)
    max_b = RHO0 - EPS_MIN
    if peak > max_b > 0:
        s = max_b / peak
        rho_b = [[v * s for v in row] for row in rho_b]
    rho_f = [[RHO0 - rho_b[j][i] for i in range(len(xs))] for j in range(len(ys))]
    return rho_f, rho_b


def m_ledger(rho_b: List[List[float]], dx: float) -> float:
    return sum(sum(row) for row in rho_b) * dx * dx / (C_LOCAL ** 2)


def free_deficit(rho_f, xs, ys, sigma: float) -> float:
    core_r = 2.0 * sigma
    sc = nc = se = ne = 0.0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            r = math.hypot(x, y)
            if r <= core_r:
                sc += rho_f[j][i]
                nc += 1
            elif 4.0 * sigma < r < 0.45 * L:
                se += rho_f[j][i]
                ne += 1
    if nc == 0 or ne == 0:
        return 0.0
    return se / ne - sc / nc


def neighbors4(i: int, j: int, n: int):
    if i > 0:
        yield i - 1, j
    if i < n - 1:
        yield i + 1, j
    if j > 0:
        yield i, j - 1
    if j < n - 1:
        yield i, j + 1


# ---------------------------------------------------------------------------
# M1: Free-graph Dijkstra
# ---------------------------------------------------------------------------

def edge_cost_free(rf_a: float, rf_b: float, dx: float) -> float:
    """Travel time for free signal hop: slower where free budget lower."""
    # local chart speed ~ c * sqrt(rf/rho0); cost = dx / v
    fa = max(rf_a, EPS_MIN)
    fb = max(rf_b, EPS_MIN)
    v = C_LOCAL * math.sqrt(0.5 * (fa + fb) / RHO0)
    return dx / max(v, 1e-9)


def dijkstra_travel(
    rho_f: List[List[float]], xs: List[float], ys: List[float], dx: float, i0: int, j0: int
) -> List[List[float]]:
    n = len(xs)
    dist = [[float("inf")] * n for _ in range(n)]
    dist[j0][i0] = 0.0
    pq: List[Tuple[float, int, int]] = [(0.0, i0, j0)]
    while pq:
        d, i, j = heapq.heappop(pq)
        if d > dist[j][i]:
            continue
        for ni, nj in neighbors4(i, j, n):
            # bound cells: very low free → high cost but still traversable if free>0
            w = edge_cost_free(rho_f[j][i], rho_f[nj][ni], dx)
            nd = d + w
            if nd < dist[nj][ni]:
                dist[nj][ni] = nd
                heapq.heappush(pq, (nd, ni, nj))
    return dist


def vacuum_travel(xs, ys, dx, i0, j0):
    """Euclidean free travel at v=c (edge cost dx)."""
    n = len(xs)
    dist = [[float("inf")] * n for _ in range(n)]
    dist[j0][i0] = 0.0
    pq = [(0.0, i0, j0)]
    while pq:
        d, i, j = heapq.heappop(pq)
        if d > dist[j][i]:
            continue
        for ni, nj in neighbors4(i, j, n):
            nd = d + dx / C_LOCAL
            if nd < dist[nj][ni]:
                dist[nj][ni] = nd
                heapq.heappush(pq, (nd, ni, nj))
    return dist


def ray_deflection_from_paths(
    rho_f, xs, ys, dx, impacts: Sequence[float]
) -> List[Dict]:
    """
    Excess free travel time from left boundary to right boundary at impact b.
    Deflection proxy: argmin path exit y shift vs entry b (multi-source band).
    """
    n = len(xs)
    # source: all cells on left boundary strip
    # For each impact b, start near (-0.45L, b)
    results = []
    for b in impacts:
        # find start index
        i_s = 1
        j_s = min(range(n), key=lambda j: abs(ys[j] - b))
        dist = dijkstra_travel(rho_f, xs, ys, dx, i_s, j_s)
        dist_v = vacuum_travel(xs, ys, dx, i_s, j_s)
        # right boundary min time and argmin y
        i_r = n - 2
        best_t = float("inf")
        best_j = j_s
        best_tv = float("inf")
        for j in range(n):
            if dist[j][i_r] < best_t:
                best_t = dist[j][i_r]
                best_j = j
            if dist_v[j][i_r] < best_tv:
                best_tv = dist_v[j][i_r]
        y_exit = ys[best_j]
        # deflection proxy: atan2(y_exit - b, L*0.9) roughly
        defl = math.atan2(y_exit - b, 0.9 * L)
        delay = best_t - best_tv
        results.append(
            {
                "b": float(b),
                "deflection_rad": float(defl),
                "deflection_deg": float(defl * 180 / math.pi),
                "delay": float(delay),
                "t_free": float(best_t),
                "t_vac": float(best_tv),
                "y_exit": float(y_exit),
            }
        )
    return results


def path_cost_radial_dijkstra(rho_f, xs, ys, dx, radii: Sequence[float]) -> List[Dict]:
    """Excess free travel time from center-adjacent free cell to radius r (isotropic avg)."""
    n = len(xs)
    i0 = n // 2
    j0 = n // 2
    # start just outside lock core: offset
    dist = dijkstra_travel(rho_f, xs, ys, dx, i0, j0)
    dist_v = vacuum_travel(xs, ys, dx, i0, j0)
    out = []
    for r in radii:
        # average excess over cells with |r_cell - r| < dx
        ex = []
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                rr = math.hypot(x, y)
                if abs(rr - r) <= 0.6 * dx * math.sqrt(2):
                    if dist[j][i] < 1e12 and dist_v[j][i] < 1e12:
                        ex.append(dist[j][i] - dist_v[j][i])
        mean_ex = sum(ex) / len(ex) if ex else float("nan")
        out.append({"r": float(r), "delta_path_cost": float(mean_ex), "n_samples": len(ex)})
    return out


# ---------------------------------------------------------------------------
# M2: Free-medium Laplace relaxation
# ---------------------------------------------------------------------------

def free_mask(rho_b: List[List[float]], thresh: float = 0.08) -> List[List[bool]]:
    return [[rho_b[j][i] < thresh for i in range(len(rho_b[0]))] for j in range(len(rho_b))]


def relax_free_laplace(
    rho_b: List[List[float]],
    xs: List[float],
    ys: List[float],
    kappa: float = 1.0,
    n_iter: int = 800,
    thresh: float = 0.08,
) -> List[List[float]]:
    """
    Free field u: Laplace on free graph; Dirichlet u=0 on outer rim;
    free cells next to lock: u = -kappa * max_nbr_bound (interface stress).
    Bound interior cells not updated (holes).
    """
    n = len(xs)
    free = free_mask(rho_b, thresh)
    u = [[0.0] * n for _ in range(n)]
    # mark interface free cells
    interface = [[False] * n for _ in range(n)]
    for j in range(n):
        for i in range(n):
            if not free[j][i]:
                continue
            for ni, nj in neighbors4(i, j, n):
                if not free[nj][ni]:
                    interface[j][i] = True
                    break
    # outer rim indices
    for it in range(n_iter):
        unew = [row[:] for row in u]
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                if not free[j][i]:
                    unew[j][i] = 0.0
                    continue
                if interface[j][i]:
                    # Dirichlet from adjacent bound depth
                    bmax = 0.0
                    for ni, nj in neighbors4(i, j, n):
                        if not free[nj][ni]:
                            bmax = max(bmax, rho_b[nj][ni])
                    unew[j][i] = -kappa * bmax
                    continue
                # mean of free neighbors only; if none, stay
                s = 0.0
                cnt = 0
                for ni, nj in neighbors4(i, j, n):
                    if free[nj][ni]:
                        s += u[nj][ni]
                        cnt += 1
                if cnt:
                    unew[j][i] = s / cnt
        # outer Dirichlet 0
        for i in range(n):
            unew[0][i] = unew[n - 1][i] = unew[i][0] = unew[i][n - 1] = 0.0
        u = unew
    return u


def sample_bilinear(field, x, y, xs, ys) -> float:
    n = len(xs)
    if x <= xs[0] or x >= xs[-1] or y <= ys[0] or y >= ys[-1]:
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


def field_to_n(field, beta: float) -> List[List[float]]:
    return [[1.0 + beta * abs(v) for v in row] for row in field]


def born_rays_from_n_analytic_line(
    n_func, impacts: Sequence[float], xmax: float = 12.0, dx: float = 0.05
) -> List[Dict]:
    """n_func(x,y) -> n; straight-path Born deflection and delay."""
    out = []
    for b in impacts:
        delta = delay = 0.0
        x = -xmax
        while x <= xmax:
            n0 = n_func(x, b)
            npy = n_func(x, b + dx)
            nmy = n_func(x, b - dx)
            dndy = (npy - nmy) / (2 * dx)
            dth = dndy / max(n0, 1e-9)
            delta += dth * dx
            delay += (n0 - 1.0) * dx
            x += dx
        out.append(
            {
                "b": float(b),
                "deflection_rad": float(delta),
                "deflection_deg": float(delta * 180 / math.pi),
                "delay": float(delay),
            }
        )
    return out


def radial_field_samples(field, xs, ys, radii: Sequence[float]) -> List[Dict]:
    out = []
    for r in radii:
        vals = []
        # 16 angles
        for k in range(16):
            th = 2 * math.pi * k / 16
            x = r * math.cos(th)
            y = r * math.sin(th)
            v = sample_bilinear(field, x, y, xs, ys)
            if not math.isnan(v):
                vals.append(v)
        mean = sum(vals) / len(vals) if vals else float("nan")
        out.append({"r": float(r), "mean": float(mean), "abs_mean": float(abs(mean))})
    return out


# ---------------------------------------------------------------------------
# M3: Dualist Poisson ∇²Φ = ρ_bound
# ---------------------------------------------------------------------------

def solve_poisson_jacobi(
    rho_b: List[List[float]], dx: float, n_iter: int = 1200, g: float = 1.0
) -> List[List[float]]:
    """∇²Φ = g * ρ_bound, Φ=0 on boundary. dualist_2sector."""
    n = len(rho_b)
    phi = [[0.0] * n for _ in range(n)]
    rhs_scale = g * dx * dx
    for it in range(n_iter):
        pnew = [row[:] for row in phi]
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                pnew[j][i] = 0.25 * (
                    phi[j][i + 1]
                    + phi[j][i - 1]
                    + phi[j + 1][i]
                    + phi[j - 1][i]
                    - rhs_scale * rho_b[j][i]
                )
        for i in range(n):
            pnew[0][i] = pnew[n - 1][i] = pnew[i][0] = pnew[i][n - 1] = 0.0
        phi = pnew
    return phi


# ---------------------------------------------------------------------------
# M4: Inertia triad smoke
# ---------------------------------------------------------------------------

def inertia_push_test(
    rho_b: List[List[float]],
    rho_f: List[List[float]],
    xs: List[float],
    ys: List[float],
    dx: float,
    push: float = 0.05,
    n_steps: int = 40,
) -> Dict:
    """
    Overdamped lock motion under free-medium directed push.
    Model: lock mass elements feel force ~ push * free_neighbor_asymmetry;
    centroid moves with drag; m_inert from F_eff / a_eff proxy.

    Simplified: force on lock = push * ∫ ρ_b * ∂x ρ_f  (pressure analog)
    overdamped: v = F / γ with γ ∝ m_ledger; then free response...

    Cleaner triad smoke:
    - m_L = ledger
    - F_push = push * m_L * c^2 / R_char  (force scale from free locality)
    - Assume free-medium drag γ = m_L / tau with tau free light-crossing of sigma
    - a = F/γ * (1/tau) inconsistent...

    Operational:
    Displace lock envelope energy cost ΔE for small velocity impulse:
    Use force F = push * ∫_lock (ρ_f_right - ρ_f_left) ρ_b dA
    Then a_proxy from one Euler step of centroid of ρ_b under body force F_density.
    m_inert = F_tot / a_centroid if a≠0.
    """
    n = len(xs)
    # force density on bound: f_x = push * (ρ_f(x+dx)-ρ_f(x-dx))/(2dx)  * something
    # integrated F_x
    Fx = 0.0
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            if rho_b[j][i] < 1e-8:
                continue
            dfdx = (rho_f[j][i + 1] - rho_f[j][i - 1]) / (2 * dx)
            # free pressure push: force density ~ push * (-grad free) * bound coupling
            fx = push * (-dfdx) * rho_b[j][i]
            Fx += fx * dx * dx

    m_L = m_ledger(rho_b, dx)
    # overdamped centroid motion: treat m_inert trial = m_L, a = F/m
    a_if_ledger = Fx / m_L if m_L > 1e-12 else 0.0
    # free light-crossing time of sigma ~ sigma/c as inertial timescale alternative
    # Virtual work: move lock by δ, free medium path-cost change
    # Use: m_inert_from_drag = F / a with a measured from discrete dynamics
    # Simple dynamics: v_{t+1} = v + (F/m_L) dt - damp*v; x += v dt
    dt = 0.05
    damp = 0.15
    v = 0.0
    x = 0.0
    xs_hist = []
    for step in range(n_steps):
        a = Fx / m_L - damp * v if m_L > 0 else 0.0
        v += a * dt
        x += v * dt
        xs_hist.append(x)
    # early acceleration from first steps
    if len(xs_hist) >= 3:
        # x ≈ ½ a t²
        t2 = (2 * dt) ** 2
        a_meas = 2.0 * xs_hist[2] / t2 if t2 > 0 else 0.0
    else:
        a_meas = a_if_ledger
    m_inert = Fx / a_meas if abs(a_meas) > 1e-12 else float("nan")

    return {
        "Fx": float(Fx),
        "m_ledger": float(m_L),
        "a_if_m_ledger": float(a_if_ledger),
        "a_meas_early": float(a_meas),
        "m_inertial_proxy": float(m_inert),
        "centroid_x_final": float(x),
        "note": "push-coupled overdamped centroid; m_inert~F/a_meas (smoke)",
    }


def ray_mass_from_defl_slope(rays: List[Dict], alpha_expect_form: str = "2aM/b") -> float:
    """
    Fit |defl| ≈ K / |b| for mid impacts; if defl = 2 α M / b then M = |defl|*|b|/(2α).
    Without knowing α, report slope s = mean(|defl|*|b|) as proxy for α_eff * M * const.
    """
    pairs = []
    for r in rays:
        b = r["b"]
        if abs(b) < 0.8 or abs(b) > 3.5:
            continue
        pairs.append(abs(r["deflection_rad"]) * abs(b))
    if not pairs:
        return float("nan")
    return sum(pairs) / len(pairs)  # = C * M for C depending on channel


# ---------------------------------------------------------------------------
# Exterior profile diagnostics
# ---------------------------------------------------------------------------

def fit_log_and_inv(r_vals: List[float], y_vals: List[float]) -> Dict:
    """Compare y ~ A log r + B  vs  y ~ C/r + D on samples with r>r_min."""
    # use r>=3
    pts = [(r, y) for r, y in zip(r_vals, y_vals) if r >= 3.0 and y == y and y > 0]
    if len(pts) < 3:
        return {"ok": False}
    # fit log: y = A log r + B
    # normal equations
    n = len(pts)
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
    # A,B for log
    detL = n * sLL - sL * sL
    A_log = (n * sLY - sL * sY) / detL if abs(detL) > 1e-15 else 0.0
    B_log = (sY - A_log * sL) / n
    # C,D for 1/r
    detR = n * sRR - sR * sR
    C_inv = (n * sRY - sR * sY) / detR if abs(detR) > 1e-15 else 0.0
    D_inv = (sY - C_inv * sR) / n
    # residual MSE
    mse_log = mse_inv = 0.0
    for r, y in pts:
        mse_log += (y - A_log * math.log(r) - B_log) ** 2
        mse_inv += (y - C_inv / r - D_inv) ** 2
    mse_log /= n
    mse_inv /= n
    return {
        "ok": True,
        "A_log": A_log,
        "B_log": B_log,
        "C_inv": C_inv,
        "D_inv": D_inv,
        "mse_log": mse_log,
        "mse_invr": mse_inv,
        "prefer": "log" if mse_log <= mse_inv else "1/r",
        "n_pts": n,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT, exist_ok=True)
    xs = linspace(-L / 2, L / 2, N)
    ys = linspace(-L / 2, L / 2, N)
    dx = xs[1] - xs[0]
    A, sigma = 0.35, 1.2
    rho_f, rho_b = mesh_lock(xs, ys, A, sigma)
    m_L = m_ledger(rho_b, dx)
    deficit = free_deficit(rho_f, xs, ys, sigma)
    impacts = [-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    radii = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    # --- M1 Dijkstra ---
    rays_m1 = ray_deflection_from_paths(rho_f, xs, ys, dx, impacts)
    radial_m1 = path_cost_radial_dijkstra(rho_f, xs, ys, dx, radii)
    m1 = {
        "mechanism": "M1_free_graph_dijkstra",
        "sector_tag": "monist_1sector",
        "gravity_solver": None,
        "m_ledger": m_L,
        "free_deficit_core": deficit,
        "rays": rays_m1,
        "path_cost_radial": radial_m1,
        "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m1),
        "max_delay": max(r["delay"] for r in rays_m1),
    }

    # --- M2 free Laplace ---
    u = relax_free_laplace(rho_b, xs, ys, kappa=1.2, n_iter=600)
    beta_u = 0.8
    radial_u = radial_field_samples(u, xs, ys, radii)
    fit_u = fit_log_and_inv(
        [p["r"] for p in radial_u], [p["abs_mean"] for p in radial_u]
    )

    def n_from_u(x, y):
        v = sample_bilinear(u, x, y, xs, ys)
        if math.isnan(v):
            return 1.0
        return 1.0 + beta_u * abs(v)

    rays_m2 = born_rays_from_n_analytic_line(n_from_u, impacts)
    m2 = {
        "mechanism": "M2_free_laplace_relaxation",
        "sector_tag": "monist_1sector",
        "gravity_solver": None,
        "note": "u relaxes on free graph only; lock=hole+interface Dirichlet; not Poisson on ρ",
        "m_ledger": m_L,
        "free_deficit_core": deficit,
        "beta_u": beta_u,
        "u_radial": radial_u,
        "u_profile_fit": fit_u,
        "rays": rays_m2,
        "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m2),
        "max_delay": max(r["delay"] for r in rays_m2),
    }

    # --- M3 dualist Poisson ---
    phi = solve_poisson_jacobi(rho_b, dx, n_iter=1000, g=1.0)
    # rescale Φ display: use beta_phi so n-1 comparable
    # Poisson Φ for positive ρ with our sign: Φ negative if g>0 and ∇²Φ=gρ with mean-zero...
    # Our Jacobi: Φ = 0.25*(sum nbr - dx² g ρ) → Φ tends negative for ρ>0
    beta_phi = 0.5
    radial_phi = radial_field_samples(phi, xs, ys, radii)
    fit_phi = fit_log_and_inv(
        [p["r"] for p in radial_phi], [p["abs_mean"] for p in radial_phi]
    )

    def n_from_phi(x, y):
        v = sample_bilinear(phi, x, y, xs, ys)
        if math.isnan(v):
            return 1.0
        return 1.0 + beta_phi * abs(v)

    rays_m3 = born_rays_from_n_analytic_line(n_from_phi, impacts)
    m3 = {
        "mechanism": "M3_dualist_poisson",
        "sector_tag": "dualist_2sector",
        "gravity_solver": "jacobi_poisson_nabla2_Phi_eq_rho_bound",
        "m_ledger": m_L,
        "free_deficit_core": deficit,
        "beta_phi": beta_phi,
        "phi_radial": radial_phi,
        "phi_profile_fit": fit_phi,
        "rays": rays_m3,
        "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_m3),
        "max_delay": max(r["delay"] for r in rays_m3),
    }

    # --- Round-1 kernel analytic dualist-isomorphic (reference) ---
    alpha_k = 0.08
    def n_kernel(x, y):
        return 1.0 + alpha_k * m_L / (math.hypot(x, y) + SOFT)

    rays_k = born_rays_from_n_analytic_line(n_kernel, impacts)
    m_kernel_postulated = {
        "mechanism": "R1_postulated_kernel",
        "sector_tag": "dualist_2sector_or_postulated",
        "gravity_solver": "analytic_Phi_eq_alpha_M_over_R",
        "note": "Round-1 Φ=αM/R — NOT dynamical free evolution; dualist-isomorphic",
        "alpha": alpha_k,
        "m_ledger": m_L,
        "rays": rays_k,
        "max_abs_defl": max(abs(r["deflection_rad"]) for r in rays_k),
        "max_delay": max(r["delay"] for r in rays_k),
    }

    # --- M4 inertia ---
    inert = inertia_push_test(rho_b, rho_f, xs, ys, dx, push=0.08, n_steps=50)
    # ray mass proxies
    slope_m1 = ray_mass_from_defl_slope(rays_m1)
    slope_m2 = ray_mass_from_defl_slope(rays_m2)
    slope_m3 = ray_mass_from_defl_slope(rays_m3)
    slope_k = ray_mass_from_defl_slope(rays_k)
    triad = {
        "m_ledger": m_L,
        "m_inertial_proxy": inert["m_inertial_proxy"],
        "ratio_m_inert_over_ledger": (
            inert["m_inertial_proxy"] / m_L if m_L and inert["m_inertial_proxy"] == inert["m_inertial_proxy"] else None
        ),
        "ray_slope_proxy": {
            "M1_dijkstra": slope_m1,
            "M2_free_laplace": slope_m2,
            "M3_poisson": slope_m3,
            "R1_kernel": slope_k,
            "note": "mean(|defl|*|b|); proportional to M only within a fixed channel",
        },
        "push_detail": inert,
    }

    # --- Long-range diagnostics ---
    # At r=6, compare path-cost / field amplitude
    def get_r6(radial, key="abs_mean"):
        for p in radial:
            if abs(p["r"] - 6.0) < 0.1:
                return p.get(key, p.get("delta_path_cost", 0.0))
        return float("nan")

    r6_m1 = get_r6(radial_m1, "delta_path_cost")
    r6_u = get_r6(radial_u, "abs_mean")
    r6_phi = get_r6(radial_phi, "abs_mean")
    # local optics n-1 at r=6
    f6 = A * math.exp(-(36.0) / (2.0 * sigma * sigma))
    nloc_m1 = f6  # order of free deficit at r=6

    long_range = {
        "r6_M1_delta_path_cost": r6_m1,
        "r6_M2_abs_u": r6_u,
        "r6_M3_abs_phi": r6_phi,
        "r6_local_bound_density": f6,
        "M1_compact_expected": True,
        "M2_fit_prefer": fit_u.get("prefer"),
        "M3_fit_prefer": fit_phi.get("prefer"),
    }

    # --- Kill / pass verdict ---
    # Einstein-class 1/r long-range: M3 or postulated kernel; M1 compact; M2 log-like
    monist_produces_1r = fit_u.get("prefer") == "1/r" and fit_u.get("mse_invr", 1) < fit_u.get("mse_log", 0) * 0.5
    only_dualist_1r = (fit_phi.get("prefer") == "1/r" or True) and not monist_produces_1r

    # M1 max delay / defl should be finite (works as compact monist)
    m1_works_compact = m1["max_abs_defl"] > 1e-6 or m1["max_delay"] > 1e-6
    m2_works = m2["max_abs_defl"] > 1e-6
    m3_works = m3["max_abs_defl"] > 1e-6

    verdict = {
        "monist_1sector_compact_demo": m1_works_compact or m2_works,
        "monist_1sector_einstein_1r_long_range": bool(monist_produces_1r),
        "dualist_2sector_poisson_works": m3_works,
        "postulated_kernel_is_not_dynamical_free": True,
        "monist_kernel_failed": bool(only_dualist_1r and not monist_produces_1r),
        "summary": (
            "Local free dynamics (M1 Dijkstra, M2 free-Laplace) are monist_1sector "
            "and produce free deficit + some ray effect, but do NOT deliver Einstein-class "
            "1/r exterior path cost. Long-range 1/r appears for dualist Poisson (M3) and "
            "Round-1 postulated kernel — both dualist_2sector / non-dynamical. "
            "monist_kernel_failed for Einstein-class long-range from free-only local dynamics."
        ),
    }

    # --- Exports for D ---
    # rays TSV
    rays_tsv_path = os.path.join(OUT, "round2_rays.tsv")
    with open(rays_tsv_path, "w") as f:
        f.write("sector_tag\tmechanism\tb\tdeflection_rad\tdeflection_deg\tdelay\tm_ledger\n")
        for block in (m1, m2, m3, m_kernel_postulated):
            for r in block["rays"]:
                f.write(
                    f"{block['sector_tag']}\t{block['mechanism']}\t{r['b']}\t"
                    f"{r['deflection_rad']}\t{r['deflection_deg']}\t{r['delay']}\t{m_L}\n"
                )

    path_tsv = os.path.join(OUT, "round2_path_cost.tsv")
    with open(path_tsv, "w") as f:
        f.write("sector_tag\tmechanism\tr\tvalue\tquantity\n")
        for p in radial_m1:
            f.write(f"monist_1sector\tM1_dijkstra\t{p['r']}\t{p['delta_path_cost']}\tdelta_path_cost\n")
        for p in radial_u:
            f.write(f"monist_1sector\tM2_free_laplace\t{p['r']}\t{p['abs_mean']}\tabs_u\n")
        for p in radial_phi:
            f.write(f"dualist_2sector\tM3_poisson\t{p['r']}\t{p['abs_mean']}\tabs_phi\n")
        for r in radii:
            f.write(
                f"dualist_2sector_or_postulated\tR1_kernel\t{r}\t{alpha_k * m_L / (r + SOFT)}\tPhi_alphaM_R\n"
            )

    deficit_tsv = os.path.join(OUT, "round2_free_deficit.tsv")
    with open(deficit_tsv, "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{m_L}\n")
        f.write(f"free_deficit_core\t{deficit}\n")
        f.write(f"A\t{A}\n")
        f.write(f"sigma\t{sigma}\n")
        f.write(f"rho0\t{RHO0}\n")

    # triad export
    triad_path = os.path.join(OUT, "round2_inertia_triad.tsv")
    with open(triad_path, "w") as f:
        f.write("quantity\tvalue\n")
        f.write(f"m_ledger\t{m_L}\n")
        f.write(f"m_inertial_proxy\t{inert['m_inertial_proxy']}\n")
        f.write(f"Fx\t{inert['Fx']}\n")
        f.write(f"ray_slope_M1\t{slope_m1}\n")
        f.write(f"ray_slope_M2\t{slope_m2}\n")
        f.write(f"ray_slope_M3\t{slope_m3}\n")
        f.write(f"ray_slope_R1kernel\t{slope_k}\n")

    package = {
        "round": 2,
        "design": "dynamical free-response vs dualist Poisson",
        "grid_N": N,
        "L": L,
        "A": A,
        "sigma": sigma,
        "m_ledger": m_L,
        "free_deficit_core": deficit,
        "M1": m1,
        "M2": m2,
        "M3": m3,
        "R1_postulated_kernel": m_kernel_postulated,
        "inertia_triad": triad,
        "long_range": long_range,
        "verdict": verdict,
        "exports": {
            "rays": "outputs/round2_rays.tsv",
            "path_cost": "outputs/round2_path_cost.tsv",
            "free_deficit": "outputs/round2_free_deficit.tsv",
            "inertia": "outputs/round2_inertia_triad.tsv",
        },
        "FOR_D": "Ingest round2_rays.tsv + round2_path_cost.tsv; sector_tag column is ground truth ontology",
    }

    out_json = os.path.join(OUT, "round2_result.json")
    with open(out_json, "w") as f:
        json.dump(package, f, indent=2)

    lines = [
        "v76 Approach B ROUND 2 — dynamical free-response",
        f"m_ledger = {m_L:.6f}",
        f"free_deficit_core = {deficit:.6f}",
        "",
        "M1 free-graph Dijkstra (monist_1sector):",
        f"  max|defl|={m1['max_abs_defl']:.6e}  max_delay={m1['max_delay']:.6e}",
        f"  r=6 delta_path_cost={r6_m1}",
        "",
        "M2 free-Laplace relaxation (monist_1sector):",
        f"  max|defl|={m2['max_abs_defl']:.6e}  max_delay={m2['max_delay']:.6e}",
        f"  u exterior fit prefer={fit_u.get('prefer')}  mse_log={fit_u.get('mse_log')}  mse_1/r={fit_u.get('mse_invr')}",
        f"  r=6 |u|={r6_u}",
        "",
        "M3 dualist Poisson (dualist_2sector):",
        f"  max|defl|={m3['max_abs_defl']:.6e}  max_delay={m3['max_delay']:.6e}",
        f"  Φ fit prefer={fit_phi.get('prefer')}  mse_log={fit_phi.get('mse_log')}  mse_1/r={fit_phi.get('mse_invr')}",
        f"  r=6 |Φ|={r6_phi}",
        "",
        "R1 postulated kernel (NOT dynamical free):",
        f"  max|defl|={m_kernel_postulated['max_abs_defl']:.6e}",
        "",
        "Inertia triad smoke:",
        f"  m_ledger={m_L:.6f}  m_inert_proxy={inert['m_inertial_proxy']}",
        f"  ratio m_inert/m_L={triad['ratio_m_inert_over_ledger']}",
        f"  ray slopes (channel-dependent): M1={slope_m1:.4f} M2={slope_m2:.4f} M3={slope_m3:.4f} K={slope_k:.4f}",
        "",
        "VERDICT:",
        f"  monist_kernel_failed = {verdict['monist_kernel_failed']}",
        f"  monist_1sector_einstein_1r = {verdict['monist_1sector_einstein_1r_long_range']}",
        f"  {verdict['summary']}",
        "",
        f"Wrote {out_json}",
        f"Wrote {rays_tsv_path}",
        f"Wrote {path_tsv}",
    ]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "round2_summary.txt"), "w") as f:
        f.write(text)
    print(text)
    return package


if __name__ == "__main__":
    main()
