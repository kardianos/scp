#!/usr/bin/env python3
"""
v76 Approach B — B2-lite optical monism sandbox (Round 1).

One continuum budget: rho_free + rho_bound = rho0.
Lock (hand-placed Gaussian) depletes free density.
Rays use refractive index n = rho0 / rho_free — NO Poisson/Einstein solver.

Outputs under work/B/outputs/
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, asdict
from typing import List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "outputs")

RHO0 = 1.0
C_LOCAL = 1.0  # free-field locality in code units
EPS_MIN = 0.05  # floor on rho_free (horizon analog deferred)
NX = NY = 256
LBOX = 20.0  # domain half-width ~ [-LBOX/2, LBOX/2] each axis; full width LBOX
# use full width L = LBOX
L = 20.0


@dataclass
class MediumStats:
    m_bound: float
    rho_free_core: float
    rho_free_exterior: float
    free_deficit_core: float  # exterior - core (positive => depleted at lock)
    budget_residual_max: float
    A_lock: float
    sigma_lock: float


# ---------------------------------------------------------------------------
# Medium
# ---------------------------------------------------------------------------

def make_grid(n: int = NX, L: float = L) -> Tuple[np.ndarray, np.ndarray, float]:
    x = np.linspace(-L / 2, L / 2, n)
    y = np.linspace(-L / 2, L / 2, n)
    dx = x[1] - x[0]
    return x, y, dx


def place_gaussian_lock(
    X: np.ndarray,
    Y: np.ndarray,
    A: float,
    sigma: float,
    cx: float = 0.0,
    cy: float = 0.0,
    rho0: float = RHO0,
    eps_min: float = EPS_MIN,
) -> Tuple[np.ndarray, np.ndarray]:
    r2 = (X - cx) ** 2 + (Y - cy) ** 2
    rho_b = A * np.exp(-r2 / (2.0 * sigma**2))
    # enforce floor: rho_b <= rho0 - eps_min
    max_b = rho0 - eps_min
    if rho_b.max() > max_b:
        rho_b *= max_b / rho_b.max()
    rho_f = rho0 - rho_b
    return rho_f, rho_b


def refractive_index(rho_f: np.ndarray, rho0: float = RHO0) -> np.ndarray:
    return rho0 / np.maximum(rho_f, 1e-15)


def medium_stats(
    rho_f: np.ndarray,
    rho_b: np.ndarray,
    dx: float,
    A: float,
    sigma: float,
    rho0: float = RHO0,
    core_radius: float = None,
) -> MediumStats:
    if core_radius is None:
        core_radius = 2.0 * sigma
    n = rho_f.shape[0]
    L = dx * (n - 1)
    x = np.linspace(-L / 2, L / 2, n)
    X, Y = np.meshgrid(x, x, indexing="xy")
    r = np.sqrt(X**2 + Y**2)
    core = r <= core_radius
    ext = (r > 4.0 * sigma) & (r < 0.45 * L)
    m_bound = float(np.sum(rho_b) * dx * dx / (C_LOCAL**2))
    rf_core = float(rho_f[core].mean()) if core.any() else float("nan")
    rf_ext = float(rho_f[ext].mean()) if ext.any() else float(rho0)
    residual = float(np.max(np.abs(rho_f + rho_b - rho0)))
    return MediumStats(
        m_bound=m_bound,
        rho_free_core=rf_core,
        rho_free_exterior=rf_ext,
        free_deficit_core=rf_ext - rf_core,
        budget_residual_max=residual,
        A_lock=A,
        sigma_lock=sigma,
    )


# ---------------------------------------------------------------------------
# Bilinear sample of n and grad n
# ---------------------------------------------------------------------------

def sample_field(field: np.ndarray, x: float, y: float, L: float, n: int) -> float:
    """Bilinear sample; field shape (n,n) with y index first (row)."""
    # map to index: x,y in [-L/2, L/2]
    fx = (x + L / 2) / L * (n - 1)
    fy = (y + L / 2) / L * (n - 1)
    if fx < 0 or fy < 0 or fx > n - 1 or fy > n - 1:
        return float("nan")
    i0 = int(math.floor(fx))
    j0 = int(math.floor(fy))
    i1 = min(i0 + 1, n - 1)
    j1 = min(j0 + 1, n - 1)
    tx = fx - i0
    ty = fy - j0
    # field[j, i] with j~y, i~x
    v00 = field[j0, i0]
    v10 = field[j0, i1]
    v01 = field[j1, i0]
    v11 = field[j1, i1]
    return float(
        (1 - ty) * ((1 - tx) * v00 + tx * v10)
        + ty * ((1 - tx) * v01 + tx * v11)
    )


def sample_grad_n(
    n_field: np.ndarray, x: float, y: float, L: float, ngrid: int, dx: float
) -> Tuple[float, float]:
    """Central differences of n on the grid, bilinear-sampled at (x,y)."""
    # build grad arrays once would be better; for simplicity sample neighbors
    h = dx
    npx = sample_field(n_field, x + h, y, L, ngrid)
    nmx = sample_field(n_field, x - h, y, L, ngrid)
    npy = sample_field(n_field, x, y + h, L, ngrid)
    nmy = sample_field(n_field, x, y - h, L, ngrid)
    if any(math.isnan(v) for v in (npx, nmx, npy, nmy)):
        return 0.0, 0.0
    return (npx - nmx) / (2 * h), (npy - nmy) / (2 * h)


# ---------------------------------------------------------------------------
# Ray integration (eikonal / geometric optics)
# State: (x, y, theta) with chart speed c/n; dtheta from grad n.
#
# Hamiltonian form: p is direction * n (optical), H = |p| - n = 0
# dx/ds_eucl = p_hat, d p / ds = grad n  (s = Euclidean arc length)
# With angle: d theta / ds = (1/n) * ( -sin theta * dn/dx + cos theta * dn/dy )
# Actually: for ds optical S = int n ds, ray curvature:
#   dT/ds = (grad n - (grad n · T) T) / n   with T unit tangent in chart
# ---------------------------------------------------------------------------

def integrate_ray(
    n_field: np.ndarray,
    x0: float,
    y0: float,
    theta0: float,
    L: float,
    ngrid: int,
    dx: float,
    s_max: float = 40.0,
    ds: float = 0.02,
) -> Tuple[np.ndarray, float, float]:
    """
    Integrate ray. Returns path (N,2), total optical time (int n ds / c),
    and final angle theta.
    """
    x, y, th = x0, y0, theta0
    path = [(x, y)]
    t_opt = 0.0
    half = 0.48 * L
    steps = int(s_max / ds)
    for _ in range(steps):
        nval = sample_field(n_field, x, y, L, ngrid)
        if math.isnan(nval) or abs(x) > half or abs(y) > half:
            break
        dndx, dndy = sample_grad_n(n_field, x, y, L, ngrid, dx)
        # unit tangent
        tx, ty = math.cos(th), math.sin(th)
        # curvature: dT/ds = (grad_perp n) / n
        gn_par = dndx * tx + dndy * ty
        # d theta / ds from rotating T toward grad n
        # dtx/ds = (dndx - gn_par*tx)/n ;  dty similar
        # dth/ds = -sin*dtx + cos*dty  (since dtx = -sin dth, dty = cos dth)
        dtx = (dndx - gn_par * tx) / nval
        dty = (dndy - gn_par * ty) / nval
        dth = -ty * dtx + tx * dty  # = ( -ty*dndx + tx*dndy ) / n
        # RK2
        th_mid = th + 0.5 * ds * dth
        x_mid = x + 0.5 * ds * math.cos(th)
        y_mid = y + 0.5 * ds * math.sin(th)
        n_mid = sample_field(n_field, x_mid, y_mid, L, ngrid)
        if math.isnan(n_mid):
            break
        dndx_m, dndy_m = sample_grad_n(n_field, x_mid, y_mid, L, ngrid, dx)
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
    return np.array(path), t_opt, th


def deflection_and_delay(
    n_field: np.ndarray,
    L: float,
    ngrid: int,
    dx: float,
    impact_params: np.ndarray,
    x_start: float = None,
) -> List[dict]:
    if x_start is None:
        x_start = -0.45 * L
    results = []
    for b in impact_params:
        # vacuum reference: straight line, n=1
        path_v, t_v, th_v = integrate_ray(
            np.ones_like(n_field), x_start, b, 0.0, L, ngrid, dx
        )
        path, t_opt, th = integrate_ray(
            n_field, x_start, b, 0.0, L, ngrid, dx
        )
        # asymptotic deflection: final angle (incoming theta=0)
        defl = th  # radians; positive = counterclockwise
        # excess delay vs vacuum path with same s if both exited
        delay = t_opt - t_v
        # exit y for rough bend measure
        y_exit = float(path[-1, 1]) if len(path) else float("nan")
        results.append(
            {
                "b": float(b),
                "deflection_rad": float(defl),
                "deflection_deg": float(defl * 180.0 / math.pi),
                "t_optical": float(t_opt),
                "t_vacuum": float(t_v),
                "delay": float(delay),
                "y_exit": y_exit,
                "n_steps": int(len(path)),
            }
        )
    return results


# ---------------------------------------------------------------------------
# Optional: discrete lock-formation step (smoke test)
# ---------------------------------------------------------------------------

def formation_step(
    rho_f: np.ndarray,
    rho_b: np.ndarray,
    gamma: float = 0.12,
    mu: float = 0.02,
    rho0: float = RHO0,
    eps_min: float = EPS_MIN,
) -> Tuple[np.ndarray, np.ndarray]:
    """Autocatalytic lock growth preserving budget: dB = gamma*B*F - mu*B.

    Grows only where bound already present (seed). Exterior B=0 stays free.
    """
    d_bound = gamma * rho_b * rho_f - mu * rho_b
    rho_b2 = np.clip(rho_b + d_bound, 0.0, rho0 - eps_min)
    rho_f2 = rho0 - rho_b2
    return rho_f2, rho_b2


# ---------------------------------------------------------------------------
# Run suite
# ---------------------------------------------------------------------------

def run_case(A: float, sigma: float, tag: str, impact_params: np.ndarray) -> dict:
    x, y, dx = make_grid()
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho_f, rho_b = place_gaussian_lock(X, Y, A=A, sigma=sigma)
    n_field = refractive_index(rho_f)
    stats = medium_stats(rho_f, rho_b, dx, A=A, sigma=sigma)
    rays = deflection_and_delay(n_field, L, NX, dx, impact_params)
    # max |deflection|
    max_defl = max(abs(r["deflection_rad"]) for r in rays) if rays else 0.0
    max_delay = max(r["delay"] for r in rays) if rays else 0.0
    return {
        "tag": tag,
        "A": A,
        "sigma": sigma,
        "stats": asdict(stats),
        "rays": rays,
        "max_abs_deflection_rad": max_defl,
        "max_delay": max_delay,
        "rho_free": rho_f,
        "rho_bound": rho_b,
        "n_field": n_field,
        "grid_x": x,
        "grid_y": y,
        "dx": dx,
    }


def save_plots(cases: List[dict]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots")
        return

    os.makedirs(OUT, exist_ok=True)

    # pick lock case for field plots
    lock = next(c for c in cases if c["A"] > 0)
    vac = next(c for c in cases if c["A"] == 0)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    extent = [-L / 2, L / 2, -L / 2, L / 2]
    im0 = axes[0].imshow(
        lock["rho_bound"], origin="lower", extent=extent, cmap="magma"
    )
    axes[0].set_title(r"$\rho_{\mathrm{bound}}$")
    fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(
        lock["rho_free"], origin="lower", extent=extent, cmap="viridis"
    )
    axes[1].set_title(r"$\rho_{\mathrm{free}}$ (deficit at lock)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    im2 = axes[2].imshow(
        lock["n_field"], origin="lower", extent=extent, cmap="cividis"
    )
    axes[2].set_title(r"$n=\rho_0/\rho_{\mathrm{free}}$")
    fig.colorbar(im2, ax=axes[2], fraction=0.046)
    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fields_lock.png"), dpi=140)
    plt.close(fig)

    # ray paths on n background
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(lock["n_field"], origin="lower", extent=extent, cmap="cividis", alpha=0.9)
    xg, yg, dx = make_grid()
    n_field = lock["n_field"]
    for b in [-3.0, -1.5, -0.5, 0.5, 1.5, 3.0]:
        path, _, _ = integrate_ray(n_field, -0.45 * L, b, 0.0, L, NX, dx)
        ax.plot(path[:, 0], path[:, 1], "w-", lw=0.9, alpha=0.9)
    ax.set_title("rays through free-density medium (no gravity solver)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-L / 2, L / 2)
    ax.set_ylim(-L / 2, L / 2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "rays_lock.png"), dpi=140)
    plt.close(fig)

    # deflection & delay vs b
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    for c, style in ((lock, "o-"), (vac, "s--")):
        bs = [r["b"] for r in c["rays"]]
        de = [r["deflection_deg"] for r in c["rays"]]
        dl = [r["delay"] for r in c["rays"]]
        axes[0].plot(bs, de, style, label=c["tag"])
        axes[1].plot(bs, dl, style, label=c["tag"])
    axes[0].set_xlabel("impact parameter b")
    axes[0].set_ylabel("deflection (deg)")
    axes[0].legend()
    axes[0].set_title("deflection vs b")
    axes[1].set_xlabel("impact parameter b")
    axes[1].set_ylabel("excess optical delay")
    axes[1].legend()
    axes[1].set_title("Shapiro-like delay vs b")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "deflection_delay.png"), dpi=140)
    plt.close(fig)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    impact = np.array([-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])

    cases_meta = []
    cases_full = []

    # vacuum control
    vac = run_case(A=0.0, sigma=1.0, tag="vacuum", impact_params=impact)
    cases_full.append(vac)
    cases_meta.append({k: v for k, v in vac.items() if k not in ("rho_free", "rho_bound", "n_field", "grid_x", "grid_y")})

    # lock (moderate A so eikonal rays do not trap)
    lock = run_case(A=0.35, sigma=1.2, tag="lock_A0.35_s1.2", impact_params=impact)
    cases_full.append(lock)
    cases_meta.append({k: v for k, v in lock.items() if k not in ("rho_free", "rho_bound", "n_field", "grid_x", "grid_y")})

    # weaker lock
    weak = run_case(A=0.15, sigma=1.2, tag="lock_A0.15_s1.2", impact_params=impact)
    cases_full.append(weak)
    cases_meta.append({k: v for k, v in weak.items() if k not in ("rho_free", "rho_bound", "n_field", "grid_x", "grid_y")})

    # formation smoke: seed small noise above threshold near center
    x, y, dx = make_grid()
    X, Y = np.meshgrid(x, y, indexing="xy")
    rho_f = np.full_like(X, RHO0)
    rho_b = np.zeros_like(X)
    # seed a free overdensity (by temporarily reducing bound elsewhere is wrong);
    # instead start with small bound seed and run formation from free above rho_star
    # Use a free "bump" by allowing temporary non-uniform total? No — stay monist:
    # start near-uniform free, seed tiny bound Gaussian, run exchange with high gamma
    # and lowered rho_star so free can lock.
    seed_A, seed_s = 0.12, 0.9
    rho_f, rho_b = place_gaussian_lock(X, Y, A=seed_A, sigma=seed_s)
    hist = []
    for step in range(60):
        rho_f, rho_b = formation_step(rho_f, rho_b, gamma=0.12, mu=0.02)
        if step % 15 == 0 or step == 59:
            st = medium_stats(rho_f, rho_b, dx, A=seed_A, sigma=seed_s)
            hist.append(
                {
                    "step": step,
                    "m_bound": st.m_bound,
                    "free_deficit_core": st.free_deficit_core,
                    "rho_b_max": float(rho_b.max()),
                }
            )

    form_stats = medium_stats(rho_f, rho_b, dx, A=float(rho_b.max()), sigma=seed_s)
    n_form = refractive_index(rho_f)
    form_rays = deflection_and_delay(n_form, L, NX, dx, impact)
    formation = {
        "tag": "formation_smoke",
        "history": hist,
        "stats": asdict(form_stats),
        "max_abs_deflection_rad": max(abs(r["deflection_rad"]) for r in form_rays),
        "max_delay": max(r["delay"] for r in form_rays),
        "rays_sample": form_rays[::2],
    }

    # Pass/fail gates
    lock_s = lock["stats"]
    vac_s = vac["stats"]
    gates = {
        "budget_residual_ok": lock_s["budget_residual_max"] < 1e-12,
        "free_deficit_positive": lock_s["free_deficit_core"] > 0.01,
        "lock_deflection_nonzero": lock["max_abs_deflection_rad"] > 1e-4,
        "lock_delay_positive": lock["max_delay"] > 1e-4,
        "vacuum_deflection_near_zero": vac["max_abs_deflection_rad"] < 1e-3,
        "vacuum_delay_near_zero": abs(vac["max_delay"]) < 1e-3,
        "weaker_smaller_defl": weak["max_abs_deflection_rad"] < lock["max_abs_deflection_rad"],
    }
    gates["round1_pass"] = all(gates.values())

    results = {
        "design": "B2-lite optical monism",
        "rho0": RHO0,
        "c_local": C_LOCAL,
        "eps_min": EPS_MIN,
        "n_index_law": "n = rho0 / rho_free",
        "gravity_solver": None,
        "cases": cases_meta,
        "formation_smoke": formation,
        "gates": gates,
    }

    out_json = os.path.join(OUT, "results.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

    # also a short human summary
    summary_path = os.path.join(OUT, "summary.txt")
    lines = [
        "v76 B2-lite sandbox results",
        f"round1_pass = {gates['round1_pass']}",
        f"lock m_bound = {lock_s['m_bound']:.6f}",
        f"lock free_deficit_core = {lock_s['free_deficit_core']:.6f}",
        f"lock max |defl| rad = {lock['max_abs_deflection_rad']:.6e}",
        f"lock max delay = {lock['max_delay']:.6e}",
        f"vac max |defl| rad = {vac['max_abs_deflection_rad']:.6e}",
        f"budget residual max = {lock_s['budget_residual_max']:.3e}",
        f"formation final m = {form_stats.m_bound:.6f}, deficit = {form_stats.free_deficit_core:.6f}",
        "gates:",
    ]
    for k, v in gates.items():
        lines.append(f"  {k}: {v}")
    # sample rays
    lines.append("sample lock rays (b, defl_deg, delay):")
    for r in lock["rays"][::2]:
        lines.append(
            f"  b={r['b']:+.1f}  defl={r['deflection_deg']:+.4f} deg  delay={r['delay']:+.5f}"
        )
    text = "\n".join(lines) + "\n"
    with open(summary_path, "w") as f:
        f.write(text)
    print(text)

    save_plots(cases_full)
    print(f"Wrote {out_json}")
    print(f"Wrote {summary_path}")
    print(f"Plots under {OUT}/")


if __name__ == "__main__":
    main()
