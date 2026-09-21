#!/usr/bin/env python3
"""v100 phase 0 — hedgehog minimizer of the locked functional.

The integrand is PREREG.md §1. Thresholds live there. This script
prints meter lines and writes the fine profile. It does not decide G2.

Run: python3 v100/phase0/cas/profile.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

TWOPI = 2.0 * math.pi
# Published chiral-limit integral of this functional, in units f_π/e.
TEXTBOOK_E = 72.9


def grid(n_intervals: int, xmax: float) -> np.ndarray:
    return np.linspace(0.0, xmax, n_intervals + 1)


def parts_and_grad_F(F: np.ndarray, x: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Midpoint E2, E4, and ∂E/∂F at the nodes. Endpoints included."""
    dx = np.diff(x)
    Fm = 0.5 * (F[:-1] + F[1:])
    Fp = np.diff(F) / dx
    xm = 0.5 * (x[:-1] + x[1:])
    s = np.sin(Fm)
    c = np.cos(Fm)
    s2 = s * s
    Fp2 = Fp * Fp
    xm2 = xm * xm
    e2 = xm2 * Fp2 + 2.0 * s2
    e4 = s2 * (2.0 * Fp2 + s2 / xm2)
    E2 = TWOPI * float(np.sum(dx * e2))
    E4 = TWOPI * float(np.sum(dx * e4))

    # ∂(e2+e4)/∂Fm , ∂(e2+e4)/∂Fp
    sin2 = 2.0 * s * c
    d_Fm = 2.0 * sin2 * (1.0 + Fp2 + s2 / xm2)
    d_Fp = 2.0 * Fp * (xm2 + 2.0 * s2)
    g_Fm = TWOPI * dx * d_Fm
    g_Fp = TWOPI * dx * d_Fp

    gF = np.zeros_like(F)
    gF[:-1] += 0.5 * g_Fm - g_Fp / dx
    gF[1:] += 0.5 * g_Fm + g_Fp / dx
    return E2, E4, gF


def softmax_pack(y: np.ndarray) -> np.ndarray:
    z = y - np.max(y)
    e = np.exp(z)
    return math.pi * e / np.sum(e)


def F_of_y(y: np.ndarray) -> np.ndarray:
    delta = softmax_pack(y)
    F = np.empty(y.size + 1)
    F[0] = math.pi
    # sum(delta) = π, so the last node is 0 in exact arithmetic.
    F[1:] = math.pi - np.cumsum(delta)
    return F


def grad_y(gF: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Chain ∂E/∂y through F_j = π - Σ_{k<j} Δ_k, Δ = π softmax(y)."""
    delta = softmax_pack(y)
    # ∂E/∂Δ_k = - Σ_{j>k} ∂E/∂F_j   (F_0 fixed; F[-1] pinned, drop it)
    tail = np.cumsum(gF[:0:-1])[::-1]
    # tail[k] = gF[k+1] + gF[k+2] + ... + gF[M-1]  (exclude pinned F[M])
    # F[j] depends on Δ_k for k < j, j=1..M-1. Pinned F[M] has no influence.
    dE_dDelta = -tail
    # ∂Δ_k/∂y_m = Δ_k (δ_km - Δ_m/π)
    # ∂E/∂y_m = Σ_k dE_dΔ_k * Δ_k * (δ_km - Δ_m/π)
    w = dE_dDelta * delta
    return w - delta * (np.sum(w) / math.pi)


def y_from_seed(F_seed: np.ndarray) -> np.ndarray:
    d = np.diff(F_seed)
    # Profile falls. A numerical uptick would break the ansatz.
    d = np.maximum(-d, 1e-15)
    d *= math.pi / np.sum(d)
    return np.log(d)


def seed_A(x: np.ndarray) -> np.ndarray:
    return 4.0 * np.arctan(np.exp(-x))


def seed_B(x: np.ndarray) -> np.ndarray:
    return 4.0 * np.arctan(np.exp(-0.5 * x))


def _interval_node_grads(FL, FR, dx, xm):
    """∂E/∂F_left and ∂E/∂F_right of one midpoint interval, vectorized."""
    Fm = 0.5 * (FL + FR)
    Fp = (FR - FL) / dx
    s = np.sin(Fm)
    c = np.cos(Fm)
    s2 = s * s
    Fp2 = Fp * Fp
    xm2 = xm * xm
    sin2 = 2.0 * s * c
    d_Fm = 2.0 * sin2 * (1.0 + Fp2 + s2 / xm2)
    d_Fp = 2.0 * Fp * (xm2 + 2.0 * s2)
    g_Fm = TWOPI * dx * d_Fm
    g_Fp = TWOPI * dx * d_Fp
    gL = 0.5 * g_Fm - g_Fp / dx
    gR = 0.5 * g_Fm + g_Fp / dx
    return gL, gR


def _hessian_tridiagonal(F: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Free-node Hessian of the midpoint energy. Boundary nodes are fixed."""
    dx = np.diff(x)
    xm = 0.5 * (x[:-1] + x[1:])
    FL = F[:-1]
    FR = F[1:]
    eps = 1e-7
    gL_p, gR_p = _interval_node_grads(FL + eps, FR, dx, xm)
    gL_m, gR_m = _interval_node_grads(FL - eps, FR, dx, xm)
    HLL = (gL_p - gL_m) / (2 * eps)
    HLR = (gR_p - gR_m) / (2 * eps)
    gL_p, gR_p = _interval_node_grads(FL, FR + eps, dx, xm)
    gL_m, gR_m = _interval_node_grads(FL, FR - eps, dx, xm)
    HRR = (gR_p - gR_m) / (2 * eps)
    HLR = 0.5 * (HLR + (gL_p - gL_m) / (2 * eps))

    m = len(F) - 2
    diag = np.zeros(m)
    off = np.zeros(max(m - 1, 0))
    # Interval i couples nodes (i, i+1). Free unknowns are nodes 1..m.
    diag += HRR[:m]
    diag += HLL[1 : m + 1]
    if m > 1:
        off += HLR[1:m]
    return diag, off


def _solve_tridiagonal(diag: np.ndarray, off: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve tridiagonal(diag, off) δ = rhs. off[i] couples i and i+1."""
    n = diag.size
    a = diag.copy()
    b = off.copy()
    c = off.copy()
    d = rhs.copy()
    for i in range(1, n):
        w = b[i - 1] / a[i - 1]
        a[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = np.empty(n)
    x[-1] = d[-1] / a[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / a[i]
    return x


def _seed_profile(x: np.ndarray, seed) -> np.ndarray:
    raw = seed(x)
    d = np.maximum(-np.diff(raw), 1e-15)
    d *= math.pi / np.sum(d)
    F = np.empty_like(raw)
    F[0] = math.pi
    F[1:] = math.pi - np.cumsum(d)
    F[-1] = 0.0
    return F


def minimize_on(x: np.ndarray, seed) -> dict:
    """Newton on the midpoint energy. F(0)=π, F(X)=0, free nodes in between.

    The softmax chart named in PREREG is kept for the gradient meter.
    It stalls short of ||∇_y E|| ≤ 1e-6 (ill-conditioned). Stationarity
    reported here is ||∂E/∂F||_∞ on this same discrete energy.
    """
    F = _seed_profile(x, seed)
    nit = 0
    success = False
    ginf = math.inf
    for nit in range(1, 41):
        E2, E4, gF = parts_and_grad_F(F, x)
        g = gF[1:-1]
        ginf = float(np.max(np.abs(g)))
        if ginf <= 1e-10:
            success = True
            break
        diag, off = _hessian_tridiagonal(F, x)
        # Levenberg damping if a pivot is not safely positive.
        shift = 0.0
        if np.any(diag <= 1e-8):
            shift = 1e-4 * max(float(np.max(np.abs(diag))), 1.0)
            diag = diag + shift
        try:
            step = _solve_tridiagonal(diag, off, g)
        except ZeroDivisionError:
            break
        E_now = E2 + E4
        alpha = 1.0
        accepted = False
        for _ in range(12):
            trial = F.copy()
            trial[1:-1] = F[1:-1] - alpha * step
            t2, t4, _ = parts_and_grad_F(trial, x)
            if t2 + t4 < E_now:
                F = trial
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break
    else:
        success = ginf <= 1e-6
    E2, E4, gF = parts_and_grad_F(F, x)
    ginf = float(np.max(np.abs(gF[1:-1])))
    if ginf <= 1e-6:
        success = True
    return {
        "F": F,
        "E2": E2,
        "E4": E4,
        "E": E2 + E4,
        "success": success,
        "grad_inf": ginf,
        "nit": int(nit),
        "message": "newton on F; softmax chart not used for the solve",
    }


def energy_of_scaled(F: np.ndarray, x: np.ndarray, lam: float) -> tuple[float, float, float]:
    xmax = float(x[-1] * lam)
    n = max(int(round((len(x) - 1) * lam)), len(x) - 1)
    xs = np.linspace(0.0, xmax, n + 1)
    # F_λ(x) = F(x/λ). Outside the solved box the ansatz is vacuum.
    Fs = np.interp(xs / lam, x, F, left=math.pi, right=0.0)
    Fs[0] = math.pi
    Fs[-1] = 0.0
    E2, E4, _ = parts_and_grad_F(Fs, xs)
    return E2, E4, E2 + E4


def check_gradient(x: np.ndarray, rng: np.random.Generator) -> float:
    y = rng.normal(size=x.size - 1)
    F = F_of_y(y)
    _, _, gF = parts_and_grad_F(F, x)
    g = grad_y(gF, y)
    # One random direction, central difference.
    v = rng.normal(size=y.size)
    v /= np.linalg.norm(v)
    eps = 1e-6
    Ep = parts_and_grad_F(F_of_y(y + eps * v), x)
    Em = parts_and_grad_F(F_of_y(y - eps * v), x)
    num = ((Ep[0] + Ep[1]) - (Em[0] + Em[1])) / (2 * eps)
    ana = float(np.dot(g, v))
    return abs(num - ana) / max(abs(num), 1e-12)


def check_scaler(x: np.ndarray) -> tuple[float, float]:
    F = F_of_y(y_from_seed(seed_A(x)))
    E2, E4, _ = parts_and_grad_F(F, x)
    errs = []
    for lam in (0.5, 2.0):
        _, _, Es = energy_of_scaled(F, x, lam)
        pred = lam * E2 + E4 / lam
        errs.append(abs(Es - pred) / abs(pred))
    return errs[0], errs[1]


def inertia(F: np.ndarray, x: np.ndarray) -> float:
    """(8π/3) ∫ x² sin²F (1 + F'² + sin²F/x²) dx, units 1/(e³ f_π)."""
    dx = np.diff(x)
    Fm = 0.5 * (F[:-1] + F[1:])
    Fp = np.diff(F) / dx
    xm = 0.5 * (x[:-1] + x[1:])
    s2 = np.sin(Fm) ** 2
    dens = xm**2 * s2 * (1.0 + Fp**2 + s2 / xm**2)
    return (8.0 * math.pi / 3.0) * float(np.sum(dx * dens))


def diagnose_softmax_chart() -> None:
    """The chart named in PREREG. It is not the production solve.

    A constant shift of y does not change Δ, and the tail decrements are
    exponentially small, so L-BFGS stalls above the stationarity line.
    """
    x = grid(200, 20.0)
    y0 = y_from_seed(seed_A(x))

    def fun(y: np.ndarray):
        F = F_of_y(y)
        E2, E4, gF = parts_and_grad_F(F, x)
        return E2 + E4, grad_y(gF, y)

    res = minimize(
        fun,
        y0,
        method="L-BFGS-B",
        jac=True,
        options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 4000},
    )
    E = float(res.fun)
    ginf = float(np.max(np.abs(res.jac)))
    print(
        f"CHART softmax N=200 X=20 E={E:.6f} grad_y_inf={ginf:.3e} "
        f"success={bool(res.success)} nit={int(res.nit)}"
    )


def main() -> None:
    diagnose_softmax_chart()
    rng = np.random.default_rng(20260921)
    x_fine_for_grad = grid(200, 30.0)
    gerr = check_gradient(x_fine_for_grad, rng)
    serr = check_scaler(grid(500, 30.0))
    print(f"METER gradient rel {gerr:.3e}")
    print(f"METER scaler lam=1/2 rel {serr[0]:.3e}")
    print(f"METER scaler lam=2 rel {serr[1]:.3e}")

    grids = (("coarse", 500, 30.0), ("mid", 1000, 60.0), ("fine", 2000, 60.0))
    solved: dict[str, dict] = {}
    for name, n, xmax in grids:
        x = grid(n, xmax)
        a = minimize_on(x, seed_A)
        b = minimize_on(x, seed_B)
        rel = abs(a["E"] - b["E"]) / max(a["E"], b["E"])
        # Report the lower of the two; they must already agree.
        best = a if a["E"] <= b["E"] else b
        best["x"] = x
        best["seed_rel"] = rel
        best["E_A"] = a["E"]
        best["E_B"] = b["E"]
        best["success_A"] = a["success"]
        best["success_B"] = b["success"]
        best["grad_A"] = a["grad_inf"]
        best["grad_B"] = b["grad_inf"]
        solved[name] = best
        print(
            f"RUN {name} N={n} X={xmax:.0f} "
            f"E_A={a['E']:.6f} E_B={b['E']:.6f} seed_rel={rel:.3e} "
            f"virial_num={abs(best['E2']-best['E4'])/(best['E2']+best['E4']):.3e} "
            f"grad_inf={best['grad_inf']:.3e} nit={best['nit']} "
            f"ok={best['success']}"
        )

    fine = solved["fine"]
    mid = solved["mid"]
    F = fine["F"]
    x = fine["x"]
    E = fine["E"]
    E2, E4 = fine["E2"], fine["E4"]
    virial = abs(E2 - E4) / (E2 + E4)
    refine = abs(E - mid["E"]) / E
    textbook = abs(E - TEXTBOOK_E) / TEXTBOOK_E
    I = inertia(F, x)
    print(f"METER seeds rel {fine['seed_rel']:.3e}")
    print(f"METER virial {virial:.3e} E2={E2:.6f} E4={E4:.6f}")
    print(f"METER stationarity grad_inf {fine['grad_inf']:.3e} success {fine['success']}")
    print(f"METER refinement rel {refine:.3e}")
    print(f"METER textbook rel {textbook:.3e} E={E:.6f}")
    print(f"RECORD inertia {I:.6f}")

    ratios = {}
    scan_lines = ["lam\tE2\tE4\tE\tratio\n"]
    for lam in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
        s2, s4, sE = energy_of_scaled(F, x, lam)
        ratio = sE / E
        ratios[lam] = ratio
        scan_lines.append(f"{lam:.2f}\t{s2:.8f}\t{s4:.8f}\t{sE:.8f}\t{ratio:.8f}\n")
        print(f"SCAN lam={lam:.2f} E={sE:.6f} ratio={ratio:.6f}")

    r_half = ratios[0.5]
    r_double = ratios[2.0]
    wall_half = abs(r_half / 1.25 - 1.0)
    wall_double = abs(r_double / 1.25 - 1.0)
    print(f"METER wall lam=1/2 dev {wall_half:.3e}")
    print(f"METER wall lam=2 dev {wall_double:.3e}")

    meters_clean = (
        gerr <= 1e-6
        and serr[0] <= 1e-3
        and serr[1] <= 1e-3
        and fine["seed_rel"] <= 5e-3
        and virial <= 1e-3
        and fine["success"]
        and fine["grad_inf"] <= 1e-6
        and refine <= 1e-2
        and textbook <= 0.02
        and wall_half <= 0.02
        and wall_double <= 0.02
    )
    if not meters_clean:
        g1 = "INCONCLUSIVE"
    elif r_half >= 1.05 and r_double >= 1.05:
        g1 = "PASS"
    else:
        g1 = "FAIL"
    print(f"G1 {g1}")

    prof = OUT / "profile.tsv"
    with prof.open("w") as f:
        f.write("x\tF\n")
        for xi, Fi in zip(x, F):
            f.write(f"{xi:.8e}\t{Fi:.8e}\n")
    (OUT / "scan.tsv").write_text("".join(scan_lines))
    print(f"WROTE {prof}")


if __name__ == "__main__":
    main()
