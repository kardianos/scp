#!/usr/bin/env python3
"""v101 experiment 2 — octave return, then the 2:1 interval chain.

Thresholds live in ../PREREG.md. The unison fractional ring is printed
as a control and is not a score.

Run: python3 v101/exp2/cas/octave.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

P_GATE = 8
A_COUP = 1.0
N_SITES = 400
E_CORE = 7.432761
TWOPI = 2.0 * math.pi


def monodromy_11(n: int) -> float:
    return TWOPI * n


def monodromy_oct(n1: int, n2: int) -> float:
    """(2:1) then (1:2), alpha = 0. c := (p*c + 2*pi*n)/q."""
    c = (2.0 * 0.0 + TWOPI * n1) / 1.0
    c = (1.0 * c + TWOPI * n2) / 2.0
    return c


def v_gate(th: np.ndarray) -> np.ndarray:
    base = 0.5 * (1.0 + np.cos(th))
    return 1.0 - base**P_GATE


def dv_gate(th: np.ndarray) -> np.ndarray:
    base = 0.5 * (1.0 + np.cos(th))
    return 0.5 * P_GATE * np.sin(th) * base ** (P_GATE - 1)


def unison_energy(th: np.ndarray, periodic: bool, twist: float) -> float:
    e = float(np.sum(v_gate(th)))
    d = np.diff(th)
    e += 0.5 * A_COUP * float(np.dot(d, d))
    if periodic:
        seam = th[0] + twist - th[-1]
        e += 0.5 * A_COUP * seam * seam
    return e


def unison_grad(th: np.ndarray, periodic: bool, twist: float) -> np.ndarray:
    g = dv_gate(th)
    g[1:] += A_COUP * (th[1:] - th[:-1])
    g[:-1] += A_COUP * (th[:-1] - th[1:])
    if periodic:
        seam = th[0] + twist - th[-1]
        g[0] += A_COUP * seam
        g[-1] -= A_COUP * seam
    return g


def interval_phi(th: np.ndarray) -> np.ndarray:
    """Φ_i = θ_i - 2 θ_{i-1}, length N-1, aligned to sites 1..N-1."""
    return th[1:] - 2.0 * th[:-1]


def interval_energy(th: np.ndarray) -> float:
    phi = interval_phi(th)
    e = float(np.sum(v_gate(phi)))
    d = np.diff(th)
    e += 0.5 * A_COUP * float(np.dot(d, d))
    return e


def interval_grad(th: np.ndarray) -> np.ndarray:
    phi = interval_phi(th)
    dv = dv_gate(phi)
    g = np.zeros_like(th)
    # Φ_i (index i-1 in phi) depends on θ_i with +1 and θ_{i-1} with -2.
    g[1:] += dv
    g[:-1] += -2.0 * dv
    g[1:] += A_COUP * (th[1:] - th[:-1])
    g[:-1] += A_COUP * (th[:-1] - th[1:])
    return g


def descend(energy, grad, th: np.ndarray, free: slice, iters: int = 400_000) -> tuple[np.ndarray, float]:
    """Overdamped descent. Step starts at 0.02 and halves when E rises."""
    th = th.copy()
    lr = 0.02
    e = energy(th)
    last_g = math.inf
    for _ in range(iters):
        g = grad(th)
        ginf = float(np.max(np.abs(g[free])))
        last_g = ginf
        if ginf <= 1e-10:
            break
        step = lr
        moved = False
        for _try in range(20):
            trial = th.copy()
            trial[free] = th[free] - step * g[free]
            et = energy(trial)
            # A 1e-18 rise is float noise on E ~ 7, not a real increase.
            if et <= e + 1e-12:
                th = trial
                e = et
                moved = True
                break
            step *= 0.5
        if not moved:
            break
    ginf = float(np.max(np.abs(grad(th)[free])))
    return th, ginf


def seed_ramp(n: int, twist: float, periodic: bool) -> np.ndarray:
    th = np.zeros(n)
    denom = float(n if periodic else n - 1)
    for i in range(n):
        th[i] = twist * i / denom
    if not periodic:
        th[0] = 0.0
        th[-1] = twist
    return th


def seed_steps(n: int, k: int, shift: int, periodic: bool, twist: float) -> np.ndarray:
    th = np.zeros(n)
    for i in range(n):
        units = 0
        for m in range(k):
            centre = int((m + 0.5) * n / k) + shift
            if i > centre:
                units += 1
        th[i] = TWOPI * units
        for m in range(k):
            centre = int((m + 0.5) * n / k) + shift
            if i == centre:
                th[i] += math.pi
    if not periodic:
        th[0] = 0.0
        th[-1] = twist
    return th


def relax_unison(periodic: bool, twist: float, k_seed: int) -> tuple[float, float, np.ndarray]:
    n = N_SITES
    seeds = [seed_ramp(n, twist, periodic)]
    if k_seed > 0:
        seeds.append(seed_steps(n, k_seed, 0, periodic, twist))
        seeds.append(seed_steps(n, k_seed, 1, periodic, twist))
    free = slice(0, n) if periodic else slice(1, n - 1)
    best_e = math.inf
    best_g = math.inf
    best = seeds[0]
    for s in seeds:
        def energy(th, _p=periodic, _t=twist):
            return unison_energy(th, _p, _t)

        def grad(th, _p=periodic, _t=twist):
            return unison_grad(th, _p, _t)

        th, ginf = descend(energy, grad, s, free)
        e = energy(th)
        if e < best_e:
            best_e, best_g, best = e, ginf, th
    return best_e, best_g, best


def relax_interval(twist: float) -> tuple[float, float, np.ndarray]:
    n = N_SITES
    seeds = [
        seed_ramp(n, twist, False),
        seed_steps(n, 1, 0, False, twist),
        seed_steps(n, 1, 1, False, twist),
    ]
    free = slice(1, n - 1)
    best_e = math.inf
    best_g = math.inf
    best = seeds[0]
    for s in seeds:
        th, ginf = descend(interval_energy, interval_grad, s, free)
        e = interval_energy(th)
        if e < best_e:
            best_e, best_g, best = e, ginf, th
    return best_e, best_g, best


def localized(th: np.ndarray) -> tuple[bool, float, float]:
    phi = interval_phi(th)
    hot = float(np.mean(v_gate(phi) > 0.01))
    d = np.abs(np.diff(th))
    mean = float(np.mean(d))
    mx = float(np.max(d))
    ratio = mx / mean if mean > 0 else math.inf
    return hot <= 0.05 and ratio >= 3.0, hot, ratio


def grad_check(rng: np.random.Generator) -> float:
    th = rng.normal(size=N_SITES)
    th[0] = 0.0
    th[-1] = TWOPI
    g = interval_grad(th)
    v = rng.normal(size=N_SITES - 2)
    v /= np.linalg.norm(v)
    eps = 1e-6
    plus = th.copy()
    minus = th.copy()
    plus[1:-1] += eps * v
    minus[1:-1] -= eps * v
    num = (interval_energy(plus) - interval_energy(minus)) / (2 * eps)
    ana = float(np.dot(g[1:-1], v))
    return abs(num - ana) / max(abs(num), 1e-12)


def nonvacuum_unison(th: np.ndarray) -> int:
    return int(np.sum(v_gate(th) > 0.01))


def main() -> None:
    # --- G2A, no chain ---
    m11 = monodromy_11(1)
    m_odd = monodromy_oct(1, 0)
    m_even = monodromy_oct(1, 1)
    psi_odd = 2.0 * m_odd
    psi_even = 2.0 * m_even
    g2a = (
        abs(m11 - TWOPI) < 1e-12
        and abs(m_odd - math.pi) < 1e-12
        and abs(m_even - TWOPI) < 1e-12
        and abs(psi_even - 4.0 * math.pi) < 1e-12
        and abs(psi_odd - TWOPI) < 1e-12
    )
    print(f"G2A 1:1 {m11/math.pi:.6f} pi")
    print(f"G2A octave_odd {m_odd/math.pi:.6f} pi")
    print(f"G2A octave_even {m_even/math.pi:.6f} pi")
    print(f"G2A psi_odd {psi_odd/math.pi:.6f} pi")
    print(f"G2A psi_even {psi_even/math.pi:.6f} pi")
    print(f"G2A {'PASS' if g2a else 'FAIL'}")

    rng = np.random.default_rng(20260921)
    gerr = grad_check(rng)
    print(f"METER gradient rel {gerr:.3e}")

    # --- G2B unison instrument ---
    e_u, g_u, _ = relax_unison(False, TWOPI, 1)
    g2b = abs(e_u / E_CORE - 1.0) <= 1e-3 and g_u <= 1e-8
    print(f"G2B E {e_u:.6f} grad {g_u:.3e} rel {abs(e_u/E_CORE-1):.3e}")
    print(f"G2B {'CLEAN' if g2b else 'DIRTY'}")

    # --- negative control, not a score ---
    e_c, g_c, th_c = relax_unison(True, TWOPI / 3.0, 0)
    nv = nonvacuum_unison(th_c)
    control_clean = nv <= 10
    print(f"CONTROL w=1/3 E {e_c:.6f} nonvacuum {nv} grad {g_c:.3e}")
    print(f"CONTROL {'CLEAN' if control_clean else 'DIRTY'}")

    # --- G2C interval chain ---
    rows = []
    states = {}
    for name, twist in (("vacuum", 0.0), ("odd", math.pi), ("closure", TWOPI)):
        e, g, th = relax_interval(twist)
        loc, hot, ratio = localized(th)
        states[name] = (e, g, loc, hot, ratio)
        print(
            f"CHAIN {name} E {e:.6f} grad {g:.3e} "
            f"hot {hot:.4f} step_ratio {ratio:.3f} loc {loc}"
        )
        rows.append(f"{name}\t{e:.8f}\t{g:.3e}\t{hot:.6f}\t{ratio:.6f}\t{int(loc)}\n")
    (OUT / "chain.tsv").write_text("boundary\tE\tgrad\thot\tfrac\tlocalized\n" + "".join(rows))

    e0, g0, _, _, _ = states["vacuum"]
    e_odd, g_odd, loc_odd, _, _ = states["odd"]
    e_cl, g_cl, loc_cl, _, _ = states["closure"]
    chain_stationary = g0 <= 1e-6 and g_odd <= 1e-6 and g_cl <= 1e-6
    vacuum_ok = e0 <= 1e-6
    meter_ok = gerr <= 1e-6 and g2b and control_clean and chain_stationary and vacuum_ok

    half = e_cl / 2.0 if e_cl != 0 else math.inf
    odd_is_half = loc_odd and abs(e_odd / half - 1.0) <= 0.10
    cost_ok = abs(e_cl / E_CORE - 1.0) <= 0.05
    odd_ok = (not loc_odd) or (not odd_is_half)

    if not g2a:
        verdict = "FAIL"
    elif not meter_ok:
        verdict = "INCONCLUSIVE"
    elif loc_cl and cost_ok and odd_ok:
        verdict = "PASS"
    else:
        verdict = "FAIL"
    print(f"G2C loc {loc_cl} cost_rel {abs(e_cl/E_CORE-1):.4f} odd_half {odd_is_half}")
    print(f"EXP2 {verdict}")


if __name__ == "__main__":
    main()
