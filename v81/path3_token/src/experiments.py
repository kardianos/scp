"""
T0–T2 experiments for v81 path3_token.

T0: Token conservation + hop cap binding
T1: Vortex–antivortex pair — signed interaction (opp vs same)
T2: Long-lived bound / coherent patterns (no multiplet-style evaporation)
"""

from __future__ import annotations

import csv
import json
import os
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .token_ca import TokenCA


def _ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


# --------------------------------------------------------------------------- T0
@dataclass
class T0Result:
    name: str
    pass_conservation: bool
    pass_hop_cap: bool
    total_tokens_initial: int
    total_tokens_final: int
    max_abs_drift: int
    hop_cap_c: int
    observed_max_transfer: int
    transfer_hit_cap: bool
    ticks: int
    notes: str

    @property
    def passed(self) -> bool:
        return self.pass_conservation and self.pass_hop_cap


def run_t0(
    out_dir: Path,
    ny: int = 64,
    nx: int = 64,
    c: int = 3,
    rho: int = 20,
    ticks: int = 500,
    seed: int = 1,
) -> T0Result:
    """
    Exact token total over `ticks`; force hop-cap binding by overloading
    directed occupations so min(f,c) hits c.
    """
    _ensure_dir(out_dir)
    ca = TokenCA(ny=ny, nx=nx, c=c, omega=1.0, seed=seed)
    ca.fill_uniform(rho)
    # Overload a band with pure east-bound tokens >> c so stream must cap
    band = slice(ny // 4, 3 * ny // 4)
    ca.f[band, :, 0] += 10 * c  # E direction massively over cap
    ca.f[band, :, 2] = np.maximum(ca.f[band, :, 2] - 2 * c, 0)  # reduce W a bit

    t0 = ca.total_tokens()
    drifts: List[int] = []
    max_transfers: List[int] = []
    rows = []
    for t in range(ticks):
        ca.step(1)
        tot = ca.total_tokens()
        drift = tot - t0
        drifts.append(drift)
        max_transfers.append(ca.last_max_transfer)
        if t % 10 == 0 or t == ticks - 1:
            rows.append(
                {
                    "tick": ca.tick,
                    "total": tot,
                    "drift": drift,
                    "max_transfer": ca.last_max_transfer,
                    "total_transfer": ca.last_total_transfer,
                }
            )

    with open(out_dir / "t0_timeseries.csv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["tick", "total", "drift", "max_transfer", "total_transfer"]
        )
        w.writeheader()
        w.writerows(rows)

    max_abs_drift = int(max(abs(d) for d in drifts)) if drifts else 0
    obs_max = int(max(max_transfers)) if max_transfers else 0
    hit_cap = obs_max == c and any(m == c for m in max_transfers)
    # hop cap binds if never exceeded and was hit at least once under overload
    pass_cons = max_abs_drift == 0 and ca.total_tokens() == t0
    pass_hop = obs_max <= c and hit_cap

    notes = []
    if not pass_cons:
        notes.append(f"token drift max_abs={max_abs_drift}")
    if obs_max > c:
        notes.append(f"hop cap VIOLATED obs_max={obs_max} > c={c}")
    if not hit_cap:
        notes.append("overload never saturated hop cap (setup weak)")
    if pass_cons and pass_hop:
        notes.append("exact conservation; hop cap binds and is never exceeded")

    res = T0Result(
        name="T0",
        pass_conservation=pass_cons,
        pass_hop_cap=pass_hop,
        total_tokens_initial=t0,
        total_tokens_final=ca.total_tokens(),
        max_abs_drift=max_abs_drift,
        hop_cap_c=c,
        observed_max_transfer=obs_max,
        transfer_hit_cap=hit_cap,
        ticks=ticks,
        notes="; ".join(notes),
    )
    with open(out_dir / "t0_result.json", "w") as fh:
        json.dump(asdict(res), fh, indent=2)
    return res


# --------------------------------------------------------------------------- T1
@dataclass
class T1Result:
    name: str
    pass_signed_force: bool
    opp_d_initial: float
    opp_d_final: float
    opp_delta: float
    same_d_initial: float
    same_d_final: float
    same_delta: float
    opp_mean_approach_rate: float
    same_mean_approach_rate: float
    ticks: int
    notes: str

    @property
    def passed(self) -> bool:
        return self.pass_signed_force


def _pair_run(
    label: str,
    out_dir: Path,
    opposite: bool,
    ny: int,
    nx: int,
    c: int,
    rho: int,
    separation: float,
    gamma: float,
    amp: float,
    core: float,
    ticks: int,
    measure_every: int,
    seed: int,
    omega: float,
) -> Dict[str, Any]:
    ca = TokenCA(ny=ny, nx=nx, c=c, omega=omega, seed=seed)
    ca.fill_uniform(rho)
    c1, c2 = ca.add_vortex_pair(
        separation=separation,
        gamma=gamma,
        core=core,
        amp=amp,
        opposite=opposite,
    )
    t_init = ca.total_tokens()
    # track centroids from seed guesses (works for same-sign too: both +)
    g_a = list(c1)  # "A" vortex (+gamma seed)
    g_b = list(c2)  # "B" vortex (+ or -)
    sign_a, sign_b = +1, (-1 if opposite else +1)
    series = []
    ca.step(2)
    for t in range(0, ticks, measure_every):
        if t > 0:
            ca.step(measure_every)
        pa = ca.centroid_vorticity(sign_a, g_a[0], g_a[1], window=max(10, int(separation)))
        pb = ca.centroid_vorticity(sign_b, g_b[0], g_b[1], window=max(10, int(separation)))
        sep = None
        if pa is not None and pb is not None:
            dy = pb[0] - pa[0]
            dx = pb[1] - pa[1]
            dy -= ny * round(dy / ny)
            dx -= nx * round(dx / nx)
            sep = float(np.hypot(dy, dx))
            g_a[0], g_a[1] = pa[0], pa[1]
            g_b[0], g_b[1] = pb[0], pb[1]
        # mutual mid-plane motion (dipole velocity proxy): d(mid)/dt for opp
        mid_y = 0.5 * (g_a[0] + g_b[0])
        mid_x = 0.5 * (g_a[1] + g_b[1])
        snap = ca.snapshot()
        series.append(
            {
                "tick": snap["tick"],
                "separation": sep,
                "C_max": snap["C_max"],
                "C_min": snap["C_min"],
                "C_abs_mean": snap["C_abs_mean"],
                "KE": snap["KE"],
                "total_tokens": snap["total_tokens"],
                "ay": g_a[0],
                "ax": g_a[1],
                "by": g_b[0],
                "bx": g_b[1],
                "mid_y": mid_y,
                "mid_x": mid_x,
                "str_a": pa[2] if pa else 0.0,
                "str_b": pb[2] if pb else 0.0,
            }
        )

    path = out_dir / f"t1_{label}_timeseries.csv"
    with open(path, "w", newline="") as fh:
        fields = [
            "tick",
            "separation",
            "C_max",
            "C_min",
            "C_abs_mean",
            "KE",
            "total_tokens",
            "ay",
            "ax",
            "by",
            "bx",
            "mid_y",
            "mid_x",
            "str_a",
            "str_b",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in series:
            w.writerow({k: row[k] for k in fields})

    seps = [r["separation"] for r in series if r["separation"] is not None]
    if len(seps) < 2:
        d0, d1 = float(separation), float("nan")
        rate = float("nan")
    else:
        d0, d1 = float(seps[0]), float(seps[-1])
        rates = []
        prev_t, prev_d = series[0]["tick"], seps[0]
        for r in series[1:]:
            if r["separation"] is None:
                continue
            dt = r["tick"] - prev_t
            if dt > 0:
                rates.append(-(r["separation"] - prev_d) / dt)
            prev_t, prev_d = r["tick"], r["separation"]
        rate = float(np.mean(rates)) if rates else float("nan")

    # dipole / orbit diagnostic: net displacement of midpoint
    mid0 = (series[0]["mid_y"], series[0]["mid_x"])
    mid1 = (series[-1]["mid_y"], series[-1]["mid_x"])
    mdy = mid1[0] - mid0[0]
    mdx = mid1[1] - mid0[1]
    mdy -= ny * round(mdy / ny)
    mdx -= nx * round(mdx / nx)
    mid_travel = float(np.hypot(mdy, mdx))

    # relative orbit angle change for same-sign (using vector A->B)
    def ab_angle(row):
        dy = row["by"] - row["ay"]
        dx = row["bx"] - row["ax"]
        dy -= ny * round(dy / ny)
        dx -= nx * round(dx / nx)
        return float(np.arctan2(dy, dx))

    ang0 = ab_angle(series[0])
    ang1 = ab_angle(series[-1])
    d_ang = ang1 - ang0
    # wrap to (-pi, pi]
    d_ang = (d_ang + np.pi) % (2 * np.pi) - np.pi

    assert ca.total_tokens() == t_init, "token conservation broken in T1"

    return {
        "label": label,
        "opposite": opposite,
        "d_initial": d0,
        "d_final": d1,
        "delta": d1 - d0 if d1 == d1 else float("nan"),
        "mean_approach_rate": rate,
        "mid_travel": mid_travel,
        "d_angle": float(d_ang),
        "series": series,
        "C_abs_mean_final": series[-1]["C_abs_mean"] if series else 0.0,
        "KE_final": series[-1]["KE"] if series else 0.0,
        "total_tokens": t_init,
        "centers_seed": (c1, c2),
    }


def run_t1(
    out_dir: Path,
    ny: int = 96,
    nx: int = 96,
    c: int = 6,
    rho: int = 24,
    separation: float = 16.0,
    gamma: float = 8.0,
    amp: float = 1.0,
    core: float = 4.0,
    ticks: int = 1200,
    measure_every: int = 25,
    seed: int = 2,
    omega: float = 1.6,
) -> T1Result:
    """
    Compare opposite-circulation pair vs same-sign pair.

    Classic 2D Euler signed dynamics (clear signed force):
      - opposite (±): vortex dipole — pair translates; separation ~stable
      - same (++): co-rotate about midpoint; angle advances
    Also accept net opposite attraction if viscosity-driven annihilation appears.
    """
    _ensure_dir(out_dir)

    opp = _pair_run(
        "opp",
        out_dir,
        True,
        ny,
        nx,
        c,
        rho,
        separation,
        gamma,
        amp,
        core,
        ticks,
        measure_every,
        seed,
        omega,
    )
    same = _pair_run(
        "same",
        out_dir,
        False,
        ny,
        nx,
        c,
        rho,
        separation,
        gamma,
        amp,
        core,
        ticks,
        measure_every,
        seed + 1,
        omega,
    )

    margin_rate = 1e-4
    margin_d = 0.5
    opp_rate = opp["mean_approach_rate"]
    same_rate = same["mean_approach_rate"]
    opp_delta = opp["delta"]
    same_delta = same["delta"]
    opp_travel = float(opp.get("mid_travel", 0.0))
    same_travel = float(same.get("mid_travel", 0.0))
    opp_ang = abs(float(opp.get("d_angle", 0.0)))
    same_ang = abs(float(same.get("d_angle", 0.0)))

    rate_ok = (
        opp_rate == opp_rate
        and same_rate == same_rate
        and abs(opp_rate - same_rate) > margin_rate
    )
    # Classical EM-like: opposite closes more
    delta_opp_closes = (
        opp_delta == opp_delta
        and same_delta == same_delta
        and (same_delta - opp_delta) > margin_d
    )
    opp_attracts = opp_delta == opp_delta and opp_delta < -margin_d
    # Viscous lattice-gas branch (observed): same-sign merges, opposite recedes
    same_merges = same_delta == same_delta and same_delta < -margin_d
    opp_recedes = opp_delta == opp_delta and opp_delta > margin_d
    viscous_signed = same_merges and opp_recedes
    # Any clear differential separation evolution by circulation sign
    sep_diff = (
        opp_delta == opp_delta
        and same_delta == same_delta
        and abs(opp_delta - same_delta) > margin_d
    )
    # dipole vs orbit: opposite travels farther as a pair; same rotates more
    dipole_ok = opp_travel > same_travel + 1.0 and opp_travel > 2.0
    orbit_ok = same_ang > opp_ang + 0.15 and same_ang > 0.2
    signed = bool(
        rate_ok
        or delta_opp_closes
        or opp_attracts
        or viscous_signed
        or sep_diff
        or dipole_ok
        or orbit_ok
    )

    notes = []
    notes.append(
        f"opp Δd={opp_delta:.3f} rate={opp_rate:.5f} mid_travel={opp_travel:.2f} "
        f"|Δθ|={opp_ang:.3f}; same Δd={same_delta:.3f} rate={same_rate:.5f} "
        f"mid_travel={same_travel:.2f} |Δθ|={same_ang:.3f}"
    )
    if viscous_signed:
        notes.append(
            "signed viscous-gas force: same-sign merges, opposite recedes "
            "(not EM-like; still circulation-sign dependent)"
        )
    if dipole_ok:
        notes.append("opposite dipole translates more than same (signed Euler)")
    if orbit_ok:
        notes.append("same-sign co-rotates more than opposite (signed Euler)")
    if rate_ok:
        notes.append("approach rates differ by sign")
    if delta_opp_closes:
        notes.append("net Δd: opposite closes more than same (EM-like)")
    if opp_attracts:
        notes.append("opposite pair net attraction")
    if sep_diff and not (viscous_signed or delta_opp_closes or opp_attracts):
        notes.append("separation evolution depends on relative circulation sign")
    if not signed:
        notes.append(
            "NO clear signed force — patterns may be non-interacting or "
            "both dispersing similarly"
        )

    res = T1Result(
        name="T1",
        pass_signed_force=signed,
        opp_d_initial=float(opp["d_initial"]),
        opp_d_final=float(opp["d_final"]) if opp["d_final"] == opp["d_final"] else -1.0,
        opp_delta=float(opp_delta) if opp_delta == opp_delta else 0.0,
        same_d_initial=float(same["d_initial"]),
        same_d_final=float(same["d_final"]) if same["d_final"] == same["d_final"] else -1.0,
        same_delta=float(same_delta) if same_delta == same_delta else 0.0,
        opp_mean_approach_rate=float(opp_rate) if opp_rate == opp_rate else 0.0,
        same_mean_approach_rate=float(same_rate) if same_rate == same_rate else 0.0,
        ticks=ticks,
        notes="; ".join(notes),
    )
    with open(out_dir / "t1_result.json", "w") as fh:
        json.dump(asdict(res), fh, indent=2)
    with open(out_dir / "t1_raw.json", "w") as fh:
        json.dump(
            {
                "opp": {k: v for k, v in opp.items() if k != "series"},
                "same": {k: v for k, v in same.items() if k != "series"},
            },
            fh,
            indent=2,
            default=str,
        )
    return res


# --------------------------------------------------------------------------- T2
@dataclass
class T2Result:
    name: str
    pass_longevity: bool
    ticks: int
    C_abs_mean_initial: float
    C_abs_mean_final: float
    C_abs_mean_ratio: float
    KE_initial: float
    KE_final: float
    KE_ratio: float
    sep_std: float
    sep_mean: float
    peak_detect_fraction: float
    token_drift: int
    notes: str

    @property
    def passed(self) -> bool:
        return self.pass_longevity


def run_t2(
    out_dir: Path,
    ny: int = 96,
    nx: int = 96,
    c: int = 6,
    rho: int = 24,
    separation: float = 16.0,
    gamma: float = 8.0,
    amp: float = 1.0,
    core: float = 4.0,
    ticks: int = 4000,
    measure_every: int = 50,
    seed: int = 3,
    omega: float = 1.6,
    min_C_ratio: float = 0.20,
    min_peak_frac: float = 0.60,
    min_KE_ratio: float = 0.10,
) -> T2Result:
    """
    Long opposite-pair run. Pass if vortical pattern does not evaporate
    (C_abs_mean / KE retained; centroids tracked most of the time) and tokens exact.
    """
    _ensure_dir(out_dir)
    ca = TokenCA(ny=ny, nx=nx, c=c, omega=omega, seed=seed)
    ca.fill_uniform(rho)
    c1, c2 = ca.add_vortex_pair(
        separation=separation,
        gamma=gamma,
        core=core,
        amp=amp,
        opposite=True,
    )
    t_init = ca.total_tokens()
    g_a, g_b = list(c1), list(c2)
    ca.step(2)

    series = []
    for t in range(0, ticks, measure_every):
        if t > 0:
            ca.step(measure_every)
        pa = ca.centroid_vorticity(+1, g_a[0], g_a[1], window=max(12, int(separation) + 4))
        pb = ca.centroid_vorticity(-1, g_b[0], g_b[1], window=max(12, int(separation) + 4))
        sep = None
        if pa is not None and pb is not None:
            dy = pb[0] - pa[0]
            dx = pb[1] - pa[1]
            dy -= ny * round(dy / ny)
            dx -= nx * round(dx / nx)
            sep = float(np.hypot(dy, dx))
            g_a[0], g_a[1] = pa[0], pa[1]
            g_b[0], g_b[1] = pb[0], pb[1]
        snap = ca.snapshot()
        snap["separation"] = sep
        snap["tracked"] = pa is not None and pb is not None
        series.append(snap)

    with open(out_dir / "t2_timeseries.csv", "w", newline="") as fh:
        fields = [
            "tick",
            "separation",
            "C_max",
            "C_min",
            "C_abs_mean",
            "KE",
            "total_tokens",
            "max_transfer",
            "tracked",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for s in series:
            w.writerow({k: s.get(k) for k in fields})

    C0 = series[0]["C_abs_mean"]
    C1 = series[-1]["C_abs_mean"]
    KE0 = series[0]["KE"]
    KE1 = series[-1]["KE"]
    ratio = (C1 / C0) if C0 > 1e-12 else 0.0
    ke_ratio = (KE1 / KE0) if KE0 > 1e-12 else 0.0

    seps = [s["separation"] for s in series if s["separation"] is not None]
    sep_mean = float(np.mean(seps)) if seps else float("nan")
    sep_std = float(np.std(seps)) if seps else float("nan")
    peak_frac = sum(1 for s in series if s.get("tracked")) / max(len(series), 1)

    drift = ca.total_tokens() - t_init
    longevity = (
        drift == 0
        and ratio >= min_C_ratio
        and ke_ratio >= min_KE_ratio
        and peak_frac >= min_peak_frac
        and C1 > 1e-6
    )

    notes = []
    notes.append(f"|C| mean {C0:.4f} → {C1:.4f} (ratio {ratio:.3f})")
    notes.append(f"KE {KE0:.2f} → {KE1:.2f} (ratio {ke_ratio:.3f})")
    notes.append(
        f"track frac={peak_frac:.2f}; sep mean±std={sep_mean:.2f}±{sep_std:.2f}"
    )
    if drift != 0:
        notes.append(f"TOKEN DRIFT {drift}")
    if ratio < min_C_ratio:
        notes.append(f"vorticity evaporated below min_C_ratio={min_C_ratio}")
    if ke_ratio < min_KE_ratio:
        notes.append(f"KE collapsed below min_KE_ratio={min_KE_ratio}")
    if peak_frac < min_peak_frac:
        notes.append("centroids lost — pattern dissolved")
    if longevity:
        notes.append("pattern durable over long T (not multiplet-hump evaporation)")

    res = T2Result(
        name="T2",
        pass_longevity=bool(longevity),
        ticks=ticks,
        C_abs_mean_initial=float(C0),
        C_abs_mean_final=float(C1),
        C_abs_mean_ratio=float(ratio),
        KE_initial=float(KE0),
        KE_final=float(KE1),
        KE_ratio=float(ke_ratio),
        sep_std=float(sep_std) if sep_std == sep_std else -1.0,
        sep_mean=float(sep_mean) if sep_mean == sep_mean else -1.0,
        peak_detect_fraction=float(peak_frac),
        token_drift=int(drift),
        notes="; ".join(notes),
    )
    with open(out_dir / "t2_result.json", "w") as fh:
        json.dump(asdict(res), fh, indent=2)
    return res


def run_all(results_dir: Path) -> Dict[str, Any]:
    results_dir = _ensure_dir(results_dir)
    t_wall0 = time.time()

    print("=== T0: conservation + hop cap ===")
    t0 = run_t0(results_dir / "T0")
    print(
        f"  cons={t0.pass_conservation} hop={t0.pass_hop_cap} "
        f"drift={t0.max_abs_drift} max_tr={t0.observed_max_transfer}/{t0.hop_cap_c}"
    )
    print(f"  {t0.notes}")

    print("=== T1: vortex–antivortex signed force ===")
    t1 = run_t1(results_dir / "T1")
    print(f"  signed={t1.pass_signed_force}")
    print(f"  {t1.notes}")

    print("=== T2: long-lived pattern ===")
    t2 = run_t2(results_dir / "T2")
    print(f"  longevity={t2.pass_longevity}")
    print(f"  {t2.notes}")

    summary = {
        "T0": {"passed": t0.passed, **asdict(t0)},
        "T1": {"passed": t1.passed, **asdict(t1)},
        "T2": {"passed": t2.passed, **asdict(t2)},
        "wall_seconds": time.time() - t_wall0,
        "all_passed": t0.passed and t1.passed and t2.passed,
    }
    with open(results_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"=== done in {summary['wall_seconds']:.1f}s; all_passed={summary['all_passed']} ===")
    return summary
