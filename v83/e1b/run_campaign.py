#!/usr/bin/env python3
"""E1b parallel campaign: B1 box, B2 seed grid, B3 mass ratio, B4 anti-lock bag surface, B5 soft-core variants."""
from __future__ import annotations

import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/home/d/code/scp")
BIN = ROOT / "bin/scp_sim"
GEN = ROOT / "bin/gen_qball_boost"
PROF = ROOT / "v74/profiles/f_w142_g005.txt"
OUT = ROOT / "v83/e1b"
SEED_DIR = OUT / "seeds"
CFG_DIR = OUT / "cfg"
RES = OUT / "results"
OMP = os.environ.get("E1B_OMP", "4")
MAX_WORKERS = int(os.environ.get("E1B_WORKERS", "4"))
G = 0.05
Q_APPROX = 315.0  # seed Q scale for continuum F


def continuum_F(r: float) -> float:
    """Attractive |F| for q=-1, continuum F = g^2 Q / (4 pi r^2)."""
    return (G * G) * Q_APPROX / (4.0 * math.pi * r * r)


def v_circ(r: float, m: float) -> float:
    return math.sqrt(max(continuum_F(r) * r / m, 1e-18))


def ensure_seed(N: int, L: float) -> Path:
    p = SEED_DIR / f"ball_N{N}_L{int(L)}.sfa"
    if p.exists():
        return p
    SEED_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [str(GEN), str(N), str(L), str(PROF), "1.42", "0", str(p)]
    subprocess.check_call(cmd)
    return p


def write_cfg(
    tag: str,
    *,
    N: int,
    L: float,
    T: float,
    seed: Path,
    lock_csv: str | None = None,
    soft_r: float = 0.0,
    soft_k: float = 0.0,
    bag_mode: int = 0,
    bag_r: float = 0.0,
    bag_k: float = 0.0,
    n_locks: int = 1,
    locks_file: str | None = None,
    diag_dt: float = 2.0,
) -> Path:
    CFG_DIR.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    cfg = CFG_DIR / f"{tag}.cfg"
    track = RES / f"track_{tag}.tsv"
    diag = RES / f"diag_{tag}.tsv"
    out = RES / f"out_{tag}.sfa"
    lines = [
        f"N={N}",
        f"L={L}",
        f"T={T}",
        "dt_factor=0.05",
        "complex_phi=1",
        "complex_gauge=1",
        f"g_gauge={G}",
        "eta=0",
        "m=1.5",
        "m_theta=1.6",
        "mu=-41.345",
        "kappa=50",
        "init=sfa",
        f"init_sfa={seed}",
        "bc_type=0",
        "damp_width=2.5",
        "damp_rate=0.05",
        f"n_locks={n_locks}",
        f"locks_track={track}",
        "locks_medium_only=0",
        f"lock_soft_r={soft_r}",
        f"lock_soft_k={soft_k}",
        f"lock_bag_mode={bag_mode}",
        f"lock_bag_r={bag_r}",
        f"lock_bag_k={bag_k}",
        f"output={out}",
        f"diag_file={diag}",
        "snap_dt=10000",
        f"diag_dt={diag_dt}",
        "precision=f32",
        f"qdiag_radius={min(L * 0.5, 12)}",
    ]
    if locks_file:
        lines.append(f"locks_file={locks_file}")
    elif lock_csv:
        lines.append(f"lock0={lock_csv}")
    else:
        raise ValueError(f"{tag}: need lock_csv or locks_file")
    cfg.write_text("\n".join(lines) + "\n")
    return cfg


def run_one(tag: str, cfg: Path) -> dict:
    log = RES / f"run_{tag}.log"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = OMP
    t0 = time.time()
    try:
        with open(log, "w") as fp:
            r = subprocess.run(
                [str(BIN), str(cfg)],
                stdout=fp,
                stderr=subprocess.STDOUT,
                env=env,
                timeout=int(os.environ.get("E1B_TIMEOUT", "7200")),
            )
        rc = r.returncode
    except subprocess.TimeoutExpired:
        rc = -9
    wall = time.time() - t0
    metrics = analyze_track(tag)
    metrics.update({"tag": tag, "rc": rc, "wall_s": wall})
    return metrics


def analyze_track(tag: str) -> dict:
    track = RES / f"track_{tag}.tsv"
    out = {
        "revs": float("nan"),
        "rmin": float("nan"),
        "rmax": float("nan"),
        "rmean": float("nan"),
        "r0": float("nan"),
        "band": 0,
        "alive": 0,
    }
    if not track.exists():
        return out
    import numpy as np

    t = np.loadtxt(track, skiprows=1)
    if t.ndim == 1:
        t = t.reshape(1, -1)
    # t id type x y z ... — keep charge locks only (type==0) when multi-lock
    if t.shape[1] >= 14:
        charge = t[t[:, 2] == 0]
        if len(charge) > 0:
            t = charge
    # t id type x y z ...
    x, y, z = t[:, 3], t[:, 4], t[:, 5]
    r = np.sqrt(x * x + y * y + z * z)
    ang = np.unwrap(np.arctan2(y, x))
    revs = abs(ang[-1] - ang[0]) / (2 * math.pi)
    r0 = float(r[0])
    rmin, rmax, rmean = float(r.min()), float(r.max()), float(r.mean())
    # soft PASS: >=0.75 rev, r in [0.55,1.7] r0 for most of late window
    late = r[len(r) // 4 :]
    band = int(
        revs >= 0.75
        and late.min() > 0.55 * r0
        and late.max() < 1.7 * r0
        and rmin > 0.4 * r0
    )
    out.update(
        {
            "revs": revs,
            "rmin": rmin,
            "rmax": rmax,
            "rmean": rmean,
            "r0": r0,
            "band": band,
            "alive": int(t[-1, -1]),
        }
    )
    return out


def jobs() -> list[tuple[str, Path]]:
    js: list[tuple[str, Path]] = []
    # Fast matrix: N=48 L=16 (hybrid Cosserat cost ~N^3). T covers ~1–2 T_orb at r~10.
    seed48 = ensure_seed(48, 16.0)
    T_scan = 280.0  # ~0.3–0.5 T_orb; enough for band triage

    # --- B2 seed grid ---
    for r in (8.0, 10.0, 12.0):
        for vf in (0.75, 0.90, 1.00, 1.10, 1.25):
            m = 2.0
            v = vf * v_circ(r, m)
            tag = f"B2_r{r:.0f}_vf{vf:.2f}".replace(".", "p")
            csv = f"-1,{m},{r},0,0,0,{v},0,0"
            cfg = write_cfg(tag, N=48, L=16, T=T_scan, seed=seed48, lock_csv=csv, soft_r=0, soft_k=0)
            js.append((tag, cfg))

    # --- B3 mass ratio at r=10 ---
    for m in (1.0, 2.0, 4.0, 8.0, 16.0):
        r = 10.0
        v = v_circ(r, m)
        tag = f"B3_m{m:.0f}_r10"
        csv = f"-1,{m},{r},0,0,0,{v},0,0"
        cfg = write_cfg(tag, N=48, L=16, T=T_scan, seed=seed48, lock_csv=csv, soft_r=0, soft_k=0)
        js.append((tag, cfg))

    # --- B5 soft core variants ---
    r, m = 10.0, 2.0
    v = v_circ(r, m)
    for soft_r, soft_k, name in (
        (0.0, 0.0, "off"),
        (1.5, 0.10, "soft15"),
        (2.5, 0.15, "soft25"),
        (3.5, 0.20, "soft35"),
    ):
        tag = f"B5_{name}_r10"
        csv = f"-1,{m},{r},0,0,0,{v},0,0"
        cfg = write_cfg(
            tag, N=48, L=16, T=T_scan, seed=seed48, lock_csv=csv, soft_r=soft_r, soft_k=soft_k
        )
        js.append((tag, cfg))

    # --- B4 surface: pinned anti-lock bag + free opposite ---
    for bag_k in (0.05, 0.15, 0.30):
        tag = f"B4_bagk{bag_k:.2f}".replace(".", "p")
        locks_path = CFG_DIR / f"locks_{tag}.tsv"
        v = v_circ(10.0, 2.0)
        locks_path.write_text(
            "# id type q m Es x y z ux uy uz pinned\n"
            "0 1 0 8 8 0 0 0 0 0 0 1\n"
            f"1 0 -1 2 2 10 0 0 0 {v} 0 0\n"
        )
        cfg = write_cfg(
            tag,
            N=48,
            L=16,
            T=T_scan,
            seed=seed48,
            n_locks=2,
            locks_file=str(locks_path),
            soft_r=0,
            soft_k=0,
            bag_mode=1,
            bag_r=5.0,
            bag_k=bag_k,
        )
        js.append((tag, cfg))

    # --- B1 larger box (N=64 L=24, still tractable) ---
    seed64b = ensure_seed(64, 24.0)
    for r, vf in ((10.0, 1.0), (12.0, 1.0), (10.0, 0.9), (12.0, 1.1)):
        m = 2.0
        v = vf * v_circ(r, m)
        tag = f"B1_r{r:.0f}_vf{vf:.1f}".replace(".", "p")
        csv = f"-1,{m},{r},0,0,0,{v},0,0"
        cfg = write_cfg(
            tag,
            N=64,
            L=24,
            T=400.0,
            seed=seed64b,
            lock_csv=csv,
            soft_r=0,
            soft_k=0,
            diag_dt=4.0,
        )
        js.append((tag, cfg))

    return js


def main() -> int:
    RES.mkdir(parents=True, exist_ok=True)
    js = jobs()
    # de-dupe tags
    seen = set()
    uniq = []
    for tag, cfg in js:
        if tag in seen:
            continue
        seen.add(tag)
        uniq.append((tag, cfg))
    print(f"E1b campaign: {len(uniq)} jobs, workers={MAX_WORKERS}, OMP={OMP}", flush=True)
    summary = RES / "summary_e1b.tsv"
    with open(summary, "w") as fp:
        fp.write("tag\trc\twall_s\trevs\tr0\trmin\trmax\trmean\tband\talive\n")
    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(run_one, tag, cfg): tag for tag, cfg in uniq}
        for fut in as_completed(futs):
            m = fut.result()
            results.append(m)
            line = (
                f"{m['tag']}\t{m['rc']}\t{m['wall_s']:.1f}\t{m['revs']:.4f}\t"
                f"{m['r0']:.4f}\t{m['rmin']:.4f}\t{m['rmax']:.4f}\t{m['rmean']:.4f}\t"
                f"{m['band']}\t{m['alive']}\n"
            )
            with open(summary, "a") as fp:
                fp.write(line)
            print(
                f"  done {m['tag']}: revs={m['revs']:.3f} r=[{m['rmin']:.2f},{m['rmax']:.2f}] "
                f"band={m['band']} wall={m['wall_s']:.0f}s",
                flush=True,
            )
    bands = [m for m in results if m.get("band")]
    print(f"\nE1b band PASS count: {len(bands)} / {len(results)}", flush=True)
    for m in sorted(results, key=lambda x: -x.get("revs", 0))[:12]:
        print(
            f"  top revs {m['tag']}: {m['revs']:.3f} band={m['band']} "
            f"r=[{m['rmin']:.2f},{m['rmax']:.2f}]",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
