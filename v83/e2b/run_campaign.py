#!/usr/bin/env python3
"""E2 parallel: C1 a_rel/sep map, C2 flavored interlock, C3 cold large-D infall. Skip C4."""
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
GENF = ROOT / "bin/gen_qball_flavored"
TRACK = ROOT / "bin/sfa_qball_track"
PROF = ROOT / "v74/profiles/f_w142_g005.txt"
INTER = ROOT / "v71/results/interlock_profile.txt"
OUT = ROOT / "v83/e2b"
SEED = OUT / "seeds"
CFG = OUT / "cfg"
RES = OUT / "results"
OMP = os.environ.get("E2_OMP", "4")
MAX_WORKERS = int(os.environ.get("E2_WORKERS", "4"))
G = 0.05
PI = math.pi


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def gen_pair(path: Path, D: float, phases_a, phases_b, profile: Path, ws=None, L: float | None = None) -> None:
    """Two balls at ±D/2 on x; phases are (d0,d1,d2) each."""
    SEED.mkdir(parents=True, exist_ok=True)
    if path.exists():
        return
    if ws is None:
        ws = (1.42, 1.42, 1.42)
    if L is None:
        L = max(14.0, D + 4.0)
    # gen_qball_flavored N L out profile w0 w1 w2 x y z d0 d1 d2 vx [more]
    half = 0.5 * D
    args = [
        str(GENF),
        "48",
        str(L),
        str(path),
        str(profile),
        str(ws[0]),
        str(ws[1]),
        str(ws[2]),
        str(-half),
        "0",
        "0",
        str(phases_a[0]),
        str(phases_a[1]),
        str(phases_a[2]),
        "0",
        str(profile),
        str(ws[0]),
        str(ws[1]),
        str(ws[2]),
        str(half),
        "0",
        "0",
        str(phases_b[0]),
        str(phases_b[1]),
        str(phases_b[2]),
        "0",
    ]
    run(args)


def write_cfg(tag: str, seed: Path, T: float, L: float = 14.0, N: int = 48) -> Path:
    CFG.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    cfg = CFG / f"{tag}.cfg"
    cfg.write_text(
        "\n".join(
            [
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
                f"output={RES / f'out_{tag}.sfa'}",
                f"diag_file={RES / f'diag_{tag}.tsv'}",
                "snap_dt=5",
                "diag_dt=1",
                "precision=f32",
                "qdiag_radius=10",
            ]
        )
        + "\n"
    )
    return cfg


def analyze_clusters(tag: str) -> dict:
    sfa = RES / f"out_{tag}.sfa"
    tsv = RES / f"clusters_{tag}.tsv"
    out = {
        "tag": tag,
        "n0": -1,
        "nend": -1,
        "sep0": float("nan"),
        "sepend": float("nan"),
        "park": 0,
        "merge": 0,
        "separate": 0,
    }
    if not sfa.exists():
        return out
    try:
        subprocess.run(
            [str(TRACK), str(sfa), "0.25", "--tsv", str(tsv)],
            check=True,
            capture_output=True,
            timeout=120,
        )
    except Exception:
        return out
    if not tsv.exists():
        return out
    import numpy as np

    d = np.genfromtxt(tsv, names=True)
    if d.ndim == 0 or len(d) == 0:
        return out
    times = np.unique(d["t"])
    t0, t1 = times[0], times[-1]
    r0 = d[d["t"] == t0]
    r1 = d[d["t"] == t1]
    out["n0"] = int(r0["n_clusters"][0]) if "n_clusters" in d.dtype.names else len(r0)
    out["nend"] = int(r1["n_clusters"][0]) if "n_clusters" in d.dtype.names else len(r1)

    def sep(rows):
        if len(rows) < 2:
            return 0.0
        # two largest Q
        if "Q" in rows.dtype.names:
            rows = np.sort(rows, order="Q")
        a, b = rows[-2], rows[-1]
        return float(math.sqrt((a["cx"] - b["cx"]) ** 2 + (a["cy"] - b["cy"]) ** 2 + (a["cz"] - b["cz"]) ** 2))

    out["sep0"] = sep(r0) if out["n0"] >= 2 else 0.0
    out["sepend"] = sep(r1) if out["nend"] >= 2 else 0.0
    if out["nend"] == 1 and out["n0"] >= 2:
        out["merge"] = 1
    elif out["nend"] >= 2 and out["sepend"] > 1.2 * max(out["sep0"], 1e-6):
        out["separate"] = 1
    elif out["nend"] >= 2 and abs(out["sepend"] - out["sep0"]) < 0.25 * out["sep0"] and out["sepend"] > 4.0:
        out["park"] = 1  # crude static-ish
    return out


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
                timeout=int(os.environ.get("E2_TIMEOUT", "3600")),
            )
        rc = r.returncode
    except subprocess.TimeoutExpired:
        rc = -9
    wall = time.time() - t0
    m = analyze_clusters(tag)
    m.update({"rc": rc, "wall_s": wall})
    return m


def jobs() -> list[tuple[str, Path]]:
    js = []
    phases = {
        "co": ((0, 0, 0), (0, 0, 0)),
        "anti": ((0, 0, 0), (PI, PI, PI)),
        "mix": ((PI, 0, 0), (0, 0, 0)),  # 1-anti on left ball only
    }

    # --- C1: force/sep map D x phase, short T ---
    for D in (6, 8, 10, 12, 14, 16):
        for name, (pa, pb) in phases.items():
            tag = f"C1_{name}_D{D}"
            Lbox = max(14.0, float(D) + 4.0)
            seed = SEED / f"pair48_{name}_D{D}.sfa"
            gen_pair(seed, float(D), pa, pb, PROF, L=Lbox)
            cfg = write_cfg(tag, seed, T=40.0, L=Lbox, N=48)
            js.append((tag, cfg))

    # --- C3: cold large-D, longer T ---
    for D in (14, 16):
        for name in ("co", "mix", "anti"):
            pa, pb = phases[name]
            tag = f"C3_{name}_D{D}"
            Lbox = max(14.0, float(D) + 4.0)
            seed = SEED / f"pair48_{name}_D{D}.sfa"
            gen_pair(seed, float(D), pa, pb, PROF, L=Lbox)
            cfg = write_cfg(tag, seed, T=100.0, L=Lbox, N=48)
            js.append((tag, cfg))

    # --- C2: flavored interlock profile two-low/one-high ---
    prof = INTER if INTER.exists() else PROF
    ws = (1.38, 1.38, 1.42) if INTER.exists() else (1.42, 1.42, 1.42)
    for D in (8, 10, 12, 14):
        for name, (pa, pb) in {
            "im1": ((PI, PI, 0), (0, 0, 0)),
            "im2": ((PI, 0, 0), (0, 0, 0)),
            "co": ((0, 0, 0), (0, 0, 0)),
        }.items():
            tag = f"C2_{name}_D{D}"
            Lbox = max(14.0, float(D) + 4.0)
            seed = SEED / f"c2_48_{name}_D{D}.sfa"
            gen_pair(seed, float(D), pa, pb, prof, ws=ws, L=Lbox)
            cfg = write_cfg(tag, seed, T=80.0, L=Lbox, N=48)
            js.append((tag, cfg))

    return js


def main() -> int:
    RES.mkdir(parents=True, exist_ok=True)
    js = jobs()
    seen = set()
    uniq = []
    for tag, cfg in js:
        if tag in seen:
            continue
        seen.add(tag)
        uniq.append((tag, cfg))
    print(f"E2 campaign: {len(uniq)} jobs, workers={MAX_WORKERS}, OMP={OMP}", flush=True)
    summary = RES / "summary_e2.tsv"
    with open(summary, "w") as fp:
        fp.write("tag\trc\twall_s\tn0\tnend\tsep0\tsepend\tpark\tmerge\tseparate\n")
    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(run_one, tag, cfg): tag for tag, cfg in uniq}
        for fut in as_completed(futs):
            m = fut.result()
            results.append(m)
            line = (
                f"{m['tag']}\t{m['rc']}\t{m['wall_s']:.1f}\t{m['n0']}\t{m['nend']}\t"
                f"{m['sep0']:.4f}\t{m['sepend']:.4f}\t{m['park']}\t{m['merge']}\t{m['separate']}\n"
            )
            with open(summary, "a") as fp:
                fp.write(line)
            print(
                f"  done {m['tag']}: n={m['n0']}->{m['nend']} sep={m['sep0']:.2f}->{m['sepend']:.2f} "
                f"park={m['park']} merge={m['merge']} sepflag={m['separate']} wall={m['wall_s']:.0f}s",
                flush=True,
            )
    parks = [m for m in results if m.get("park")]
    merges = [m for m in results if m.get("merge")]
    seps = [m for m in results if m.get("separate")]
    print(
        f"\nE2: park={len(parks)} merge={len(merges)} separate={len(seps)} / {len(results)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
