#!/usr/bin/env python3
"""
v75 Stage-1 self-tune controller — B1 soft θ, frozen action per trial,
predictive backtracking, JSONL ledger. NO kernel edits.

See v75/SELF_TUNE_C.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PROFILES = ROOT / "v74" / "profiles"
GEN_MULTI = ROOT / "bin" / "gen_qball_multi"
GEN_PAIR = ROOT / "bin" / "gen_mf_pair_boost"
TRACK = ROOT / "bin" / "mf_pair_track"
SIM_CPU = ROOT / "bin" / "scp_sim"

# Profile map: omega -> filename
PROFILE_MAP = {
    1.42: "f_w142_g005.txt",
    1.46: "f_w146_g005.txt",
    1.485: "f_w1485_g005.txt",
}


@dataclass
class Theta:
    n_C: int
    omega_C: float
    n_L: int
    omega_L: float
    R_shell: float
    D_nuc: float = 10.0
    g_gauge: float = 0.05
    N: int = 192
    L: float = 48.0

    def key(self) -> str:
        return (
            f"nC{self.n_C}_wC{self.omega_C:g}_nL{self.n_L}_"
            f"wL{self.omega_L:g}_R{self.R_shell:g}"
        )


@dataclass
class Node:
    trial_id: str
    parent_id: Optional[str]
    theta: Theta
    note: str = ""
    children: list[str] = field(default_factory=list)


def octahedron_centers(edge: float) -> list[tuple[float, float, float]]:
    a = edge / math.sqrt(2.0)
    return [
        (a, 0.0, 0.0),
        (-a, 0.0, 0.0),
        (0.0, a, 0.0),
        (0.0, -a, 0.0),
        (0.0, 0.0, a),
        (0.0, 0.0, -a),
    ]


def tetrahedron_centers(edge: float) -> list[tuple[float, float, float]]:
    """4 vertices of regular tetrahedron with edge `edge`, centered at origin."""
    # Regular tetrahedron before scaling
    raw = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    md = math.dist(raw[0], raw[1])
    s = edge / md
    pts = [(s * x, s * y, s * z) for x, y, z in raw]
    cx = sum(p[0] for p in pts) / 4
    cy = sum(p[1] for p in pts) / 4
    cz = sum(p[2] for p in pts) / 4
    return [(x - cx, y - cy, z - cz) for x, y, z in pts]


def pair_centers(sep: float) -> list[tuple[float, float, float]]:
    return [(sep / 2, 0.0, 0.0), (-sep / 2, 0.0, 0.0)]


def shell_centers(n: int, R: float) -> list[tuple[float, float, float]]:
    """n points on a sphere of radius R (octa / tetra / pair / poles+eq)."""
    if n == 1:
        return [(R, 0.0, 0.0)]
    if n == 2:
        return pair_centers(2 * R)  # separation 2R → each at distance R from origin? pair at ±R
    if n == 4:
        # tetra on sphere radius R
        raw = tetrahedron_centers(1.0)
        norms = [math.sqrt(x * x + y * y + z * z) for x, y, z in raw]
        return [(R * x / nrm, R * y / nrm, R * z / nrm) for (x, y, z), nrm in zip(raw, norms)]
    if n == 6:
        return [
            (R, 0.0, 0.0),
            (-R, 0.0, 0.0),
            (0.0, R, 0.0),
            (0.0, -R, 0.0),
            (0.0, 0.0, R),
            (0.0, 0.0, -R),
        ]
    raise ValueError(f"unsupported n_L/n_C={n}")


def nuclear_centers(n: int, D_nuc: float) -> list[tuple[float, float, float]]:
    if n == 1:
        return [(0.0, 0.0, 0.0)]
    if n == 2:
        return pair_centers(D_nuc)
    if n == 4:
        return tetrahedron_centers(D_nuc)
    if n == 6:
        return octahedron_centers(D_nuc)
    raise ValueError(f"unsupported n_C={n}")


def profile_path(omega: float, profiles_dir: Path) -> Path:
    if omega not in PROFILE_MAP:
        raise KeyError(f"no profile for omega={omega}; have {list(PROFILE_MAP)}")
    p = profiles_dir / PROFILE_MAP[omega]
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def write_multi_seed(
    out: Path,
    N: int,
    L: float,
    omega: float,
    centers: list[tuple[float, float, float]],
    profiles_dir: Path,
) -> None:
    prof = profile_path(omega, profiles_dir)
    cmd = [str(GEN_MULTI), str(N), str(L), str(out)]
    for x, y, z in centers:
        cmd += [str(prof), f"{omega:.6f}", "0", f"{x:.10f}", f"{y:.10f}", f"{z:.10f}"]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_diag(path: Path) -> list[dict[str, float]]:
    """Parse scp_sim diag TSV (tab-separated)."""
    if not path.is_file():
        return []
    rows: list[dict[str, float]] = []
    with open(path) as f:
        header: Optional[list[str]] = None
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if header is None:
                if parts[0] == "t" or not _is_float(parts[0]):
                    header = parts
                    continue
                header = []  # no header row
            if not _is_float(parts[0]):
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) < 18:
                continue
            rows.append({
                "t": vals[0],
                "E_total": vals[9],
                "Q_phi": vals[14],
                "s_max": vals[17],
                "gauss_max": vals[25] if len(vals) > 25 else 0.0,
            })
    return rows


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def parse_track(path: Path) -> list[dict[str, float]]:
    if not path.is_file():
        return []
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if parts[0] == "frame":
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            # frame t massC massL Qc Ql cxC cyC czC cxL cyL czL D
            if len(vals) < 7:
                continue
            rows.append({
                "frame": vals[0],
                "t": vals[1],
                "massC": vals[2],
                "massL": vals[3],
                "Qc": vals[4],
                "Ql": vals[5],
                "D": vals[12] if len(vals) > 12 else float("nan"),
            })
    return rows


def cost_from_metrics(m: dict[str, float]) -> tuple[float, dict[str, float]]:
    Qc0 = max(m.get("Qc0", 0) or 0, 1e-12)
    mL0 = max(m.get("massL0", 0) or 0, 1e-12)
    E0 = max(abs(m.get("E0", 0) or 0), 1e-12)
    c_Q = min(max((Qc0 - m.get("Qc", Qc0)) / Qc0, 0.0), 1.0)
    c_L = min(max((mL0 - m.get("massL", mL0)) / mL0, 0.0), 1.0)
    c_E = min(max(abs(m.get("E", E0) - E0) / E0 / 0.20, 0.0), 1.0)
    c_G = 1.0 if m.get("gauss_max", 0) > 1e-10 else 0.0
    c = 0.35 * c_Q + 0.35 * c_L + 0.20 * c_E + 0.10 * c_G
    return c, {"c_Q": c_Q, "c_L": c_L, "c_E": c_E, "c_G": c_G, "cost": c}


def predict_from_diag(rows: list[dict[str, float]]) -> Optional[str]:
    if not rows:
        return "PRED_NODIAG"
    r0, r1 = rows[0], rows[-1]
    Q0 = abs(r0.get("Q_phi", 0)) or 1e-12
    E0 = abs(r0.get("E_total", 0)) or 1e-12
    if abs(r1.get("Q_phi", 0)) / Q0 < 0.80:
        return "PRED_Q"
    if abs(r1.get("E_total", 0)) / E0 < 0.88:
        return "PRED_E"
    s_peak = max(r.get("s_max", 0) for r in rows)
    if s_peak > 1.5:
        return "PRED_S"
    g_peak = max(abs(r.get("gauss_max", 0)) for r in rows)
    if g_peak > 1e-10:
        return "PRED_GAUSS"
    return None


def predict_from_track(rows: list[dict[str, float]]) -> Optional[str]:
    if not rows:
        return "PRED_NOTRACK"
    m0, m1 = rows[0], rows[-1]
    mL0 = m0.get("massL", 0) or 1e-12
    if m1.get("massL", 0) / mL0 < 0.05:
        return "FAIL_L_ABSORB"
    if (mL0 - m1.get("massL", 0)) / mL0 > 0.25:
        return "PRED_L"
    return None


def build_tree(N: int, L: float) -> list[Node]:
    """Ordered depth-first-friendly node list with parent links."""
    def th(**kw):
        base = dict(N=N, L=L, D_nuc=10.0, g_gauge=0.05)
        base.update(kw)
        return Theta(**base)

    nodes = [
        Node("B1", None, th(n_C=6, omega_C=1.46, n_L=6, omega_L=1.46, R_shell=18),
             "parked c6_light + L6 R18"),
        Node("B1a", "B1", th(n_C=6, omega_C=1.46, n_L=6, omega_L=1.46, R_shell=22),
             "larger shell"),
        Node("B1b", "B1", th(n_C=6, omega_C=1.46, n_L=6, omega_L=1.485, R_shell=22),
             "soft L + large shell"),
        Node("B2", None, th(n_C=4, omega_C=1.46, n_L=4, omega_L=1.46, R_shell=18),
             "Z4 + L4"),
        Node("B2a", "B2", th(n_C=4, omega_C=1.46, n_L=4, omega_L=1.46, R_shell=22),
             "Z4 larger shell"),
        Node("B3", None, th(n_C=2, omega_C=1.46, n_L=2, omega_L=1.46, R_shell=18),
             "Z2 + L2"),
        Node("B3a", "B3", th(n_C=2, omega_C=1.46, n_L=2, omega_L=1.46, R_shell=22),
             "Z2 larger shell"),
        Node("B4", None, th(n_C=1, omega_C=1.42, n_L=6, omega_L=1.46, R_shell=18),
             "single heavy + L6"),
        Node("B4a", "B4", th(n_C=1, omega_C=1.42, n_L=6, omega_L=1.485, R_shell=22),
             "single heavy + soft L6 large R"),
        Node("B5", None, th(n_C=6, omega_C=1.46, n_L=2, omega_L=1.46, R_shell=20),
             "c6 + few light"),
    ]
    return nodes


class Controller:
    def __init__(
        self,
        work: Path,
        sim: Path,
        profiles_dir: Path,
        T_screen: float,
        T_full: float,
        N: int,
        L: float,
        snap_dt_screen: float,
        snap_dt_full: float,
        dry_run: bool = False,
    ):
        self.work = work
        self.sim = sim
        self.profiles_dir = profiles_dir
        self.T_screen = T_screen
        self.T_full = T_full
        self.N = N
        self.L = L
        self.snap_dt_screen = snap_dt_screen
        self.snap_dt_full = snap_dt_full
        self.dry_run = dry_run
        self.work.mkdir(parents=True, exist_ok=True)
        self.ledger_path = work / "ledger.jsonl"
        self.summary_path = work / "ledger_summary.tsv"
        self.results: list[dict[str, Any]] = []
        if not self.summary_path.is_file():
            self.summary_path.write_text(
                "trial_id\tphase\toutcome\tcost\tc_Q\tc_L\tpred\taction\t"
                "n_C\tomega_C\tn_L\tomega_L\tR_shell\tmassL0\tmassL\tQc0\tQc\n"
            )

    def log(self, rec: dict[str, Any]) -> None:
        with open(self.ledger_path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        self.results.append(rec)
        th = rec.get("theta", {})
        m = rec.get("metrics", {})
        parts = rec.get("cost_parts", {})
        line = (
            f"{rec.get('trial_id')}\t{rec.get('phase')}\t{rec.get('outcome')}\t"
            f"{m.get('cost', float('nan')):.4f}\t{parts.get('c_Q', float('nan')):.4f}\t"
            f"{parts.get('c_L', float('nan')):.4f}\t{rec.get('pred')}\t{rec.get('action')}\t"
            f"{th.get('n_C')}\t{th.get('omega_C')}\t{th.get('n_L')}\t{th.get('omega_L')}\t"
            f"{th.get('R_shell')}\t{m.get('massL0', float('nan')):.4g}\t"
            f"{m.get('massL', float('nan')):.4g}\t{m.get('Qc0', float('nan')):.4g}\t"
            f"{m.get('Qc', float('nan')):.4g}\n"
        )
        with open(self.summary_path, "a") as f:
            f.write(line)
        print(
            f"  LEDGER {rec.get('trial_id')} {rec.get('phase')} "
            f"{rec.get('outcome')} cost={m.get('cost', float('nan')):.3f} "
            f"pred={rec.get('pred')} action={rec.get('action')}",
            flush=True,
        )

    def make_seeds(self, tag: str, th: Theta) -> tuple[Path, Path]:
        seed_C = self.work / f"{tag}_C.sfa"
        seed_L = self.work / f"{tag}_L.sfa"
        c_cent = nuclear_centers(th.n_C, th.D_nuc)
        l_cent = shell_centers(th.n_L, th.R_shell)
        write_multi_seed(seed_C, th.N, th.L, th.omega_C, c_cent, self.profiles_dir)
        write_multi_seed(seed_L, th.N, th.L, th.omega_L, l_cent, self.profiles_dir)
        return seed_C, seed_L

    def write_cfg(
        self,
        tag: str,
        th: Theta,
        seed_C: Path,
        seed_L: Path,
        T: float,
        snap_dt: float,
    ) -> Path:
        cfg = self.work / f"{tag}.cfg"
        out = self.work / f"{tag}.sfa"
        diag = self.work / f"{tag}_diag.tsv"
        text = f"""# auto self-tune {tag}
N = {th.N}
L = {th.L}
T = {T}
dt_factor = 0.025
m = 1.5
m_theta = 1.6
eta = 0
complex_phi = 1
complex_gauge = 1
g_gauge = {th.g_gauge}
n_fabrics = 3
mf_lock_CQ = 1
mf_stage = 1
q_C = 0
q_Q = 1
q_L = -1
bc_type = 0
damp_width = 4.0
damp_rate = 0.01
init = sfa
init_sfa = {seed_C.name}
init_sfa_L = {seed_L.name}
output = {out.name}
diag_file = {diag.name}
precision = 0
snap_dt = {snap_dt}
diag_dt = 1.0
"""
        cfg.write_text(text)
        return cfg

    def run_sim(self, cfg: Path) -> int:
        if self.dry_run:
            print(f"  [dry-run] would run {self.sim} {cfg}")
            return 0
        env = dict(os.environ)
        env.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 4))
        # Run with cwd=work so relative seed/output paths resolve.
        # IMPORTANT: do not use capture_output=True — long GPU runs fill the
        # pipe buffer and deadlock the child (no diag progress).
        log_path = cfg.with_suffix(".simlog")
        with open(log_path, "w") as logf:
            r = subprocess.run(
                [str(self.sim), str(cfg.name)],
                cwd=str(self.work),
                env=env,
                stdout=logf,
                stderr=subprocess.STDOUT,
            )
        if r.returncode != 0:
            tail = ""
            try:
                tail = log_path.read_text()[-2000:]
            except OSError:
                pass
            print(f"  SIM FAIL rc={r.returncode}\n{tail}", file=sys.stderr)
        return r.returncode

    def track(self, sfa: Path, out_tsv: Path) -> None:
        if self.dry_run:
            return
        if not sfa.is_file():
            print(f"  WARN: missing SFA {sfa}", file=sys.stderr)
            return
        subprocess.run(
            [str(TRACK), str(sfa), str(out_tsv)],
            check=False,
            capture_output=True,
        )

    def evaluate_trial(
        self,
        node: Node,
        phase: str,
        T: float,
        snap_dt: float,
    ) -> dict[str, Any]:
        tag = f"{node.trial_id}_{phase}"
        th = node.theta
        th.N, th.L = self.N, self.L
        t0 = time.time()
        print(f"\n=== {tag}  {node.note}  θ={th.key()}  T={T} ===", flush=True)

        seed_C, seed_L = self.make_seeds(tag, th)
        cfg = self.write_cfg(tag, th, seed_C, seed_L, T, snap_dt)
        cfg_hash = hashlib.sha256(cfg.read_bytes()).hexdigest()[:12]

        rc = self.run_sim(cfg)
        wall = time.time() - t0
        diag_path = self.work / f"{tag}_diag.tsv"
        sfa_path = self.work / f"{tag}.sfa"
        track_path = self.work / f"{tag}_track.tsv"

        if rc != 0 and not self.dry_run:
            rec = {
                "trial_id": tag,
                "parent_id": node.parent_id,
                "phase": phase,
                "t_sim": T,
                "theta": asdict(th),
                "metrics": {},
                "cost_parts": {},
                "pred": "SIM_CRASH",
                "outcome": "FAIL_SIM",
                "action": "reject",
                "cfg_hash": cfg_hash,
                "seed_note": node.note,
                "wall_s": wall,
            }
            self.log(rec)
            return rec

        diag_rows = parse_diag(diag_path)
        pred_d = predict_from_diag(diag_rows)

        self.track(sfa_path, track_path)
        track_rows = parse_track(track_path)
        pred_t = predict_from_track(track_rows)

        metrics: dict[str, float] = {}
        if diag_rows:
            metrics["E0"] = diag_rows[0]["E_total"]
            metrics["E"] = diag_rows[-1]["E_total"]
            metrics["Qc0_diag"] = diag_rows[0]["Q_phi"]
            metrics["Qc_diag"] = diag_rows[-1]["Q_phi"]
            metrics["s_max_peak"] = max(r["s_max"] for r in diag_rows)
            metrics["gauss_max"] = max(abs(r["gauss_max"]) for r in diag_rows)
        if track_rows:
            metrics["massC0"] = track_rows[0]["massC"]
            metrics["massC"] = track_rows[-1]["massC"]
            metrics["massL0"] = track_rows[0]["massL"]
            metrics["massL"] = track_rows[-1]["massL"]
            metrics["Qc0"] = track_rows[0]["Qc"]
            metrics["Qc"] = track_rows[-1]["Qc"]
            metrics["Ql0"] = track_rows[0]["Ql"]
            metrics["Ql"] = track_rows[-1]["Ql"]
        else:
            # fallback diag-only Q
            metrics["Qc0"] = metrics.get("Qc0_diag", 0)
            metrics["Qc"] = metrics.get("Qc_diag", 0)
            metrics["massL0"] = 1.0
            metrics["massL"] = 1.0

        cost, parts = cost_from_metrics(metrics)
        metrics["cost"] = cost

        pred = pred_d or pred_t or "none"
        # Outcome classification
        if pred == "PRED_GAUSS" or parts["c_G"] > 0:
            outcome = "FAIL_GAUSS"
        elif pred in ("PRED_Q",) or parts["c_Q"] > 0.20:
            outcome = "FAIL_Q"
        elif pred in ("PRED_L", "FAIL_L_ABSORB") or parts["c_L"] > 0.25:
            outcome = "FAIL_L"
        elif pred in ("PRED_E",) or parts["c_E"] > 0.9:
            outcome = "FAIL_E"
        elif pred == "PRED_S":
            outcome = "FAIL_S"
        elif (
            phase == "full"
            and cost <= 0.15
            and parts["c_Q"] <= 0.15
            and parts["c_L"] <= 0.15
            and parts["c_G"] == 0
        ):
            outcome = "PASS"
        elif phase == "screen" and cost <= 0.35 and pred == "none":
            outcome = "PARTIAL"  # promote
        elif cost <= 0.25:
            outcome = "PARTIAL"
        else:
            outcome = "FAIL_COST"

        if outcome == "PASS":
            action = "stop_success"
        elif phase == "screen" and outcome == "PARTIAL":
            action = "promote_full"
        elif outcome.startswith("FAIL"):
            action = "backtrack"
        else:
            action = "reject"

        rec = {
            "trial_id": tag,
            "parent_id": node.parent_id,
            "phase": phase,
            "t_sim": T,
            "theta": asdict(th),
            "metrics": metrics,
            "cost_parts": parts,
            "pred": pred,
            "outcome": outcome,
            "action": action,
            "cfg_hash": cfg_hash,
            "seed_note": node.note,
            "wall_s": wall,
        }
        self.log(rec)

        # Free large SFA after tracking to save disk (keep diag/track)
        if sfa_path.is_file() and sfa_path.stat().st_size > 50_000_000:
            try:
                sfa_path.unlink()
                print(f"  deleted large SFA {sfa_path.name}", flush=True)
            except OSError:
                pass
        return rec

    def should_skip_child(self, parent_rec: Optional[dict], child: Node) -> bool:
        """Predictive prune: skip children that cannot fix the parent's wall."""
        if not parent_rec:
            return False
        oc = parent_rec.get("outcome", "")
        th = child.theta
        pth = parent_rec.get("theta", {})
        # If parent nuclear-failed and child has same or higher nuclear Q load, skip
        if oc in ("FAIL_Q", "FAIL_S", "FAIL_E"):
            # child must reduce n_C or raise omega_C or (same n_C but we still allow R variants for L)
            if th.n_C > pth.get("n_C", 99):
                return True
            if th.n_C == pth.get("n_C") and th.omega_C < pth.get("omega_C", 0):
                return True
        if oc == "FAIL_L":
            # need larger R or softer L
            if th.R_shell <= pth.get("R_shell", 0) and th.omega_L <= pth.get("omega_L", 0):
                if th.n_L >= pth.get("n_L", 0) and th.n_C >= pth.get("n_C", 0):
                    return True
        return False

    def run_campaign(self) -> str:
        nodes = build_tree(self.N, self.L)
        by_id = {n.trial_id: n for n in nodes}
        parent_screen: dict[str, dict] = {}
        best: Optional[dict] = None
        any_pass = False
        evaluated = 0

        print("=" * 64)
        print("v75 SELF_TUNE Stage 1 — B1 soft θ + predictive backtrack")
        print(f"  work={self.work}")
        print(f"  sim={self.sim}")
        print(f"  N={self.N} L={self.L} T_screen={self.T_screen} T_full={self.T_full}")
        print(f"  tree nodes={len(nodes)}")
        print("=" * 64)

        for node in nodes:
            # prune against parent screen outcome
            pref = parent_screen.get(node.parent_id) if node.parent_id else None
            if node.parent_id and pref and self.should_skip_child(pref, node):
                rec = {
                    "trial_id": f"{node.trial_id}_screen",
                    "parent_id": node.parent_id,
                    "phase": "screen",
                    "t_sim": 0,
                    "theta": asdict(node.theta),
                    "metrics": {},
                    "cost_parts": {},
                    "pred": "PRUNED",
                    "outcome": "SKIP",
                    "action": "backtrack",
                    "cfg_hash": "",
                    "seed_note": f"pruned vs parent {node.parent_id}",
                    "wall_s": 0,
                }
                self.log(rec)
                continue

            srec = self.evaluate_trial(node, "screen", self.T_screen, self.snap_dt_screen)
            evaluated += 1
            parent_screen[node.trial_id] = srec
            if best is None or (srec.get("metrics", {}).get("cost", 9) < best.get("metrics", {}).get("cost", 9)):
                if srec.get("metrics"):
                    best = srec

            if srec["action"] == "promote_full":
                frec = self.evaluate_trial(node, "full", self.T_full, self.snap_dt_full)
                evaluated += 1
                if best is None or (frec.get("metrics", {}).get("cost", 9) < best.get("metrics", {}).get("cost", 9)):
                    if frec.get("metrics"):
                        best = frec
                if frec["outcome"] == "PASS":
                    any_pass = True
                    self._write_verdict("DEFINITIVE_SUCCESS", frec, evaluated)
                    return "DEFINITIVE_SUCCESS"

            # free seed files for this node (regen if needed)
            for p in self.work.glob(f"{node.trial_id}_*"):
                if p.suffix == ".sfa" and p.stat().st_size > 10_000_000:
                    try:
                        p.unlink()
                    except OSError:
                        pass

        if any_pass:
            self._write_verdict("DEFINITIVE_SUCCESS", best, evaluated)
            return "DEFINITIVE_SUCCESS"

        # Tree exhausted
        self._write_verdict("DEFINITIVE_FAIL_SOFT_THETA", best, evaluated)
        return "DEFINITIVE_FAIL_SOFT_THETA"

    def _write_verdict(self, verdict: str, best: Optional[dict], n: int) -> None:
        path = self.work / "VERDICT.md"
        lines = [
            f"# Self-tune Stage 1 verdict: **{verdict}**\n",
            f"\nEvaluated trials (screen+full): {n}\n",
            f"Ledger: `{self.ledger_path}`\n",
            f"Summary: `{self.summary_path}`\n",
        ]
        if best:
            lines.append("\n## Best trial\n")
            lines.append("```json\n")
            lines.append(json.dumps(best, indent=2)[:4000])
            lines.append("\n```\n")
        if verdict == "DEFINITIVE_FAIL_SOFT_THETA":
            lines.append(
                "\n## Interpretation\n\n"
                "All soft-θ branches failed packaging criteria. Isolation may still "
                "hold (nonzero massL) but atom packaging is not achieved by B1 seed "
                "knobs alone. Next design moves per SELF_TUNE_C.md §8: sub-Q_max / "
                "lower-g nuclear redesign, BO relax, or Stage 2 Option C if the wall "
                "is engagement-shaped.\n"
            )
        else:
            lines.append(
                "\n## Interpretation\n\n"
                "Found soft θ that meets PASS at T_full under frozen B1 action. "
                "Freeze this θ; Option C not required for this scale claim.\n"
            )
        path.write_text("".join(lines))
        print("\n" + "=" * 64)
        print(f"VERDICT: {verdict}")
        print(f"  wrote {path}")
        print("=" * 64)


def main() -> int:
    ap = argparse.ArgumentParser(description="v75 B1 soft-θ self-tune controller")
    ap.add_argument("--mode", choices=("local", "gpu"), default="local")
    ap.add_argument("--sim", type=Path, default=None, help="path to scp_sim binary")
    ap.add_argument("--work", type=Path, default=Path("/space/scp/v75/self_tune"))
    ap.add_argument("--profiles", type=Path, default=DEFAULT_PROFILES)
    ap.add_argument("--gen-multi", type=Path, default=None, help="gen_qball_multi binary")
    ap.add_argument("--track", type=Path, default=None, help="mf_pair_track binary")
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--L", type=float, default=None)
    ap.add_argument("--T_screen", type=float, default=None)
    ap.add_argument("--T_full", type=float, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", type=str, default=None, help="comma trial ids to run only")
    args = ap.parse_args()

    global GEN_MULTI, TRACK
    if args.gen_multi:
        GEN_MULTI = args.gen_multi.resolve()
    if args.track:
        TRACK = args.track.resolve()

    if args.mode == "local":
        N = args.N or 64
        L = args.L or 20.0
        T_screen = args.T_screen or 20.0
        T_full = args.T_full or 40.0
        snap_s, snap_f = max(T_screen / 4, 5.0), max(T_full / 4, 5.0)
        sim = args.sim or SIM_CPU
        work = args.work if args.work != Path("/space/scp/v75/self_tune") else Path("/tmp/scp_self_tune_smoke")
    else:
        N = args.N or 192
        L = args.L or 48.0
        T_screen = args.T_screen or 150.0
        T_full = args.T_full or 400.0
        # Sparse snaps: track only needs a few frames; saves disk/IO on N=192 MF
        snap_s, snap_f = max(T_screen / 3.0, 50.0), max(T_full / 8.0, 50.0)
        sim = args.sim or Path("./scp_sim_mf_cuda")
        work = args.work

    for req in (GEN_MULTI, TRACK):
        if not req.is_file():
            print(f"FATAL: missing {req}", file=sys.stderr)
            return 1
    if not args.dry_run and not Path(sim).is_file():
        print(f"FATAL: sim binary not found: {sim}", file=sys.stderr)
        return 1

    ctl = Controller(
        work=work,
        sim=Path(sim).resolve() if Path(sim).is_file() else Path(sim),
        profiles_dir=args.profiles,
        T_screen=T_screen,
        T_full=T_full,
        N=N,
        L=L,
        snap_dt_screen=snap_s,
        snap_dt_full=snap_f,
        dry_run=args.dry_run,
    )

    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        # Patch module-level builder used by Controller.run_campaign
        global build_tree
        _orig = build_tree

        def _filtered(N: int, L: float):
            return [n for n in _orig(N, L) if n.trial_id in want]

        build_tree = _filtered  # type: ignore

    verdict = ctl.run_campaign()
    return 0 if verdict == "DEFINITIVE_SUCCESS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
