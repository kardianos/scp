#!/usr/bin/env python3
"""Generate B4 orbit seeds + multi-Z seeds/cfgs for GPU campaign.

Does NOT run sims — writes seeds and .cfg under work dir for scp-runner upload.
"""
from __future__ import annotations
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "bin" / "gen_mf_shell_orbit"
GEN_MULTI = ROOT / "bin" / "gen_qball_multi"
PROF142 = ROOT / "v74" / "profiles" / "f_w142_g005.txt"
PROF146 = ROOT / "v74" / "profiles" / "f_w146_g005.txt"

N, LBOX = 192, 48.0
R = 18.0
# v_circ ≈ 0.071 for single L around C
VC = 0.071


def mf_cfg(tag: str, seed_c: str, seed_l: str, T: float = 400.0, snap: float = 20.0) -> str:
    return f"""# {tag}
N = {N}
L = {LBOX}
T = {T}
dt_factor = 0.025
m = 1.5
m_theta = 1.6
eta = 0
complex_phi = 1
complex_gauge = 1
g_gauge = 0.05
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
init_sfa = {seed_c}
init_sfa_L = {seed_l}
output = {tag}.sfa
diag_file = {tag}_diag.tsv
precision = 0
snap_dt = {snap}
diag_dt = 1.0
"""


def shell_orbit(work: Path, tag: str, n_L: int, vt: float, oC=1.42, oL=1.46, R=18.0):
    sc = work / f"{tag}_C.sfa"
    sl = work / f"{tag}_L.sfa"
    cmd = [
        str(GEN), str(N), str(LBOX), str(R),
        str(PROF142), str(oC), str(PROF146), str(oL),
        str(n_L), str(vt), str(sc), str(sl),
    ]
    print("RUN", " ".join(cmd[:8]), f"... nL={n_L} vt={vt}")
    subprocess.check_call(cmd)
    cfg = work / f"{tag}.cfg"
    cfg.write_text(mf_cfg(tag, sc.name, sl.name))
    return cfg


def multi_z6(work: Path, tag: str, R_shell: float = 22.0):
    """Z6 octa D=10 ω=1.46 + L6 shell R rest."""
    sc = work / f"{tag}_C.sfa"
    sl = work / f"{tag}_L.sfa"
    a = 10.0 / math.sqrt(2.0)
    # C octa
    cmd = [str(GEN_MULTI), str(N), str(LBOX), str(sc)]
    for p in [(a,0,0),(-a,0,0),(0,a,0),(0,-a,0),(0,0,a),(0,0,-a)]:
        cmd += [str(PROF146), "1.46", "0", f"{p[0]:.10f}", f"{p[1]:.10f}", f"{p[2]:.10f}"]
    print("RUN gen_qball_multi C octa...")
    subprocess.check_call(cmd)
    # L shell via gen_mf_shell_orbit with n_C dummy — use shell_orbit but only L;
    # easier: gen_mf_shell_orbit writes both; for multi-Z C we already have sc
    # generate L only via shell with throwaway C
    tmp = work / f"{tag}_tmp"
    tmp_c = work / f"{tag}_tmp_C.sfa"
    cmd = [
        str(GEN), str(N), str(LBOX), str(R_shell),
        str(PROF146), "1.46", str(PROF146), "1.46",
        "6", "0", str(tmp_c), str(sl),
    ]
    subprocess.check_call(cmd)
    tmp_c.unlink(missing_ok=True)
    cfg = work / f"{tag}.cfg"
    cfg.write_text(mf_cfg(tag, sc.name, sl.name, T=400, snap=50))
    return cfg


def main():
    work = Path(sys.argv[1] if len(sys.argv) > 1 else "/space/scp/v75/b4o_mz")
    work.mkdir(parents=True, exist_ok=True)
    if not GEN.is_file():
        print("FATAL: build gen_mf_shell_orbit first", file=sys.stderr)
        return 1
    # Orbit: single L pair-like for clean D(t) + shell co-rot
    for name, nL, frac in [
        ("b4o_pair_sub", 1, 0.7),
        ("b4o_pair_vc", 1, 1.0),
        ("b4o_pair_super", 1, 1.3),
        ("b4o_shell_vc", 6, 1.0),
    ]:
        shell_orbit(work, name, nL, VC * frac)
    # Multi-Z parked packaging
    multi_z6(work, "mz2_z6_L6_R22", R_shell=22.0)
    print(f"Wrote seeds+cfg under {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
