#!/usr/bin/env python3
"""
v65 X1 — homeostatic self-tuning controller (mechanism B), NO kernel edit.

The control law IS the dynamics: an outer loop wraps the unmodified CPU kernel and
self-tunes kappa toward the stability edge kappa_crit via SOC feedback:
  - run a short relaxation from a fixed seed at the current kappa,
  - sense collapse from the core density P_max (kernel diag),
  - COLLAPSE  (P_max NaN or > P_collapse) -> push kappa UP   (x up_factor)   [stability feedback]
  - STABLE                                -> drift kappa DOWN (x down_factor) [action pull, dE/dkappa>0]
The fixed point is the balance = the self-organized critical edge. Run from a HIGH and a
LOW start; both should converge to a common band around kappa_crit (E3 measured ~1),
basin-independent. That is the SOC attractor.

Autonomous: this script runs the whole loop; no hand-tuning, no GPU round-trips.
Refs: v65/SELF_TUNING.md (mechanism B, X1), v65/FINDINGS_E13.md (kappa_crit~1).
"""
import subprocess, os, math, sys

BIN   = "/home/d/code/scp/bin/scp_sim"
GEN   = "/home/d/code/scp/bin/gen_oscillon"
WORK  = "/tmp/x1"
N, L, T = 48, 10.0, 1.5
P_COLLAPSE = 1.5          # P_max above this (or NaN) => collapsing
UP, DOWN   = 1.30, 0.78   # SOC feedback factors (push up on collapse, drift down when stable)
STEPS      = 13

os.makedirs(WORK, exist_ok=True)
SEED = f"{WORK}/seed.sfa"
subprocess.run([GEN,"-o",SEED,"-N",str(N),"-L",str(L),"-A","1.0","-sigma","1.5",
                "-precision","f32"], check=True, capture_output=True)

def run_step(kappa, tag):
    """Run one relaxation at kappa; return P_max (inf if collapsed/NaN)."""
    cfg = f"{WORK}/{tag}.cfg"; diag = f"{WORK}/{tag}_diag.tsv"; out = f"{WORK}/{tag}.sfa"
    with open(cfg,"w") as f:
        f.write(f"""init=template
init_sfa={SEED}
T={T}
dt_factor=0.025
m=1.5
m_theta=0.0
eta=0.5
mu=-41.345
kappa={kappa}
mode=0
bc_type=0
damp_width=2.0
damp_rate=0.02
A_bg=0.0
delta=0,3.0005,4.4325
output={out}
precision=f16
snap_dt={T}
diag_dt=0.5
diag_file={diag}
""")
    env = dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count()))
    subprocess.run([BIN,cfg], check=False, capture_output=True, env=env)
    # parse last diag row; P_max is column 12 (1-based) = index 11
    pmax = float("inf")
    try:
        with open(diag) as f:
            rows = [r for r in f if r.strip() and not r.startswith("t\t") and not r[0].isalpha()]
        if rows:
            toks = rows[-1].split()
            v = toks[11]
            pmax = float(v)
            if not math.isfinite(pmax): pmax = float("inf")
    except Exception:
        pmax = float("inf")
    return pmax

def controller(kappa0, name):
    kappa = kappa0
    traj = []
    print(f"\n--- SOC controller, start kappa={kappa0} ({name}) ---")
    print(f"{'step':>4} {'kappa':>9} {'P_max':>10} {'state':>9} {'-> next kappa':>14}")
    for n in range(STEPS):
        pmax = run_step(kappa, f"{name}_{n}")
        collapsed = (not math.isfinite(pmax)) or (pmax > P_COLLAPSE)
        nxt = kappa * (UP if collapsed else DOWN)
        state = "COLLAPSE" if collapsed else "stable"
        pm_str = "inf/NaN" if not math.isfinite(pmax) else f"{pmax:.4f}"
        print(f"{n:>4} {kappa:>9.4f} {pm_str:>10} {state:>9} {nxt:>14.4f}")
        traj.append((kappa, pmax, collapsed))
        kappa = nxt
    # self-tuned value = geometric mean of last 6 kappas (the parked band)
    tail = [t[0] for t in traj[-6:]]
    kstar = math.exp(sum(math.log(k) for k in tail)/len(tail))
    lo, hi = min(tail), max(tail)
    print(f"  => parked band kappa* = {kstar:.3f}  (range {lo:.3f}-{hi:.3f})")
    return kstar, (lo, hi), traj

if __name__ == "__main__":
    print("="*64)
    print("v65 X1 — homeostatic self-tuning of kappa to the stability edge")
    print("="*64)
    k_hi, band_hi, _ = controller(8.0, "hi")    # high start (stable, drifts down to edge)
    k_lo, band_lo, _ = controller(0.7,  "lo")   # low start (collapse, pushed up)
    print("\n" + "="*64)
    print("RESULT")
    print("-"*6)
    print(f"  from HIGH start (kappa0=50):  kappa* = {k_hi:.3f}  band {band_hi[0]:.2f}-{band_hi[1]:.2f}")
    print(f"  from LOW  start (kappa0=0.7): kappa* = {k_lo:.3f}  band {band_lo[0]:.2f}-{band_lo[1]:.2f}")
    agree = abs(math.log(k_hi/k_lo)) < math.log(2.0)   # within a factor of 2
    if agree:
        print(f"  => CONVERGED: both starts self-tune to a COMMON edge kappa* ~ "
              f"{math.exp((math.log(k_hi)+math.log(k_lo))/2):.2f}, basin-independent.")
        print("     The SOC attractor at the stability edge is REAL (mechanism B works).")
    else:
        print("  => starts did NOT converge to a common value; SOC fixed point not robust.")
    print(f"  (E3-measured stability threshold kappa_crit ~ 1 for reference.)")
    print("="*64)
