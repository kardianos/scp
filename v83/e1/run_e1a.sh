#!/usr/bin/env bash
# E1a: hybrid force scan — gauged Q-ball + free heavy opposite lock
set -euo pipefail
ROOT=/home/d/code/scp
BIN=$ROOT/bin/scp_sim
SEED=$ROOT/v83/e1/seeds/ball_N48.sfa
OUT=$ROOT/v83/e1/results
mkdir -p "$OUT"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

run_one() {
  local tag=$1 D=$2 medium_only=$3
  local cfg=$OUT/e1a_${tag}.cfg
  local track=$OUT/track_${tag}.tsv
  local diag=$OUT/diag_${tag}.tsv
  local log=$OUT/run_${tag}.log
  # lock at +D on x; q=-1 opposite to positive Q-ball; heavy free (pinned skips EM f)
  cat >"$cfg" <<EOF
N=48
L=12
T=8
dt_factor=0.05
complex_phi=1
complex_gauge=1
g_gauge=0.05
eta=0
m=1.5
m_theta=1.6
mu=-41.345
kappa=50
A=0
A_bg=0
init=sfa
init_sfa=$SEED
bc_type=0
damp_width=2.0
damp_rate=0.08
n_locks=1
lock0=-1,500,$D,0,0,0,0,0,0
locks_track=$track
locks_medium_only=$medium_only
lock_soft_r=0
lock_soft_k=0
lock_bag_mode=0
output=$OUT/out_${tag}.sfa
diag_file=$diag
snap_dt=1000
diag_dt=0.5
precision=f32
qdiag_radius=8
EOF
  echo "=== E1a $tag D=$D medium_only=$medium_only ==="
  "$BIN" "$cfg" >"$log" 2>&1 || { echo "RUN FAIL $tag"; tail -30 "$log"; return 1; }
  # parse mean Fx over mid window (skip t=0)
  python3 - <<PY
import numpy as np
t = np.loadtxt("$track", skiprows=1)
# cols: t id type x y z ux uy uz fx fy fz pinned alive
if t.ndim==1: t=t.reshape(1,-1)
tt, x, fx = t[:,0], t[:,3], t[:,9]
mask = (tt>1.0) & (tt<7.0)
if mask.sum()<2:
    mask = tt>0
fxm = fx[mask].mean()
uxm = t[mask,6].mean()
xm = x[mask].mean()
# attract toward origin: for lock at +D, F_x should be negative if attract
print(f"tag=$tag D=$D  <x>={xm:.4f}  <Fx>={fxm:.6e}  <ux>={uxm:.6e}  n={mask.sum()}")
open("$OUT/summary_${tag}.txt","w").write(f"$tag\t$D\t{fxm:.8e}\t{xm:.6f}\t{uxm:.8e}\n")
PY
  # gauss from log/diag
  rg -n "gauss|Gauss|Q_phi|Q_locks|E_total" "$log" | head -20 || true
  if [[ -f "$diag" ]]; then
    head -1 "$diag"
    tail -3 "$diag"
  fi
}

# Vacuum control: no multiplet force, empty phi still from seed... use medium_only
# Hybrid scan
for D in 5 6 8 10 12; do
  run_one "hyb_D${D}" "$D" 0
done
# medium_only control at D=8 (Φ frozen off — vacuum multiplet path + seed still loaded?)
# locks_medium_only skips multiplet force; seed still deposits matter ρ if not zeroed
# For pure lock self-force: use init=oscillon? Better: hyb vs note

# Self-force baseline: same seed but lock q=-1 far — actually same D with locks_medium_only
# doesn't remove matter ρ from Gauss. Need n_locks with vacuum init.
run_vac() {
  local D=$1
  local tag=vac_D${D}
  local cfg=$OUT/e1a_${tag}.cfg
  local track=$OUT/track_${tag}.tsv
  cat >"$cfg" <<EOF
N=48
L=12
T=8
dt_factor=0.05
complex_phi=1
complex_gauge=1
g_gauge=0.05
eta=0
m=1.5
m_theta=1.6
mu=-41.345
kappa=50
init=oscillon
A=0
A_bg=0
bc_type=0
damp_width=2.0
damp_rate=0.08
n_locks=1
lock0=-1,500,$D,0,0,0,0,0,0
locks_track=$track
locks_medium_only=1
lock_soft_r=0
lock_soft_k=0
lock_bag_mode=0
output=$OUT/out_${tag}.sfa
diag_file=$OUT/diag_${tag}.tsv
snap_dt=1000
diag_dt=0.5
precision=f32
EOF
  echo "=== E1a $tag vacuum self-force ==="
  "$BIN" "$cfg" >"$OUT/run_${tag}.log" 2>&1 || { tail -20 "$OUT/run_${tag}.log"; return 1; }
  python3 - <<PY
import numpy as np
t = np.loadtxt("$track", skiprows=1)
if t.ndim==1: t=t.reshape(1,-1)
tt, fx = t[:,0], t[:,9]
mask = (tt>1.0)&(tt<7.0)
fxm = fx[mask].mean() if mask.sum() else fx.mean()
print(f"tag=$tag D=$D  <Fx_self>={fxm:.6e}")
open("$OUT/summary_${tag}.txt","w").write(f"$tag\t$D\t{fxm:.8e}\n")
PY
}
run_vac 8

echo "==== E1a summaries ===="
cat "$OUT"/summary_*.txt 2>/dev/null | sort
