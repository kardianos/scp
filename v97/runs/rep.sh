#!/bin/bash
# v97 PLAN §5b replication queue (user: bigger box, smaller timestep).
# L=24 dt=0.01 snap_every=400 keeps frame cadence at 4 t.u.
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=./fabrun
C24=(L=24 dt=0.01 T=240 amp_nat=2 amp_logate=1 hop_order=1
     snap_every=400 diag_every=4000)
RING=(exp=ring ring_n=6 seed_mw=2 bath=1 wf_on=1 k_rad=0.05 promo=2
      amp_door=1 noise_amp=0.4)

run() {
  local name="$1"; shift
  "$BIN" "${C24[@]}" "$@" snap_file=runs/rep/${name}.fcs \
    > "runs/rep/${name}.log" 2>&1
  echo "$name done drift=$(awk '/RESULT drift_rel/{print $2}' runs/rep/${name}.log | tail -1)"
}
mkdir -p runs/rep
case "${1:?arm}" in
  rc)  for s in 20260802 111 271828; do run rc_s$s  "${RING[@]}" kx=0   seed=$s & done; wait ;;
  rf)  for s in 20260802 111 271828; do run rf_s$s  "${RING[@]}" kx=1.1 seed=$s & done; wait ;;
  qp)  run quench "${RING[@]}" kx=0 noise_amp=0.1 T=800 seed=20260802 ;;
  p3)  run p3ring "${RING[@]}" kx=0 promo=3 seed=20260802 ;;
  *) echo "unknown arm $1"; exit 2 ;;
esac
