#!/bin/bash
# v97 promo leads — Go-only, H1 scored. Default promo=0 is the control.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p runs/leads
BIN=./fabrun
COMMON=(exp=blob bath=1 L=16 dt=0.02 T=240 amp_nat=2 amp_logate=1 hop_order=1
        seed=20260802 snap_every=200 diag_every=2000)

run() {
  local name="$1"; shift
  echo "=== $name ==="
  "$BIN" "${COMMON[@]}" "$@" snap_file=runs/leads/${name}.fcs \
    > "runs/leads/${name}.log" 2>&1
  echo "$name done drift=$(awk '/RESULT drift_rel/{print $3}' runs/leads/${name}.log | tail -1)"
}

run p0_ctl   amp=0.8 sigma=2.5 kx=1.1 wf_on=1 k_rad=0.05 promo=0 &
run p1_slot  amp=0.8 sigma=2.5 kx=1.1 wf_on=1 k_rad=0.05 promo=1 &
run p2_space amp=0.8 sigma=2.5 kx=1.1 wf_on=1 k_rad=0.05 promo=2 &
run p3_close amp=0.8 sigma=2.5 kx=1.1 wf_on=1 k_rad=0.05 promo=3 &
wait
run p2_hot   amp=2.0 sigma=2.5 kx=1.1 wf_on=0 k_rad=0.05 promo=2 &
run p2_noise amp=0.8 sigma=2.5 kx=1.1 wf_on=1 k_rad=0.05 noise_amp=0.4 promo=2 &
run p2_ring  exp=ring ring_n=6 bath=0 seed_mw=2 kx=0 wf_on=0 k_rad=0.05 promo=2 T=240 &
wait
echo ALL_LEADS_DONE
