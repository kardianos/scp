#!/bin/bash
# v97 PLAN step 2 — the registered A.5 arm, unbundled (one knob per arm).
# Usage: runs/a5.sh <arm>   (arms: rg2a rg2b rg2c rg2d)
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=./fabrun
COMMON=(L=16 dt=0.02 T=240 amp_nat=2 amp_logate=1 hop_order=1
        seed=20260802 snap_every=200 diag_every=2000)
RING=(exp=ring ring_n=6 seed_mw=2 bath=1 wf_on=1 k_rad=0.05 promo=2)

run() {
  local name="$1"; shift
  echo "=== $name ==="
  "$BIN" "${COMMON[@]}" "$@" snap_file=runs/a5/${name}.fcs \
    > "runs/a5/${name}.log" 2>&1
  echo "$name done drift=$(awk '/RESULT drift_rel/{print $3}' runs/a5/${name}.log | tail -1)"
  awk '/RESULT conv/{print}' "runs/a5/${name}.log" | tail -1
}

case "${1:?arm}" in
  rg2a) run rg2a "${RING[@]}" kx=0 ;;                 # 2a: A.5 + promo=2, nothing else
  rg2b) run rg2b "${RING[@]}" kx=0 amp_door=1 ;;      # 2b: + amp_door (A.5 untested refill)
  rg2c) run rg2c "${RING[@]}" kx=0 amp_door=1 noise_amp=0.4 ;;  # 2a-lit: light the field honestly
  rg2d) run rg2d "${RING[@]}" amp_door=1 noise_amp=0.4 kx=1.1 ;; # 2c: helix tilt alone on best-of (NULL: seed_mw overrides kx)
  rg2dp) run rg2dp exp=ring ring_n=6 seed_mw=0 bath=1 wf_on=1 k_rad=0.05 promo=2 amp_door=1 noise_amp=0.4 kx=1.1 ;; # tilted unwound ring
  rg2e) run rg2e "${RING[@]}" kx=0 amp_door=1 noise_amp=0.4 T=1000 ;;   # park probe: longer T
  rg2f) run rg2f "${RING[@]}" amp_door=1 noise_amp=0.4 kx=1.1 ;;    # composed helix
  *) echo "unknown arm $1"; exit 2 ;;
esac
# post-rg2d discovery: seed_mw ASSIGNS th2 after the kx block, so
# kx+seed_mw ≡ seed_mw (rg2d was byte-identical to rg2c). The composed
# helix needs a fab/run.go seed edit (user sign-off). Config-reachable:
# the tilted UNWOUND ring — H1 does not require winding (H3 does).
# rg2f: the composed helix (seed_mw=2 + kx tilt), enabled by the
# user-authorized fab/run.go seed-order fix (battery ALL GREEN 87 after).
