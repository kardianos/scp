#!/bin/bash
# v93 FACE B: lit-bath refill — does the arg(psi) door (amp_door=1) retain the
# ring winding when the door actually fires? (RING_DNLS §A.5: under V3a the
# dark bath gave cond=0, so amp_door was untested. noise_amp lights the medium.)
# Prediction (RESUME §7-B): amp_door=0 scrambles on refill; amp_door=1 is the
# only candidate hold.
cd "$(dirname "$0")/.."
mkdir -p runs/faceB
# Lit medium: noise_amp seeds field everywhere -> door fires (cond>0). V3a law.
BASE="exp=ring ring_n=6 bath=1 L=16 dt=0.02 T=80 seed=20260802
  k_rad=0.05 p_rad=4 wf_on=1 noise_amp=2.0 seed_mw=2
  diag_every=1000 snap_every=8"
echo "### FACE B: lit-bath refill (noise_amp=2 -> cond>0; V3a) ###"
run () {
  local name=$1; shift
  ./freecell $BASE "$@" snap_file=runs/faceB/$name.fcs > runs/faceB/$name.log 2>&1
  echo "--- $name ($*) ---"
  awk '/RESULT conv/{printf "  cond=%s evap=%s rad=%s\n",$4,$6,$10}' runs/faceB/$name.log
  ./fcsdump -mode cells runs/faceB/$name.fcs 2>/dev/null \
    | python3 report/analyze_ring.py /dev/stdin 8 8 2 2>/dev/null \
    | grep -E "t n_ring|0.0 |8.0 |16.0 |24.0 |32.0 |40.0 |56.0 |80.0 " | sed 's/^/    /'
}
# 1. additive baseline (no unitary channel; door default, th2 not written by door)
run lit_add    amp_nat=0
# 2. unitary (Strang) + amp_door=0 (magnitude-only refill; field phase discarded)
run lit_uni_d0 amp_nat=2 amp_logate=1 hop_order=1 amp_door=0
# 3. unitary (Strang) + amp_door=1 (coherent arg(psi) merge — the candidate hold)
run lit_uni_d1 amp_nat=2 amp_logate=1 hop_order=1 amp_door=1
# 4. dark control (noise_amp=0) — confirms A.5: cond=0, amp_door inert
run dark_uni_d1 amp_nat=2 amp_logate=1 hop_order=1 amp_door=0 noise_amp=0

echo "=== faceB done ==="
