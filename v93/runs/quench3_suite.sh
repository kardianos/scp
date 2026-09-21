#!/bin/bash
cd "$(dirname "$0")/.."
mkdir -p runs/quench
# QUENCH-3 standing-vortex retention suite (v91 §7.5): spin_m=2 field winding,
# door carries phase (qp_phase=1), measure matter winding R2d retention.
# v93 adds amp_nat>0 (unitary coherent dense transport) — does it retain?
S="exp=slit slit_mask=3 L=64 T=300 dt=0.02 sigma=8 slit_sy=8 amp=2 slit_srcx=32 kx=0 spin_m=2 seed=20260802 diag_every=2000 snap_every=1000 bath=1"
run () {
  local name=$1; shift
  echo "=== $name : $* ==="
  ./freecell $S "$@" snap_file=runs/quench/$name.fcs > runs/quench/$name.log 2>&1
  echo "--- drift + census (final) ---"
  grep -E "^# RESULT (drift_rel|births|geometry)" runs/quench/$name.log
  echo "--- winding R2d retention curve ---"
  ./fcsdump -mode cells runs/quench/$name.fcs 2>/dev/null | python3 report/analyze_winding.py /dev/stdin 32 32 0 2>/dev/null | awk 'NR==1||/^ *[0-9]/'
  echo
}
run q3_g0_ctl  qp_phase=0 amp_nat=0     # control: door does NOT carry phase
run q3_g1_v91  qp_phase=1 amp_nat=0     # v91: door carries, additive transport
run q3_g2_a1   qp_phase=1 amp_nat=1     # v93: + unitary transport (amp_nat=1)
run q3_g3_a2   qp_phase=1 amp_nat=2     # v93: + unitary transport (amp_nat=2)
echo "=== QUENCH-3 suite done ==="
