#!/bin/bash
cd "$(dirname "$0")/.."
mkdir -p runs/quench
# QUENCH-3 arg(psi) door suite: does the coherent-amplitude door retain winding?
# (v91 G1 baseline = qp_phase=1 amp_nat=0 gave R2d~0.029.)
S="exp=slit slit_mask=3 L=64 T=300 dt=0.02 sigma=8 slit_sy=8 amp=2 slit_srcx=32 kx=0 spin_m=2 seed=20260802 diag_every=2000 snap_every=1000 bath=1"
run () {
  local name=$1; shift
  echo "=== $name : $* ==="
  ./freecell $S "$@" snap_file=runs/quench/$name.fcs > runs/quench/$name.log 2>&1
  grep -E "^# RESULT (drift_rel|births)" runs/quench/$name.log | sed 's/^# /  /'
  echo "  R2d retention curve:"
  ./fcsdump -mode cells runs/quench/$name.fcs 2>/dev/null | python3 report/analyze_winding.py /dev/stdin 32 32 0 2>/dev/null | awk '/^ *[0-9]/{printf "    t=%6.1f R2d=%s\n",$1,$(NF-3)}'
  echo
}
run h0_nodoor  amp_door=0 amp_nat=0                  # no door (control)
run h1_door    amp_door=1 amp_nat=0                  # arg(psi) door, additive transport
run h2_door_uni amp_door=1 amp_nat=2                 # arg(psi) door + unitary transport
run h3_uni_clk  amp_door=1 amp_nat=2 q_detune=0      # + uniform clock (shear hypothesis)
echo "=== arg(psi) door suite done ==="
