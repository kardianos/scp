#!/bin/bash
# v93 R1 ENDURANCE (RECONCILIATION R1): is t=300 a plateau or a slope? Deep +
# law corner, Strang, T=3000. (R2 isolation redefined: bath=0 gives NC=0 in
# this substrate, so the leak/keep question is folded in here via melt PR/rms
# + a final keep census.) snap_every=2500 -> ~13 frames for the time series.
cd "$(dirname "$0")/.."
mkdir -p runs/reckon
BASE="exp=blob bath=1 L=16 dt=0.02 sigma=2.5 kx=0 seed_mw=0 amp_logate=1
  hop_order=1 k_rad=0 wf_on=0 diag_every=10000 snap_every=2500 seed=20260802"
run () {
  local name=$1; shift
  ./freecell $BASE T=3000 "$@" snap_file=runs/reckon/$name.fcs > runs/reckon/$name.log 2>&1
  ./fcsdump -mode cells runs/reckon/$name.fcs 2>/dev/null \
    | python3 report/analyze_melt.py /dev/stdin 16 > runs/reckon/$name.melt 2>/dev/null
  echo "--- $name ($*) ---"
  awk '/RESULT conv/{printf "  cond=%s evap=%s rad=%s\n",$5,$6,$8}' runs/reckon/$name.log
  awk 'NR>1{t=$1; if(t==0||t==300||t==1000||t==2000||t==3000) print "    t="t" Emax="$4" PR="$5" rms="$6}' runs/reckon/$name.melt
}
run r1_deep amp=2   amp_nat=2 q_detune=3.6
run r1_law  amp=0.5 amp_nat=2 q_detune=1.2
echo "=== r1 done ==="
