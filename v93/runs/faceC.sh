#!/bin/bash
# v93 FACE C: condensation-lane characterization (RESUME §7-C / RING_DNLS B).
# The unitary channel + the law's q_detune detuning spontaneously condense a
# spread packet into long-lived dense hoards (B.1/B.2, T=80). Characterize as
# MASS FORMATION: (1) long-T stability to T=300; (2) radiance interaction
# (does the V3a tax select/equilibrate hoard sizes?); (3) the hoard spectrum.
# Strang (hop_order=1, the Face A winner); RING_DNLS route B used sequential.
cd "$(dirname "$0")/.."
mkdir -p runs/faceC
echo "### FACE C: condensation lane (Strang, T=300) ###"
BASE="exp=blob bath=1 L=16 dt=0.02 T=300 sigma=2.5 kx=0 seed_mw=0
  amp_logate=1 hop_order=1 diag_every=2000 snap_every=500 seed=20260802"
run () {
  local name=$1; shift
  ./freecell $BASE "$@" snap_file=runs/faceC/$name.fcs > runs/faceC/$name.log 2>&1
  ./fcsdump -mode cells runs/faceC/$name.fcs 2>/dev/null \
    | python3 report/analyze_melt.py /dev/stdin 16 > runs/faceC/$name.melt 2>/dev/null
  echo "--- $name ($*) ---"
  awk '/RESULT conv/{printf "  cond=%s evap=%s rad=%s\n",$5,$6,$8}' runs/faceC/$name.log
  awk 'NR==1{print "    "$0} NR>1{t=$1; if(t==0||t==8||t==40||t==100||t==200||t==300) print "    "$0}' runs/faceC/$name.melt
}
# --- C.1 long-T stability (no radiance, no door; pure conservative condensation) ---
run c_law05   amp=0.5 amp_nat=2 q_detune=1.2 k_rad=0 wf_on=0
run c_deep05  amp=0.5 amp_nat=2 q_detune=3.6 k_rad=0 wf_on=0
run c_deep2   amp=2   amp_nat=2 q_detune=3.6 k_rad=0 wf_on=0
# --- C.2 radiance interaction (V3a tax on the hoards) ---
run c_law05_V3a  amp=0.5 amp_nat=2 q_detune=1.2 k_rad=0.05 p_rad=4 wf_on=1
run c_deep05_V3a amp=0.5 amp_nat=2 q_detune=3.6 k_rad=0.05 p_rad=4 wf_on=1
run c_deep2_V3a  amp=2   amp_nat=2 q_detune=3.6 k_rad=0.05 p_rad=4 wf_on=1
# --- C.3 melt control (qd=0: pure linear, no nonlinearity -> should disperse) ---
run c_melt05  amp=0.5 amp_nat=2 q_detune=0 k_rad=0 wf_on=0

echo "=== hoard spectra (final frame) ==="
for name in c_law05 c_deep05 c_deep2 c_law05_V3a c_deep2_V3a c_melt05; do
  echo "--- $name @t=300 ---"
  ./fcsdump -mode cells runs/faceC/$name.fcs 2>/dev/null \
    | python3 report/analyze_hoard.py /dev/stdin 2>/dev/null | sed 's/^/    /'
done
echo "=== faceC done ==="
