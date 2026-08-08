#!/bin/bash
# v93 DNLS SELF-TRAP MAP (RESUME §7 route B): the load-detuning
# w2e = w2/(1+q_detune*x) already in the kernel is a discrete-NLS
# nonlinearity (Em-dependent phase precession + Em-dependent res
# detuning of the hop). Map the existence region of self-trapped
# (non-melting) packets under the LINEAR unitary channel: does the
# packet's own load bind it against the diffraction that killed item-4?
# Axes: amp_nat (hop), amp (packet depth = nonlinearity strength),
# q_detune (off-law diagnostic: 0 = pure linear baseline; battery
# precedent — the xs bars themselves run q_detune=0), kx (band edge).
# Meters: Em_max / PR / rms via report/analyze_melt.py (all-cell, not
# tagged). NOTE dt fixed at 0.02 — absolute rates carry the known
# dt-Trotter artifact; comparisons across arms at fixed dt are clean.
cd "$(dirname "$0")/.."
mkdir -p runs/dnls
S="exp=blob bath=1 L=16 dt=0.02 T=80 sigma=2.5 kx=0 seed_mw=0 k_rad=0 wf_on=0 diag_every=500 snap_every=200 seed=20260802 amp_logate=1"
run () {
  local name=$1; shift
  ./freecell $S "$@" snap_file=runs/dnls/$name.fcs > runs/dnls/$name.log 2>&1
  ./fcsdump -mode cells runs/dnls/$name.fcs 2>/dev/null \
    | python3 report/analyze_melt.py /dev/stdin 16 > runs/dnls/$name.melt 2>/dev/null
  echo "--- $name ($*) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s\n",$4}' runs/dnls/$name.log
  # compact: t=0/4/16/40/80 rows + t_half of Em_max
  awk 'NR==1{print "    "$0} NR>1{t=$1; if(t==0||t==4||t==16||t==40||t==80) print "    "$0}' runs/dnls/$name.melt
  awk 'NR==2{m0=$4} NR>2 && $4<0.5*m0 && !h{printf "  t_half(Em_max)=%s\n",$1; h=1} END{if(!h)printf "  t_half(Em_max)=>80\n"}' runs/dnls/$name.melt
}
# --- axis 1: hop strength (law q_detune=1.2, amp=0.5) ---
run an05 amp=0.5 amp_nat=0.5
run an1  amp=0.5 amp_nat=1
run an2  amp=0.5 amp_nat=2
run an4  amp=0.5 amp_nat=4
# --- axis 2: packet depth (nonlinearity; amp_nat=2) ---
run am025 amp=0.25 amp_nat=2
run am1   amp=1    amp_nat=2
run am2   amp=2    amp_nat=2
# --- axis 3: q_detune (off-law diagnostic; amp_nat=2 amp=0.5) ---
run qd0  amp=0.5 amp_nat=2 q_detune=0
run qd06 amp=0.5 amp_nat=2 q_detune=0.6
run qd36 amp=0.5 amp_nat=2 q_detune=3.6
run qd12 amp=0.5 amp_nat=2 q_detune=12
# --- axis 3b: depth x strong detune (the self-trap corner probe) ---
run qd36am2 amp=2 amp_nat=2 q_detune=3.6
run qd12am2 amp=2 amp_nat=2 q_detune=12
# --- axis 4: carrier toward the band edge (amp_nat=2 amp=0.5 law qd) ---
run kx11 amp=0.5 amp_nat=2 kx=1.1
run kx20 amp=0.5 amp_nat=2 kx=2.0
run kx26 amp=0.5 amp_nat=2 kx=2.6
# --- additive controls (melt under the old law, same meters) ---
run add_am05 amp=0.5 amp_nat=0
run add_am2  amp=2   amp_nat=0
echo "=== dnls_map done ==="
