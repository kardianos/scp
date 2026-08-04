#!/bin/bash
# Pad-campaign master wave 1. Serves pads 03,05,06,07,10,11,12,13,15,16,
# 17,18,19,20,22,23,24,25,30,31,33 (round 1 data).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
J=runs/pad/jobs1.txt
: > "$J"
# M1: instrumented glow baths (full QATOM stream + snapshots)
echo "./freecell exp=bath T=2000 noise_amp=0.5 diag_every=200 qatom_every=1 snap_every=500 snap_file=runs/pad/m1_glow_k005.fcs $K > runs/pad/m1_glow_k005.log 2>&1" >> "$J"
echo "./freecell exp=bath T=2000 noise_amp=0.5 diag_every=200 qatom_every=1 snap_every=500 snap_file=runs/pad/m1_glow_k0.fcs k_rad=0 > runs/pad/m1_glow_k0.log 2>&1" >> "$J"
# M2: instrumented cavity (QATOM every fire)
echo "./freecell exp=ring ring_n=6 bath=0 T=20000 diag_every=100 qatom_every=1 snap_every=10000 snap_file=runs/pad/m2_cav.fcs $K > runs/pad/m2_cav.log 2>&1" >> "$J"
# M3: noise_amp x k sweep (nucleation/ripening/tracker)
for a in 0.15 0.3 0.5 0.8; do
  for kk in 0 0.05; do
    kn=$(echo $kk | tr -d '.')
    if [ "$kk" = "0" ]; then KA="k_rad=0"; else KA="$K"; fi
    echo "./freecell exp=bath T=600 noise_amp=$a diag_every=200 snap_every=500 snap_file=runs/pad/m3_a${a}_k${kn}.fcs $KA > runs/pad/m3_a${a}_k${kn}.log 2>&1" >> "$J"
  done
done
# M5: nv sweep at x* (minimum-mass pad)
for n in 6 12 24 48; do
  echo "./freecell exp=ring ring_n=$n ring_x=0.62 bath=0 T=20000 diag_every=500 $K > runs/pad/m5_nv${n}.log 2>&1" >> "$J"
done
# M6: at-cap corpse N-scaling (Hawking pad)
for n in 6 12 24; do
  echo "./freecell exp=ring ring_n=$n ring_x=1.0 bath=0 T=5000 diag_every=100 $K > runs/pad/m6_corpse${n}.log 2>&1" >> "$J"
done
# M7id: identity decay — lone vs bonded vs ring (confinement pad)
DOFFP=$(awk 'BEGIN{det=1+1.2*0.65; printf "%.10f", 7.0-3.14159265358979*det/2.9}')
echo "./freecell exp=pair pair_x0=0.65 pair_doff=$DOFFP bath=1 noise_amp=0.5 convtag=1 T=600 diag_every=200 $K > runs/pad/m7_lone.log 2>&1" >> "$J"
echo "./freecell exp=pair pair_x0=0.65 pair_doff=0.05 bath=1 noise_amp=0.5 convtag=1 T=600 diag_every=200 $K > runs/pad/m7_bonded.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.65 bath=1 noise_amp=0.5 convtag=1 T=600 diag_every=200 $K > runs/pad/m7_ring.log 2>&1" >> "$J"
# M7g4: maintained pair, space profile via snapshots (far-field pad)
DOFF62=$(awk 'BEGIN{det=1+1.2*0.62; printf "%.10f", 7.0-3.14159265358979*det/2.9}')
echo "./freecell exp=pair pair_x0=0.62 pair_doff=$DOFF62 bath=1 noise_amp=0.5 convtag=1 T=1000 diag_every=200 snap_every=1000 snap_file=runs/pad/m7_g4_k005.fcs $K > runs/pad/m7_g4_k005.log 2>&1" >> "$J"
echo "./freecell exp=pair pair_x0=0.62 pair_doff=$DOFF62 bath=1 noise_amp=0.5 convtag=1 T=1000 diag_every=200 snap_every=1000 snap_file=runs/pad/m7_g4_k0.fcs k_rad=0 > runs/pad/m7_g4_k0.log 2>&1" >> "$J"
# M7and: glow point-injection into a dark bath, disorder dial (Anderson pad)
for rj in 0.01 0.06 0.15; do
  echo "./freecell exp=pair pair_x0=1.0 pair_doff=$DOFF62 bath=1 noise_amp=0 convtag=1 T=300 diag_every=100 snap_every=250 snap_file=runs/pad/m7_and_rj${rj}.fcs rjit=$rj $K > runs/pad/m7_and_rj${rj}.log 2>&1" >> "$J"
done
wc -l "$J"
xargs -P 10 -I CMD bash -c CMD < "$J"
echo WAVE1-DONE
