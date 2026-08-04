#!/bin/bash
# R1b refinement — RADIANCE.md §1b. Arms k000/k005/k010, T=480,
# fine x bracket, 3 seeds. Windowed CONVTAG every 2 t.u.
set -u
cd "$(dirname "$0")/../.."
XS="0.52 0.56 0.60 0.64 0.68 0.72 0.76 0.80"
SEEDS="20260802 111 314159"
ARMS="k000:0 k005:0.05 k010:0.10"
JOBS=runs/radiance/jobs_r1b.txt
: > "$JOBS"
for arm in $ARMS; do
  IFS=: read -r name k <<< "$arm"
  for s in $SEEDS; do
    for x in $XS; do
      doff=$(awk -v x="$x" 'BEGIN{det=1+1.2*x; printf "%.10f", 7.0-3.14159265358979323846*det/2.9}')
      xtag=$(awk -v x="$x" 'BEGIN{printf "%05d", x*10000+0.5}')
      echo "./freecell exp=pair bath=1 freeze_geo=1 convtag=1 noise_amp=0.5 T=480 diag_every=100 seed=$s pair_x0=$x pair_doff=$doff k_rad=$k p_rad=4 rad_clock=0 > runs/radiance/r1b_${name}_s${s}_x${xtag}.log 2>&1" >> "$JOBS"
    done
  done
done
wc -l "$JOBS"
xargs -P 12 -I CMD bash -c CMD < "$JOBS"
echo "R1b done: $(ls runs/radiance/r1b_*.log | wc -l) logs"
