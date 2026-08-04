#!/bin/bash
# R1 coarse sweep driver — RADIANCE.md §1. 9 arms x 17 x-points, T=120.
# doff = 7.0 - d*(x) exact => picker target constant => one probed pair
# (475/1565) across the entire grid.
set -u
cd "$(dirname "$0")/../.."
XS="0 0.10 0.21 0.28 0.35 0.4167 0.48 0.52 0.56 0.60 0.65 0.70 0.75 0.8333 0.92 1.00 1.05"
ARMS="k000:0:4:0 k001:0.01:4:0 k002:0.02:4:0 k005:0.05:4:0 k010:0.10:4:0 k030:0.30:4:0 c1k005:0.05:4:1 p2k005:0.05:2:0 p6k005:0.05:6:0"
JOBS=runs/radiance/jobs.txt
: > "$JOBS"
for arm in $ARMS; do
  IFS=: read -r name k p c <<< "$arm"
  for x in $XS; do
    doff=$(awk -v x="$x" 'BEGIN{det=1+1.2*x; printf "%.10f", 7.0-2*3.14159265358979323846/(2.9/det+2.9)}')
    xtag=$(awk -v x="$x" 'BEGIN{printf "%05d", x*10000+0.5}')
    echo "./freecell exp=pair bath=1 freeze_geo=1 convtag=1 noise_amp=0.5 T=120 diag_every=200 seed=20260802 pair_x0=$x pair_x1=0 pair_doff=$doff k_rad=$k p_rad=$p rad_clock=$c > runs/radiance/r1_${name}_x${xtag}.log 2>&1" >> "$JOBS"
  done
done
wc -l "$JOBS"
xargs -P 12 -I CMD bash -c CMD < "$JOBS"
echo "R1 sweep done: $(ls runs/radiance/r1_*.log | wc -l) logs"
