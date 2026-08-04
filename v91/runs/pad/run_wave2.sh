#!/bin/bash
# Pad wave 2 (bath subset): aged-beam (09), 5:3 scan (34), blob-snap (33),
# comp-dim (29), blob2 floors (30).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
BASE="exp=slit slit_mask=3 L=64 sigma=4 slit_sy=3 kx=2.0 amp=4 slit_srcx=14 slit_screenx=40 sect_meter=1 sect_r0=7 sect_r1=11 diag_every=200 convtag=1 tag_r=6"
J=runs/pad/jobs2.txt
: > "$J"
echo "./freecell $BASE T=660 slit_t0=600 slit_t1=630 sect_t0=600 sect_t1=630 seed=20260802 $K > runs/pad/w2_aged_k005.log 2>&1" >> "$J"
echo "./freecell $BASE T=60 slit_t0=0 slit_t1=30 sect_t0=0 sect_t1=30 seed=20260802 $K > runs/pad/w2_fresh_k005.log 2>&1" >> "$J"
for x in 0.53 0.55 0.555 0.56 0.565 0.57 0.59; do
  doff=$(awk -v x="$x" 'BEGIN{det=1+1.2*x; printf "%.10f", 7.0-3.14159265358979*det/2.9}')
  xt=$(awk -v x="$x" 'BEGIN{printf "%05d", x*10000+0.5}')
  echo "./freecell exp=pair pair_x0=$x pair_doff=$doff bath=1 noise_amp=0.5 convtag=1 T=480 diag_every=100 seed=20260802 $K > runs/pad/w2_53scan_x$xt.log 2>&1" >> "$J"
done
echo "./freecell exp=blob convtag=1 noise_amp=0.5 T=2000 diag_every=200 qatom_every=20 snap_every=500 snap_file=runs/pad/w2_blob.fcs $K > runs/pad/w2_blob.log 2>&1" >> "$J"
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.2 convtag=1 T=600 diag_every=200 snap_every=500 snap_file=runs/pad/w2_compdim.fcs $K > runs/pad/w2_compdim.log 2>&1" >> "$J"
echo "./freecell exp=blob2 blob2_sep=10 convtag=1 noise_amp=0.5 T=1500 diag_every=200 snap_every=1000 snap_file=runs/pad/w2_blob2.fcs $K > runs/pad/w2_blob2.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 8 -I CMD bash -c CMD < "$J"
echo WAVE2-DONE
