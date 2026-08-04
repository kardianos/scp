#!/bin/bash
# R2 — stability of real objects at the selected point (k_rad=0.05 p_rad=4
# rad_clock=0) vs k000 controls, in bath. Vacuum cavity arms already run
# (ring6_vac_k005, uud_vac_k005, + k0 controls). comp6_bath_k0 already
# running separately. Ledgers on (convtag=1). Seed 20260802.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
JOBS=runs/r2/jobs_r2.txt
: > "$JOBS"
echo "./freecell exp=ring ring_n=6 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=500 $K > runs/r2/ring6_bath_k005.log 2>&1" >> "$JOBS"
echo "./freecell exp=ring ring_n=6 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=500 k_rad=0 > runs/r2/ring6_bath_k0.log 2>&1" >> "$JOBS"
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=500 $K > runs/r2/comp6_bath_k005.log 2>&1" >> "$JOBS"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=500 $K > runs/r2/uud_bath_k005.log 2>&1" >> "$JOBS"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=500 k_rad=0 > runs/r2/uud_bath_k0.log 2>&1" >> "$JOBS"
echo "./freecell exp=blob convtag=1 noise_amp=0.5 T=2000 diag_every=500 $K > runs/r2/blob_bath_k005.log 2>&1" >> "$JOBS"
echo "./freecell exp=blob convtag=1 noise_amp=0.5 T=2000 diag_every=500 k_rad=0 > runs/r2/blob_bath_k0.log 2>&1" >> "$JOBS"
wc -l "$JOBS"
xargs -P 7 -I CMD bash -c CMD < "$JOBS"
echo "R2 done"
