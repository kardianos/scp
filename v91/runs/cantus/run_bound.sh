#!/bin/bash
# The cantus UPPER BOUND (registered in CANTUS.md §4.3 before running):
# cant_tau=1e18 freezes sgg at seed (object bonds amp=1 forever, bath 0)
# and cxl at the birth chord — the strongest expressible lock+memory.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
B="k_cant=1 k_tune=0.2 cant_seed=1 cant_grow=0 cant_tau=1e18"
D=runs/cantus
J=$D/jobs_bound.txt
: > "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $B > $D/hb_e.log 2>&1" >> "$J"
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $B > $D/hb_c.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=2000 seed=20260802 $K $B > $D/hb_d.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K $B > $D/hb_a.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K $B > $D/hb_3.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 5 -I CMD bash -c CMD < "$J"
echo "BOUND done"
