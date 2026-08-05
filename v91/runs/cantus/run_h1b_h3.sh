#!/bin/bash
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
J=runs/cantus/jobs_h1b3.txt
: > "$J"
for kc in 0 0.3 1 3; do
  for kt in 0 0.2; do
    tagc="kc${kc/./p}_kt${kt/./p}"
    echo "./freecell exp=pair bath=0 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.28 pair_x1=0.8367 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K k_cant=$kc k_tune=$kt > runs/cantus/h1b_${tagc}.log 2>&1" >> "$J"
  done
done
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K k_cant=1 k_tune=0.2 > runs/cantus/h3_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K > runs/cantus/h3_ctl.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 5 -I CMD bash -c CMD < "$J"
echo "H1B/H3 done"
