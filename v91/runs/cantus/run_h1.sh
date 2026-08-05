#!/bin/bash
# H1 — the lock instrument at pair scale (CANTUS.md §3.1/§3.2).
# H1a: pair OF the bath at x=0.47, seeded INSIDE the lens (doff -0.15;
#      the rung 1.694 sits outside s_pull-deflated contact 1.594 — the
#      squeezed-bond regime IS the x-ceiling physics), ambient jitter
#      real, T=1200. Control = k_cant=0, same seed.
# H1b: fifth pair (3:2, m=2) vacuum+noise, x_U=0.28 x_D=0.8367 (chart),
#      T=1200, same grid. The kr=1 negative-template face.
# All arms: k005 radiance harness. Seed 20260802.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
JOBS=runs/cantus/jobs_h1.txt
: > "$JOBS"
for kc in 0 0.3 1 3; do
  for kt in 0 0.2; do
    tagc="kc${kc/./p}_kt${kt/./p}"
    echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K k_cant=$kc k_tune=$kt > runs/cantus/h1a_${tagc}.log 2>&1" >> "$JOBS"
    echo "./freecell exp=pair bath=0 pair_pp=3 pair_qq=2 pair_m=2 pair_x0=0.28 pair_x1=0.8367 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K k_cant=$kc k_tune=$kt > runs/cantus/h1b_${tagc}.log 2>&1" >> "$JOBS"
  done
done
wc -l "$JOBS"
xargs -P 8 -I CMD bash -c CMD < "$JOBS"
echo "H1 done"
