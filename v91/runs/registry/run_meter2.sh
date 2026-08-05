#!/bin/bash
# Clean rerun of the 6 arms whose logs got double-written after a
# false kill diagnosis (instrument lesson in REGISTRY.md §3.4: a
# missing RESULT line means UNFINISHED — verify by PID). New binary:
# q90 column + f0 forms present but INERT here (reg_gate=0).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
D=runs/registry
J=$D/jobs_meter2.txt
: > "$J"
echo "./freecell exp=bath T=600 noise_amp=0.5 diag_every=100 seed=20260802 $K reg_tau=100 > $D/m1_bath_t100.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=100 > $D/m2_uud_t100.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=30 > $D/m3_pair_t30.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=100 > $D/m3_pair_t100.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=30 > $D/m4_i5_t30.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=100 > $D/m4_i5_t100.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 6 -I CMD bash -c CMD < "$J"
echo "METER2 done"
