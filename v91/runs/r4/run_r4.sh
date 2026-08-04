#!/bin/bash
# R4 — forging closure: FORGE E3 (v_law big-void protocol, verbatim from
# v90/runs/forge/v_law_*.log headers) at k0 vs the selected radiance point.
# lump = void grad_r0=6 grad_frac=0.15; ctl = uniform bath. tag_r=6 region
# ledger. Law regime (V2g defaults), beam via slit apparatus mask=3.
set -u
cd "$(dirname "$0")/../.."
BASE="exp=slit slit_mask=3 L=64 T=32 sigma=4 slit_sy=3 kx=2.0 amp=4 slit_srcx=14 slit_screenx=40 slit_t0=0 slit_t1=30 sect_meter=1 sect_r0=7 sect_r1=11 sect_t0=0 sect_t1=30 diag_every=100 tag_r=6 convtag=1"
K005="k_rad=0.05 p_rad=4 rad_clock=0"
JOBS=runs/r4/jobs_r4.txt
: > "$JOBS"
for s in 20260802 111; do
  echo "./freecell $BASE seed=$s grad_r0=6 grad_frac=0.15 k_rad=0 > runs/r4/r4_lump_k0_s$s.log 2>&1" >> "$JOBS"
  echo "./freecell $BASE seed=$s k_rad=0 > runs/r4/r4_ctl_k0_s$s.log 2>&1" >> "$JOBS"
  echo "./freecell $BASE seed=$s grad_r0=6 grad_frac=0.15 $K005 > runs/r4/r4_lump_k005_s$s.log 2>&1" >> "$JOBS"
  echo "./freecell $BASE seed=$s $K005 > runs/r4/r4_ctl_k005_s$s.log 2>&1" >> "$JOBS"
done
wc -l "$JOBS"
xargs -P 8 -I CMD bash -c CMD < "$JOBS"
echo "R4 done"
