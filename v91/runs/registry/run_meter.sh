#!/bin/bash
# REGISTRY meter arms M1-M4 (REGISTRY.md §3.1): reg_gate=0, no cantus
# force — pure ledger measurement of the bond-vs-churn identity gap.
# 3 taus x 4 arms + tau-free warmup sanity. k005 radiance ambient,
# seed carried from the CANTUS harness.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
D=runs/registry
J=$D/jobs_meter.txt
: > "$J"
for TAU in 10 30 100; do
  echo "./freecell exp=bath T=600 noise_amp=0.5 diag_every=100 seed=20260802 $K reg_tau=$TAU > $D/m1_bath_t$TAU.log 2>&1" >> "$J"
  echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=$TAU > $D/m2_uud_t$TAU.log 2>&1" >> "$J"
  echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=$TAU > $D/m3_pair_t$TAU.log 2>&1" >> "$J"
  echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=600 diag_every=100 seed=20260802 $K reg_tau=$TAU > $D/m4_i5_t$TAU.log 2>&1" >> "$J"
done
wc -l "$J"
xargs -P 6 -I CMD bash -c CMD < "$J"
echo "METER done"
