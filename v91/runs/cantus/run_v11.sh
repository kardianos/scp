#!/bin/bash
# v1.1 campaign relaunch (CANTUS.md §3.3). k_cant=0 arms are NOT rerun:
# the v1.1 change is byte-inert at defaults, so every v1 control log
# remains valid. v1 kc>0 logs are renamed *_v1.log by the driver before
# relaunch (the v1 record is kept — it measured the bath transition).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
SEL="k_cant=1 k_tune=0.2"
LOW="k_cant=0.3 k_tune=0.2"
D=runs/cantus
# preserve the v1 record
for f in h1a_kc0_kt0p2 h1b_kc0_kt0p2 h1a_kc0p3_kt0 h1a_kc0p3_kt0p2 h1a_kc1_kt0 h1a_kc1_kt0p2 \
         h1a_kc3_kt0 h1a_kc3_kt0p2 h1b_kc0p3_kt0 h1b_kc0p3_kt0p2 \
         h1b_kc1_kt0 h1b_kc1_kt0p2 h1b_kc3_kt0 h1b_kc3_kt0p2 \
         h2a_sel h2a_low h2a_tau150 h2a_seed1 h4_sel; do
  [ -f "$D/$f.log" ] && mv "$D/$f.log" "$D/${f}_v1.log"
done
J=$D/jobs_v11.txt
: > "$J"
for kc in 0 0.3 1 3; do
  for kt in 0 0.2; do
    if [ "$kc" = "0" ] && [ "$kt" = "0" ]; then continue; fi
    tagc="kc${kc/./p}_kt${kt/./p}"
    echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K k_cant=$kc k_tune=$kt > $D/h1a_${tagc}.log 2>&1" >> "$J"
    echo "./freecell exp=pair bath=0 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.28 pair_x1=0.8367 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K k_cant=$kc k_tune=$kt > $D/h1b_${tagc}.log 2>&1" >> "$J"
  done
done
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K $SEL > $D/h2a_sel.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K $LOW > $D/h2a_low.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K $SEL cant_seed=1 > $D/h2a_seed1.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.62 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $SEL > $D/h2b_sel.log 2>&1" >> "$J"
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $SEL > $D/h2c_sel.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=2000 seed=20260802 $K $SEL > $D/h2d_sel.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $SEL > $D/h2e_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K $SEL > $D/h3_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K > $D/h3_ctl.log 2>&1" >> "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 seed=20260802 $K $SEL > $D/h4_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=0 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.28 pair_x1=0.8367 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K $SEL cant_seed=1 > $D/h1b_seed1.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K $SEL cant_seed=1 > $D/h3_seed1.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 14 -I CMD bash -c CMD < "$J"
true
echo "V11 done"
