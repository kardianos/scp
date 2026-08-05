#!/bin/bash
# v1.1i — the INSTRUMENT arms (CANTUS.md §3.4): seeded object coherence,
# cant_grow=0 => the bath stays exactly V2g+radiance. Answers P-Hi1..5.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
I="cant_seed=1 cant_grow=0"
SEL="k_cant=1 k_tune=0.2"
LOWT="k_cant=0.3 k_tune=0.2"
D=runs/cantus
J=$D/jobs_v11i.txt
: > "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $SEL $I > $D/hi_e_sel.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $LOWT $I > $D/hi_e_low.log 2>&1" >> "$J"
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $SEL $I > $D/hi_c_sel.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=2000 seed=20260802 $K $SEL $I > $D/hi_d_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K $LOWT $I > $D/hi_a_low.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=1000 seed=20260802 $K $SEL $I > $D/hi_a_sel.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_p=3 pair_q=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K $SEL $I > $D/hi_3_sel.log 2>&1" >> "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 seed=20260802 $K $SEL $I > $D/hi_4_sel.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 8 -I CMD bash -c CMD < "$J"
echo "V11I done"
