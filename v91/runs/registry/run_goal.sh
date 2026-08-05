#!/bin/bash
# P-REG5/P-REG6 — THE GOAL ARMS (REGISTRY.md §3.2, form/f0 selected
# by the §3.4 guard). Primary: self-arming UUD chord in live bath at
# the CANTUS coupled point. Secondaries: fast-tau (arming race),
# seeded-instrument ceiling, pair bar, i5 measurement. Controls
# verifiably run at the same T with the same binary.
set -u
cd "$(dirname "$0")/../.."
GATE="${GATE:?set GATE=1|2}"
F0="${F0:?set F0=0.0002|0.0005}"
K="k_rad=0.05 p_rad=4 rad_clock=0"
R="reg_tau=10 reg_gate=$GATE reg_f0=$F0"
S="k_cant=1 k_tune=0.2 cant_grow=1"
D=runs/registry
J=$D/jobs_goal.txt
: > "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $R $S cant_tau=50 > $D/goal_uud.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $R $S cant_tau=25 > $D/goal_uud_fast.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $R $S cant_tau=50 cant_seed=1 > $D/goal_uud_seed.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K > $D/goal_uud_ctl.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=100 seed=20260802 $K $R $S cant_tau=50 > $D/sec_pair.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=100 seed=20260802 $K > $D/sec_pair_ctl.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $R $S cant_tau=50 > $D/sec_i5.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K > $D/sec_i5_ctl.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 8 -I CMD bash -c CMD < "$J"
echo "GOAL done"
