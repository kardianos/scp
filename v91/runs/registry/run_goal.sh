#!/bin/bash
# P-REG5'/P-REG6' — THE GOAL ARMS in the honest-instrument frame
# (REGISTRY.md §3.5, registered after the §4.2 guard rejected all
# self-growing combos). cant_grow=0 => the bath can never arm
# (medium honest by construction); seeded gauges on the object's own
# bonds are registry-MAINTAINED (F-D, f0=2e-4, tau=10). Controls at
# the same T with the same binary.
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
R="reg_tau=10 reg_gate=2 reg_f0=0.0002"
S="k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1"
D=runs/registry
J=$D/jobs_goal.txt
: > "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K > $D/uud_ctl.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $S cant_tau=50 $R > $D/uud_seed_reg.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $S cant_tau=50 reg_tau=10 > $D/uud_seed_plain.log 2>&1" >> "$J"
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 $K $S cant_tau=1e18 reg_tau=10 > $D/uud_frozen.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=100 seed=20260802 $K $S cant_tau=50 $R > $D/pair_seed_reg.log 2>&1" >> "$J"
echo "./freecell exp=pair bath=1 pair_x0=0.47 pair_x1=0.47 pair_doff=-0.15 noise_amp=0.5 convtag=1 T=1200 diag_every=100 seed=20260802 $K > $D/pair_ctl.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $S cant_tau=50 $R > $D/i5_seed_reg.log 2>&1" >> "$J"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K > $D/i5_ctl.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 8 -I CMD bash -c CMD < "$J"
echo "GOAL done"
