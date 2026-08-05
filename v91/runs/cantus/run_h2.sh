#!/bin/bash
# H2/H3 — the R-campaign killers + the chord at the balance, at the
# EXPECTED selection point (k_cant=1, k_tune=0.2, cant_tau=50; CANTUS.md
# §2). Launched in parallel with H1 for wall-time; the selection is
# confirmed (or re-run) from H1's measured curves BEFORE any claim.
# Controls: committed logs (comp6/uud/i5/x62-bath) + in-arm k_cant=0
# where no committed control exists (h2a vacuum 0.47, h3 fifth pair).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
SEL="k_cant=1 k_tune=0.2"
LOW="k_cant=0.3 k_tune=0.2"
JOBS=runs/cantus/jobs_h2.txt
: > "$JOBS"
# H2a ring6@0.47 VACUUM (squeezed-bond regime; quantitative gg
# prediction: locked ring parks at lens edge with gg ~ 0.7-0.8;
# unlocked control's differential error drifts free)
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K $SEL > runs/cantus/h2a_sel.log 2>&1" >> "$JOBS"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K $LOW > runs/cantus/h2a_low.log 2>&1" >> "$JOBS"
echo "./freecell exp=ring ring_n=6 ring_x=0.47 ring_doff=-0.15 bath=0 noise_amp=0.5 convtag=1 T=20000 diag_every=2000 seed=20260802 $K > runs/cantus/h2a_ctl.log 2>&1" >> "$JOBS"
# H2b ring6@0.62 bath (SECONDARY: skin-lock observation)
echo "./freecell exp=ring ring_n=6 ring_x=0.62 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $SEL > runs/cantus/h2b_sel.log 2>&1" >> "$JOBS"
# H2c composite nv=6 (UUD rings; boundary-fifth face)
echo "./freecell exp=rings rings_kind=1 rings_nv=6 rings_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=2000 seed=20260802 $K $SEL > runs/cantus/h2c_sel.log 2>&1" >> "$JOBS"
# H2d UUD triad (FCQ face: ggm, D bloat)
echo "./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=2000 diag_every=2000 seed=20260802 $K $SEL > runs/cantus/h2d_sel.log 2>&1" >> "$JOBS"
# H2e i5 rerun (small matter at the sweet spot — PRIMARY 10x bar)
echo "./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 $K $SEL > runs/cantus/h2e_sel.log 2>&1" >> "$JOBS"
# H3 the chord at the balance: fifth pair (0.35, 0.94), mean 0.645 ~ x*
echo "./freecell exp=pair bath=1 pair_pp=3 pair_qq=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K $SEL > runs/cantus/h3_sel.log 2>&1" >> "$JOBS"
echo "./freecell exp=pair bath=1 pair_pp=3 pair_qq=2 pair_m=2 pair_x0=0.35 pair_x1=0.94 pair_doff=-0.1 noise_amp=0.5 convtag=1 T=5000 diag_every=1000 seed=20260802 $K > runs/cantus/h3_ctl.log 2>&1" >> "$JOBS"
wc -l "$JOBS"
xargs -P 9 -I CMD bash -c CMD < "$JOBS"
echo "H2/H3 done"
