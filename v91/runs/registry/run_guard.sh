#!/bin/bash
# P-REG3 vacuum guard as FORM/f0 SELECTOR (REGISTRY.md §3.4): bare
# bath, FULL stack live, self-arming, registry-gated. Bars: cond/rad
# within 10% of the k005 control (2203.972100 / 2016.222021),
# tune_total <= 1000, nlock <= 2x null (3546), bath a_max < 0.5
# throughout. Control reference: runs/cantus/h4_ctl.log (k005, no
# cantus) and runs/cantus/h4_null.log (zero-force gauges).
set -u
cd "$(dirname "$0")/../.."
K="k_rad=0.05 p_rad=4 rad_clock=0"
S="k_cant=1 k_tune=0.2 cant_tau=50 cant_grow=1 reg_tau=10"
D=runs/registry
J=$D/jobs_guard.txt
: > "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 diag_every=40 seed=20260802 $K $S reg_gate=1 reg_f0=0.0002 > $D/g_B_f2.log 2>&1" >> "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 diag_every=40 seed=20260802 $K $S reg_gate=1 reg_f0=0.0005 > $D/g_B_f5.log 2>&1" >> "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 diag_every=40 seed=20260802 $K $S reg_gate=2 reg_f0=0.0002 > $D/g_D_f2.log 2>&1" >> "$J"
echo "./freecell exp=bath T=480 noise_amp=0.5 diag_every=40 seed=20260802 $K $S reg_gate=2 reg_f0=0.0005 > $D/g_D_f5.log 2>&1" >> "$J"
wc -l "$J"
xargs -P 4 -I CMD bash -c CMD < "$J"
echo "GUARD done"
