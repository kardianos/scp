#!/bin/bash
# K3 registered runs (../v102/kernel/K3_PREREG.md); binary = the K2e commit's kernel
S=/tmp/claude-1000/-home-d-code-scp/f4d23169-e1d3-47f4-925b-9679d844f070/scratchpad/k3
B=$S/freecell_k3
LAW="k_rad=0.05 p_rad=4 rad_clock=0 wf_on=1 wf_floor=0.01 wf_far=99"
COM="T=2000 diag_every=5000 lockq=1 lockq_p11=2 lockq_proj=1 lockq_every=10"
F32="lockq_M=2 lockq_n=10 fmodes=2"; FW="lockq_M=4 lockq_n=7 fmodes=4"
cd /home/d/code/scp/v97
MAXJ=${MAXJ:-6}
job() { while [ $(jobs -rp | wc -l) -ge $MAXJ ]; do sleep 5; done; "$@" & }
out=$S/runs
# part B first (seeded, no light)
job sh -c "$B $LAW exp=bath $COM $F32 lockq_rate=0 lockq_force=4 lockq_fmode=0 > $out/B1mid_M2.out 2>&1"
job sh -c "$B $LAW exp=bath $COM $FW  lockq_rate=0 lockq_force=3 lockq_fmode=0 > $out/B1mid_M4.out 2>&1"
job sh -c "$B $LAW exp=bath $COM $F32 lockq_rate=0 lockq_force=1 lockq_fmode=0 > $out/B1end_M2.out 2>&1"
job sh -c "$B $LAW exp=bath $COM $FW  lockq_rate=0 lockq_force=1 lockq_fmode=0 > $out/B1end_M4.out 2>&1"
job sh -c "$B $LAW exp=bath $COM $F32 lockq_rate=0 lockq_force=3 lockq_fmode=0 lockq_force2=6 lockq_fmode2=1 > $out/B2_M2.out 2>&1"
job sh -c "$B $LAW exp=bath $COM $FW  lockq_rate=0 lockq_force=2 lockq_fmode=0 lockq_force2=4 lockq_fmode2=3 > $out/B2_M4.out 2>&1"
# part A (light, mid placement)
for sd in 1 2 3 4 5; do
  job sh -c "$B $LAW exp=pulse $COM $FW  fseed_mode=-1 lockq_dx=-5.0 amp=4.0 sigma=3.0 lockq_seed=$sd > $out/A_M4_b_s$sd.out 2>&1"
  job sh -c "$B $LAW exp=pulse $COM $F32 fseed_mode=-1 lockq_dx=-7.5 amp=4.0 sigma=3.0 lockq_seed=$sd > $out/A_M2_b_s$sd.out 2>&1"
done
for sd in 1 2 3 4 5; do
  job sh -c "$B $LAW exp=pulse $COM $F32 fseed_mode=-1 lockq_dx=-7.5 amp=1.0 sigma=2.2 lockq_seed=$sd > $out/A_M2_d_s$sd.out 2>&1"
  job sh -c "$B $LAW exp=pulse $COM $FW  fseed_mode=-1 lockq_dx=-5.0 amp=1.0 sigma=2.2 lockq_seed=$sd > $out/A_M4_d_s$sd.out 2>&1"
done
wait
touch $S/done
