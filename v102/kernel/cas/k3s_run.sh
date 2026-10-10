#!/bin/bash
# K3 supplementary (labelled, not registered): part A again with high-entropy seeds splitmix64(1..5).
# Same commands as run.sh otherwise. Reason: xorshift64 seeded with 1..5 gives a first draw ~6e-11*seed
# and a second draw seed/16, so the registered runs' first event was not random.
S=/tmp/claude-1000/-home-d-code-scp/f4d23169-e1d3-47f4-925b-9679d844f070/scratchpad/k3
B=$S/freecell_k3
LAW="k_rad=0.05 p_rad=4 rad_clock=0 wf_on=1 wf_floor=0.01 wf_far=99"
COM="T=2000 diag_every=5000 lockq=1 lockq_p11=2 lockq_proj=1 lockq_every=10"
F32="lockq_M=2 lockq_n=10 fmodes=2"; FW="lockq_M=4 lockq_n=7 fmodes=4"
SEEDS=(0 16294208416658607535 7960286522194355700 487617019471545679 17909611376780542444 1961750202426094747)
cd /home/d/code/scp/v97
nrun() { ps -eo args | awk '$1 ~ /freecell_k3$/' | wc -l; }
job() { local f=$1; shift; [ -e $f ] && return; while [ $(nrun) -ge $(cat $S/supp/maxj) ]; do sleep 10; done; sh -c "$* > $f 2>&1" & sleep 2; }
out=$S/supp
for k in 1 2 3 4 5; do sd=${SEEDS[$k]}
  job $out/A_M4_d_s$k.out "$B $LAW exp=pulse $COM $FW  fseed_mode=-1 lockq_dx=-5.0 amp=1.0 sigma=2.2 lockq_seed=$sd"
  job $out/A_M2_d_s$k.out "$B $LAW exp=pulse $COM $F32 fseed_mode=-1 lockq_dx=-7.5 amp=1.0 sigma=2.2 lockq_seed=$sd"
  job $out/A_M4_b_s$k.out "$B $LAW exp=pulse $COM $FW  fseed_mode=-1 lockq_dx=-5.0 amp=4.0 sigma=3.0 lockq_seed=$sd"
  job $out/A_M2_b_s$k.out "$B $LAW exp=pulse $COM $F32 fseed_mode=-1 lockq_dx=-7.5 amp=4.0 sigma=3.0 lockq_seed=$sd"
done
while [ $(nrun) -gt 0 ]; do sleep 30; done
touch $S/supp.done
