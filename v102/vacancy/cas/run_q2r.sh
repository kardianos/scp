#!/bin/bash
# Q2 ring scale (../PREREG.md §3) and Q4's vacancy overlaps (§5): one vacancy
# on the largest ring each AO cell runs (3+2 n = 8, W n = 6), static (t = 0,
# vacancy at cell 0) and mobile (t = 0.1 J0), LO and AO.  The background's
# correlations at the same n come from the Q3 ground-state runs.
cd "$(dirname "$0")"
X=./vx
run() { local name=$1; shift; [ -s out/$name.out ] && return; $X "$@" > out/$name.out 2> out/$name.err; echo "$name $(date +%T)"; }
for model in lo ao; do
  run q2r_${model}_M2_n8_static  model=$model M=2 n=8 geom=ring vac=static vsite=0 sz=1 k=2 nvec=1 corr=1 cs=1
  run q2r_${model}_M2_n8_mobile  model=$model M=2 n=8 geom=ring vac=mobile t=0.1 sz=1 k=2 nvec=1 corr=1
  run q2r_${model}_M4_n6_static  model=$model M=4 n=6 geom=ring vac=static vsite=0 sz=1 iz=1 k=2 nvec=1 corr=1 cs=1
  run q2r_${model}_M4_n6_mobile  model=$model M=4 n=6 geom=ring vac=mobile t=0.1 sz=1 iz=1 k=2 nvec=1 corr=1
done
# Q3's "ripple beside a particle": the largest LO W ring with one static vacancy,
# the lowest levels of the vacancy's own sector, labelled
run q3v_lo_M4_n8_static model=lo M=4 n=8 geom=ring vac=static vsite=0 sz=1 iz=1 k=12 nvec=12
run q3v_lo_M4_n8_static_rip model=lo M=4 n=8 geom=ring vac=static vsite=0 sz=3 iz=3 k=2 nvec=1
touch out/q2r.done
