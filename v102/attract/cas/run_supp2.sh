#!/bin/bash
# supplementary (not registered): do pairs still bind where the ladder appears?
cd /home/d/code/scp/v102/attract/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
TP=../../perm/cas/tp_chain.py
run() { f=out/tp_M$1_$2_$3_$4_t$5_chi$6.out; [ -s $f ] && return
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py/' | wc -l) -ge 8 ]; do sleep 10; done
  $PY -I $TP $1 $2 $3 $4 $5 $6 > $f 2> ${f%.out}.err & sleep 2; }
for t in 3.0 5.0 7.0; do for L in 24 40; do run 2 mob2 $L 0 $t 64; run 2 mob2 $L 2 $t 64; done; run 2 mob 25 1 $t 64; run 2 mob 41 1 $t 64; done
run 4 mob2 16 0,0 7.0 128; run 4 mob2 16 2,0 7.0 128; run 4 mob 17 1,1 7.0 128
wait
touch out/supp2.done
