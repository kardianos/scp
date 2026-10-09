#!/bin/bash
# attract batch; TeNPy chains via ../../perm/cas/tp_chain.py, all-order 3+2 via a2_allorder.py
cd /home/d/code/scp/v102/attract/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
TP=../../perm/cas/tp_chain.py
tp() {  # M mode Lb target t chi
  f=out/tp_M$1_$2_$3_$4_t$5_chi$6.out; [ -s $f ] && return
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py|a2_allorder.py/' | wc -l) -ge 8 ]; do sleep 15; done
  $PY -I $TP $1 $2 $3 $4 $5 $6 > $f 2> ${f%.out}.err &
  sleep 2
}
a2() {  # gr Lb twoSz mode
  f=out/a2_g$1_$4_$2_$3.out; [ -s $f ] && return
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py|a2_allorder.py/' | wc -l) -ge 8 ]; do sleep 15; done
  python3 -I a2_allorder.py $1 $2 $3 64 $4 > $f 2> ${f%.out}.err &
  sleep 2
}
# B1 / A1: the moving W hole
for t in 0.1 1.0 3.0; do for L in 16 24; do for tg in 1,1 3,3 3,1; do tp 4 mob $L $tg $t 128; done; done; done
# A3: 3+2 two holes, static pairs (relative to d = 10) and mobile pairs (with E1, E0 at matching size)
for d in 2 3 4 5 6 10; do for tg in 0 2; do tp 2 pair$d 24 $tg 0 64; done; done
for t in 0.1 1.0; do for L in 24 40; do for tg in 0 2; do tp 2 mob2 $L $tg $t 64; done; done; done
for t in 0.1 1.0; do tp 2 mob 25 1 $t 64; tp 2 mob 41 1 $t 64; done
tp 2 nohole 26 0 0 64; tp 2 nohole 42 0 0 64
# A1: 3+2 at larger t
for t in 5.0 7.0; do for L in 24 40; do for tg in 1 3; do tp 2 mob $L $tg $t 64; done; done; done
# A3: W static pairs
for d in 2 3 4 6 10; do for tg in 0,0 2,0 2,2; do tp 4 pair$d 24 $tg 0 128; done; done
# A2: 3+2 all-order at larger g
for g in 0.2 0.3; do for L in 24 40; do a2 $g $L 1 hole; a2 $g $L 3 hole; a2 $g $L 0 nohole; a2 $g $L 2 nohole; done; done
a2 0.1 40 0 nohole; a2 0.1 40 2 nohole; a2 0.1 24 1 hole; a2 0.1 24 3 hole; a2 0.1 24 0 nohole; a2 0.1 24 2 nohole
# A4: W (3,1)
tp 4 hole 24 3,1 0 128; tp 4 hole 32 3,1 0 128
# A3: W mobile pairs at t = 0.1 J0 (heaviest last)
for L in 16 24; do for tg in 0,0 2,0 2,2; do tp 4 mob2 $L $tg 0.1 128; done; done
tp 4 mob 17 1,1 0.1 128; tp 4 mob 25 1,1 0.1 128; tp 4 nohole 18 0,0 0 128; tp 4 nohole 26 0,0 0 128
wait
touch out/attract.done
