#!/bin/bash
# bind M2: cluster energetics (TeNPy chains via ../../perm/cas/tp_chain.py)
cd /home/d/code/scp/v102/bind/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
TP=../../perm/cas/tp_chain.py
MAXJ=${MAXJ:-7}
tp() {  # M nbulk Nh target chi p11
  local M=$1 nb=$2 nh=$3 tg=$4 chi=$5 p=$6 mode Lb
  if [ $nh -eq 0 ]; then mode=nohole; Lb=$nb; elif [ $nh -eq 1 ]; then mode=mob; Lb=$((nb-1)); else mode=mob$nh; Lb=$((nb-nh)); fi
  local f=out/m2_M${M}_n${nb}_h${nh}_${tg}_${p}.out; [ -s $f ] && return
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py/' | wc -l) -ge $MAXJ ]; do sleep 10; done
  $PY -I $TP $M $mode $Lb $tg 0.1 $chi 20 $p > $f 2> ${f%.out}.err &
  sleep 1
}
sec() {  # M Nh -> minimal sector
  if [ $1 -eq 4 ]; then [ $(($2 % 2)) -eq 0 ] && echo 0,0 || echo 1,1; else [ $(($2 % 2)) -eq 0 ] && echo 0 || echo 1; fi
}
# 3+2 first (fast)
for nb in 26 42; do tp 2 $nb 0 0 64 zero; for p in zero core lower upper; do for nh in 1 2 3 4 5 6; do tp 2 $nb $nh $(sec 2 $nh) 64 $p; done; done; done
for nh in 1 2 3 4 5 6; do tp 2 26 $nh $(sec 2 $nh) 64 score; done
for p in zero upper; do for nh in 2 3 4 5 6; do [ $((nh % 2)) -eq 0 ] && tp 2 26 $nh 2 64 $p || tp 2 26 $nh 3 64 $p; done; done
# W
for nb in 18 26; do tp 4 $nb 0 0,0 128 zero; for p in zero core lower upper; do for nh in 1 2 3 4 5; do tp 4 $nb $nh $(sec 4 $nh) 128 $p; done; done; done
for nh in 1 2 3 4 5; do tp 4 18 $nh $(sec 4 $nh) 128 score; done
for p in zero core lower upper; do tp 4 18 6 0,0 128 $p; done
for p in zero upper; do for nh in 2 3 4 5; do
  if [ $((nh % 2)) -eq 0 ]; then tp 4 18 $nh 2,0 128 $p; tp 4 18 $nh 0,2 128 $p; else tp 4 18 $nh 3,1 128 $p; tp 4 18 $nh 1,3 128 $p; fi
done; done
wait
touch out/m2.done
