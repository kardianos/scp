#!/bin/bash
# bind, supplementary (after R12 ruled singlet core): score at the longer lengths
cd /home/d/code/scp/v102/bind/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
TP=../../perm/cas/tp_chain.py
tp() {  # M nbulk Nh target chi p11
  local M=$1 nb=$2 nh=$3 tg=$4 chi=$5 p=$6 mode Lb
  if [ $nh -eq 1 ]; then mode=mob; Lb=$((nb-1)); else mode=mob$nh; Lb=$((nb-nh)); fi
  local f=out/m2_M${M}_n${nb}_h${nh}_${tg}_${p}.out; [ -s $f ] && return
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py/' | wc -l) -ge ${MAXJ:-3} ]; do sleep 10; done
  $PY -I $TP $M $mode $Lb $tg 0.1 $chi 20 $p > $f 2> ${f%.out}.err &
  sleep 1
}
for nh in 1 2 3 4 5 6; do [ $((nh % 2)) -eq 0 ] && tg=0 || tg=1; tp 2 42 $nh $tg 64 score; done
wait
until [ -e out/m2.done ]; do sleep 30; done
MAXJ=7
for nh in 1 2 3 4 5; do [ $((nh % 2)) -eq 0 ] && tg=0,0 || tg=1,1; tp 4 26 $nh $tg 128 score; done
wait
touch out/supp.done
