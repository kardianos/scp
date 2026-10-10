#!/bin/bash
# Q3 supplement (labelled, not registered): the next size in two W cells whose
# verdict sits at R ~ 1 — AO ring n = 8 (the cell that reads "line" at n = 4,
# 6) and LO open L_bulk = 7.  Run after the registered grid.
cd "$(dirname "$0")"
until [ -e out/q3.done ]; do sleep 30; done
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12
run() { local name=$1; shift; [ -s out/$name.out ] && return; ./vx "$@" > out/$name.out 2> out/$name.err; echo "$name done $(date +%T)"; }
for spec in "ao ring 8" "lo open 9"; do
  set -- $spec
  b=q3s_$1_M4_$2_n$3
  run ${b}_e0   model=$1 M=4 n=$3 geom=$2 sz=0 iz=0 k=2 nvec=1 m=200
  run ${b}_rip  model=$1 M=4 n=$3 geom=$2 sz=2 iz=2 k=2 nvec=1 m=200
  run ${b}_sing model=$1 M=4 n=$3 geom=$2 sz=0 iz=0 k=3 nvec=2 pen=1 m=200
done
touch out/q3s.done
