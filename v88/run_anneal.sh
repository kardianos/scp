#!/bin/bash
# v88 annealing sweep: 3 seeds x {quench, thermal, detune-control}
cd /home/d/code/scp/v88
for sched in quench thermal detune; do
  for seed in 11111 22222 33333; do
    echo "### $sched seed=$seed"
    timeout 3000 ./fabric_harmonic $sched 3 $seed 2>&1 | grep -E "^  \[|DONE" | tail -6
  done
done
