#!/bin/bash
cd /home/d/code/scp/v102/options/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
printf "0.05 F0\n0.1 F0\n0.2 F0\n0.05 F1\n0.1 F1\n0.2 F1\n0.1 X0\n" | \
  xargs -P 7 -n 2 sh -c 'python3 -I e5_cool.py $0 $1 > out/e5_g$0_$1.out 2> out/e5_g$0_$1.err'
touch out/e5.done
