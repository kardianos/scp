#!/bin/bash
cd /home/d/code/scp/v102/options/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
( for gr in 0.025 0.05 0.1 0.15 0.2 0.3 0.4 0.5; do echo "static $gr"; echo "margin $gr"; done; ) | \
  xargs -P 6 -n 2 sh -c 'python3 -I e8_g.py $0 $1 > out/e8_$0_$1.out 2> out/e8_$0_$1.err'
touch out/e8.done
