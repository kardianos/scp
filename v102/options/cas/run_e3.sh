#!/bin/bash
cd /home/d/code/scp/v102/options/cas
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for rc in "T0 lock" "T0 geo" "T0 held" "T4 lock" "T4 held"; do for x in X0 X1; do for e in cold res; do echo "$rc $x $e"; done; done; done | \
  xargs -P 8 -n 4 sh -c 'python3 -I e23_boundary.py E3 $0 $1 $2 $3 > out/e3_$0_$1_$2_$3.out 2> out/e3_$0_$1_$2_$3.err'
touch out/e3.done
