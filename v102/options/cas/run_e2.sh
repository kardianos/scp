#!/bin/bash
cd /home/d/code/scp/v102/options/cas
for r in T0 T2 T3 T4; do for c in lock geo held; do echo "$r $c"; done; done | \
  xargs -P 12 -n 2 sh -c 'python3 -I e23_boundary.py E2 $0 $1 > out/e2_$0_$1.out 2> out/e2_$0_$1.err'
touch out/e2.done
