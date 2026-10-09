#!/bin/bash
cd /home/d/code/scp/v102/perm/cas
for a in "2 6" "2 8" "2 10" "4 4" "4 6"; do set -- $a; OPENBLAS_NUM_THREADS=1 python3 -I x3_rings.py $1 $2 > out/x3_rings_M$1_n$2.out 2> out/x3_rings_M$1_n$2.err & done; wait
touch out/x3r.done
