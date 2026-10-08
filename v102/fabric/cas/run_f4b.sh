#!/bin/bash
cd /home/d/code/scp/v102/fabric/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
printf "24 1\n24 3\n40 1\n40 3\n56 1\n56 3\n" | xargs -P 6 -n 2 sh -c 'python3 -I f4b_dmrg.py $0 $1 64 hole > out/f4b_$0_$1.out 2> out/f4b_$0_$1.err'
touch out/f4b.done
