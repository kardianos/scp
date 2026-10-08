#!/bin/bash
cd /home/d/code/scp/v102/options/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
printf "a\nb\nc\n" | xargs -P 3 -n 1 sh -c 'python3 -I e9_size.py $0 > out/e9_$0.out 2> out/e9_$0.err'
touch out/e9.done
