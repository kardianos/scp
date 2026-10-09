#!/bin/bash
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
{
for L in 40 24; do
  echo "nohole $L 0 0"; echo "nohole $L 2 0"
  echo "hole $L 1 0"; echo "hole $L 3 0"
  for t in 0.1 1.0 3.0; do echo "mob $L 1 $t"; echo "mob $L 3 $t"; done
done
} | xargs -P 4 -L 1 sh -c 'nice -n 10 python3 -I lo_dmrg.py 2 $0 $1 $2 $3 48 > out/x3c_M2_$0_$1_$2_t$3.out 2> out/x3c_M2_$0_$1_$2_t$3.err'
touch out/x3c.done
