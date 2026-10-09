#!/bin/bash
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
{
for L in 16 12 8; do
  echo "nohole $L 0,0"; echo "nohole $L 2,2"; echo "hole $L 1,1"; echo "hole $L 3,3"
done
} | xargs -P 6 -L 1 sh -c 'python3 -I lo_dmrg.py 4 $0 $1 $2 0 64 > out/x12_W_$0_$1_${2}_chi64.out 2> out/x12_W_$0_$1_${2}_chi64.err'
touch out/x12d.done
