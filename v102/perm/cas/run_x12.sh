#!/bin/bash
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
{
for L in 28 20 12; do
  for tg in 0,0 2,0 2,2; do echo "nohole $L $tg"; done
  for tg in 1,1 3,1 3,3; do echo "hole $L $tg"; done
done
echo "hole 20 1,3"
echo "nohole 20 0,0 48"
echo "nohole 20 2,0 48"
} | xargs -P 8 -L 1 sh -c 'chi=${3:-32}; python3 -I lo_dmrg.py 4 $0 $1 $2 0 $chi > out/x12_W_$0_$1_${2}_chi$chi.out 2> out/x12_W_$0_$1_${2}_chi$chi.err'
touch out/x12.done
