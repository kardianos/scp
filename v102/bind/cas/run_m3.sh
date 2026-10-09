#!/bin/bash
# bind M1 (reported) and M3a: exact rings
cd /home/d/code/scp/v102/bind/cas
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
python3 -I m1_link.py 4 4 2 > out/m1_W_n4.out 2>&1
python3 -I m1_link.py 2 8 3 > out/m1_32_n8.out 2>&1
for p in zero core score lower upper; do
  python3 -I m3_yield.py 2 10 $p > out/m3a_M2_n10_$p.out 2>&1
  python3 -I m3_yield.py 4 6 $p > out/m3a_M4_n6_$p.out 2>&1
done
touch out/m3.done
