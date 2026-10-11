#!/bin/bash
# The four 3+2 c-on jobs of run_grid.sh, rerun with a larger rate cache
# (cache=2e8 entries).  The cache changes speed only: a tiny and a large
# cache give byte-identical output (3+2 ac exp=hold, checked before this
# rerun).  Same source as gx; binary gx_c adds only the cache= key.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for r in c ac bc abc; do
  ( ./gx_c M=2 n=6 k=4 rules=$r exp=all ntraj=200 cache=200000000 > out/g_M2_$r.out 2> out/g_M2_$r.err; echo "g_M2_$r done $(date +%T)" >> out/grid.log ) &
done
wait
touch out/grid_c32.done
