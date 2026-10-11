#!/bin/bash
# T0 (../PREREG.md §4): gates for every rule set in both fabrics, and the
# leading-order lock energies of ../../vacancy/cas/vx.c for gate 3.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8
for r in none a b c ab ac bc abc; do
  ./gx M=2 n=6 k=4 rules=$r exp=t0 > out/t0_M2_$r.out
  ./gx M=4 n=5 k=3 rules=$r exp=t0 > out/t0_M4_$r.out
done
../../vacancy/cas/vx model=lo M=2 n=4 geom=open sz=0 k=2 nvec=1 > out/t0_vx_M2.out 2>/dev/null
../../vacancy/cas/vx model=lo M=4 n=3 geom=open sz=0 iz=0 k=2 nvec=1 > out/t0_vx_M4.out 2>/dev/null
touch out/t0.done
