#!/bin/bash
# Q2 (../PREREG.md §3) long chains: the W background's correlation (L_bulk 32,
# and 24 as a report), the mobile W vacancy (L_bulk 24, the 14-state site),
# and the 3+2 all-order static vacancy at g/Delta = 0.05 and 0.2 (L_bulk 40).
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=3 OMP_NUM_THREADS=3
TP=~/.venvs/scp-tenpy/bin/python
P3=../../perm/cas
A2=../../attract/cas
( $TP -I vq2_tp.py 4 32 128 > out/q2_W_nohole_32.out 2> out/q2_W_nohole_32.err ) &
( $TP -I vq2_tp.py 4 24 128 > out/q2_W_nohole_24.out 2> out/q2_W_nohole_24.err ) &
( $TP -I $P3/tp_chain.py 4 mob 24 1,1 0.1 128 > out/q2_W_mob_24_t0.1.out 2> out/q2_W_mob_24_t0.1.err ) &
( python3 -I $A2/a2_allorder.py 0.05 40 1 64 hole > out/q2_32_ao_hole_40_g0.05.out 2> out/q2_32_ao_hole_40_g0.05.err ) &
( python3 -I $A2/a2_allorder.py 0.2 40 1 64 hole > out/q2_32_ao_hole_40_g0.2.out 2> out/q2_32_ao_hole_40_g0.2.err ) &
wait
touch out/q2.done
