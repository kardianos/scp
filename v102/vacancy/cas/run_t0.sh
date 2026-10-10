#!/bin/bash
# T0 (../PREREG.md §1): vx against six recorded cases.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=4
X=./vx
$X model=lo M=4 n=6 geom=ring vac=mobile t=0   sz=1 iz=1 k=8 > out/t0_lo_W6_t0.out 2>&1
$X model=lo M=4 n=6 geom=ring vac=mobile t=0   sz=3 iz=1 k=2 > out/t0_lo_W6_t0_31.out 2>&1
$X model=lo M=4 n=6 geom=ring vac=mobile t=0.1 sz=1 iz=1 k=8 > out/t0_lo_W6_t0.1.out 2>&1
$X model=lo M=2 n=8 geom=ring vac=mobile t=0   sz=1 k=9 > out/t0_lo_32n8_t0.out 2>&1
$X model=lo M=2 n=8 geom=ring vac=mobile t=0.1 sz=1 k=9 > out/t0_lo_32n8_t0.1.out 2>&1
$X model=ao M=2 n=4 geom=ring sz=0 k=4 > out/t0_ao_32n4.out 2>&1
$X model=ao M=4 n=4 geom=ring sz=0 iz=0 k=4 > out/t0_ao_W4.out 2>&1
$X model=ao M=2 n=8 geom=ring sz=0 cmax=3 k=2 > out/t0_ao_32n8_e0.out 2>&1
$X model=ao M=2 n=8 geom=ring sz=2 cmax=3 k=2 > out/t0_ao_32n8_rip.out 2>&1
$X model=ao M=2 n=8 geom=ring sz=0 cmax=3 k=3 nvec=2 pen=1 > out/t0_ao_32n8_sing.out 2>&1
touch out/t0.done
