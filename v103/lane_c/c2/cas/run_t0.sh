#!/bin/sh
# v103 c2 tool gates (PREREG §4)
cd "$(dirname "$0")"
gcc -O3 -march=native -fopenmp -Wall -Wextra -o gap1d gap1d.c -l:liblapack.so.3 -l:libopenblas.so.0 -lm || exit 1
OMP_NUM_THREADS=8 ./gap1d mode=t0 > out/t0.out 2> out/t0.err
tail -1 out/t0.out
