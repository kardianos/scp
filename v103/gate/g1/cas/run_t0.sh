#!/bin/sh
# v103 gate g1 tool gates (PREREG §6). Must end "T0 PASS" before any registered run.
cd "$(dirname "$0")"
mkdir -p out
gcc -O3 -march=native -Wall -Wextra -o geo3d geo3d.c -l:liblapack.so.3 -l:libopenblas.so.0 -lm || exit 1
OPENBLAS_NUM_THREADS=1 ./geo3d mode=t0 > out/t0.out 2> out/t0.err
tail -1 out/t0.out
