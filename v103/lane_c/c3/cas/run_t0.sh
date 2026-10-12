#!/bin/sh
# v103 lane c step c3 tool gates (PREREG §6). Must end "T0 PASS" before any registered run.
cd "$(dirname "$0")"
mkdir -p out
gcc -O3 -march=native -Wall -Wextra -o c3pair c3pair.c -lm || exit 1
./c3pair mode=t0 > out/t0.out 2> out/t0.err
tail -1 out/t0.out
