#!/bin/sh
# v103 gate g1 registered grid: 8 sets x L = 6, 8, 10, 12, 16 (PREREG §4).
# One single-threaded process per (set, L); heavy runs (L = 16, link rule on) first.
# Writes out/g1_a?b?c?_L?.out/.err, appends to out/run.log, touches out/run.done.
# Run only after out/t0.out ends "T0 PASS".
cd "$(dirname "$0")"
tail -1 out/t0.out 2>/dev/null | grep -q "T0 PASS" || { echo "T0 has not passed: run ./run_t0.sh first"; exit 1; }
mkdir -p out
: > out/run.log
{
  for L in 16 12 10 8 6; do
    for b in 1 0; do
      for a in 1 0; do
        for c in 1 0; do
          echo "g1_a${a}b${b}c${c}_L${L} L=$L a=$a b=$b c=$c"
        done
      done
    done
  done
} | OPENBLAS_NUM_THREADS=1 xargs -P 12 -L 1 sh -c 'n=$0; ./geo3d mode=run "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/run.log'
touch out/run.done
