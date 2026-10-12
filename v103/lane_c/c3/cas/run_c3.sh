#!/bin/sh
# v103 lane c step c3 registered grid: 8 sets x L = 8, 12, 16 x placements k = 0..3 (PREREG §3), 96 processes.
# One single-threaded process per (set, L, k); heavy runs (L = 16, link rule on) first, 12 at a time.
# Writes out/c3_a?b?c?_L?_k?.out/.err, appends to out/run.log, touches out/run.done.
# Run only after out/t0.out ends "T0 PASS".
cd "$(dirname "$0")"
tail -1 out/t0.out 2>/dev/null | grep -q "T0 PASS" || { echo "T0 has not passed: run ./run_t0.sh first"; exit 1; }
mkdir -p out
: > out/run.log
{
  for L in 16 12 8; do
    for b in 1 0; do
      for a in 1 0; do
        for c in 1 0; do
          for k in 0 1 2 3; do
            echo "c3_a${a}b${b}c${c}_L${L}_k${k} L=$L a=$a b=$b c=$c k=$k"
          done
        done
      done
    done
  done
} | xargs -P 12 -L 1 sh -c 'n=$0; ./c3pair mode=run "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/run.log'
touch out/run.done
