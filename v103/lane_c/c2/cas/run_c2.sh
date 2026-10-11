#!/bin/sh
# v103 c2 registered grid: 16 sets x L = 12, 16 (PREREG §3). Heavy runs (L = 16, a off, c on) first.
cd "$(dirname "$0")"
mkdir -p out
: > out/run.log
{
  for L in 16 12; do
    for c in 1 0; do
      for a in 0 1; do
        for b in 0 1; do
          for d in 0 1; do
            echo "c2_a${a}b${b}c${c}d${d}_L${L} L=$L a=$a b=$b c=$c d=$d"
          done
        done
      done
    done
  done
} | OMP_NUM_THREADS=4 xargs -P 3 -L 1 sh -c 'n=$0; ./gap1d mode=run "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/run.log'
touch out/run.done
