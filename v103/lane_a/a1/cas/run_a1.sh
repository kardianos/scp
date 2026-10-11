#!/bin/bash
# The registered a1 grid (../PREREG.md §2): eight rule sets x L = 6, 8; n̄ = 2.
#   a: cap 2 (on) / 3 (off);  b: M = 2 (on) / 1;  c: Jx = 5 (on) / 0.
cd "$(dirname "$0")"
mkdir -p out
BX=../../../tools/bx
jobs=$(mktemp)
for L in 8 6; do for r in none a b c ab ac bc abc; do
  C=3; M=1; J=0
  case $r in *a*) C=2;; esac; case $r in *b*) M=2;; esac; case $r in *c*) J=5;; esac
  echo "a1_${r}_L$L L=$L M=$M C=$C nb=2 Jx=$J" >> $jobs
done; done
echo "# $(wc -l < $jobs) jobs $(date +%T)" >> out/run.log
OMP_NUM_THREADS=3 xargs -P 4 -L 1 sh -c 'n=$0; '"$BX"' mode=a1 "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/run.log' < $jobs
rm -f $jobs
touch out/run.done
