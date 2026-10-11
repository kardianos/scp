#!/bin/bash
# The registered b1 grid (../PREREG.md §2): eight rule sets, ladders in n̄.
#   a: cap n̄ (on) / none (off);  b: U = 4 (on) / 0;  c: M = 2 (on) / 1.
cd "$(dirname "$0")"
mkdir -p out
BX=../../../tools/bx
jobs=$(mktemp)
add() { local L=$1 r=$2 nb=$3; C=0; U=0; M=1
  case $r in *a*) C=$nb;; esac; case $r in *b*) U=4;; esac; case $r in *c*) M=2;; esac
  echo "b1_${r}_L${L}_n$nb L=$L M=$M C=$C nb=$nb U=$U" >> $jobs; }
# heavy runs first
for r in none b; do add 6 $r 8; done
for r in c bc; do add 6 $r 3; done
for nb in 8 6 4 3 2 1; do for r in none a b c ab ac bc abc; do add 4 $r $nb; done; done
for nb in 6 4 3 2 1; do for r in none b; do add 6 $r $nb; done; done
for nb in 8 6 4 3 2 1; do for r in a ab; do add 6 $r $nb; done; done
for nb in 2 1; do for r in c bc; do add 6 $r $nb; done; done
for nb in 3 2 1; do for r in ac abc; do add 6 $r $nb; done; done
echo "# $(wc -l < $jobs) jobs $(date +%T)" >> out/run.log
OMP_NUM_THREADS=2 xargs -P 5 -L 1 sh -c 'n=$0; '"$BX"' mode=b1 "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/run.log' < $jobs
rm -f $jobs
touch out/run.done
