#!/bin/bash
# The registered grid (../PREREG.md §3): every rule set x {2-D, 3-D} x {E2, E4},
# plus the 1-D ring for the c-on sets.  One process per run, single core.
cd "$(dirname "$0")"
mkdir -p out
jobs=$(mktemp)
S="seeds=8 sweeps=1000000 warm=10000"
for r in none a b c ab ac bc abc; do
  case $r in *c*) l2="lat=square L=32"; l3="lat=cubic L=12";; *) l2="lat=kagome L=18"; l3="lat=pyro L=5";; esac
  for D in 2 4; do
    echo "c1_3d_${r}_E$D $l3 rules=$r D=$D $S" >> $jobs
    echo "c1_2d_${r}_E$D $l2 rules=$r D=$D $S" >> $jobs
    case $r in *c*) echo "c1_1d_${r}_E$D lat=ring L=128 rules=$r D=$D $S" >> $jobs;; esac
  done
done
echo "# $(wc -l < $jobs) jobs $(date +%T)" >> out/grid.log
xargs -P 12 -L 1 sh -c 'n=$0; ./chain mode=run "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/grid.log' < $jobs
rm -f $jobs
touch out/grid.done
