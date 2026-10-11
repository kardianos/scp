#!/bin/bash
# The registered grid (../PREREG.md §3): one process per (fabric, rule set)
# running X4, X1, X2, X3, X5 in order (each sector is diagonalized once per
# process), plus the X5 side variant for the c-off rule sets.  W processes
# start first (they are the long pole).  A job whose output ends in CHECK
# is skipped.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
jobs=$(mktemp)
add() { local name=$1; shift; grep -q '^CHECK' out/$name.out 2>/dev/null && return; echo "$name $*" >> $jobs; }
for r in none a b c ab ac bc abc; do add g_M4_$r M=4 n=5 k=3 rules=$r exp=all ntraj=100; done
for r in none a b ab; do add gs_M4_$r M=4 n=5 k=3 rules=$r exp=vac side=1 ntraj=100; done
for r in none a b c ab ac bc abc; do add g_M2_$r M=2 n=6 k=4 rules=$r exp=all ntraj=200; done
for r in none a b ab; do add gs_M2_$r M=2 n=6 k=4 rules=$r exp=vac side=1 ntraj=200; done
echo "# $(wc -l < $jobs) jobs  $(date +%T)"
xargs -P 14 -L 1 sh -c 'n=$0; ./gx "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)"' < $jobs
rm -f $jobs
touch out/grid.done
