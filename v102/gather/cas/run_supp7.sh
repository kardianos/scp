#!/bin/bash
# The registered 3+2 supplement (../PREREG.md §2: n = 7, Λ = 0 … 4, run
# only once the main grid has finished), labelled and not scored.
# Its run counts were not registered.  It uses the main grid's 3+2
# dilute count (200) and, for dense runs, the cost rule §6 registered for
# W (12 trajectories at T = 1000, averaged over [500, 1000]): the largest
# n = 7 sector (6405) is W's scale.  T0 runs first, for all eight rule
# sets, with gate 3 against vx.c's LO open chain of 5 cells.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
for r in none a b c ab ac bc abc; do
  s=0; case $r in *c*) ;; *) s=1;; esac
  echo "$r $s"
done | xargs -P 8 -L 1 sh -c './gx M=2 n=7 k=5 rules=$0 exp=t0 side=$1 > out/t0_s7_$0.out 2> out/t0_s7_$0.err'
../../vacancy/cas/vx model=lo M=2 n=5 geom=open sz=0 k=2 nvec=1 > out/t0_vx_s7.out 2>/dev/null
touch out/t0_s7.done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
jobs=$(mktemp)
for r in none a b c ab ac bc abc; do echo "s7_M2_$r M=2 n=7 k=5 rules=$r exp=all ntraj=200 ntd=12 Td=1000" >> $jobs; done
for r in none a b ab; do echo "s7s_M2_$r M=2 n=7 k=5 rules=$r exp=vac side=1 ntraj=200 ntd=12 Td=1000" >> $jobs; done
echo "# supplement $(wc -l < $jobs) jobs  $(date +%T)" >> out/supp7.log
xargs -P 12 -L 1 sh -c 'n=$0; ./gx "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)" >> out/supp7.log' < $jobs
rm -f $jobs
touch out/supp7.done
