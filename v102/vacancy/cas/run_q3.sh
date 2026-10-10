#!/bin/bash
# Q3 (../PREREG.md §4) and the Q2 ring-scale runs (run_q2r.sh's list): every
# job is one vx run.  Small jobs run one thread each, 12 at once; the W LO
# ring n = 8 jobs (D = 3.2e6) run afterwards with 12 threads.  A job whose
# output exists and is non-empty is skipped.
cd "$(dirname "$0")"
jobs=$(mktemp); big=$(mktemp)
add() { local name=$1; shift; [ -s out/$name.out ] && return; echo "$name $*" >> $jobs; }
addbig() { local name=$1; shift; [ -s out/$name.out ] && return; echo "$name $*" >> $big; }
for model in lo ao; do
  for M in 2 4; do
    for geom in ring open; do
      case "$model$M$geom" in
        lo2ring) sizes="6 8 10 12 14";; lo2open) sizes="8 10 12 14";;
        lo4ring) sizes="4 6 8";;        lo4open) sizes="5 6 7 8";;
        ao2ring) sizes="4 6 8";;        ao2open) sizes="6 8";;
        ao4ring) sizes="4 6";;          ao4open) sizes="5 6";;
      esac
      for n in $sizes; do
        if [ $M = 2 ]; then s0="sz=0"; s1="sz=2"; else s0="sz=0 iz=0"; s1="sz=2 iz=2"; fi
        b=q3_${model}_M${M}_${geom}_n$n
        f=add; [ "$model$M$geom$n" = lo4ring8 ] && f=addbig
        $f ${b}_e0 model=$model M=$M n=$n geom=$geom $s0 k=2 nvec=1 corr=1 cs=1
        $f ${b}_rip model=$model M=$M n=$n geom=$geom $s1 k=2 nvec=1
        $f ${b}_sing model=$model M=$M n=$n geom=$geom $s0 k=3 nvec=2 pen=1
      done
    done
  done
done
for model in lo ao; do
  add q2r_${model}_M2_n8_static  model=$model M=2 n=8 geom=ring vac=static vsite=0 sz=1 k=2 nvec=1 corr=1 cs=1
  add q2r_${model}_M2_n8_mobile  model=$model M=2 n=8 geom=ring vac=mobile t=0.1 sz=1 k=2 nvec=1 corr=1
  add q2r_${model}_M4_n6_static  model=$model M=4 n=6 geom=ring vac=static vsite=0 sz=1 iz=1 k=2 nvec=1 corr=1 cs=1
  add q2r_${model}_M4_n6_mobile  model=$model M=4 n=6 geom=ring vac=mobile t=0.1 sz=1 iz=1 k=2 nvec=1 corr=1
done
addbig q3v_lo_M4_n8_static model=lo M=4 n=8 geom=ring vac=static vsite=0 sz=1 iz=1 k=12 nvec=12
add q3v_lo_M4_n8_static_rip model=lo M=4 n=8 geom=ring vac=static vsite=0 sz=3 iz=3 k=2 nvec=1
echo "# $(wc -l < $jobs) small jobs, $(wc -l < $big) big jobs  $(date +%T)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
xargs -P 12 -L 1 sh -c 'n=$0; ./vx "$@" > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)"' < $jobs
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12
while read -r n args; do ./vx $args > out/$n.out 2> out/$n.err; echo "$n done $(date +%T)"; done < $big
rm -f $jobs $big
touch out/q3.done
