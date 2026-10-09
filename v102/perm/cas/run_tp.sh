#!/bin/bash
# TeNPy long chains (charge-conserving DMRG), the registered lengths; at most 8 at once
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
items=(
"4 hole 32 1,1 0 128" "4 nohole 32 0,0 0 128" "4 nohole 32 2,2 0 128"
"4 hole 24 1,1 0 128" "4 nohole 24 0,0 0 128" "4 nohole 24 2,2 0 128"
"4 hole 16 1,1 0 128" "4 nohole 16 0,0 0 128" "4 nohole 16 2,2 0 128"
"4 hole 32 3,3 0 128" "4 hole 24 3,3 0 128" "4 nohole 24 2,0 0 128"
"2 mob 40 1 1.0 64" "2 mob 40 3 1.0 64" "2 mob 24 1 1.0 64" "2 mob 24 3 1.0 64"
"2 mob 40 1 3.0 64" "2 mob 40 3 3.0 64" "2 mob 24 1 3.0 64" "2 mob 24 3 3.0 64"
"2 mob 40 1 0.1 64" "2 mob 24 1 0.1 64" "2 mob 24 3 0.1 64"
"2 nohole 24 0 0 64" "2 nohole 24 2 0 64" "2 nohole 40 0 0 64" "2 nohole 40 2 0 64"
)
for it in "${items[@]}"; do
  set -- $it
  f=out/tp_M$1_$2_$3_$4_t$5_chi$6.out
  [ -s $f ] && continue
  while [ $(ps -eo args | awk '$3=="tp_chain.py"' | wc -l) -ge 8 ]; do sleep 20; done
  $PY -I tp_chain.py $1 $2 $3 $4 $5 $6 > $f 2> ${f%.out}.err &
  sleep 3
done
wait
touch out/tp.done
