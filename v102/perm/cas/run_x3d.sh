#!/bin/bash
# pruned remainder of the 3+2 mobile-hole chains; at most 4 at once
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
items=("mob 40 1 3.0" "mob 40 3 3.0" "nohole 24 0 0" "nohole 24 2 0" "mob 24 1 1.0" "mob 24 1 3.0" "mob 24 3 3.0")
for it in "${items[@]}"; do
  set -- $it
  f=out/x3c_M2_$1_$2_$3_t$4.out
  [ -s $f ] && continue
  while [ $(ps -eo args | awk '$1=="python3" && $3=="lo_dmrg.py" && $4=="2"' | wc -l) -ge 4 ]; do sleep 20; done
  nice -n 10 python3 -I lo_dmrg.py 2 $1 $2 $3 $4 48 > $f 2> out/x3c_M2_$1_$2_$3_t$4.err &
  sleep 5
done
wait
touch out/x3c.done
