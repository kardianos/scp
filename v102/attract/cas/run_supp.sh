#!/bin/bash
# supplementary (not registered): W moving hole at t = 7 J0, where 3+2 binds its S = 3/2 level
cd /home/d/code/scp/v102/attract/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=~/.venvs/scp-tenpy/bin/python
for it in "16 1,1" "16 3,3" "24 1,1" "24 3,3"; do
  set -- $it
  f=out/tp_M4_mob_$1_$2_t7.0_chi128.out; [ -s $f ] && continue
  while [ $(ps -eo args | awk '$3 ~ /tp_chain.py|a2_allorder.py/' | wc -l) -ge 9 ]; do sleep 15; done
  $PY -I ../../perm/cas/tp_chain.py 4 mob $1 $2 7.0 128 > $f 2> ${f%.out}.err &
  sleep 2
done
wait
touch out/supp.done
