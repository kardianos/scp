#!/bin/bash
# W items still needed; penalty margins: (2,2) and (3,3) targets cost 2 J0 to miss.
cd /home/d/code/scp/v102/perm/cas
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
items=("nohole 20 2,2 32" "hole 20 1,1 32" "hole 20 3,3 32" "nohole 12 0,0 32" "nohole 12 2,2 32" "hole 12 1,1 32" "hole 12 3,3 32" "nohole 12 0,0 48" "nohole 12 2,2 48")
for it in "${items[@]}"; do
  set -- $it
  f=out/x12_W_$1_$2_$3_chi$4.out
  if [ -s $f ] || pgrep -f "^python3 -I lo_dmrg.py 4 $1 $2 $3 0 $4\$" >/dev/null; then continue; fi
  while [ $(pgrep -fc "^python3 -I lo_dmrg.py 4 ") -ge 8 ]; do sleep 20; done
  python3 -I lo_dmrg.py 4 $1 $2 $3 0 $4 > $f 2> out/x12_W_$1_$2_$3_chi$4.err &
  sleep 5
done
wait
touch out/x12b.done
