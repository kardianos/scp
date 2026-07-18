#!/bin/sh
cd /home/d/code/scp/v76/work/B
# Prefer smaller N first for quick confirmation, then production
python3 sandbox_r3_3d_free.py --N 28 --L 14 --iters 400 --A 0.4 --sigma 1.0 --kappa 1.0 --gamma 0.5
