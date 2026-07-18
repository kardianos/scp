#!/bin/bash
# v77 NM Round 1 runner
set -euo pipefail
cd "$(dirname "$0")"
python3 offline_r1_multilock.py
# Optional full SOR (slower):
# python3 sandbox_r1_multilock.py --N 24 --iters 300 --charge_mode same_sign
echo "NM offline r1 complete"
