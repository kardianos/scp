#!/bin/bash
cd "$(dirname "$0")/.."
echo "=== LINEAR fd current, seed-robustness ==="
for s in 20260802 314159; do
  r=$(./freecell exp=blob bath=1 L=16 dt=0.02 T=80 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 seed=$s amp_nat=2 amp_logate=1 p1_meter=1 2>/dev/null | grep "RESULT p1x")
  fd=$(echo "$r" | sed -E 's/.*fd=([+0-9.e-]+).*/\1/')
  echo "seed=$s  p1x: $r"
done
echo "=== done ==="
