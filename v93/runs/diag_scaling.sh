#!/bin/bash
cd "$(dirname "$0")/.."
echo "=== amp_nat scaling (seed 111, bath=1) ==="
for a in 1 3 10 30; do
  r=$(./freecell exp=blob bath=1 L=16 dt=0.02 T=80 seed=111 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=$a 2>/dev/null | grep "RESULT blob_drift")
  sp=$(echo "$r" | sed -n 's/.*speed=\([0-9.e+-]*\).*cos_to_kdir=\([0-9.e+-]*\).*/\1  cos=\2/p')
  echo "bath=1 amp_nat=$a  $sp"
done
echo "=== bath=0 (vacuum blob) seed 111 ==="
for a in 1 10 30; do
  r=$(./freecell exp=blob bath=0 L=16 dt=0.02 T=80 seed=111 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=$a 2>/dev/null | grep "RESULT blob_drift")
  sp=$(echo "$r" | sed -n 's/.*speed=\([0-9.e+-]*\).*cos_to_kdir=\([0-9.e+-]*\).*/\1  cos=\2/p')
  echo "bath=0 amp_nat=$a  $sp"
done
echo "=== done ==="
