#!/bin/bash
cd "$(dirname "$0")/.."
echo "=== L1-A face 3: FINE sweep around the amp_nat=2 resonance, 5 seeds ==="
printf "%-8s %-10s %-12s %-10s\n" amp_nat seed speed cos
for a in 1.5 1.75 2.0 2.25; do
  for s in 111 20260802 314159 137035 161803; do
    r=$(./freecell exp=blob bath=1 L=16 dt=0.02 T=80 seed=$s diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=$a 2>/dev/null | grep "RESULT blob_drift")
    sp=$(echo "$r" | sed -n 's/.*speed=\([0-9.e+-]*\).*/\1/p')
    co=$(echo "$r" | sed -n 's/.*cos_to_kdir=\([0-9.e+-]*\).*/\1/p')
    printf "%-8s %-10s %-12s %-10s\n" "$a" "$s" "${sp:-NA}" "${co:-NA}"
  done
done
echo "=== done ==="
