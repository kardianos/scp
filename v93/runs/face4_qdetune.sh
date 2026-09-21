#!/bin/bash
cd "$(dirname "$0")/.."
echo "=== face 4 probe: does load-detuning (q_detune) shear the tilt? ==="
echo "=== amp_nat=2.0, q_detune sweep, 5 seeds (wf_on k_rad=0 baseline) ==="
printf "%-10s %-10s %-12s %-10s\n" q_detune seed speed cos
for q in 1.2 0.6 0.3 0.0; do
  for s in 111 20260802 314159 137035 161803; do
    r=$(./freecell exp=blob bath=1 L=16 dt=0.02 T=80 seed=$s diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=2.0 q_detune=$q 2>/dev/null | grep "RESULT blob_drift")
    sp=$(echo "$r" | sed -n 's/.*speed=\([0-9.e+-]*\).*/\1/p')
    co=$(echo "$r" | sed -n 's/.*cos_to_kdir=\([0-9.e+-]*\).*/\1/p')
    printf "%-10s %-10s %-12s %-10s\n" "$q" "$s" "${sp:-NA}" "${co:-NA}"
  done
done
echo "=== done ==="
