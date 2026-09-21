#!/bin/bash
# v93 L1 physics sweep: L1-B conservation (drift floor) + L1-A e3b (blob translation).
# Law: V2g compiled defaults + wf_on=1 (the documented e3b "wf-on, no radiance" baseline).
cd "$(dirname "$0")/.."
mkdir -p runs/l1
echo "=== L1-B: conservation drift floor vs amp_nat (bath T=40, V3a law) ==="
for a in 0 0.5 1.0 2.0; do
  d=$(./freecell exp=bath T=40 k_rad=0.05 p_rad=4 wf_on=1 wf_floor=0.01 wf_far=99 amp_nat=$a \
        2>/dev/null | grep -E "RESULT drift_rel|drift" | head -1)
  echo "amp_nat=$a  $d"
done
echo
echo "=== L1-A: e3b blob translation (exp=blob amp=0.5 sigma=2.5 kx=1.1 L=16 T=80 wf_on=1) ==="
printf "%-8s %-10s %-12s %-8s\n" amp_nat seed speed cos_to_kdir
for a in 0 0.25 0.5 1.0 2.0; do
  for s in 111 20260802 314159; do
    r=$(./freecell exp=blob bath=1 L=16 dt=0.02 T=80 seed=$s diag_every=200 \
          amp=0.5 sigma=2.5 kx=1.1 wf_on=1 wf_floor=0.01 wf_far=99 amp_nat=$a \
          2>/dev/null | grep "RESULT blob_drift")
    sp=$(echo "$r" | sed -n 's/.*speed=\([0-9.e+-]*\).*/\1/p')
    co=$(echo "$r" | sed -n 's/.*cos_to_kdir=\([0-9.e+-]*\).*/\1/p')
    printf "%-8s %-10s %-12s %-8s\n" "$a" "$s" "${sp:-NA}" "${co:-NA}"
  done
done
echo "=== sweep done ==="
