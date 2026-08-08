#!/bin/bash
# v93 L1-A RE-RUN: the adoption gate. Does Strang (hop_order=1) preserve the
# coherent seed-robust +x dense translation (fd) that sequential (0) gave
# (items 1/2: fd ~+287/+297/+277)? If yes -> adopt Strang as default.
cd "$(dirname "$0")/.."
mkdir -p runs/l1a_rerun
echo "### L1-A re-run: e3b blob amp_nat=2 amp_logate=1 p1_meter=1 (fd=current) ###"
printf "%-8s %-10s %16s %14s\n" hop_order seed fd drift_rel
for ho in 0 1; do
  for s in 111 20260802 314159; do
    ./freecell exp=blob bath=1 L=16 dt=0.02 T=80 seed=$s diag_every=200 \
      amp=0.5 sigma=2.5 kx=1.1 wf_on=1 amp_nat=2 amp_logate=1 p1_meter=1 \
      hop_order=$ho > runs/l1a_rerun/ho${ho}_s${s}.log 2>&1
    fd=$(grep -oE "fd=[+0-9.eE-]+" runs/l1a_rerun/ho${ho}_s${s}.log | tail -1 | cut -d= -f2)
    dr=$(awk '/RESULT drift_rel/{print $4}' runs/l1a_rerun/ho${ho}_s${s}.log | tail -1)
    printf "%-8s %-10s %16s %14s\n" $ho $s "$fd" "$dr"
  done
done
echo "=== l1a rerun done ==="
