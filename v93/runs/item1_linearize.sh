#!/bin/bash
# v93 ITEM 1: linearize-tau experiment (reviewer item 1).
# Does dropping the phase-dependent gate (tau = amp_nat*base, linear Schrodinger
# hop) recover coherent propagation post-empty-cell-fix? And does v_g(kx) follow
# the tight-binding sin(kx*d) dispersion, or is it artifact?
cd "$(dirname "$0")/.."
BLOB="exp=blob bath=1 L=16 dt=0.02 T=80 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1"
sc () { grep "RESULT blob_drift" | sed -E 's/.*speed=([0-9.e+-]+).*cos_to_kdir=([0-9.e+-]+).*/\1 cos=\2/'; }

echo "=== 1a) GATE-ON (amp_logate=0) vs GATE-OFF/LINEAR (amp_logate=1), amp_nat=2 ==="
printf "%-10s %-8s %-10s %-18s\n" gate seed amp_nat "speed/cos"
for s in 111 20260802 314159; do
  for gl in 0 1; do
    r=$(./freecell $BLOB seed=$s amp_nat=2 amp_logate=$gl 2>/dev/null | sc)
    printf "%-10s %-8s %-10s %-18s\n" "$( [ $gl = 0 ] && echo ON || echo OFF )" "$s" "2" "$r"
  done
done

echo
echo "=== 1b) kx SWEEP (LINEAR, amp_nat=2 amp_logate=1, seed 111) -- v_g ~ sin(kx*d)? ==="
echo "    (d~1.4; sin(kx*d) peaks near kx~1.1, zero near kx~2.2)"
printf "%-8s %-18s\n" kx "speed/cos"
for kx in 0.5 0.8 1.1 1.5 2.0; do
  r=$(./freecell $BLOB kx=$kx seed=111 amp_nat=2 amp_logate=1 2>/dev/null | sc)
  printf "%-8s %-18s\n" "$kx" "$r"
done

echo
echo "=== 1c) dt-INVARIANCE (LINEAR, amp_nat*dt product fixed; real rate invariant, artifacts move) ==="
printf "%-8s %-10s %-18s\n" dt amp_nat "speed/cos"
# amp_nat*dt fixed at 0.04 (amp_nat=2,dt=0.02) and (amp_nat=4,dt=0.01)
for cfg in "0.02 2" "0.01 4"; do set -- $cfg; dt=$1; an=$2;
  r=$(./freecell exp=blob bath=1 L=16 dt=$dt T=80 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 seed=111 amp_nat=$an amp_logate=1 2>/dev/null | sc)
  printf "%-8s %-10s %-18s\n" "$dt" "$an" "$r"
done
echo "=== item 1 done ==="
