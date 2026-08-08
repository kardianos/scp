#!/bin/bash
# v93 corrected re-measurement (post-review fixes): live-matter conservation,
# C==Go at strength, and QUENCH-3 WITH the empty-cell-reset + symmetric-reverse-door fixes.
cd "$(dirname "$0")/.."
mkdir -p runs/quench
LAW="k_rad=0.05 p_rad=4 wf_on=1 wf_floor=0.01 wf_far=99"
BLOB="exp=blob bath=1 L=16 dt=0.02 T=80 diag_every=200 amp=0.5 sigma=2.5 kx=1.1 wf_on=1"

echo "========================================================="
echo "1) L1-B CONSERVATION ON LIVE MATTER (e3b blob, has Em), not the matterless bath"
echo "========================================================="
printf "%-9s %-9s %-12s\n" amp_nat amp_door drift_rel
for a in 0 2; do
  for d in 0 1; do
    dr=$(./freecell $BLOB amp_nat=$a amp_door=$d 2>/dev/null | awk '/RESULT drift_rel/{print $4}')
    printf "%-9s %-9s %-12s\n" "$a" "$d" "${dr:-NA}"
  done
done
echo "(matterless-bath baseline for reference:)"
for a in 0 2; do
  dr=$(./freecell exp=bath T=40 $LAW amp_nat=$a 2>/dev/null | awk '/RESULT drift_rel/{print $4}')
  echo "  bath amp_nat=$a  drift=$dr"
done

echo
echo "========================================================="
echo "2) C==Go AT STRENGTH (e3b amp_nat=2, 3 seeds) -- tolerance, not byte"
echo "========================================================="
printf "%-10s %-22s %-22s\n" seed "C speed/cos" "Go speed/cos"
for s in 111 20260802 314159; do
  c=$(./freecell $BLOB seed=$s amp_nat=2 2>/dev/null | grep "RESULT blob_drift" | sed -E 's/.*speed=([0-9.e+-]+).*cos_to_kdir=([0-9.e+-]+).*/\1\/\2/')
  g=$(./fabrun  $BLOB seed=$s amp_nat=2 2>/dev/null | grep "RESULT blob_drift" | sed -E 's/.*speed=([0-9.e+-]+).*cos_to_kdir=([0-9.e+-]+).*/\1\/\2/')
  printf "%-10s %-22s %-22s\n" "$s" "$c" "$g"
done

echo
echo "========================================================="
echo "3) QUENCH-3 WITH THE FIXES (empty-cell reset + symmetric reverse door)"
echo "    Does retention (R2d) improve vs the pre-fix ~0.02 floor?"
echo "========================================================="
Q="exp=slit slit_mask=3 L=64 T=300 dt=0.02 sigma=8 slit_sy=8 amp=2 slit_srcx=32 kx=0 spin_m=2 seed=20260802 diag_every=2000 snap_every=1000 bath=1"
qrun () {
  local name=$1; shift
  ./freecell $Q "$@" snap_file=runs/quench/$name.fcs > runs/quench/$name.log 2>&1
  echo "--- $name ($*) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s  ",$4} /RESULT births/{printf "births=%s\n",$4}' runs/quench/$name.log
  echo "  R2d(t):"
  ./fcsdump -mode cells runs/quench/$name.fcs 2>/dev/null | python3 report/analyze_winding.py /dev/stdin 32 32 0 2>/dev/null \
    | grep -oE "^ +[0-9]+\.?[0-9]*.*R2d=[0-9.-]+" | sed -E 's/\|.*R2d=//' | awk '{printf "    t=%6.1f R2d=%s\n",$1,$NF}' | awk 'NR%2==1'
}
qrun fix_g0  amp_door=0 amp_nat=0                  # control (no door)
qrun fix_g1  amp_door=1 amp_nat=0                  # arg-door (coherent condensation + now symmetric evap)
qrun fix_g2  amp_door=1 amp_nat=2                  # arg-door + unitary
qrun fix_g3  amp_door=1 amp_nat=2 q_detune=0       # + uniform clock
echo "=== corrected re-measure done ==="
