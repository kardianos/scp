#!/bin/bash
# v93 ITEM 4: matter-winding retention (hand-seeded m=+2 vortex, no field/door).
# Isolates hold-vs-imprint: does unitary transport preserve a matter winding
# better than the additive path? (removes the field-decoherence confound)
cd "$(dirname "$0")/.."
mkdir -p runs/quench
# dense matter vortex, no radiance (k_rad=0), no workfn/door (wf_on=0 -> no condensation)
S="exp=blob bath=1 L=16 dt=0.02 T=80 sigma=2.5 amp=0.5 kx=0 k_rad=0 wf_on=0 diag_every=500 snap_every=100 seed_mw=2"
run () {
  local name=$1; shift
  ./freecell $S "$@" snap_file=runs/quench/$name.fcs > runs/quench/$name.log 2>&1
  echo "--- $name ($*) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s\n",$4}' runs/quench/$name.log
  echo "  R2d(t) [matter m=+2 winding coherence]:"
  ./fcsdump -mode cells runs/quench/$name.fcs 2>/dev/null | python3 report/analyze_winding.py /dev/stdin 8 8 0 2>/dev/null \
    | grep -oE "^ +[0-9]+\.?[0-9]*.*R2d=[0-9.-]+" | sed -E 's/\|.*R2d=//' | awk '{printf "    t=%5.1f R2d=%s\n",$1,$NF}' | awk 'NR%2==1'
}
run mw_add  amp_nat=0                 # additive transport (control)
run mw_uni  amp_nat=2 amp_logate=1    # unitary linear transport
echo "=== item 4 done ==="
