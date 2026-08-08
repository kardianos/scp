#!/bin/bash
# v93 RING RETENTION (RESUME §7 route A): can a winding survive on a
# topologically-supported CLOSED DENSE RING (|psi|>0 everywhere on the
# cycle) where the localized blob failed (item 4: unitary scrambled a
# hand-seeded m=+2 blob vortex to R2d~0.2 in 4 t.u.)?
# Substrate: exp=ring + seed_mw (th2 = m*phi on the cycle). No field, no
# door (k_rad=0 wf_on=0), kx=0 — pure hold test, item-4 comparable.
cd "$(dirname "$0")/.."
mkdir -p runs/quench
S="L=16 dt=0.02 T=80 k_rad=0 wf_on=0 kx=0 seed=20260802 diag_every=500 snap_every=100"
run () {
  local name=$1 m=$2; shift 2
  ./freecell $S "$@" snap_file=runs/quench/$name.fcs > runs/quench/$name.log 2>&1
  echo "--- $name ($*) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s\n",$4}' runs/quench/$name.log
  ./fcsdump -mode cells runs/quench/$name.fcs 2>/dev/null \
    | python3 report/analyze_ring.py /dev/stdin 8 8 $m 2>/dev/null \
    | awk 'NR<=1 || NR%4==2' | sed 's/^/    /'
}
# --- vacuum nv=6 (exact C6 cycle graph; no leak sites) ---
run rv6_add_m2 2 exp=ring ring_n=6 bath=0 seed_mw=2 amp_nat=0
run rv6_uni_m2 2 exp=ring ring_n=6 bath=0 seed_mw=2 amp_nat=2 amp_logate=1
run rv6_gat_m2 2 exp=ring ring_n=6 bath=0 seed_mw=2 amp_nat=2
run rv6_uni_m1 1 exp=ring ring_n=6 bath=0 seed_mw=1 amp_nat=2 amp_logate=1
# --- bath-embedded nv=6 (leak sites available: the real test) ---
run rb6_add_m2 2 exp=ring ring_n=6 bath=1 seed_mw=2 amp_nat=0
run rb6_uni_m2 2 exp=ring ring_n=6 bath=1 seed_mw=2 amp_nat=2 amp_logate=1
# --- finer cycle nv=12, bath ---
run rb12_uni_m2 2 exp=ring ring_n=12 bath=1 seed_mw=2 amp_nat=2 amp_logate=1
# --- no-winding control (chain lock left in place) ---
run rb6_uni_m0 2 exp=ring ring_n=6 bath=1 seed_mw=0 amp_nat=2 amp_logate=1
echo "=== ringret done ==="
