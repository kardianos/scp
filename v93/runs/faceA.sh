#!/bin/bash
# v93 FACE A: symmetrized hop schedule (hop_order) — RESUME §7-A / RING_DNLS §A.4.
# Two acceptance bars:
#  (1) frozen-scaffold vacuum C6 ring holds W=+2 to t=1000 (the sweep-order
#      ceiling: sequential slips by t~50-100 with all noise off; Strang/random
#      should hold).  + a live-geometry companion.
#  (2) e3b blob fd current is more dt-invariant under Strang than sequential
#      (item 1's Trotter artifact): fix amp_nat, vary dt.
cd "$(dirname "$0")/.."
mkdir -p runs/faceA
echo "### FACE A: symmetrized hop schedule (hop_order 0=seq 1=strang 2=rand) ###"

# --- (1) vacuum frozen-scaffold C6 ring retention, T=1000 ---
echo "=== (1) vacuum C6 ring, T=1000 (seed_mw=2, amp_nat=2 amp_logate=1) ==="
ring_run () {
  local name=$1 ho=$2; shift 2
  ./freecell L=16 dt=0.02 T=1000 k_rad=0 wf_on=0 kx=0 seed=20260802 \
    diag_every=2000 snap_every=200 \
    exp=ring ring_n=6 bath=0 seed_mw=2 amp_nat=2 amp_logate=1 \
    freeze_geo=1 sigma_tumble=0 hop_order=$ho "$@" \
    snap_file=runs/faceA/$name.fcs > runs/faceA/$name.log 2>&1
  echo "--- $name (hop_order=$ho $3) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s\n",$4}' runs/faceA/$name.log
  ./fcsdump -mode cells runs/faceA/$name.fcs 2>/dev/null \
    | python3 report/analyze_ring.py /dev/stdin 8 8 2 2>/dev/null \
    | sed 's/^/    /'
}
ring_run ringA_seq    0
ring_run ringA_strang  1
ring_run ringA_rand   2

# live-geometry companion (freeze_geo=0, tumble on) — the physical case
echo "=== (1b) vacuum C6 ring LIVE geometry, T=1000 ==="
ring_run_live () {
  local name=$1 ho=$2
  ./freecell L=16 dt=0.02 T=1000 k_rad=0 wf_on=0 kx=0 seed=20260802 \
    diag_every=2000 snap_every=200 \
    exp=ring ring_n=6 bath=0 seed_mw=2 amp_nat=2 amp_logate=1 \
    freeze_geo=0 hop_order=$ho \
    snap_file=runs/faceA/$name.fcs > runs/faceA/$name.log 2>&1
  echo "--- $name (hop_order=$ho live) ---"
  awk '/RESULT drift_rel/{printf "  drift=%s\n",$4}' runs/faceA/$name.log
  ./fcsdump -mode cells runs/faceA/$name.fcs 2>/dev/null \
    | python3 report/analyze_ring.py /dev/stdin 8 8 2 2>/dev/null \
    | sed 's/^/    /'
}
ring_run_live ringA_seq_live   0
ring_run_live ringA_strang_live 1
ring_run_live ringA_rand_live  2

# --- (2) e3b blob dt-invariance: fix amp_nat=2, vary dt; read fd ---
echo "=== (2) e3b blob dt-invariance (fd current), amp_nat=2 fixed ==="
e3b_run () {
  local name=$1; shift
  ./freecell exp=blob bath=1 L=16 T=40 amp=0.5 sigma=2.5 kx=1.1 wf_on=1 \
    amp_logate=1 p1_meter=1 seed=20260802 diag_every=1000 "$@" \
    > runs/faceA/$name.log 2>&1
  local fd=$(grep -oE "fd=[+0-9.eE-]+" runs/faceA/$name.log | tail -1)
  local dr=$(awk '/RESULT drift_rel/{print $4}' runs/faceA/$name.log | tail -1)
  echo "  $name : $fd  drift=$dr"
}
for ho in 0 1 2; do
  e3b_run e3b_ho${ho}_dt02 dt=0.02 amp_nat=2 hop_order=$ho
  e3b_run e3b_ho${ho}_dt01 dt=0.01 amp_nat=2 hop_order=$ho
done

echo "=== faceA done ==="
