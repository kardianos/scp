#!/bin/bash
# Batch 2: clean grain, mode door, composition H0. Fork only.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
BIN=chartexp/chartrun
OUT=chartexp/runs
mkdir -p "$OUT"
go test ./chartexp/fab/
go build -o "$BIN" ./chartexp/cmd/chartrun

BASE=(exp=ring ring_n=6 seed_mw=2 bath=0 L=16 dt=0.02
      amp_nat=2 amp_logate=1 hop_order=1
      seed=20260802 wf_on=0 k_rad=0.05)

run() {
  local name="$1"; shift
  echo "=== $name ==="
  "$BIN" "${BASE[@]}" "$@" > "$OUT/${name}.log" 2>&1
  echo "$name  $(grep -E 'RESULT (drift|conv |chart)' "$OUT/${name}.log" | tr '\n' ' ')"
}

# (1) clean grain — vacuum identity + lit vs off
run b2_v_off    T=80 diag_every=2000 chart=0
run b2_v_grain  T=80 diag_every=2000 chart=2 chart_scope=star chart_acct=energy
run b2_v_gmot   T=80 diag_every=2000 chart=2 chart_scope=motif chart_acct=energy
run b2_l_off    T=80 diag_every=2000 noise_amp=0.4 chart=0
run b2_l_grain  T=80 diag_every=2000 noise_amp=0.4 chart=2 chart_scope=star chart_acct=energy
run b2_l_gmot   T=80 diag_every=2000 noise_amp=0.4 chart=2 chart_scope=motif chart_acct=energy

# (2) mode door
run b2_l_mode   T=80 diag_every=2000 noise_amp=0.4 chart=4 chart_scope=motif chart_acct=energy
run b2_l_mstar  T=80 diag_every=2000 noise_amp=0.4 chart=4 chart_scope=star chart_acct=energy
run b2_v_mode   T=80 diag_every=2000 chart=4 chart_scope=motif chart_acct=energy

# (3) composition H0 time series
run b2_h0vac    T=200 diag_every=400 chart=1 chart_scope=motif chart_acct=energy
run b2_h0lit    T=80  diag_every=200 noise_amp=0.4 chart=1 chart_scope=motif chart_acct=energy

echo ALL_BATCH2_DONE
