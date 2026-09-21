#!/bin/bash
# Tiny then medium inner-chart experiments. Fork only; does not touch v97/fab.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
BIN=chartexp/chartrun
OUT=chartexp/runs
mkdir -p "$OUT"
go test ./chartexp/fab/
go build -o "$BIN" ./chartexp/cmd/chartrun

COMMON=(exp=ring ring_n=6 seed_mw=2 bath=0 L=16 dt=0.02
        amp_nat=2 amp_logate=1 hop_order=1
        seed=20260802 wf_on=0 k_rad=0.05)

run() {
  local name="$1"; shift
  echo "=== $name ==="
  "$BIN" "${COMMON[@]}" "$@" > "$OUT/${name}.log" 2>&1
  echo "$name  $(grep -E 'RESULT (drift|conv |chart)' "$OUT/${name}.log" | tr '\n' ' ')"
}

# --- tiny: T=40 vacuum C6 ---
TINY=(T=40 diag_every=2000)
run t_off    "${TINY[@]}" chart=0
run t_proj   "${TINY[@]}" chart=1 chart_scope=star chart_acct=energy
run t_grain  "${TINY[@]}" chart=2 chart_scope=star chart_acct=energy
run t_mix    "${TINY[@]}" chart=3 chart_scope=star chart_acct=energy
run t_mass   "${TINY[@]}" chart=1 chart_scope=star chart_acct=mass
run t_area   "${TINY[@]}" chart=1 chart_scope=star chart_acct=area
run t_motif  "${TINY[@]}" chart=1 chart_scope=motif chart_acct=energy
run t_mix_m  "${TINY[@]}" chart=3 chart_scope=star chart_acct=mass

# --- medium: T=200 vacuum C6 (retention window) ---
MED=(T=200 diag_every=2000)
run m_off    "${MED[@]}" chart=0
run m_proj   "${MED[@]}" chart=1 chart_scope=star chart_acct=energy
run m_grain  "${MED[@]}" chart=2 chart_scope=star chart_acct=energy
run m_mix    "${MED[@]}" chart=3 chart_scope=star chart_acct=energy
run m_mix_m  "${MED[@]}" chart=3 chart_scope=star chart_acct=mass
run m_motif  "${MED[@]}" chart=3 chart_scope=motif chart_acct=energy

echo ALL_CHARTEXP_DONE
