#!/usr/bin/env bash
# Wait for sustain campaign done marker, then print SUSTAIN_REPORT head.
# Usage: bash theory/runs/watch_and_score.sh theory/runs/sustain_T3000_...
set -euo pipefail
DIR="${1:?dir}"
while [[ ! -f "${DIR}.done" ]]; do
  n=$(pgrep -c fabrun 2>/dev/null || echo 0)
  # live progress from any log that has t=
  # rg -h is --help; use --no-filename. Logs: "t=   20.00 ..."
  last=$(rg --no-filename '^t=' "$DIR"/*.log 2>/dev/null | tail -1 || true)
  # fallback: partial FCS last CELL line via fcsdump if no log yet
  if [[ -z "$last" ]]; then
    last="logs=0 fcs=$(ls "$DIR"/*.fcs 2>/dev/null | wc -l)"
  fi
  echo "$(date -Is) fabrun=$n $last"
  sleep 120
done
echo "DONE $(date -Is)"
[[ -f "$DIR/SUSTAIN_REPORT.md" ]] && cat "$DIR/SUSTAIN_REPORT.md"
[[ -f "$DIR/results.csv" ]] && column -t -s, "$DIR/results.csv" | head -40
