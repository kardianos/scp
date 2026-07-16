#!/bin/bash
# Launch F19 Z6 campaign in background. Poll with sim_exec using ONLY the
# instant status snippet (v75/analysis/remote_status_snippet.sh) — no sleep.
set -euo pipefail
cd /root/pn/z6
BIN=${BIN:-/root/scp_sim_pn_cuda}
LOG=campaign.log
echo "F19 Z6 campaign start $(date -u)" | tee "$LOG"
for cfg in z6_n0_nuc z6_n6_nuc z6_a_n0 z6_a_n6; do
  echo "===== $cfg $(date -u) =====" | tee -a "$LOG"
  "$BIN" "${cfg}.cfg" 2>&1 | tee "${cfg}.runlog" | tail -30
  echo "===== done $cfg rc=$? $(date -u) =====" | tee -a "$LOG"
done
echo "ALL DONE $(date -u)" | tee -a "$LOG"
