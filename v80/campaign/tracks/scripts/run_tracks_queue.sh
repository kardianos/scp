#!/usr/bin/env bash
# Force + low-vt orbit queue. KEEP all SFAs (no prune).
set -u
BIN=${BIN:-/root/v80/bin/scp_sim_cuda}
WORK=${WORK:-/root/v80/work_tracks}
STATUS=$WORK/STATUS.txt
LOG=$WORK/queue.log
mkdir -p "$WORK/done"
cd "$WORK"
echo "QUEUE_START $(date -Is)" | tee -a "$LOG" "$STATUS"

run_one() {
  local step="$1" cfg="$2"
  local name; name=$(basename "$cfg" .cfg)
  echo "=== RUN $step $name $(date -Is) ===" | tee -a "$LOG" "$STATUS"
  if ! "$BIN" "$cfg" > "${name}.runlog" 2>&1; then
    echo "FAIL $step $name" | tee -a "$LOG" "$STATUS"
    echo FAIL > "done/${name}.fail"
    return 1
  fi
  echo "OK $step $name $(date -Is)" | tee -a "$LOG" "$STATUS"
  echo OK > "done/${name}.ok"
  ls -lh "${name}.sfa" 2>/dev/null | tee -a "$LOG" || true
  df -h / | tee -a "$LOG"
  return 0
}

# Phase 2 force opposite
for D in 8 12 16 20 24; do
  run_one F "cfg/mf_force_D${D}.cfg" || true
done
# Phase 2 same-sign control
for D in 12 16 20; do
  run_one R "cfg/elite_same_D${D}.cfg" || true
done
# Phase 3 orbit low vt
for tag in 0p03 0p05; do
  run_one O "cfg/mf_orbit_R16_vt${tag}.cfg" || true
done

echo "QUEUE_DONE $(date -Is)" | tee -a "$LOG" "$STATUS"
ls -lh *.sfa 2>/dev/null | tee -a "$LOG"
df -h / | tee -a "$LOG"
