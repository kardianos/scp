#!/usr/bin/env bash
# Gate O long-hold: hard multi-rev then soft T=4000 capture.
# KEEP SFAs; track after each job; prune only after track TSV OK.
set -u
BIN=${BIN:-/root/v80/bin/scp_sim_cuda}
TRACK=${TRACK:-/root/v80/bin/mf_pair_track}
WORK=${WORK:-/root/v80/work_orbit2}
STATUS=$WORK/STATUS.txt
LOG=$WORK/queue.log
mkdir -p "$WORK/done" "$WORK/tracks"
cd "$WORK"
echo "QUEUE_START $(date -Is)" | tee -a "$LOG" "$STATUS"

run_one() {
  local name="$1" cfg="$2"
  echo "=== RUN $name $(date -Is) ===" | tee -a "$LOG" "$STATUS"
  df -h / | tee -a "$LOG"
  if ! "$BIN" "$cfg" > "${name}.runlog" 2>&1; then
    echo "FAIL $name" | tee -a "$LOG" "$STATUS"
    echo FAIL > "done/${name}.fail"
    return 1
  fi
  echo "OK $name $(date -Is)" | tee -a "$LOG" "$STATUS"
  echo OK > "done/${name}.ok"
  ls -lh "${name}.sfa" 2>/dev/null | tee -a "$LOG" || true
  # track immediately
  if [ -f "${name}.sfa" ]; then
    echo "TRACK $name $(date -Is)" | tee -a "$LOG" "$STATUS"
    if "$TRACK" "${name}.sfa" "tracks/${name}_track.tsv" >>"$LOG" 2>&1; then
      echo "TRACK_OK $name $(wc -l < tracks/${name}_track.tsv) lines" | tee -a "$LOG" "$STATUS"
    else
      echo "TRACK_FAIL $name" | tee -a "$LOG" "$STATUS"
    fi
  fi
  df -h / | tee -a "$LOG"
  return 0
}

# 1) hard multi-rev attempt (~1.5 h)
run_one mf_orbit_R16_vt0p50_T1200 cfg/mf_orbit_R16_vt0p50_T1200.cfg || true

# free hard SFA after track if disk tight (<25G free) — keep track+diag
avail=$(df -BG / | awk 'NR==2{gsub("G",""); print $4}')
if [ "${avail:-0}" -lt 25 ] && [ -f tracks/mf_orbit_R16_vt0p50_T1200_track.tsv ]; then
  sz=$(stat -c%s tracks/mf_orbit_R16_vt0p50_T1200_track.tsv 2>/dev/null || echo 0)
  if [ "$sz" -gt 1000 ]; then
    rm -f mf_orbit_R16_vt0p50_T1200.sfa
    echo "PRUNED hard SFA after track (disk avail was ${avail}G)" | tee -a "$LOG" "$STATUS"
  fi
fi

# 2) soft long capture T=4000 (~5 h)
run_one mf_orbit_R16_vt0p03_T4000 cfg/mf_orbit_R16_vt0p03_T4000.cfg || true

echo "QUEUE_DONE $(date -Is)" | tee -a "$LOG" "$STATUS"
ls -lh *.sfa tracks/* 2>/dev/null | tee -a "$LOG"
df -h / | tee -a "$LOG"
