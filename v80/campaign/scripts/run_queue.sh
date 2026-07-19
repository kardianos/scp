#!/usr/bin/env bash
# Serial GPU campaign queue for v80 overnight.
# Run from /root/v80/work with binary at /root/v80/scp_sim_cuda
set -u
BIN=${BIN:-/root/v80/scp_sim_cuda}
WORK=${WORK:-/root/v80/work}
STATUS=$WORK/STATUS.txt
LOG=$WORK/queue.log
FLAG_FORCE_OK=$WORK/FLAG_FORCE_OK
FLAG_PAIR_OK=$WORK/FLAG_PAIR_OK

mkdir -p "$WORK/out" "$WORK/done"
cd "$WORK"
echo "QUEUE_START $(date -Is)" | tee -a "$LOG" "$STATUS"

run_one() {
  local step="$1" cfg="$2"
  local name
  name=$(basename "$cfg" .cfg)
  echo "=== RUN $step $name $(date -Is) ===" | tee -a "$LOG" "$STATUS"
  # symlink seeds if cfg references bare names
  export OMP_NUM_THREADS=1
  if ! "$BIN" "$cfg" > "${name}.runlog" 2>&1; then
    echo "FAIL $step $name rc=$?" | tee -a "$LOG" "$STATUS"
    echo "FAIL" > "$WORK/done/${name}.fail"
    return 1
  fi
  echo "OK $step $name $(date -Is)" | tee -a "$LOG" "$STATUS"
  echo "OK" > "$WORK/done/${name}.ok"
  # keep diag+runlog; drop large SFA after keep last frame optional — delete SFA to save disk
  if [[ -f "${name}.sfa" ]]; then
    local sz
    sz=$(stat -c%s "${name}.sfa" 2>/dev/null || echo 0)
    # keep SFA only if < 400MB (N=128 short); else delete after success (diags enough)
    if [[ "$sz" -gt 400000000 ]]; then
      rm -f "${name}.sfa"
      echo "pruned large SFA $name ($sz)" | tee -a "$LOG"
    fi
  fi
  df -h / | tee -a "$LOG"
  return 0
}

# --- S0 ---
run_one S0 cfg/smoke_mf.cfg || { echo "S0_FAIL abort" | tee -a "$STATUS"; exit 1; }

# --- S1 force matrix ---
FORCE_OK=0
for D in 12 16 20; do
  run_one S1 cfg/mf_force_D${D}.cfg || true
  run_one S1 cfg/elite_opp_D${D}.cfg || true
  run_one S1 cfg/elite_same_D${D}.cfg || true
done
# crude auto gate: if any mf_force runlog contains COMPLETE and gauss e-14
if grep -l "COMPLETE" mf_force_D*_runlog mf_force_D*.runlog 2>/dev/null | head -1 | grep -q .; then
  FORCE_OK=1
fi
if ls mf_force_D16.runlog >/dev/null 2>&1 && grep -q "COMPLETE" mf_force_D16.runlog; then
  FORCE_OK=1
  touch "$FLAG_FORCE_OK"
fi
echo "FORCE_OK=$FORCE_OK" | tee -a "$STATUS"

# --- S2 pair ---
PAIR_OK=0
if [[ "$FORCE_OK" -eq 1 ]]; then
  run_one S2 cfg/mf_pair_rest_D20.cfg || true
  run_one S2 cfg/mf_headon_D20.cfg || true
  if grep -q "COMPLETE" mf_pair_rest_D20.runlog 2>/dev/null; then
    PAIR_OK=1
    touch "$FLAG_PAIR_OK"
  fi
else
  echo "SKIP S2 (force not OK)" | tee -a "$STATUS"
fi
echo "PAIR_OK=$PAIR_OK" | tee -a "$STATUS"

# --- S3 orbit ---
if [[ "$FORCE_OK" -eq 1 ]]; then
  for tag in 0p03 0p05 0p08; do
    run_one S3 cfg/mf_orbit_R16_vt${tag}.cfg || true
  done
else
  echo "SKIP S3 (force not OK)" | tee -a "$STATUS"
fi

# --- S4 ---
if [[ -f "$FLAG_PAIR_OK" ]] || [[ "$FORCE_OK" -eq 1 ]]; then
  # prefer long orbit if pair ok else short N192 force
  if [[ -f "$FLAG_PAIR_OK" ]]; then
    run_one S4 cfg/mf_h1_orbit_N192.cfg || true
  else
    run_one S4 cfg/mf_force_N192_D20.cfg || true
  fi
else
  echo "SKIP S4" | tee -a "$STATUS"
fi

echo "QUEUE_DONE $(date -Is)" | tee -a "$LOG" "$STATUS"
df -h / | tee -a "$LOG"
