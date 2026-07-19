#!/usr/bin/env bash
# Local 15-min poll: status + rsync diags/runlogs; score hints.
set -u
HOST=ssh7.vast.ai
PORT=19142
RWORK=/root/v80/work
CAMP=/home/d/code/scp/v80/campaign
STAMP=$(date -Is)
LOG=$CAMP/logs/poll.log
mkdir -p "$CAMP/logs" "$CAMP/results" \
  "$CAMP/S0_smoke/results" "$CAMP/S1_force/results" \
  "$CAMP/S2_ql_pair/results" "$CAMP/S3_orbit/results" "$CAMP/S4_hydrogenoid/results"

echo "==== POLL $STAMP ====" | tee -a "$LOG"

ssh -o ConnectTimeout=20 -o StrictHostKeyChecking=no -p $PORT root@$HOST "
  echo TIME:\$(date -Is)
  if [[ -f $RWORK/queue.pid ]]; then
    pid=\$(cat $RWORK/queue.pid)
    if ps -p \$pid >/dev/null 2>&1; then echo QUEUE_PID=\$pid RUNNING; else echo QUEUE_PID=\$pid DEAD; fi
  else echo NO_PID; fi
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader 2>/dev/null || true
  echo --- STATUS ---
  tail -20 $RWORK/STATUS.txt 2>/dev/null || true
  echo --- DONE ---
  ls $RWORK/done 2>/dev/null || true
  echo --- DISK ---
  df -h / | tail -1
  echo --- TAIL RUNLOG ---
  ls -t $RWORK/*.runlog 2>/dev/null | head -1 | xargs -r tail -5
" 2>&1 | tee -a "$LOG" | tee "$CAMP/logs/last_poll.txt"

# rsync diags and runlogs
rsync -az -e "ssh -p $PORT -o StrictHostKeyChecking=no" \
  "root@$HOST:$RWORK/*_diag.tsv" "root@$HOST:$RWORK/*.runlog" \
  "root@$HOST:$RWORK/STATUS.txt" "root@$HOST:$RWORK/queue.log" \
  "root@$HOST:$RWORK/FLAG_*" \
  "$CAMP/results/" 2>>"$LOG" || true

# distribute into step folders by name prefix
shopt -s nullglob
for f in "$CAMP/results"/*; do
  b=$(basename "$f")
  case "$b" in
    smoke*) dest=$CAMP/S0_smoke/results ;;
    mf_force_D*|elite_*) dest=$CAMP/S1_force/results ;;
    mf_pair*|mf_headon*) dest=$CAMP/S2_ql_pair/results ;;
    mf_orbit*) dest=$CAMP/S3_orbit/results ;;
    mf_h1*|mf_force_N192*) dest=$CAMP/S4_hydrogenoid/results ;;
    *) dest=$CAMP/results ;;
  esac
  cp -u "$f" "$dest/" 2>/dev/null || true
done

# if queue done, note it
if grep -q QUEUE_DONE "$CAMP/results/STATUS.txt" 2>/dev/null; then
  echo "CAMPAIGN_COMPLETE $STAMP" | tee -a "$LOG"
  touch "$CAMP/logs/COMPLETE"
fi
df -h /space /home 2>/dev/null | tee -a "$LOG"
echo "poll ok $STAMP" >> "$LOG"
