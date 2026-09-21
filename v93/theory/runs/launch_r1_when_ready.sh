#!/usr/bin/env bash
set -euo pipefail
PRIMARY=/home/d/code/scp/v93/theory/runs/sustain_T3000_20260808_104548
V93=/home/d/code/scp/v93
while [[ ! -f "${PRIMARY}.done" ]]; do
  echo "$(date -Is) waiting primary.done fabrun=$(pgrep -c fabrun 2>/dev/null || echo 0)"
  sleep 120
done
echo "$(date -Is) PRIMARY DONE — launching R1"
if [[ -f /home/d/code/scp/v93/theory/runs/R1_LAUNCHED ]]; then
  echo "R1 already launched"
  exit 0
fi
touch /home/d/code/scp/v93/theory/runs/R1_LAUNCHED
cd "$V93"
OUT=theory/runs/R1_T10000_$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT"
ln -sfn "$(pwd)/$OUT" theory/runs/ACTIVE_R1
setsid env GOMAXPROCS=4 ./theory/bin/probe -suite r1 -kernel fabrun -T 10000 -j 4 \
  -seed 20260802 -snap-every 10000 -out "$OUT" \
  > "${OUT}.launch.log" 2>&1
echo EXIT=$? >> "${OUT}.launch.log"
touch "${OUT}.done"
echo "$(date -Is) R1 finished" >> /home/d/code/scp/v93/theory/HOURLOG.md
