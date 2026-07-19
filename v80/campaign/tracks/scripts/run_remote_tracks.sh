#!/usr/bin/env bash
# Run on remote: pair-track force SFAs + elite same (seed trajectories).
set -u
WORK=${WORK:-/root/v80/work_tracks}
BIN=${BIN:-/root/v80/bin}
SEEDS=${SEEDS:-/root/v80/seeds}
mkdir -p "$WORK/tracks"
cd "$WORK"
LOG=$WORK/track_jobs.log
echo "TRACK_START $(date -Is)" | tee -a "$LOG"

have_track() {
  local out="$1"
  [ -f "$out" ] && [ "$(stat -c%s "$out" 2>/dev/null || echo 0)" -gt 100 ]
}

for D in 8 12 16 20 24; do
  f=mf_force_D${D}.sfa
  out=tracks/mf_force_D${D}_track.tsv
  if [ -f "$f" ] && ! have_track "$out"; then
    rm -f "$out"
    echo "TRACK force D=$D $(date -Is)" | tee -a "$LOG"
    if "$BIN/mf_pair_track" "$f" "$out" >>"$LOG" 2>&1; then
      echo "OK force D=$D $(wc -l < "$out") lines" | tee -a "$LOG"
    else
      echo "FAIL force D=$D" | tee -a "$LOG"
    fi
  elif have_track "$out"; then
    echo "SKIP force D=$D already tracked" | tee -a "$LOG"
  else
    echo "WAIT force D=$D no SFA yet" | tee -a "$LOG"
  fi
done

for D in 12 16 20; do
  f=$SEEDS/elite_same_D${D}.sfa
  out=tracks/elite_same_D${D}_track.tsv
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$sz" -gt 100000000 ] && ! have_track "$out"; then
    rm -f "$out"
    echo "TRACK elite D=$D sz=$sz $(date -Is)" | tee -a "$LOG"
    if "$BIN/sfa_qball_track" "$f" --tsv "$out" >>"$LOG" 2>&1; then
      echo "OK elite D=$D $(wc -l < "$out") lines" | tee -a "$LOG"
    else
      echo "FAIL elite D=$D" | tee -a "$LOG"
    fi
  elif have_track "$out"; then
    echo "SKIP elite D=$D already tracked" | tee -a "$LOG"
  else
    echo "WAIT elite D=$D sz=$sz" | tee -a "$LOG"
  fi
done

# Orbits if present
for tag in 0p03 0p05; do
  f=mf_orbit_R16_vt${tag}.sfa
  out=tracks/mf_orbit_R16_vt${tag}_track.tsv
  if [ -f "$f" ] && ! have_track "$out"; then
    rm -f "$out"
    echo "TRACK orbit vt=$tag $(date -Is)" | tee -a "$LOG"
    if "$BIN/mf_pair_track" "$f" "$out" >>"$LOG" 2>&1; then
      echo "OK orbit vt=$tag $(wc -l < "$out") lines" | tee -a "$LOG"
    else
      echo "FAIL orbit vt=$tag" | tee -a "$LOG"
    fi
  fi
done

echo "TRACK_DONE $(date -Is)" | tee -a "$LOG"
ls -lh tracks/ | tee -a "$LOG"
