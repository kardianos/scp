#!/bin/bash
# Run all Gen 1 candidates sequentially, T=10 probe, 16 threads
# Usage: ./run_gen1.sh

set -e
cd /home/d/code/scp/v40
SIM=/home/d/code/scp/sfa/sim/scp_sim
ANALYZE=tools/analyze_sfa
export OMP_NUM_THREADS=16

echo "=== Gen 1 Batch Run (T=10) ==="
echo ""

for i in $(seq 1 8); do
  dir="gen_001/candidate_00${i}"
  desc=$(head -1 "$dir/config.cfg" | sed 's/^# Gen 1, Candidate [0-9]*: //')
  echo "--- C${i}: $desc ---"

  # Clean previous output
  rm -f "$dir/output.sfa" "$dir/analysis.json"

  # Run simulation
  if $SIM "$dir/config.cfg" -T 10 > "$dir/sim_log.txt" 2>&1; then
    echo "  SIM: OK ($(tail -1 "$dir/sim_log.txt" | grep -o '[0-9.]* min'))"
  else
    echo "  SIM: FAILED (exit $?)"
    continue
  fi

  # Fix SFA index if needed
  tools/fixup_all "$dir/output.sfa" > /dev/null 2>&1 || true

  # Analyze
  if $ANALYZE "$dir/output.sfa" --json "$dir/analysis.json" > "$dir/analysis_log.txt" 2>&1; then
    score=$(grep S_final "$dir/analysis.json" | head -1 | grep -o '[0-9.]*')
    epot=$(grep E_pot_final "$dir/analysis.json" | head -1 | grep -o '[-0-9.]*')
    echo "  SCORE: S=$score E_pot=$epot"
  else
    echo "  ANALYZE: FAILED"
  fi
  echo ""
done

echo "=== Summary ==="
printf "%-5s %-40s %10s %10s\n" "ID" "Description" "S_final" "E_pot"
for i in $(seq 1 8); do
  dir="gen_001/candidate_00${i}"
  desc=$(head -1 "$dir/config.cfg" | sed 's/^# Gen 1, Candidate [0-9]*: //')
  if [ -f "$dir/analysis.json" ]; then
    score=$(grep S_final "$dir/analysis.json" | head -1 | grep -o '[0-9.]*')
    epot=$(grep E_pot_final "$dir/analysis.json" | head -1 | grep -o '[-0-9.]*')
    printf "C%-4d %-40s %10s %10s\n" $i "$desc" "$score" "$epot"
  else
    printf "C%-4d %-40s %10s %10s\n" $i "$desc" "FAILED" "FAILED"
  fi
done
