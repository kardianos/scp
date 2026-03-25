#!/bin/bash
# Run all 8 Gen 1 candidates in batches of 2 (8 threads each)
# T=20 probe runs
set -e
BASE="/home/d/code/scp/v40/gen_001"
SIM="/home/d/code/scp/sfa/sim/scp_sim"
ANALYZE="/home/d/code/scp/v40/tools/analyze_sfa"

echo "=== Gen 1 Evolutionary Search (T=20 probe) ==="
echo "Started: $(date)"

run_pair() {
    local a=$1 b=$2
    echo "--- Batch: candidate_$a + candidate_$b --- $(date)"
    (cd "$BASE/candidate_$a" && OMP_NUM_THREADS=8 $SIM config.cfg && $ANALYZE output.sfa --json analysis.json && echo "candidate_$a DONE") &
    (cd "$BASE/candidate_$b" && OMP_NUM_THREADS=8 $SIM config.cfg && $ANALYZE output.sfa --json analysis.json && echo "candidate_$b DONE") &
    wait
    echo "--- Batch complete --- $(date)"
}

run_pair 001 002
run_pair 003 004
run_pair 005 006
run_pair 007 008

echo "=== All 8 candidates complete ==="
echo "Finished: $(date)"
