#!/bin/bash
# run_f2_tests.sh — Complete F2 experimental battery
#
# Test 1: E(D) — total energy vs separation (energy minimization check)
# Test 2: Asymmetry vs gradient strength
#
# Each sub-test uses 4 threads, multiple in parallel

set -e
mkdir -p data/f2_energy data/f2_gradient

echo "=========================================="
echo "F2 Complete Test Battery"
echo "=========================================="

# ---- Build ----
echo "Building..."
gcc -O3 -march=native -fopenmp -o energy_vs_D src/v33_energy_vs_D.c -lm
gcc -O3 -march=native -fopenmp -o v33_gtest src/v33_gradient_test.c -lm
gcc -O3 -march=native -o footprint src/footprint_asymmetry.c -lm
echo "Build complete."

# =============================================
# TEST 1: E(D) — Energy vs Separation
# =============================================
echo ""
echo "=== TEST 1: Energy vs Braid Separation ==="
echo "Running single-braid baseline + pair at D=8,10,12,15,18,20,25,30,40,50..."

RESULTS_E="data/f2_energy/results.tsv"
echo -e "config\tD\tE_avg\tsamples" > "$RESULTS_E"

# Single braid baseline
OMP_NUM_THREADS=4 ./energy_vs_D -single -N 128 -L 40 -T 30 >> "$RESULTS_E" &
PID_S=$!

# Two-braid pairs at various D (4 threads each, up to 5 concurrent)
PIDS=()
for D in 8 10 12 15 18 20 25 30 40 50; do
    OMP_NUM_THREADS=4 ./energy_vs_D -D $D -N 128 -L 40 -T 30 >> "$RESULTS_E" &
    PIDS+=($!)
    # Limit concurrency to 5
    if [ ${#PIDS[@]} -ge 5 ]; then
        wait ${PIDS[0]}
        PIDS=("${PIDS[@]:1}")
    fi
done

# Wait for all energy jobs
wait $PID_S
for pid in "${PIDS[@]}"; do wait $pid; done

echo "Energy sweep complete. Results in $RESULTS_E"
cat "$RESULTS_E"

# =============================================
# TEST 2: Gradient Strength Sweep
# =============================================
echo ""
echo "=== TEST 2: Footprint Asymmetry vs Gradient Strength ==="
echo "Running gradient tests at 4 different strengths..."

RESULTS_G="data/f2_gradient/results.tsv"
echo -e "A_high\tA_low\tratio\tdrift\tR_low_R_high" > "$RESULTS_G"

# Gradient tests: N=128 T=60 (enough for footprint + some drift)
# 4 gradient strengths
declare -a AH=(0.105 0.11 0.13 0.15)
declare -a AL=(0.095 0.09 0.07 0.05)
declare -a LABELS=("gentle" "mild" "moderate" "strong")

GPIDS=()
for idx in 0 1 2 3; do
    ah=${AH[$idx]}
    al=${AL[$idx]}
    lab=${LABELS[$idx]}
    outdir="data/f2_gradient/${lab}"
    mkdir -p "$outdir"
    echo "  Launching: A_high=$ah A_low=$al ($lab) → $outdir"
    OMP_NUM_THREADS=4 ./v33_gtest -N 128 -L 20 -T 60 -ah $ah -al $al -o "$outdir" > "${outdir}/log.txt" 2>&1 &
    GPIDS+=($!)
done

# Wait for all gradient tests
for pid in "${GPIDS[@]}"; do wait $pid; done
echo "Gradient tests complete."

# Run footprint analysis on each
echo ""
echo "=== Footprint Analysis ==="
for idx in 0 1 2 3; do
    ah=${AH[$idx]}
    al=${AL[$idx]}
    lab=${LABELS[$idx]}
    outdir="data/f2_gradient/${lab}"
    echo ""
    echo "--- $lab (A_high=$ah, A_low=$al, ratio=$(echo "scale=2; $ah*$ah/($al*$al)" | bc)) ---"
    # Analyze t=0 and latest snapshot
    for snap in "$outdir"/field_t*.bin; do
        echo "  Snapshot: $snap"
        ./footprint "$snap" $ah $al 2>/dev/null | grep -E "Half-width|Shift|Footprint extent|center"
        echo ""
    done
done

echo ""
echo "=========================================="
echo "F2 TEST BATTERY COMPLETE"
echo "=========================================="
echo "Results:"
echo "  Energy: $RESULTS_E"
echo "  Gradient: data/f2_gradient/{gentle,mild,moderate,strong}/"
