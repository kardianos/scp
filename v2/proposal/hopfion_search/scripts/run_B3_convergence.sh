#!/bin/bash
# B=3 Convergence Study
#
# Investigates the B=3 anomaly: 18.7% energy error at N=160, L=6
# while B=2 and B=4 have <0.5% error at the same resolution.
#
# Hypothesis: The B=3 tetrahedral rational map creates sharper angular
# features than B=2 (axial) or B=4 (cubic), requiring finer resolution.
#
# Tests: B=3 at increasing N (init + energy only, no relaxation).
# Also re-tests B=2, B=4 at a few points for comparison.
#
# Memory estimates (init-only, no force/gradient arrays):
#   N=256:  ~1.1 GB       N=384:  ~3.6 GB
#   N=320:  ~2.1 GB       N=512:  ~8.6 GB
#
# Machine: 115 GB RAM available — can comfortably run N=512.
#
# Usage: bash run_B3_convergence.sh 2>&1 | tee B3_convergence_results.txt

set -e

echo "========================================="
echo " B=3 Convergence Study"
echo " $(date)"
echo " Threads: $(nproc)"
echo " RAM: $(free -g | awk '/Mem:/{print $7}') GB available"
echo "========================================="
echo

# Build
make -j$(nproc) verify3d 2>&1 | tail -1
echo

# ============================================
# Part 1: B=3 convergence at L=6
# Full sweep from N=160 (known bad) to N=512
# ============================================
echo "========================================="
echo " Part 1: B=3 Convergence at L=6"
echo "========================================="
echo

for N in 160 192 256 320 384 512; do
    h=$(echo "scale=6; 12.0/$N" | bc)
    echo ">>> B=3 N=$N L=6 h=$h <<<"
    ./verify3d -N $N -L 6 -B 3 -profile profile_B3.dat 2>&1
    echo
done

# ============================================
# Part 2: B=3 at L=8 (separate periodic BC effect from discretization)
# ============================================
echo "========================================="
echo " Part 2: B=3 at L=8"
echo "========================================="
echo

for N in 256 384 512; do
    h=$(echo "scale=6; 16.0/$N" | bc)
    echo ">>> B=3 N=$N L=8 h=$h <<<"
    ./verify3d -N $N -L 8 -B 3 -profile profile_B3.dat 2>&1
    echo
done

# ============================================
# Part 3: B=1,2,4 comparison at higher N
# Confirms that the anomaly is B=3 specific
# ============================================
echo "========================================="
echo " Part 3: B=1,2,4 comparison at N=256 L=6"
echo "========================================="
echo

for B in 1 2 4; do
    echo ">>> B=$B N=256 L=6 <<<"
    ./verify3d -N 256 -L 6 -B $B -profile profile_B${B}.dat 2>&1 | grep -E "(E_total|E2/E4|Error|Q error|Q )"
    echo
done

echo "========================================="
echo " DONE — $(date)"
echo "========================================="
