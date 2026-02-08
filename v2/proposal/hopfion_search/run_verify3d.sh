#!/bin/bash
# 3D Verification Batch Script
#
# Runs convergence study for B=1 Skyrmion on 3D lattice.
# Part 1: Initialization accuracy (no relaxation, fast)
# Part 2: Gradient flow relaxation (tests topology preservation)
#
# Usage: bash run_verify3d.sh 2>&1 | tee verify3d_results.txt
#
# Estimated runtimes (16 threads):
#   N=64:  seconds
#   N=96:  ~10 sec per gradient step
#   N=128: ~30 sec per gradient step
#   N=160: ~90 sec per gradient step
#   N=192: ~180 sec per gradient step
#   N=256: ~600 sec per gradient step (probably overnight for 200 steps)

set -e

echo "========================================="
echo " 3D Skyrmion Verification — Batch Run"
echo " $(date)"
echo " Threads: $(nproc)"
echo "========================================="
echo

# Build
make -j$(nproc) verify3d 2>&1 | tail -1
echo

# ============================================
# Part 1: Initialization Accuracy (no relaxation)
# ============================================
echo "========================================="
echo " PART 1: Initialization Accuracy"
echo "========================================="
echo

echo "--- Varying N at L=6 ---"
for N in 64 96 128 160 192; do
    echo ">>> N=$N L=6 h=$(echo "scale=6; 12.0/$N" | bc) <<<"
    ./verify3d -N $N -L 6 -profile profile_B1.dat 2>&1 | grep -E "(E_total|E2/E4|Error|Q error)"
    echo
done

echo "--- Varying N at L=8 ---"
for N in 128 192 256; do
    echo ">>> N=$N L=8 h=$(echo "scale=6; 16.0/$N" | bc) <<<"
    ./verify3d -N $N -L 8 -profile profile_B1.dat 2>&1 | grep -E "(E_total|E2/E4|Error|Q error)"
    echo
done

# ============================================
# Part 2: Gradient Flow Relaxation
# Tests whether topology is preserved at each resolution.
# The soliton core radius is ~sqrt(c4) = 0.354 for e=4.
# Need h << 0.354 to prevent topological unwinding.
# ============================================
echo "========================================="
echo " PART 2: Gradient Flow Relaxation"
echo " (Tests topology preservation)"
echo "========================================="
echo

# N=128, L=6 (h=0.09375, ~3.8 points across core) — moderate
echo ">>> N=128 L=6 h=0.09375 (100 relaxation steps) <<<"
./verify3d -N 128 -L 6 -relax 1 -steps 100 -dt 1e-3 -profile profile_B1.dat 2>&1
echo

# N=160, L=6 (h=0.075, ~4.7 points across core) — should be better
echo ">>> N=160 L=6 h=0.075 (200 relaxation steps) <<<"
./verify3d -N 160 -L 6 -relax 1 -steps 200 -dt 1e-3 -profile profile_B1.dat 2>&1
echo

# N=192, L=6 (h=0.0625, ~5.7 points across core) — best chance
echo ">>> N=192 L=6 h=0.0625 (200 relaxation steps) <<<"
./verify3d -N 192 -L 6 -relax 1 -steps 200 -dt 1e-3 -profile profile_B1.dat 2>&1
echo

# ============================================
# Part 3: Higher-B Initialization (no relaxation)
# ============================================
echo "========================================="
echo " PART 3: Higher-B Initialization"
echo "========================================="
echo

for B in 2 3 4; do
    echo ">>> B=$B N=160 L=6 <<<"
    ./verify3d -N 160 -L 6 -B $B -profile profile_B${B}.dat 2>&1 | grep -E "(E_total|E2/E4|Error|Q error|Q )"
    echo
done

echo "========================================="
echo " DONE — $(date)"
echo "========================================="
