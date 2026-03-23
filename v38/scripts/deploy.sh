#!/bin/bash
# deploy.sh -- Upload V38 code to GPU server, compile, and run experiments
#
# Usage:
#   ./deploy.sh <user@host> [phase]
#
# Phases:
#   1  Single braid N=256 characterization (default)
#   2  Two-braid interaction N=256
#   3  (Reserved for JAX evolutionary search)
#   4  Seed validation at N=128/256
#
# Prerequisites:
#   - SSH key access to target machine
#   - NVIDIA GPU with CUDA toolkit installed
#   - libzstd-dev installed (apt install libzstd-dev)

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <user@host> [phase=1]"
    echo "  phase 1: Single braid N=256 (est. 30 min on V100)"
    echo "  phase 2: Two-braid interaction N=256 (est. 60 min)"
    echo "  phase 4: Seed validation (requires seed file)"
    exit 1
fi

HOST="$1"
PHASE="${2:-1}"
REMOTE_DIR="/root/v38"
LOCAL_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SFA_DIR="$(cd "$LOCAL_DIR/../sfa/format" && pwd)"

echo "=== V38 GPU Deployment ==="
echo "Host: $HOST"
echo "Phase: $PHASE"
echo "Local: $LOCAL_DIR"
echo ""

# ---- Upload code ----
echo "--- Uploading code ---"
ssh "$HOST" "mkdir -p $REMOTE_DIR/src $REMOTE_DIR/output $REMOTE_DIR/sfa/format"

# Upload source files
scp "$LOCAL_DIR/src/seedrun_cuda.cu" "$HOST:$REMOTE_DIR/src/"
scp "$LOCAL_DIR/src/braid_analyze_cuda.cu" "$HOST:$REMOTE_DIR/src/"
scp "$LOCAL_DIR/scripts/monitor.sh" "$HOST:$REMOTE_DIR/"

# Upload SFA header (single copy, referenced by relative path from src/)
scp "$SFA_DIR/sfa.h" "$HOST:$REMOTE_DIR/sfa/format/"

echo "Upload complete."

# ---- Compile on remote ----
echo ""
echo "--- Compiling CUDA on remote ---"

# Fix include path: source uses #include "../../sfa/format/sfa.h"
# From v38/src/, that resolves to: scp/sfa/format/sfa.h (local layout)
# On remote at /root/v38/src/, it would try /root/sfa/format/sfa.h (missing)
# Solution: create /root/sfa symlink -> /root/v38/sfa so the relative path works
ssh "$HOST" bash << 'COMPILE_EOF'
set -e
cd /root/v38

# Create symlink so ../../sfa from src/ resolves correctly
# /root/v38/src/../../sfa -> /root/sfa -> /root/v38/sfa
ln -sfn /root/v38/sfa /root/sfa

echo "Building braid_analyze..."
nvcc -O3 -arch=sm_70 -o braid_analyze src/braid_analyze_cuda.cu -lzstd -lm

echo "Building seedrun_cuda..."
nvcc -O3 -arch=sm_70 -o seedrun_cuda src/seedrun_cuda.cu -lzstd -lm

echo "Compilation complete."
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
COMPILE_EOF

echo "Compilation done."
echo ""

# ---- Start GPU monitor in background ----
echo "--- Starting GPU monitor ---"
ssh "$HOST" "chmod +x $REMOTE_DIR/monitor.sh && nohup $REMOTE_DIR/monitor.sh > $REMOTE_DIR/gpu_monitor.log 2>&1 &"
echo "Monitor started (logging to gpu_monitor.log)"
echo ""

# ---- Run experiment based on phase ----
case "$PHASE" in
    1)
        echo "=== Phase 1: Single Braid N=256 Characterization ==="
        echo "Estimated time: 30 min on V100"
        ssh "$HOST" bash << 'RUN_P1'
cd /root/v38
mkdir -p output/braid_n256

echo "Starting single braid N=256 run..."
./braid_analyze \
    -N 256 -L 20 -T 200 \
    -eta 0.5 -mt 0 -m 1.5 \
    -snap 10 -aevery 10 \
    -o output/braid_n256 \
    2>&1 | tee output/braid_n256/run.log

echo ""
echo "=== Phase 1 Complete ==="
echo "Output files:"
ls -lh output/braid_n256/
RUN_P1
        ;;
    2)
        echo "=== Phase 2: Two-Braid Interaction N=256 ==="
        echo "Estimated time: 60 min on V100"
        ssh "$HOST" bash << 'RUN_P2'
cd /root/v38
mkdir -p output/two_braid_n256

echo "Starting two-braid N=256 run (D=12)..."
./braid_analyze \
    -N 256 -L 20 -T 300 \
    -eta 0.5 -mt 0 -m 1.5 \
    -braids 2 -D 12 \
    -snap 10 -aevery 10 \
    -o output/two_braid_n256 \
    2>&1 | tee output/two_braid_n256/run.log

echo ""
echo "=== Phase 2 Complete ==="
echo "Output files:"
ls -lh output/two_braid_n256/
RUN_P2
        ;;
    4)
        echo "=== Phase 4: Seed Validation ==="
        if [ $# -lt 3 ]; then
            echo "Usage for phase 4: $0 <host> 4 <seed_file>"
            exit 1
        fi
        SEED="$3"
        echo "Uploading seed: $SEED"
        scp "$SEED" "$HOST:$REMOTE_DIR/seed.bin"

        ssh "$HOST" bash << 'RUN_P4'
cd /root/v38
mkdir -p output/seed_validate

echo "Starting seed validation run..."
./seedrun_cuda \
    -seed seed.bin \
    -N 256 -L 15 -T 500 \
    -eta 0.5 -mt 0 -m 1.5 \
    -snap 5 -diag 1 \
    -damp_width 4 -damp_rate 0.005 \
    -o output/seed_validate \
    2>&1 | tee output/seed_validate/run.log

echo ""
echo "=== Phase 4 Complete ==="
echo "Output files:"
ls -lh output/seed_validate/
RUN_P4
        ;;
    *)
        echo "Unknown phase: $PHASE"
        exit 1
        ;;
esac

echo ""
echo "--- Downloading results ---"
mkdir -p "$LOCAL_DIR/results/phase${PHASE}"
scp -r "$HOST:$REMOTE_DIR/output/" "$LOCAL_DIR/results/phase${PHASE}/"
scp "$HOST:$REMOTE_DIR/gpu_monitor.log" "$LOCAL_DIR/results/phase${PHASE}/"

echo ""
echo "=== Deployment Complete ==="
echo "Results in: $LOCAL_DIR/results/phase${PHASE}/"
echo "GPU monitor log: $LOCAL_DIR/results/phase${PHASE}/gpu_monitor.log"
