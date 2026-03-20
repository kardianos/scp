#!/bin/bash
# Phase 1: Fine m^2 scan
# Launches simulations for given m^2 values
# The -m flag squares its argument, so we pass sqrt(m^2)

BASEDIR=/home/d/code/scp/v34/G_metastability
BIN=$BASEDIR/v33_G
export OMP_NUM_THREADS=4

# m^2 values and their sqrt for the -m flag
declare -A M2_TO_M
M2_TO_M[0.25]=0.5
M2_TO_M[0.50]=0.7071067811865476
M2_TO_M[0.75]=0.8660254037844387
M2_TO_M[1.00]=1.0
M2_TO_M[1.25]=1.118033988749895
M2_TO_M[1.50]=1.2247448713915890
M2_TO_M[1.75]=1.3228756555322954
M2_TO_M[2.00]=1.4142135623730951
M2_TO_M[2.25]=1.5

for M2 in "$@"; do
    MVAL=${M2_TO_M[$M2]}
    if [ -z "$MVAL" ]; then
        echo "ERROR: Unknown m^2=$M2"
        continue
    fi
    OUTDIR=$BASEDIR/data/G1_scan/m${M2}
    mkdir -p "$OUTDIR"
    echo "Launching m^2=${M2} (m=${MVAL}) -> $OUTDIR"
    $BIN -N 128 -L 20 -T 500 -braids 1 -bg 0.1 -m $MVAL \
         -diag 5 -snap 100 -o "$OUTDIR" \
         > "$OUTDIR/log.txt" 2>&1 &
    echo "  PID: $!"
done

echo "Batch launched. Monitor with: ps aux | grep v33_G"
