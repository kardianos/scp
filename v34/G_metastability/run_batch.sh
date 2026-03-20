#!/bin/bash
# Launch a batch of Phase 1 simulations
# Usage: ./run_batch.sh m_value1 m_value2 ...
# m_value is sqrt(m^2), so for m^2=0.25 pass 0.5

BASEDIR=/home/d/code/scp/v34/G_metastability
BIN=$BASEDIR/v33_G
export OMP_NUM_THREADS=4

for MVAL in "$@"; do
    # Compute m^2 label from m value (for directory naming)
    M2_LABEL=$MVAL
    OUTDIR=$BASEDIR/data/G1_scan/m${M2_LABEL}
    mkdir -p "$OUTDIR"

    echo "Launching m^2=${M2_LABEL} (m=${MVAL})..."
    $BIN -N 128 -L 20 -T 500 -braids 1 -bg 0.1 -m $MVAL \
         -diag 5 -snap 100 -o "$OUTDIR" \
         > "$OUTDIR/log.txt" 2>&1 &
    echo "  PID: $!"
done

echo "All launched. Use 'ps aux | grep v33_G' to monitor."
