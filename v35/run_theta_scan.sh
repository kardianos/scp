#!/bin/bash
# Phase A: θ self-interaction parameter scan
# Run 3-4 concurrent jobs from 18 combinations
cd /home/d/code/scp/v35

MU_VALS="0 -1000 -10000 -100000 -500000 -1000000"
KAP_VALS="50 5000 50000"

MAX_JOBS=4

for mu in $MU_VALS; do
    for kap in $KAP_VALS; do
        # Wait if we have too many jobs
        while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
            sleep 2
        done

        label="mu${mu}_kap${kap}"
        outdir="data/theta_self/${label}"
        echo "=== Launching: mu_t=${mu} kap_t=${kap} -> ${outdir} ==="

        OMP_NUM_THREADS=4 ./v33_theta_self \
            -N 64 -L 20 -T 100 -eta 0.5 -mt 0 \
            -mu_t $mu -kap_t $kap \
            -diag 2 -o "$outdir" \
            > "${outdir}.log" 2>&1 &
    done
done

echo "All jobs launched. Waiting..."
wait
echo "=== All complete ==="
