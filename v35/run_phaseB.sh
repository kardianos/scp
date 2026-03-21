#!/bin/bash
# Phase B: θ-braid initialization tests
cd /home/d/code/scp/v35

MAX_JOBS=3

# Test configurations: (mu_t, kap_t, A_theta, label)
configs=(
    # mu=-10000, kap=50000 with different theta amplitudes
    "-10000 50000 0.05 mu1e4_At0.05"
    "-10000 50000 0.10 mu1e4_At0.10"
    "-10000 50000 0.20 mu1e4_At0.20"
    "-10000 50000 0.50 mu1e4_At0.50"
    # mu=-1000, kap=50000 (slower growth, might preserve structure longer)
    "-1000 50000 0.05 mu1e3_At0.05"
    "-1000 50000 0.20 mu1e3_At0.20"
    # Control: no theta self-interaction (does the theta-braid disperse?)
    "0 50000 0.20 control_At0.20"
    # Also test with NO phi-braid, just a theta-braid
    "-10000 50000 0.20 nophi_mu1e4_At0.20"
)

for conf in "${configs[@]}"; do
    read mu kap At label <<< "$conf"

    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 2
    done

    outdir="data/theta_self/phaseB_${label}"
    echo "=== Launching: $label (mu_t=$mu kap_t=$kap A_theta=$At) ==="

    extra_args=""
    if [[ "$label" == nophi_* ]]; then
        extra_args="-bg 0"
        # For no-phi case, we use a tiny bg to avoid exact zero phi
        # Actually we need a phi-braid for the Cosserat coupling...
        # Let's just set ETA=0 so theta is independent
        extra_args="-eta 0 -bg 0"
    fi

    OMP_NUM_THREADS=4 ./v33_theta_self \
        -N 64 -L 20 -T 200 -eta 0.5 -mt 0 \
        -mu_t $mu -kap_t $kap \
        -A_theta $At -theta_x 10 \
        -diag 2 -o "$outdir" $extra_args \
        > "${outdir}.log" 2>&1 &
done

echo "All Phase B jobs launched. Waiting..."
wait
echo "=== All Phase B complete ==="
