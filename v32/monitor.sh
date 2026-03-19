#!/bin/bash
# Monitor running E1/E2 experiments
# Usage: ./monitor.sh [interval_seconds]

INTERVAL=${1:-60}

while true; do
    clear
    echo "=== V32 Campaign Monitor ($(date)) ==="
    echo ""

    # Check running processes
    echo "--- Running Processes ---"
    ps aux | grep v32_run | grep -v grep | grep -v defunct | \
        awk '{printf "PID %s CPU=%s%% MEM=%.1fGB (%s)\n", $2, $3, $6/1024/1024, $11" "$12" "$13" "$14" "$15}'
    echo ""

    # Check E2 progress
    if [ -f data/E2/timeseries.tsv ]; then
        echo "--- E2 (Two Braids, N=512, L=40) ---"
        lines=$(wc -l < data/E2/timeseries.tsv)
        echo "Data points: $((lines - 1))"
        if [ "$lines" -gt 1 ]; then
            tail -1 data/E2/timeseries.tsv | awk -F'\t' '{
                printf "t=%.1f  E=%.3e  D=%.2f  fc=%.3f  w=%.3f  max_rho=%.2e\n", $1, $6, $9, $7, $10, $8
            }'
            # Estimate completion
            t=$(tail -1 data/E2/timeseries.tsv | awk -F'\t' '{print $1}')
            python3 -c "
t=$t
T=2000
if t > 0:
    # 51100 steps total
    frac = t/T
    print(f'Progress: {100*frac:.1f}%')
" 2>/dev/null
        fi
        echo ""
    fi

    # Check E1 progress
    if [ -f data/E1/timeseries.tsv ]; then
        echo "--- E1 (Single Braid, N=512, L=20) ---"
        lines=$(wc -l < data/E1/timeseries.tsv)
        echo "Data points: $((lines - 1))"
        if [ "$lines" -gt 1 ]; then
            tail -1 data/E1/timeseries.tsv | awk -F'\t' '{
                printf "t=%.1f  E=%.3e  fc=%.3f  w=%.3f  max_rho=%.2e\n", $1, $6, $7, $10, $8
            }'
            t=$(tail -1 data/E1/timeseries.tsv | awk -F'\t' '{print $1}')
            python3 -c "
t=$t
T=2000
if t > 0:
    frac = t/T
    print(f'Progress: {100*frac:.1f}%')
" 2>/dev/null
        fi
        echo ""
    fi

    # Memory usage
    echo "--- Memory ---"
    free -g | grep Mem

    # Disk usage
    echo ""
    echo "--- Disk Usage ---"
    du -sh data/E1/ data/E2/ 2>/dev/null

    echo ""
    echo "Press Ctrl+C to stop monitoring (processes continue)"
    sleep $INTERVAL
done
