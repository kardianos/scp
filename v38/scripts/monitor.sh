#!/bin/bash
# monitor.sh -- GPU utilization monitor
#
# Runs in background, checks GPU every 10s. Alerts if GPU is idle.
# Per CLAUDE.md: ALL remote simulations MUST run on GPU.
#
# Usage: nohup ./monitor.sh > gpu_monitor.log 2>&1 &

INTERVAL=10      # seconds between checks
IDLE_THRESH=5    # GPU util% below which we consider idle
IDLE_LIMIT=6     # consecutive idle checks before alert (60s)
LOGFILE="${1:-gpu_monitor.log}"

idle_count=0
total_checks=0
total_idle=0

echo "=== GPU Monitor Started ($(date)) ==="
echo "Checking every ${INTERVAL}s, alert after ${IDLE_LIMIT} consecutive idle"
echo ""

while true; do
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # Query GPU utilization
    gpu_info=$(nvidia-smi --query-gpu=utilization.gpu,utilization.memory,memory.used,memory.total,temperature.gpu,power.draw \
                          --format=csv,noheader,nounits 2>/dev/null)

    if [ $? -ne 0 ]; then
        echo "[$timestamp] ERROR: nvidia-smi failed -- no GPU detected?"
        sleep "$INTERVAL"
        continue
    fi

    gpu_util=$(echo "$gpu_info" | cut -d',' -f1 | tr -d ' ')
    mem_util=$(echo "$gpu_info" | cut -d',' -f2 | tr -d ' ')
    mem_used=$(echo "$gpu_info" | cut -d',' -f3 | tr -d ' ')
    mem_total=$(echo "$gpu_info" | cut -d',' -f4 | tr -d ' ')
    temp=$(echo "$gpu_info" | cut -d',' -f5 | tr -d ' ')
    power=$(echo "$gpu_info" | cut -d',' -f6 | tr -d ' ')

    total_checks=$((total_checks + 1))

    # Check for idle GPU
    if [ "$gpu_util" -lt "$IDLE_THRESH" ] 2>/dev/null; then
        idle_count=$((idle_count + 1))
        total_idle=$((total_idle + 1))

        if [ "$idle_count" -ge "$IDLE_LIMIT" ]; then
            echo "[$timestamp] *** ALERT: GPU IDLE for ${idle_count}x${INTERVAL}s! util=${gpu_util}% mem=${mem_used}/${mem_total}MiB ***"
            echo "[$timestamp] *** Check if simulation is running on CPU instead of GPU! ***"
            # Also check what processes are using the GPU
            echo "[$timestamp] GPU processes:"
            nvidia-smi --query-compute-apps=pid,name,used_memory --format=csv,noheader 2>/dev/null || true
            echo ""
        else
            echo "[$timestamp] IDLE: gpu=${gpu_util}% mem=${mem_used}/${mem_total}MiB T=${temp}C P=${power}W (idle ${idle_count}/${IDLE_LIMIT})"
        fi
    else
        if [ "$idle_count" -ge "$IDLE_LIMIT" ]; then
            echo "[$timestamp] GPU resumed: util=${gpu_util}%"
        fi
        idle_count=0
        echo "[$timestamp] OK: gpu=${gpu_util}% mem=${mem_used}/${mem_total}MiB T=${temp}C P=${power}W"
    fi

    # Periodic summary every 5 minutes (30 checks)
    if [ $((total_checks % 30)) -eq 0 ]; then
        idle_pct=$((total_idle * 100 / total_checks))
        echo "[$timestamp] --- Summary: ${total_checks} checks, ${total_idle} idle (${idle_pct}%) ---"
    fi

    sleep "$INTERVAL"
done
