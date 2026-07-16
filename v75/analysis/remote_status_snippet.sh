#!/bin/bash
# Instant remote status — NEVER use sleep in sim_exec.
# Paste body into: sim_exec(cmd=...)
set -e
echo "DATE=$(date -u)"
echo "=== CAMPAIGN ==="
cat /root/pn/z6/campaign.log 2>/dev/null || echo "(no campaign.log)"
echo "=== PROCS ==="
ps aux | grep -E 'scp_sim|run_z6' | grep -v grep || echo "(no sim proc)"
echo "=== GPU ==="
nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "(no smi)"
echo "=== DISK ==="
df -h /root 2>/dev/null | tail -1
echo "=== DIAGS ==="
ls -la /root/pn/z6/*_diag.tsv 2>/dev/null || echo "(no diags)"
echo "=== NOHUP TAIL ==="
tail -15 /root/pn/z6/nohup.out 2>/dev/null || echo "(no nohup)"
