#!/bin/bash
echo "=== Gravity run (N=400, alpha=0.1) ==="
tail -3 /home/d/code/scp/v22/data/production_grav.txt
echo ""
echo "=== Control run (N=350, alpha=0) ==="
tail -3 /home/d/code/scp/v22/data/control_nograv_n350.txt
echo ""
echo "=== Separation comparison ==="
# Extract time and separation from both
echo "Time  Grav_sep  Ctrl_sep  Delta"
paste <(grep "^  t=" /home/d/code/scp/v22/data/production_grav.txt | awk '{print $1"="$2, $3"="$4}') \
      <(grep "^  t=" /home/d/code/scp/v22/data/control_nograv_n350.txt | awk '{print $3"="$4}') 2>/dev/null | head -5
