#!/bin/bash
# Master script: deploy and run field search on Vast.ai GPU
# Usage: ./gpu_run.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

export VAST_API_KEY=$(cat /home/d/code/scp/gpu/key.txt | cut -d' ' -f2)

echo "=== Phase 0: Find and create GPU instance ==="
# Search for cheapest GPU with >=16GB VRAM
echo "Searching for GPU..."
OFFER=$(vastai search offers "num_gpus=1 disk_space>=30 gpu_ram>=14 inet_down>=100 dph_total<=0.60 reliability>0.95" --order 'dph_total' --limit 1 --raw 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(d[0]['id'])" 2>/dev/null || echo "")

if [ -z "$OFFER" ]; then
    echo "No offers found at <$0.60/hr. Widening search..."
    OFFER=$(vastai search offers "num_gpus=1 disk_space>=30 gpu_ram>=14" --order 'dph_total' --limit 1 --raw 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(d[0]['id'])" 2>/dev/null)
fi

echo "Best offer: $OFFER"
echo "Creating instance..."
RESULT=$(vastai create instance $OFFER --image pytorch/pytorch:2.4.1-cuda12.4-cudnn9-devel --disk 30 --onstart-cmd 'pip install "jax[cuda12]" optax -q' --raw 2>/dev/null)
INSTANCE=$(echo "$RESULT" | python3 -c "import json,sys; print(json.load(sys.stdin)['new_contract'])")
echo "Instance: $INSTANCE"

echo "Waiting for instance to start..."
for i in $(seq 1 60); do
    STATUS=$(vastai show instance $INSTANCE --raw 2>/dev/null | python3 -c "import json,sys; print(json.load(sys.stdin)['actual_status'])" 2>/dev/null || echo "loading")
    if [ "$STATUS" = "running" ]; then
        echo "Instance running!"
        break
    fi
    echo "  Status: $STATUS (attempt $i/60)"
    sleep 10
done

# Get SSH info
SSH_PORT=$(vastai show instance $INSTANCE --raw 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['ssh_port'])" 2>/dev/null)
SSH_HOST=$(vastai show instance $INSTANCE --raw 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['ssh_host'])" 2>/dev/null)
echo "SSH: $SSH_HOST:$SSH_PORT"

SSH_CMD="ssh -o StrictHostKeyChecking=no -o ConnectTimeout=10 -p $SSH_PORT root@$SSH_HOST"

echo ""
echo "=== Phase 1: Upload and install ==="
scp -o StrictHostKeyChecking=no -P $SSH_PORT field_search.py root@$SSH_HOST:/root/field_search.py
echo "Uploaded field_search.py"

# Wait for JAX installation
echo "Waiting for JAX install..."
for i in $(seq 1 30); do
    $SSH_CMD "python3 -c 'import jax; print(jax.devices())'" 2>/dev/null && break
    echo "  JAX not ready (attempt $i)..."
    sleep 10
done

echo ""
echo "=== Phase 2: Run validation tests ==="
$SSH_CMD "cd /root && python3 field_search.py --validate 2>&1" | tee validate.log
if grep -q "FAIL" validate.log; then
    echo "VALIDATION FAILED — check validate.log"
    echo "Instance $INSTANCE still running. Debug with:"
    echo "  ssh -p $SSH_PORT root@$SSH_HOST"
    exit 1
fi

echo ""
echo "=== Phase 3: Run field search ==="
# Run in background on the GPU with nohup
$SSH_CMD "cd /root && nohup python3 field_search.py --n_seeds 32 --N 32 --L 15 --epochs 1500 --outdir /root/results_n32 > /root/search.log 2>&1 &"
echo "Search launched in background. Monitor with:"
echo "  $SSH_CMD 'tail -20 /root/search.log'"
echo ""
echo "Instance ID: $INSTANCE"
echo "SSH: $SSH_CMD"
echo ""
echo "To download results:"
echo "  scp -P $SSH_PORT root@$SSH_HOST:/root/results_n32/best_*.npy ."
echo ""
echo "To destroy instance:"
echo "  vastai destroy instance $INSTANCE"
