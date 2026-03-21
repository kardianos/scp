#!/bin/bash
# Deploy structure search to Vast.ai GPU instance
# Usage: ./deploy.sh <instance_id>

INSTANCE_ID=$1
if [ -z "$INSTANCE_ID" ]; then
    echo "Usage: ./deploy.sh <instance_id>"
    echo "First create an instance with:"
    echo "  vastai create instance <offer_id> --image pytorch/pytorch:2.4.1-cuda12.4-cudnn9-devel --disk 30 --onstart-cmd 'pip install jax[cuda12] optax'"
    exit 1
fi

# Wait for instance to be ready
echo "Waiting for instance $INSTANCE_ID..."
vastai wait instance $INSTANCE_ID

# Get SSH info
SSH_INFO=$(vastai ssh-url $INSTANCE_ID)
echo "SSH: $SSH_INFO"

# Upload the search script
echo "Uploading structure_search.py..."
vastai copy $INSTANCE_ID structure_search.py /root/structure_search.py

# Run the search
echo "Running structure search..."
vastai execute $INSTANCE_ID "cd /root && pip install jax[cuda12] optax -q && python structure_search.py --N 32 --L 12 --T_opt 30 --n_epochs 300 --seed braid --outdir /root/results 2>&1 | tee /root/run.log"

echo "Downloading results..."
vastai copy $INSTANCE_ID /root/results/ ./results/
vastai copy $INSTANCE_ID /root/run.log ./run.log

echo "Done. Destroy instance with: vastai destroy instance $INSTANCE_ID"
