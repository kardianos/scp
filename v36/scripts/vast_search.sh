#!/bin/bash
# Search for best V100 offer on Vast.ai
# Usage: ./vast_search.sh [num_gpus]
NUM_GPUS=${1:-1}
echo "Searching for V100 × $NUM_GPUS with ≥20 GB disk, ≥200 Mbps..."
vastai search offers \
    "gpu_name=V100 num_gpus=$NUM_GPUS disk_space>=20 inet_down>=200 dph_total<=${2:-0.50}" \
    --order 'dph_total' \
    --limit 10
