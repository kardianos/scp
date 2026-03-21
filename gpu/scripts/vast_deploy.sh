#!/bin/bash
# Deploy simulation to a Vast.ai V100 instance
# Usage: ./vast_deploy.sh <instance_id>
#
# Prerequisites: vastai CLI configured, SSH key trusted
set -e

INSTANCE_ID=$1
if [ -z "$INSTANCE_ID" ]; then
    echo "Usage: $0 <instance_id>"
    echo "First run: ./vast_search.sh to find an offer"
    echo "Then: vastai create instance <offer_id> --image nvidia/cuda:12.2.0-devel-ubuntu22.04 --disk 20 --ssh"
    echo "Then: $0 <instance_id>"
    exit 1
fi

# Get SSH connection details
echo "Getting instance info..."
SSH_INFO=$(vastai ssh-url $INSTANCE_ID)
SSH_HOST=$(echo $SSH_INFO | sed 's/ssh:\/\///' | cut -d: -f1)
SSH_PORT=$(echo $SSH_INFO | cut -d: -f2 | tr -d '/')
echo "SSH: $SSH_HOST:$SSH_PORT"

# Upload source files
echo "Uploading source..."
SRC_DIR="$(dirname $0)/.."
scp -P $SSH_PORT \
    $SRC_DIR/src/*.cu $SRC_DIR/src/*.c \
    $(dirname $0)/../../sfa/format/sfa.h \
    root@$SSH_HOST:~/

# Setup + compile + run
echo "Compiling on remote..."
ssh -p $SSH_PORT root@$SSH_HOST << 'REMOTE'
    apt-get update -qq && apt-get install -y -qq libzstd-dev > /dev/null 2>&1
    nvcc -O3 -arch=sm_70 -o sim cosserat_cuda.cu -lzstd -lm -lpthread
    echo "Build successful. Binary ready."
    nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv
REMOTE

echo ""
echo "Instance ready. To run:"
echo "  ssh -p $SSH_PORT root@$SSH_HOST './sim -N 256 -T 300 -o output.sfa'"
echo ""
echo "To download results:"
echo "  scp -P $SSH_PORT root@$SSH_HOST:~/output.sfa ./"
echo ""
echo "To destroy when done:"
echo "  vastai destroy instance $INSTANCE_ID"
