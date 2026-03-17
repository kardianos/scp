#!/bin/bash
set -e

LAMBDAS="0.0 0.2 0.4 0.6 0.7 0.8 0.85 0.90 0.95"

for lam in $LAMBDAS; do
    echo "=== lambda=$lam mode=120 ==="
    ./pw120 -lambda $lam -mode 120 -o data 2>&1 | tail -10
    echo ""
    echo "=== lambda=$lam mode=0 ==="
    ./pw120 -lambda $lam -mode 0 -o data 2>&1 | tail -10
    echo ""
done
