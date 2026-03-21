#!/bin/bash
# Record SFA files for all surviving braid geometries
# Run after the crossed chirality experiments complete
cd /home/d/code/scp/v37

echo "=== Recording SFA for all surviving geometries ==="

# Build the SFA-enabled version if not already built
if [ ! -f v37_sfa ]; then
    gcc -O3 -march=native -fopenmp -o v37_sfa src/v37_sfa.c -lzstd -lm
fi

# Build the crossed version with SFA
if [ -f src/v37_crossed.c ]; then
    # Add SFA support to crossed version
    cp src/v37_crossed.c src/v37_crossed_sfa.c
    # Insert SFA include after the first #include block
    sed -i '/#include <string.h>/a\\n#define SFA_IMPLEMENTATION\n#include "sfa.h"' src/v37_crossed_sfa.c 2>/dev/null
    gcc -O3 -march=native -fopenmp -o v37_crossed_sfa src/v37_crossed_sfa.c -lzstd -lm 2>/dev/null
fi

echo ""
echo "=== 1. braid3(z) — proven winner ==="
mkdir -p data/sfa_braid3z
OMP_NUM_THREADS=8 ./v37_sfa -geom braid3 -N 128 -L 15 -T 300 -eta 0.5 -mt 0 -diag 5 -snap 3 -o data/sfa_braid3z 2>&1 | tail -10
echo "SFA: $(ls -lh data/sfa_braid3z.sfa 2>/dev/null)"

echo ""
echo "=== 2. Checking crossed results ==="
# Check which chiralities survived
for ch in UUU UUD UDD DDD; do
    ts="data/crossed_${ch}/timeseries.tsv"
    if [ -f "$ts" ] && [ $(wc -l < "$ts") -gt 10 ]; then
        lines=$(wc -l < "$ts")
        last=$(tail -1 "$ts")
        t=$(echo "$last" | awk '{print $1}')
        Ep=$(echo "$last" | awk '{printf "%.1f", $6}')
        echo "  $ch: t=$t Ep=$Ep ($lines lines)"

        # If it survived (E_pot still significantly negative), record SFA
        # Check if E_pot magnitude > 1 (alive)
        alive=$(echo "$last" | awk '{if ($6 < -1.0) print "YES"; else print "NO"}')
        if [ "$alive" = "YES" ]; then
            echo "  $ch SURVIVED — recording SFA..."
            if [ -f v37_crossed_sfa ]; then
                mkdir -p data/sfa_crossed_${ch}
                OMP_NUM_THREADS=8 ./v37_crossed_sfa -geom crossed -chiral $ch -N 128 -L 15 -T 300 -eta 0.5 -mt 0 -diag 5 -snap 3 -o data/sfa_crossed_${ch} 2>&1 | tail -5
                echo "  SFA: $(ls -lh data/sfa_crossed_${ch}.sfa 2>/dev/null)"
            else
                echo "  (SFA build failed — using timeseries only)"
            fi
        else
            echo "  $ch DEAD — skipping SFA"
        fi
    else
        echo "  $ch: no data or too few points"
    fi
done

echo ""
echo "=== 3. Truncated helix (previous winner) ==="
# Already have this from earlier run
if [ -f data/truncated_sfa.sfa ]; then
    echo "  Already recorded: $(ls -lh data/truncated_sfa.sfa)"
else
    echo "  Recording truncated helix..."
    mkdir -p data/sfa_truncated
    OMP_NUM_THREADS=8 ./v37_sfa -geom truncated -N 128 -L 15 -T 300 -eta 0.5 -mt 0 -diag 5 -snap 3 -o data/sfa_truncated 2>&1 | tail -5
    echo "  SFA: $(ls -lh data/sfa_truncated.sfa 2>/dev/null)"
fi

echo ""
echo "=== Recording complete ==="
echo "Available SFA files:"
ls -lh data/*.sfa data/sfa_*/*.sfa 2>/dev/null
echo ""
echo "View with: /home/d/code/scp/sfa/viewer/volview <file.sfa>"
