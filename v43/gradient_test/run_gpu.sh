#!/bin/bash
# GPU gradient force test — run on Vast.ai V100
# Usage: ssh into the instance, then run this script
set -e

echo "=== Gradient Force Test — V43 ==="
echo "Testing proton and braid response to density gradient"
echo ""

# --- Compile GPU kernel ---
echo "Building GPU kernel..."
cd /root
nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm
echo "Build OK"

# --- Generate seeds ---
echo ""
echo "=== Generating seeds ==="

# Build seed generators
gcc -O3 -fopenmp -o gen_composite gen_composite.c -lzstd -lm
gcc -O3 -fopenmp -o gen_braid gen_braid.c -lzstd -lm

# Proton seeds (from averaged template)
echo "Generating proton seed (gentle gradient)..."
./gen_composite -N 512 -L 100 \
    -template proton_averaged.sfa -place 0,0,0 \
    -gradient_A_high 0.12 -gradient_A_low 0.08 \
    -precision f32 -o proton_gentle_seed.sfa

echo "Generating proton seed (steep gradient)..."
./gen_composite -N 512 -L 100 \
    -template proton_averaged.sfa -place 0,0,0 \
    -gradient_A_high 0.15 -gradient_A_low 0.05 \
    -precision f32 -o proton_steep_seed.sfa

# Braid seeds (single z-aligned braid in gradient)
echo "Generating braid seed (gentle gradient)..."
./gen_braid -N 512 -L 100 -A 0.8 -R 3.0 \
    -precision f32 -o braid_gentle_seed_tmp.sfa
# gen_braid doesn't support gradients, so we use gen_composite with no template
# Actually, gen_braid produces a braid at center. We need to add gradient background.
# Use gen_composite with the braid as a "template" — but gen_braid output isn't a 64^3 template.
# Instead: generate braid seed normally, then the sim kernel handles the gradient via bc_type=1
# The gradient is maintained by the pinned BCs, not the initial condition.
# So we just need a braid in a UNIFORM background, and bc_type=1 creates the gradient.

# Actually — the gradient must be in the INITIAL condition too, for pinning to work.
# The pinned BC saves the initial state. If initial has no gradient, pinning maintains no gradient.
# So we need gradient in the seed.

# For the braid: init with gradient background + braid at center
# gen_braid only does uniform background. We need to modify or use gen_composite.
# gen_composite can place any template. A single braid IS a template if we extract it.
# But simpler: generate a uniform braid seed, then modify in the sim with bc_type=1.
# Wait — the sim with bc_type=1 pins the INITIAL state. So the initial state MUST have the gradient.

# Solution: use gen_composite with a braid template.
# First generate a small braid template:
echo "Generating braid template..."
./gen_braid -N 64 -L 8.3 -A 0.8 -R 3.0 -precision f32 -o braid_template.sfa

echo "Generating braid seed (gentle gradient) via gen_composite..."
./gen_composite -N 512 -L 100 \
    -template braid_template.sfa -place 0,0,0 \
    -gradient_A_high 0.12 -gradient_A_low 0.08 \
    -precision f32 -o braid_gentle_seed.sfa

echo "Generating braid seed (steep gradient) via gen_composite..."
./gen_composite -N 512 -L 100 \
    -template braid_template.sfa -place 0,0,0 \
    -gradient_A_high 0.15 -gradient_A_low 0.05 \
    -precision f32 -o braid_steep_seed.sfa

echo "Seeds generated."

# --- Run experiments ---
echo ""
echo "=== Running gradient tests ==="

# 1. Proton, gentle gradient, 6-field
echo "[1/8] Proton gentle eta=0.5..."
./scp_sim_cuda gradient_proton_gentle.cfg
echo "  Done."

# 2. Proton, steep gradient, 6-field
echo "[2/8] Proton steep eta=0.5..."
./scp_sim_cuda gradient_proton_steep.cfg
echo "  Done."

# 3. Braid, gentle gradient, 6-field
echo "[3/8] Braid gentle eta=0.5..."
./scp_sim_cuda gradient_braid_gentle.cfg
echo "  Done."

# 4. Braid, steep gradient, 6-field
echo "[4/8] Braid steep eta=0.5..."
./scp_sim_cuda gradient_braid_steep.cfg
echo "  Done."

# 5-8: Repeat with eta=0 (3-field, gravity only)
for obj in proton braid; do
    for grad in gentle steep; do
        echo "[eta=0] ${obj} ${grad}..."
        cfg="${obj}_${grad}_eta0.cfg"
        # Create eta=0 config by copying and overriding eta
        cp gradient_${obj}_${grad}.cfg $cfg
        sed -i 's/eta = 0.5/eta = 0.0/' $cfg
        sed -i "s/output = ${obj}_${grad}_output/output = ${obj}_${grad}_eta0_output/" $cfg
        sed -i "s/diag_file = ${obj}_${grad}_diag/diag_file = ${obj}_${grad}_eta0_diag/" $cfg
        ./scp_sim_cuda $cfg
        echo "  Done."
    done
done

echo ""
echo "=== All 8 runs complete ==="
echo ""

# --- Quick analysis ---
echo "=== Centroid tracking ==="
for f in *_diag.tsv; do
    echo "--- $f ---"
    # The diag files don't have centroid. We need to run spatial_analysis on the SFA.
    # For now, just show energy evolution
    head -1 "$f"
    tail -1 "$f"
    echo ""
done

echo "=== Done. Download diag.tsv and output.sfa files for analysis ==="
ls -lh *_diag.tsv *_output.sfa
