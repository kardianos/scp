#!/bin/bash
# test_temporal/run_test.sh — Integration test for gen_temporal_seed
#
# 1. Runs proton template simulation (N=64, T=100) with vec output to build temporal model
# 2. Extracts temporal mean seed
# 3. Runs second simulation from that seed (T=50)
# 4. Compares energy drift between original and re-seeded runs

set -e

BINDIR="/home/d/code/scp/bin"
SIM="$BINDIR/scp_sim"
SEED_TOOL="/home/d/code/scp/sfa/seed/gen_temporal_seed"
TEMPLATE="/home/d/code/scp/v43/proton_formation/proton_template.sfa"
WORKDIR="$(cd "$(dirname "$0")" && pwd)"

cd "$WORKDIR"

echo "=== Temporal Seed Extraction Test ==="
echo ""

# Check prerequisites
if [ ! -f "$SIM" ]; then
    echo "ERROR: scp_sim not found at $SIM"
    echo "  Run: make -C /home/d/code/scp/sfa install"
    exit 1
fi
if [ ! -f "$SEED_TOOL" ]; then
    echo "ERROR: gen_temporal_seed not found at $SEED_TOOL"
    echo "  Run: gcc -O3 -fopenmp -o $SEED_TOOL /home/d/code/scp/sfa/seed/gen_temporal_seed.c -lzstd -lm"
    exit 1
fi
if [ ! -f "$TEMPLATE" ]; then
    echo "ERROR: proton template not found at $TEMPLATE"
    exit 1
fi

# Clean previous run
rm -f step1_output.sfa step1_diag.tsv step1.cfg
rm -f temporal_mean_seed.sfa
rm -f step2_output.sfa step2_diag.tsv step2.cfg

# ---- Step 1: Run proton template simulation with vec output ----
echo "--- Step 1: Template simulation (N=64, T=100, vec_snap_dt=0.1) ---"
cat > step1.cfg <<'EOF'
N=64
L=8.3
T=100
dt_factor=0.025
m=1.5
eta=0.5
mu=-41.345
kappa=50
bc_type=2
A_bg=0.1
delta=0,3.0005,4.4325
init=template
init_sfa=/home/d/code/scp/v43/proton_formation/proton_template.sfa
output=step1_output.sfa
diag_file=step1_diag.tsv
precision=f32
snap_dt=50
diag_dt=1.0
vec_snap_dt=0.1
vec_iframe_interval=50
vec_block_size=8
EOF

echo "Running simulation..."
time "$SIM" step1.cfg
echo ""

# ---- Step 2: Extract temporal mean seed ----
echo "--- Step 2: Extract temporal mean seed ---"
"$SEED_TOOL" step1_output.sfa temporal_mean_seed.sfa -amp_report
echo ""

# ---- Step 3: Run from temporal mean seed ----
echo "--- Step 3: Re-seed simulation (N=64, T=50) ---"
cat > step2.cfg <<'EOF'
N=64
L=8.3
T=50
dt_factor=0.025
m=1.5
eta=0.5
mu=-41.345
kappa=50
bc_type=2
A_bg=0.1
delta=0,3.0005,4.4325
init=sfa
init_sfa=temporal_mean_seed.sfa
output=step2_output.sfa
diag_file=step2_diag.tsv
precision=f32
snap_dt=25
diag_dt=1.0
EOF

echo "Running simulation from temporal mean seed..."
time "$SIM" step2.cfg
echo ""

# ---- Step 4: Compare energy drift ----
echo "--- Step 4: Energy comparison ---"
echo ""

# Extract initial and final energies from both runs
# diag.tsv columns: t E_phi_kin E_theta_kin E_grad E_mass E_pot E_tgrad E_tmass E_coupling E_total ...

echo "Step 1 (template run, T=100):"
STEP1_E0=$(awk 'NR==2 {print $10}' step1_diag.tsv)
STEP1_EF=$(awk 'END {print $10}' step1_diag.tsv)
STEP1_T0=$(awk 'NR==2 {print $1}' step1_diag.tsv)
STEP1_TF=$(awk 'END {print $1}' step1_diag.tsv)
echo "  E_total(t=$STEP1_T0) = $STEP1_E0"
echo "  E_total(t=$STEP1_TF) = $STEP1_EF"

# Compute energy stats for step 1 (skip header)
STEP1_STATS=$(awk 'NR>1 {
    e=$10; n++; sum+=e;
    if(NR==2) {min=e; max=e; e0=e}
    if(e<min) min=e; if(e>max) max=e
}
END {
    mean=sum/n
    # Re-scan for stddev
    printf "%.6e %.6e %.6e %.6e %d\n", e0, mean, min, max, n
}' step1_diag.tsv)
echo "  E stats: $STEP1_STATS"

echo ""
echo "Step 2 (temporal mean seed, T=50):"
STEP2_E0=$(awk 'NR==2 {print $10}' step2_diag.tsv)
STEP2_EF=$(awk 'END {print $10}' step2_diag.tsv)
STEP2_T0=$(awk 'NR==2 {print $1}' step2_diag.tsv)
STEP2_TF=$(awk 'END {print $1}' step2_diag.tsv)
echo "  E_total(t=$STEP2_T0) = $STEP2_E0"
echo "  E_total(t=$STEP2_TF) = $STEP2_EF"

STEP2_STATS=$(awk 'NR>1 {
    e=$10; n++; sum+=e;
    if(NR==2) {min=e; max=e; e0=e}
    if(e<min) min=e; if(e>max) max=e
}
END {
    mean=sum/n
    printf "%.6e %.6e %.6e %.6e %d\n", e0, mean, min, max, n
}' step2_diag.tsv)
echo "  E stats: $STEP2_STATS"

echo ""
echo "--- Energy drift comparison ---"
# Compute relative energy spread (max-min)/mean for both runs
awk 'NR>1 {e=$10; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 1 (template):     E_spread/E_mean = %.6f%%\n", spread*100}' step1_diag.tsv

awk 'NR>1 {e=$10; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 2 (mean seed):    E_spread/E_mean = %.6f%%\n", spread*100}' step2_diag.tsv

echo ""

# Compare kinetic energy (indicator of breathing amplitude)
echo "--- Kinetic energy comparison (breathing indicator) ---"
awk 'NR>1 {e=$2; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 1 (template):     E_kin_spread/E_kin_mean = %.4f%%\n", spread*100}' step1_diag.tsv

awk 'NR>1 {e=$2; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 2 (mean seed):    E_kin_spread/E_kin_mean = %.4f%%\n", spread*100}' step2_diag.tsv

echo ""
echo "--- phi_max comparison (field amplitude stability) ---"
awk 'NR>1 {e=$11; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 1 (template):     phi_max spread = %.4f%%\n", spread*100}' step1_diag.tsv

awk 'NR>1 {e=$11; n++; sum+=e; if(NR==2){min=e;max=e} if(e<min)min=e; if(e>max)max=e}
END {mean=sum/n; spread=(max-min)/mean; printf "Step 2 (mean seed):    phi_max spread = %.4f%%\n", spread*100}' step2_diag.tsv

echo ""
echo "=== Test complete ==="
