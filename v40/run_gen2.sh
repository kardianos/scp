#!/bin/bash
# Gen 2: Refinement of Gen 1 winners + 3-braid configurations
# Survivors: C4 (1.5× amp), C3 (counter-braid), C7 (high ellip)
# New: 3-braid xyz with chirality variations

set -e
cd /home/d/code/scp/v40
SIM=/home/d/code/scp/sfa/sim/scp_sim
ANALYZE=tools/analyze_sfa
MODIFY=tools/modify_sfa
GEN3=tools/gen_3braid
BASELINE_SFA=/home/d/code/scp/v39/data/braid_control.sfa
export OMP_NUM_THREADS=16
T=20

echo "=== Gen 2: 12 Candidates (T=$T) ==="
echo ""

mkdir -p gen_002

run_candidate() {
    local id=$1 desc="$2" seed_cmd="$3"
    local dir="gen_002/candidate_${id}"
    mkdir -p "$dir"

    echo "--- ${id}: ${desc} ---"

    # Generate seed
    eval "$seed_cmd"

    # Write config
    cat > "$dir/config.cfg" << CFGEOF
# Gen 2, ${id}: ${desc}
N = 128
L = 15.0
T = ${T}
dt_factor = 0.025
m = 1.5
m_theta = 0.0
eta = 0.5
mu = -41.345
kappa = 50.0
mode = 0
damp_width = 3.0
damp_rate = 0.01
init = sfa
init_sfa = $(realpath "$dir/seed.sfa")
init_frame = 0
output = $(realpath "$dir/output.sfa")
precision = f32
snap_dt = 5.0
diag_dt = 2.0
diag_file = $(realpath "$dir/diag.tsv")
CFGEOF

    # Run simulation
    rm -f "$dir/output.sfa"
    if $SIM "$dir/config.cfg" > "$dir/sim_log.txt" 2>&1; then
        echo "  SIM: OK"
    else
        echo "  SIM: FAILED (exit $?)"
        return
    fi

    # Fix SFA and analyze
    tools/fixup_all "$dir/output.sfa" > /dev/null 2>&1 || true
    if $ANALYZE "$dir/output.sfa" --json "$dir/analysis.json" > "$dir/analysis_log.txt" 2>&1; then
        score=$(grep -o '"S_final": [0-9.]*' "$dir/analysis.json" | grep -o '[0-9.]*')
        epot=$(grep -o '"E_pot_final": [-0-9.]*' "$dir/analysis.json" | grep -o '[-0-9.]*')
        echo "  SCORE: S=${score} E_pot=${epot}"
    else
        echo "  ANALYZE: FAILED"
    fi
    echo ""
}

# === From C4 (amplitude scaling) ===

run_candidate "C4a" "Scale 1.3x" \
    "$MODIFY $BASELINE_SFA --scale 1.3 -o gen_002/candidate_C4a/seed.sfa"

run_candidate "C4b" "Scale 1.7x" \
    "$MODIFY $BASELINE_SFA --scale 1.7 -o gen_002/candidate_C4b/seed.sfa"

run_candidate "C4c" "Scale 2.0x" \
    "$MODIFY $BASELINE_SFA --scale 2.0 -o gen_002/candidate_C4c/seed.sfa"

# === From C3 (counter-braid) — combine with amplitude scaling ===

run_candidate "C3a" "Counter-braid at 1.5x amplitude" \
    "$MODIFY $BASELINE_SFA --scale 1.5 --add-braid 5 0 0 0 --flip-chirality -o gen_002/candidate_C3a/seed.sfa"

run_candidate "C3b" "Counter-braid tighter (shift 3)" \
    "$MODIFY $BASELINE_SFA --add-braid 3 0 0 0 --flip-chirality -o gen_002/candidate_C3b/seed.sfa"

run_candidate "C3c" "Counter-braid at 1.7x amplitude" \
    "$MODIFY $BASELINE_SFA --scale 1.7 --add-braid 5 0 0 0 --flip-chirality -o gen_002/candidate_C3c/seed.sfa"

# === From C7 (high ellipticity) — combine with amplitude ===

run_candidate "C7a" "High ellip (0.5) at 1.5x" \
    "/home/d/code/scp/sfa/seed/gen_braid -N 128 -L 15 -A 1.2 -ellip 0.5 -o gen_002/candidate_C7a/seed.sfa"

run_candidate "C7b" "Very high ellip (0.65)" \
    "/home/d/code/scp/sfa/seed/gen_braid -N 128 -L 15 -A 0.8 -ellip 0.65 -o gen_002/candidate_C7b/seed.sfa"

# === 3-braid xyz configurations (user request) ===

run_candidate "X3a" "3-braid xyz UUU, A=0.4" \
    "$GEN3 -N 128 -L 15 -A 0.4 -chirality UUU -o gen_002/candidate_X3a/seed.sfa"

run_candidate "X3b" "3-braid xyz UUD, A=0.4" \
    "$GEN3 -N 128 -L 15 -A 0.4 -chirality UUD -o gen_002/candidate_X3b/seed.sfa"

run_candidate "X3c" "3-braid xyz UDD, A=0.4" \
    "$GEN3 -N 128 -L 15 -A 0.4 -chirality UDD -o gen_002/candidate_X3c/seed.sfa"

run_candidate "X3d" "3-braid xyz UUD, A=0.5 (higher)" \
    "$GEN3 -N 128 -L 15 -A 0.5 -chirality UUD -o gen_002/candidate_X3d/seed.sfa"

# === Summary ===

echo "=== Gen 2 Summary ==="
printf "%-6s %-42s %10s %10s\n" "ID" "Description" "S_final" "E_pot"
printf "%-6s %-42s %10s %10s\n" "---" "---" "---" "---"
for dir in gen_002/candidate_*/; do
    id=$(basename "$dir" | sed 's/candidate_//')
    desc=$(head -1 "$dir/config.cfg" | sed 's/^# Gen 2, [^:]*: //')
    if [ -f "$dir/analysis.json" ]; then
        score=$(grep -o '"S_final": [0-9.]*' "$dir/analysis.json" | grep -o '[0-9.]*')
        epot=$(grep -o '"E_pot_final": [-0-9.]*' "$dir/analysis.json" | grep -o '[-0-9.]*')
        printf "%-6s %-42s %10s %10s\n" "$id" "$desc" "$score" "$epot"
    else
        printf "%-6s %-42s %10s %10s\n" "$id" "$desc" "FAILED" "FAILED"
    fi
done
