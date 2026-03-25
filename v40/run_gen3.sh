#!/bin/bash
# Gen 3: Evolve top candidates at larger grid (N=192, L=25)
# From user review: pursue X3c, C4c, X3b, C3a, C4a, C4b
# Focus: X3c refinement, UDD amplitude sweep, 3-braid with tighter R_tube

set -e
cd /home/d/code/scp/v40
SIM=/home/d/code/scp/sfa/sim/scp_sim
ANALYZE=tools/analyze_sfa
MODIFY=tools/modify_sfa
GEN3=tools/gen_3braid
export OMP_NUM_THREADS=16
T=30
N=192
L=25

echo "=== Gen 3: Larger Grid (N=$N L=$L T=$T) ==="
echo "  dx = $(echo "scale=4; 2*$L/($N-1)" | bc)"
echo ""

mkdir -p gen_003

run_candidate() {
    local id=$1 desc="$2" seed_cmd="$3"
    local dir="gen_003/candidate_${id}"
    mkdir -p "$dir"
    echo "--- ${id}: ${desc} ---"
    eval "$seed_cmd"
    cat > "$dir/config.cfg" << CFGEOF
# Gen 3, ${id}: ${desc}
N = ${N}
L = ${L}.0
T = ${T}.0
dt_factor = 0.025
m = 1.5
m_theta = 0.0
eta = 0.5
mu = -41.345
kappa = 50.0
mode = 0
damp_width = 5.0
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

    rm -f "$dir/output.sfa"
    if $SIM "$dir/config.cfg" > "$dir/sim_log.txt" 2>&1; then
        echo "  SIM: OK"
    else
        echo "  SIM: FAILED (exit $?)"
        return
    fi
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

# === X3c evolution: UDD amplitude and geometry sweep ===

run_candidate "UDD_A04" "UDD A=0.4 (baseline at larger grid)" \
    "$GEN3 -N $N -L $L -A 0.4 -chirality UDD -o gen_003/candidate_UDD_A04/seed.sfa"

run_candidate "UDD_A05" "UDD A=0.5 (higher amplitude)" \
    "$GEN3 -N $N -L $L -A 0.5 -chirality UDD -o gen_003/candidate_UDD_A05/seed.sfa"

run_candidate "UDD_A06" "UDD A=0.6 (push amplitude)" \
    "$GEN3 -N $N -L $L -A 0.6 -chirality UDD -o gen_003/candidate_UDD_A06/seed.sfa"

run_candidate "UDD_R2" "UDD A=0.5 R_tube=2.0 (tighter)" \
    "$GEN3 -N $N -L $L -A 0.5 -R 2.0 -chirality UDD -o gen_003/candidate_UDD_R2/seed.sfa"

run_candidate "UDD_R4" "UDD A=0.5 R_tube=4.0 (wider)" \
    "$GEN3 -N $N -L $L -A 0.5 -R 4.0 -chirality UDD -o gen_003/candidate_UDD_R4/seed.sfa"

# === X3b evolution: UUD at optimal amplitude ===

run_candidate "UUD_A05" "UUD A=0.5 (larger grid)" \
    "$GEN3 -N $N -L $L -A 0.5 -chirality UUD -o gen_003/candidate_UUD_A05/seed.sfa"

run_candidate "UUD_A06" "UUD A=0.6 (push)" \
    "$GEN3 -N $N -L $L -A 0.6 -chirality UUD -o gen_003/candidate_UUD_A06/seed.sfa"

# === C4c/C4a evolution: scaled braid at larger grid ===

run_candidate "S15" "Scale 1.5x (larger grid)" \
    "/home/d/code/scp/sfa/seed/gen_braid -N $N -L $L -A 1.2 -o gen_003/candidate_S15/seed.sfa"

run_candidate "S20" "Scale 2.0x (larger grid)" \
    "/home/d/code/scp/sfa/seed/gen_braid -N $N -L $L -A 1.6 -o gen_003/candidate_S20/seed.sfa"

# === C3a evolution: counter-braid at larger grid ===

run_candidate "CB15" "Counter-braid 1.5x (larger grid)" \
    "/home/d/code/scp/sfa/seed/gen_braid -N $N -L $L -A 1.2 -o gen_003/candidate_CB15/seed.sfa && \
     $MODIFY gen_003/candidate_CB15/seed.sfa --add-braid 5 0 0 0 --flip-chirality -o gen_003/candidate_CB15/seed2.sfa && \
     mv gen_003/candidate_CB15/seed2.sfa gen_003/candidate_CB15/seed.sfa"

# === C7 evolution: high ellipticity at larger grid with amplitude ===

run_candidate "HE_A12" "High ellip (0.5) A=1.2 (larger grid)" \
    "/home/d/code/scp/sfa/seed/gen_braid -N $N -L $L -A 1.2 -ellip 0.5 -o gen_003/candidate_HE_A12/seed.sfa"

# === Hybrid: UDD + counter-chirality element ===

run_candidate "UDD_A05_e04" "UDD A=0.5 ellip=0.4" \
    "$GEN3 -N $N -L $L -A 0.5 -ellip 0.4 -chirality UDD -o gen_003/candidate_UDD_A05_e04/seed.sfa"

# === Summary ===
echo "=== Gen 3 Summary ==="
printf "%-15s %-45s %10s %10s\n" "ID" "Description" "S_final" "E_pot"
printf "%-15s %-45s %10s %10s\n" "---" "---" "---" "---"
for dir in gen_003/candidate_*/; do
    id=$(basename "$dir" | sed 's/candidate_//')
    desc=$(head -1 "$dir/config.cfg" | sed 's/^# Gen 3, [^:]*: //')
    if [ -f "$dir/analysis.json" ]; then
        score=$(grep -o '"S_final": [0-9.]*' "$dir/analysis.json" | grep -o '[0-9.]*')
        epot=$(grep -o '"E_pot_final": [-0-9.]*' "$dir/analysis.json" | grep -o '[-0-9.]*')
        printf "%-15s %-45s %10s %10s\n" "$id" "$desc" "$score" "$epot"
    else
        printf "%-15s %-45s %10s %10s\n" "$id" "$desc" "FAILED" "FAILED"
    fi
done
