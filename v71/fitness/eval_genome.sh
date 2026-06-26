#!/usr/bin/env bash
# eval_genome.sh — evaluate ONE shape genome -> one scalar fitness line.
# Inverted methodology (FUTURE.md): viability-first, then QFI.
#   J = viability * (1 + WQ*min(nQFI,2)) - drift_penalty
#   viability = S_mean * alive_frac     (shape retention under kernel evolution)
# Usage: eval_genome.sh <genome> <id> <workdir> <bin_dir> <sim_bin>
#   sim_bin = path to scp_sim (CPU) or scp_sim_cuda (GPU)
set -u
GEN=$1; ID=$2; W=$3; BIN=$4; SIM=$5
N=${N:-48}; L=${L:-12.0}; T=${T:-20.0}; SNAP=${SNAP:-0.5}; ETA=${ETA:-0.25}; WQ=${WQ:-0.5}
mkdir -p $W
seed=$W/${ID}.sfa; run=$W/${ID}_run.sfa; cfg=$W/${ID}.cfg

$BIN/gen_blob_field $N $L $GEN $seed > $W/${ID}_gen.log 2>&1 || { echo "FITNESS -1e9 err gen"; exit 0; }
cat > $cfg <<EOF
N = $N
L = $L
T = $T
dt_factor = 0.025
m = 1.5
m_theta = 1.6
eta = $ETA
mu = -41.345
kappa = 50.0
mode = 0
complex_phi = 1
complex_gauge = 1
g_gauge = 0.05
init = sfa
init_sfa = $seed
bc_type = 0
damp_width = 2.5
damp_rate = 0.01
output = $run
precision = f32
snap_dt = $SNAP
diag_dt = 5.0
diag_file = $W/${ID}_diag.tsv
EOF
$SIM $cfg > $W/${ID}_run.log 2>&1 || { echo "FITNESS -1e9 err sim"; exit 0; }

an=$($BIN/analyze_sfa $run 2>/dev/null)
smean=$(echo "$an" | awk -F'=' '/S_mean/{gsub(/ /,"",$2);print $2+0}')
sfinal=$(echo "$an" | awk -F'=' '/S_final/{gsub(/ /,"",$2);print $2+0}')
alivef=$(echo "$an" | awk -F: '/Alive frames/{split($2,p,"/"); print (p[2]>0)?(p[1]/p[2]):0}')
drift=$(grep -oE 'drift -?[0-9.]+%' $W/${ID}_run.log | tail -1 | grep -oE '[0-9.]+')
nq=$($BIN/sfa_qfi $run --auto-T --q 0,0,0 --q 1,0,0 --q 2,0,0 2>/dev/null | awk '$1=="rhoQ"{if($7>m)m=$7} END{printf "%.6e", m+0}')

smean=${smean:-0}; sfinal=${sfinal:-0}; alivef=${alivef:-0}; drift=${drift:-0}; nq=${nq:-0}
# stability is a HARD multiplicative gate: a true particle conserves energy.
# g_stab = exp(-(drift/DTOL)^2) kills exploding/radiating configs (drift in %).
# viability = sqrt(S_mean * S_final) * alive_frac: rewards PERSISTENCE (the core
# must still be present at the END of the run, not just on average) -> "larger time".
DTOL=${DTOL:-1.0}
v=$(awk -v sm=$smean -v sf=$sfinal -v a=$alivef 'BEGIN{print sqrt((sm<0?0:sm)*(sf<0?0:sf))*a}')
J=$(awk -v v=$v -v d=$drift -v q=$nq -v wq=$WQ -v dt=$DTOL 'BEGIN{
  qb=(q<2)?q:2;
  gstab=exp(-(d/dt)*(d/dt));
  print v*(1+wq*qb)*gstab }')
echo "FITNESS $J viability $v smean $smean sfinal $sfinal alive $alivef drift $drift nqfi $nq"
# clean up per-eval artifacts (keep only the printed scalar) to bound disk use
rm -f "$seed" "$run" "$cfg" "$W/${ID}_diag.tsv" "$W/${ID}_run.log" "$W/${ID}_gen.log"
