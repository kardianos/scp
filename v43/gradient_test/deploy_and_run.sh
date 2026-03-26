#!/bin/bash
# Deploy and run gradient force test on V100
# Usage: ./deploy_and_run.sh <ssh_port> <ssh_host>
set -e

PORT=$1
HOST=$2
SSH="ssh -o StrictHostKeyChecking=no -p $PORT root@$HOST"
SCP="scp -o StrictHostKeyChecking=no -P $PORT"

echo "=== Uploading code and template ==="
$SCP ../../sfa/sim/scp_sim.cu root@$HOST:/root/
$SCP ../../sfa/format/sfa.h root@$HOST:/root/
$SCP ../proton_formation/proton_averaged.sfa root@$HOST:/root/

echo "=== Building on remote ==="
$SSH 'cd /root && cp sfa.h sfa_format_sfa.h && \
    sed -i "s|../format/sfa.h|sfa_format_sfa.h|" scp_sim.cu && \
    apt-get install -y -qq libzstd-dev 2>/dev/null && \
    nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm 2>&1 | grep -c error && \
    echo "Build errors above" || echo "Build OK"'

echo "=== Writing configs and running ==="
$SSH 'cd /root

# Write all 8 configs using template init
for obj in proton braid; do
  for grad in gentle steep; do
    for eta_val in 0.5 0.0; do
      suffix="${obj}_${grad}"
      if [ "$eta_val" = "0.0" ]; then suffix="${suffix}_eta0"; fi
      if [ "$grad" = "gentle" ]; then AH=0.12; AL=0.08; else AH=0.15; AL=0.05; fi

      # For braid: use built-in braid init, not template
      if [ "$obj" = "braid" ]; then
        INIT_MODE="braid"
        INIT_SFA=""
      else
        INIT_MODE="template"
        INIT_SFA="init_sfa = proton_averaged.sfa"
      fi

      cat > ${suffix}.cfg << EOCFG
N = 384
L = 100
T = 200
m = 1.5
eta = ${eta_val}
mu = -41.345
kappa = 50
init = ${INIT_MODE}
${INIT_SFA}
A = 0.8
A_bg = 0.1
R_tube = 3.0
bc_type = 1
gradient_A_high = ${AH}
gradient_A_low = ${AL}
gradient_margin = 5
snap_dt = 20
diag_dt = 50
output = ${suffix}_output.sfa
diag_file = ${suffix}_diag.tsv
precision = f16
EOCFG
      echo "Config: ${suffix}.cfg (init=$INIT_MODE eta=${eta_val} AH=${AH})"
    done
  done
done

echo ""
echo "=== Running all 8 experiments ==="
for cfg in proton_gentle proton_steep braid_gentle braid_steep \
           proton_gentle_eta0 proton_steep_eta0 braid_gentle_eta0 braid_steep_eta0; do
    echo ""
    echo "=== $cfg ==="
    ./scp_sim_cuda ${cfg}.cfg 2>&1 | grep -E "INIT:|COMPLETE|Wall:|drift|theta_rms"
done

echo ""
echo "=== ALL 8 COMPLETE ==="
ls -lh *_diag.tsv *_output.sfa 2>/dev/null'
