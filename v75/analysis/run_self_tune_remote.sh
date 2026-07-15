#!/bin/bash
set -e
cd /root/st_bundle
export PYTHONUNBUFFERED=1
rm -rf /root/self_tune
mkdir -p /root/self_tune
: > /root/self_tune_run.log
exec python3 -u self_tune_controller.py \
  --mode gpu \
  --sim /root/scp_sim_9ce10ec2f9072352 \
  --gen-multi /root/st_bundle/bin/gen_qball_multi \
  --track /root/st_bundle/bin/mf_pair_track \
  --profiles /root/st_bundle/profiles \
  --work /root/self_tune \
  --N 192 --L 48 \
  --T_screen 150 --T_full 400 \
  >> /root/self_tune_run.log 2>&1
