#!/bin/bash
# speedup-ladder benchmark: 3 repeats each box, median wall seconds.
# usage: ./bench.sh <label> [binary]   (binary default ../cellfab)
set -e
cd "$(dirname "$0")"
BIN=${2:-../cellfab}
LAWS=../battery/laws_V2g.cfg
for box in bench_blob bench_noise; do
    cat $LAWS $box.cfg > /tmp/bench_$box.cfg
    ts=()
    for r in 1 2 3; do
        s=$(date +%s.%N)
        $BIN /tmp/bench_$box.cfg > /tmp/bench_$box.log 2>&1
        e=$(date +%s.%N)
        ts+=($(echo "$e - $s" | bc))
    done
    med=$(printf '%s\n' "${ts[@]}" | sort -n | sed -n 2p)
    steps=$(awk -F= '/^T=/{t=$2} /^dt=/{d=$2} END{print t/d}' $box.cfg)
    echo "$1 $box wall_s=$med ms_per_step=$(echo "scale=3; 1000*$med/$steps" | bc)"
done
