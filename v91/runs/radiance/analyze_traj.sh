#!/bin/bash
# Trajectory extractor: per-run time series from diag + CONVTAG lines.
# Usage: analyze_traj.sh <log> -> rows: t Em_tag ct_cond ct_evap
# (ct_* cumulative; window rates = differences between rows)
awk '
  /^t=/ {
    sub("^t=", "", $1); t=$2==""?$1:$1;   # "t=  4.00" splits as "t=" "4.00"? robust: strip prefix
    n=split($0,seg,"|"); split(seg[5],a," "); em=a[1];
    tt=seg[1]; gsub(/[t= ]/,"",tt); sub(/\+.*/,"",tt);
    pend_t=tt; pend_em=em;
  }
  /^# CONVTAG/ {
    for(i=1;i<=NF;i++){ if($i~"^cond=")c=substr($i,6); if($i~"^evap=")e=substr($i,6); if($i~"^t=")ct=substr($i,3) }
    printf "%s %s %s %s\n", ct, pend_em, c, e;
  }
' "$1"
