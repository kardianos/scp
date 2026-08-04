#!/bin/bash
# R1b flux-curve: pool CONVTAG windows across runs/seeds per arm,
# per-cell flux F(x) = d(ct_cond - ct_evap)/dt / 2 binned by per-cell
# x = (Em_tag/2)/cap at window start. Output: arm xbin Fmean n
cd "$(dirname "$0")"
for arm in k000 k005 k010; do
  for f in r1b_${arm}_s*_x*.log; do
    [ -f "$f" ] || continue
    grep -q '# RESULT drift_rel' "$f" || continue
    ./analyze_traj.sh "$f"
  done | awk -v arm="$arm" '
    { t=$1+0; em=$2+0; c=$3+0; e=$4+0;
      if (t > pt && pt >= 0) {
        dt=t-pt; F=((c-pc)-(e-pe))/2.0/dt;   # per-cell
        x=(pem/2.0)/2.5;
        b=int(x/0.02)*0.02;
        S[b]+=F; N[b]++;
      }
      pt=t; pem=em; pc=c; pe=e;
    }
    $1=="0.00" { pt=0; pem=$2+0; pc=$3+0; pe=$4+0 }   # reset at each new run
    END { for (b in S) printf "%s %.2f %+.6f %d\n", arm, b, S[b]/N[b], N[b] }
  ' | sort -k2,2g
done
