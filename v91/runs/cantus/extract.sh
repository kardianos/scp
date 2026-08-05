#!/bin/bash
# extract.sh — harvest the H-campaign logs into one table per arm class.
# Refuses logs without RESULT lines (the k030 stale-extractor lesson).
cd "$(dirname "$0")"
echo "== H1a (bath pair x=0.47, squeezed): arm | gate-alive-frac(last half) | d_fin | off_live | Em_i/Em_j | a_tag | tune | drift"
for f in h1a_*.log; do
  [ -f "$f" ] || continue
  grep -q "RESULT drift_rel" "$f" || { echo "$f INCOMPLETE"; continue; }
  gf=$(grep "^t=" "$f" | awk -F'gf=' 'NF>1 {split($2,a," "); split(a[2],b,"="); n++; if (NR>0) {split($2,g," "); gfv=g[1]+0; gbv=substr(g[2],4)+0}} END {print ""}')
  # gate stats: parse gf= and gb= from pair diag columns, last half rows
  stats=$(grep "^t=" "$f" | grep -o "gf=[0-9.]* gb=[0-9.]*" | awk '{split($1,a,"="); split($2,b,"="); print a[2]*b[2]}' | awk '{v[NR]=$1} END {n=NR; alive=0; for(i=int(n/2)+1;i<=n;i++) if(v[i]>0.5) alive++; printf "%.2f", n>0 ? alive/(n-int(n/2)) : -1}')
  dead=$(grep -c "PAIR-CHANNEL-DEAD" "$f")
  rp=$(grep "# RESULT pair d_final" "$f" | sed 's/.*d_final=\([0-9.]*\).*off_live=\([+-][0-9.]*\).*Em_i=\([0-9.]*\) Em_j=\([0-9.]*\).*/\1 \2 \3 \4/')
  [ -z "$rp" ] && rp="DEAD - - -"
  ct=$(grep "# RESULT cantus" "$f" | sed 's/.*a_tag=\([0-9.]*\).*tune_total=\([0-9.]*\).*/\1 \2/')
  [ -z "$ct" ] && ct="- -"
  dr=$(grep "# RESULT drift_rel" "$f" | awk '{print $3}')
  echo "$f | gg>0.5: $stats (deadrows $dead) | $rp | $ct | $dr"
done
echo
echo "== H1b (fifth pair vacuum): same columns"
for f in h1b_*.log; do
  [ -f "$f" ] || continue
  grep -q "RESULT drift_rel" "$f" || { echo "$f INCOMPLETE"; continue; }
  stats=$(grep "^t=" "$f" | grep -o "gf=[0-9.]* gb=[0-9.]*" | awk '{split($1,a,"="); split($2,b,"="); print a[2]*b[2]}' | awk '{v[NR]=$1} END {n=NR; alive=0; for(i=int(n/2)+1;i<=n;i++) if(v[i]>0.5) alive++; printf "%.2f", n>0 ? alive/(n-int(n/2)) : -1}')
  dead=$(grep -c "PAIR-CHANNEL-DEAD" "$f")
  rp=$(grep "# RESULT pair d_final" "$f" | sed 's/.*d_final=\([0-9.]*\).*off_live=\([+-][0-9.]*\).*Em_i=\([0-9.]*\) Em_j=\([0-9.]*\).*/\1 \2 \3 \4/')
  [ -z "$rp" ] && rp="DEAD - - -"
  ct=$(grep "# RESULT cantus" "$f" | sed 's/.*a_tag=\([0-9.]*\).*tune_total=\([0-9.]*\).*/\1 \2/')
  [ -z "$ct" ] && ct="- -"
  dr=$(grep "# RESULT drift_rel" "$f" | awk '{print $3}')
  echo "$f | gg>0.5: $stats (deadrows $dead) | $rp | $ct | $dr"
done
echo
echo "== H2/H3: arm | ret trajectory (t=1000/2500/5000 or end) | edges gg>0.5 at end | ledger in/out | cantus | drift"
for f in h2*.log h3*.log; do
  [ -f "$f" ] || continue
  grep -q "RESULT drift_rel" "$f" || { echo "$f INCOMPLETE"; continue; }
  ret=$(grep "^t=" "$f" | awk '{n=split($0,a,"|"); split(a[5],b," "); split(a[1],c," "); printf "%s:%s ", c[1], b[2]}' | awk '{for(i=1;i<=NF;i++) if (i==1 || i==int(NF/3) || i==int(2*NF/3) || i==NF) printf "%s ", $i; print ""}')
  edges=$(grep "# EDGE" "$f" | grep -o "gg=[0-9.]*" | awk -F= '$2>0.5 {n++} END {print n+0}')
  edgetot=$(grep -c "# EDGE" "$f")
  deadch=$(grep -c "CHANNEL DEAD" "$f")
  cvt=$(grep "# RESULT convtag" "$f" | sed 's/.*cond=\([0-9.]*\) evap=\([0-9.]*\).*/in=\1 out=\2/')
  ct=$(grep "# RESULT cantus" "$f" | sed 's/.*a_tag=\([0-9.]*\).*nlock=\([0-9]*\).*tune_total=\([0-9.]*\).*/a=\1 nl=\2 tune=\3/')
  dr=$(grep "# RESULT drift_rel" "$f" | awk '{print $3}')
  echo "$f | $ret| edges>0.5: $edges/$edgetot (dead $deadch) | $cvt | $ct | $dr"
done
echo
echo "== H4 (bath census): arm | RESULT conv | cantus census | drift"
for f in h4_*.log; do
  [ -f "$f" ] || continue
  grep -q "RESULT drift_rel" "$f" || { echo "$f INCOMPLETE"; continue; }
  cv=$(grep "# RESULT conv " "$f" | head -1)
  ct=$(grep "# RESULT cantus" "$f" | head -1)
  dr=$(grep "# RESULT drift_rel" "$f" | awk '{print $3}')
  echo "$f | $cv | ${ct:-no-cantus-line} | $dr"
done
