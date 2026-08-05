#!/bin/bash
# score_v11.sh — score the v1.1 arms against the CANTUS.md pre-registered
# bars. Refuses incomplete logs. Controls: committed R2/pad logs (k0
# arms unchanged by the byte-inert candidate) + in-campaign k_cant=0.
cd "$(dirname "$0")"
pairx() { # file -> final tagged-pair x + Em sum + gate-alive frac (last half)
  local f=$1
  grep -q "RESULT drift_rel" "$f" || { echo "INCOMPLETE"; return; }
  local em=$(grep "RESULT pair d" "$f" | sed 's/.*Em_i=\([0-9.]*\) Em_j=\([0-9.]*\).*/\1 \2/')
  [ -z "$em" ] && { echo "PAIR-DEAD"; return; }
  local ems=$(echo $em | awk '{printf "%.3f", $1+$2}')
  local gfrac=$(grep "^t=" "$f" | grep -o "gf=[0-9.]* gb=[0-9.]*" | awk '{split($1,a,"="); split($2,b,"="); print a[2]*b[2]}' | awk '{v[NR]=$1} END {n=NR; if(n==0){print "-"; exit} a=0; h=int(n/2); for(i=h+1;i<=n;i++) if(v[i]>0.5) a++; printf "%.2f", a/(n-h)}')
  local xf=$(grep "^t=" "$f" | tail -1 | grep -o "x=[0-9.]*" | tail -1)
  local ca=$(grep "RESULT cantus" "$f" | sed 's/.*a_tag=\([0-9.]*\).*tune_total=\([0-9.]*\).*/a=\1 tune=\2/')
  local dr=$(grep "RESULT drift_rel" "$f" | awk '{print $4}')
  echo "Em_sum=$ems ${xf:-x=?} gg50=$gfrac ${ca:-} drift=$dr"
}
echo "=== H1a (bath pair 0.47; seeded Em_sum=2.35; v1 k0 control: Em_sum=0.308, x=0.12)"
for f in h1a_kc0_kt0.log h1a_kc0_kt0p2.log h1a_kc0p3_kt0.log h1a_kc0p3_kt0p2.log h1a_kc1_kt0.log h1a_kc1_kt0p2.log h1a_kc3_kt0.log h1a_kc3_kt0p2.log; do
  [ -f "$f" ] && echo "$f :: $(pairx $f)"
done
echo
echo "=== H1b (fifth pair vacuum; chart-hold = d_star_live stays ~1.45 (3:2), not ~3.5 (1:1))"
for f in h1b_kc0_kt0.log h1b_kc0_kt0p2.log h1b_kc0p3_kt0.log h1b_kc0p3_kt0p2.log h1b_kc1_kt0.log h1b_kc1_kt0p2.log h1b_kc3_kt0.log h1b_kc3_kt0p2.log h1b_seed1.log; do
  [ -f "$f" ] || continue
  dsl=$(grep "RESULT pair d" "$f" | sed 's/.*d_star_live=\([0-9.]*\).*/\1/')
  echo "$f :: d_star_live=${dsl:-DEAD} $(pairx $f)"
done
echo
echo "=== H2 objects (ret + edges; controls: comp6 t_half 85-90, uud 140, i5 230, ring6x62bath 140)"
for f in h2a_sel.log h2a_low.log h2a_seed1.log h2a_ctl.log h2b_sel.log h2c_sel.log h2d_sel.log h2e_sel.log; do
  [ -f "$f" ] || continue
  grep -q "RESULT drift_rel" "$f" || { echo "$f :: INCOMPLETE ($(grep -c '^t=' $f) rows)"; continue; }
  ret=$(grep "^t=" "$f" | tail -1 | awk '{n=split($0,a,"|"); split(a[5],b," "); print b[2]}')
  thalf=$(grep "^t=" "$f" | awk '{n=split($0,a,"|"); split(a[5],b," "); split(a[1],c," "); if (b[2]<0.5 && !done) {sub(/t=/,"",c[1]); print c[1]; done=1}}' | head -1)
  edges=$(grep "# EDGE" "$f" | grep -o "gg=[0-9.]*" | awk -F= '$2>0.5 {n++} END {print n+0}')
  etot=$(grep -c "# EDGE" "$f")
  dead=$(grep -c "CHANNEL DEAD" "$f")
  cvt=$(grep "RESULT convtag" "$f" | sed 's/.*cond=\([0-9.]*\) evap=\([0-9.]*\).*/in=\1 out=\2/')
  ca=$(grep "RESULT cantus" "$f" | sed 's/.*a_tag=\([0-9.]*\).*nlock=\([0-9]*\).*tune_total=\([0-9.]*\).*/a=\1 nlock=\2 tune=\3/')
  dr=$(grep "RESULT drift_rel" "$f" | awk '{print $4}')
  echo "$f :: ret_end=$ret t_half=${thalf:-none} edges_gg50=$edges/$etot(dead$dead) $cvt $ca drift=$dr"
done
echo
echo "=== H3 (chord at balance; hold = pitch ratio 3:2 => d_star_live ~1.5 not ~3.6)"
for f in h3_sel.log h3_seed1.log h3_ctl.log; do
  [ -f "$f" ] || continue
  dsl=$(grep "RESULT pair d" "$f" | sed 's/.*d_star_live=\([0-9.]*\).*/\1/')
  echo "$f :: d_star_live=${dsl:-DEAD} $(pairx $f)"
done
echo
echo "=== H4 (P-H8': nlock <135 (5% NC); cond/rad within ±15% of ctl 2204/2016)"
for f in h4_sel.log h4_ctl.log h4_sel_v1.log; do
  [ -f "$f" ] || continue
  cv=$(grep "RESULT conv " "$f" | sed 's/.*cond=\([0-9.]*\) evap=[0-9.]* backs=\([0-9.]*\) rad=\([0-9.]*\).*/cond=\1 backs=\2 rad=\3/')
  ca=$(grep "RESULT cantus" "$f" | sed 's/.*nlock=\([0-9]*\).*tune_total=\([0-9.]*\).*/nlock=\1 tune=\2/')
  echo "$f :: $cv ${ca:-k0}"
done
