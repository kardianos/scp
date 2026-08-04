#!/bin/bash
# R1 extractor: one TSV row per run.
# cols: arm k_rad p_rad rad_clock x pair in(ct_cond) out_evap out_rough backs
#       rad_global evap_global Em_final Em0 drift
cd "$(dirname "$0")"
printf "arm\tk\tp\tc\tx\tpair\tin\tout_evap\tout_rough\tbacks\trad_glob\tevap_glob\tcond_glob\tEm_fin\tEm0\tdrift\n"
for f in r1_*.log; do
  grep -q '# RESULT drift_rel' "$f" || continue   # never harvest incomplete logs
  arm=$(basename "$f" .log | sed 's/^r1_//; s/_x[0-9]*$//')
  awk -v arm="$arm" '
    /# v91 radiance/   { for(i=1;i<=NF;i++){ if($i~"^k_rad=")k=substr($i,7); if($i~"^p_rad=")p=substr($i,7); if($i~"^rad_clock=")c=substr($i,11)} }
    /^# SEED pair:/    { for(i=1;i<=NF;i++){ if($i~"^i=")pi=substr($i,3); if($i~"^j=")pj=substr($i,3); if($i~"^x=")x=substr($i,3) } sub("/.*","",x) }
    /^# SEED pair radii/ { em0=$NF; sub("Em/voice=","",em0) }
    /^t=/              { n=split($0,seg,"|"); split(seg[5],a," "); em=a[1] }
    /^# RESULT convtag/ { for(i=1;i<=NF;i++){ if($i~"^cond=")ci=substr($i,6); if($i~"^evap=")ce=substr($i,6); if($i~"^rough=")cr=substr($i,7); if($i~"^backs=")cb=substr($i,7) } }
    /^# RESULT conv /  { for(i=1;i<=NF;i++){ if($i~"^rad=")rg=substr($i,5); if($i~"^evap=")eg=substr($i,6); if($i~"^cond=")cg=substr($i,6) } }
    /^# RESULT drift_rel/ { dr=$4 }
    END { printf "%s\t%s\t%s\t%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
          arm,k,p,c,x,pi,pj,ci,ce,cr,cb,rg,eg,cg,em,em0,dr }
  ' "$f"
done | sort -t$'\t' -k1,1 -k5,5g
