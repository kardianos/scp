#!/bin/bash
# R3 zoo: the finite (interval, x*) particle table the comb admits.
# Usage: zoo_r3.sh <xstar>  ->  for each comb interval p:q (p*q<=6),
# the m=1..3 bond rungs d* = 2*pi*m*C/(q*w + p*w') for two voices at
# x* (w = w2/(1+q_detune*x*), unison partner) — the predicted sizes
# of radiance-balanced objects.
X=${1:?usage: zoo_r3.sh xstar}
awk -v x="$X" 'BEGIN{
  C=1; w2=2.9; qd=1.2;
  det=1+qd*x; w=w2/det;
  printf "x*=%.3f  w(x*)=%.5f  det=%.3f  (V2g: w2=2.9 qd=1.2)\n", x, w, det;
  printf "%-6s %-10s %-10s %-10s\n","p:q","d*(m=1)","d*(m=2)","d*(m=3)";
  n=split("1:1 2:1 3:1 4:1 5:1 6:1 3:2", I, " ");
  for(i=1;i<=n;i++){
    split(I[i],pq,":"); p=pq[1]; q=pq[2];
    s=q*w+p*w;   # both voices at x*: w_i=w_j=w
    printf "%-6s", I[i];
    for(m=1;m<=3;m++) printf " %-10.4f", 2*3.14159265358979*m*C/s;
    printf "\n";
  }
  printf "unison m=1 = pi*C/w = %.4f (the R3 headline rung)\n", 3.14159265358979*C/w;
}'
