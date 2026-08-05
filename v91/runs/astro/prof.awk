# ASTRO radial Es/Em profiler over `fcsdump -mode cells` output.
# Usage: fcsdump -mode cells X.fcs | awk -v L=16 -v tmin=500 -v tmax=99999 \
#          -v dr=0.5 -v rmax=8 -v cx=-1 -v cy=-1 -v cz=-1 -f prof.awk
# Center: tag-cell centroid per frame (default), or fixed (cx,cy,cz) if given.
# Output: r_mid  n_avg  Es_mean  Es_sem  Em_mean  frames
# Es_sem is over FRAME means (frames treated as independent samples; the
# 10 t.u. cadence exceeds the bath churn time — stated in ASTRO.md §1).
BEGIN { fixed = (length(cx) > 0); nb = int(rmax/dr) + 1 }
NR == 1 { next }
{
    t = $1 + 0
    if (t < tmin || t > tmax) next
    if (t != curt) { flushframe(); curt = t; nrow = 0 }
    nrow++
    X[nrow]=$3; Y[nrow]=$4; Z[nrow]=$5; ES[nrow]=$7; EM[nrow]=$8; EE[nrow]=$9
    TG[nrow]=$11
}
END { flushframe(); emit() }
function wr(d) { while (d >  L/2) d -= L; while (d < -L/2) d += L; return d }
function flushframe(   i,k,tx,ty,tz,tn,dx,dy,dz,r,b,fs,fn,fm,fp,m,p) {
    if (nrow == 0) return
    if (fixed) { tx=cx; ty=cy; tz=cz }
    else {
        # centroid of tagged cells relative to the first tagged cell
        # (wrap-safe for compact objects)
        tn=0
        for (i=1; i<=nrow; i++) if (TG[i]) {
            if (tn==0) { rx=X[i]; ry=Y[i]; rz=Z[i]; tx=0; ty=0; tz=0 }
            tx += wr(X[i]-rx); ty += wr(Y[i]-ry); tz += wr(Z[i]-rz); tn++
        }
        if (tn==0) return
        tx = rx + tx/tn; ty = ry + ty/tn; tz = rz + tz/tn
    }
    for (b=0; b<nb; b++) { fs[b]=0; fn[b]=0; fm[b]=0 }
    for (i=1; i<=nrow; i++) {
        if (TG[i]) continue                # object cells excluded: medium only
        dx=wr(X[i]-tx); dy=wr(Y[i]-ty); dz=wr(Z[i]-tz)
        if (ringR > 0) {                    # torus distance to the ring hoop
            rxy = sqrt(dx*dx+dy*dy)
            r = sqrt((rxy-ringR)*(rxy-ringR) + dz*dz)
        } else r = sqrt(dx*dx+dy*dy+dz*dz)
        b=int(r/dr); if (b>=nb) continue
        fs[b]+=ES[i]; fn[b]++; fm[b]+=EM[i]
        fp[b]+=ES[i] + 0.3*(EM[i]+EE[i])   # pi: pass-S pressure, s_disp=0.3
    }
    for (b=0; b<nb; b++) if (fn[b]>0) {
        m = fs[b]/fn[b]; p = fp[b]/fn[b]
        S1[b]+=m; S2[b]+=m*m; SN[b]++; SC[b]+=fn[b]; SM[b]+=fm[b]/fn[b]
        SP[b]+=p; SP2[b]+=p*p
    }
    nframes++
}
function emit(   b,m,v) {
    printf "# frames=%d window=[%g,%g] center=%s\n", nframes, tmin, tmax, \
           fixed ? "fixed" : "tag-centroid"
    print "# r_mid n_avg Es_mean Es_sem Em_mean pi_mean pi_sem frames_bin"
    for (b=0; b<nb; b++) {
        if (SN[b] < 2) continue
        m = S1[b]/SN[b]; v = S2[b]/SN[b] - m*m; if (v<0) v=0
        p = SP[b]/SN[b]; pv = SP2[b]/SN[b] - p*p; if (pv<0) pv=0
        printf "%.3f %.1f %.6f %.6f %.6f %.6f %.6f %d\n", (b+0.5)*dr, \
               SC[b]/SN[b], m, sqrt(v/SN[b]), SM[b]/SN[b], p, \
               sqrt(pv/SN[b]), SN[b]
    }
}
