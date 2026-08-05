# ASTRO species spectrometer over a freecell LOG (QATOM + diag rows).
# Usage: awk -v i0=2700 -v tmin=200 -f spec.awk run.log
#   i0     = first object cell id (tri: NC-3; cells i0,i0+1 are U, i0+2 is D)
#   tmin   = ignore events before this t (settle)
# Reads the diag rows' xUDD=(xU1,xU2,xD) column as the per-cell load
# meter; every object QATOM event gets a parameter-free prediction
#   DF (emission):  w_pred = 2.9  / (1 + 1.2*x_cell(t))
#   FD (absorption): w_pred = 1.65 / (1 + 1.2*x_cell(t))
# Output sections:
#   EV  t i dir w_meas w_pred resid        (object events)
#   HB  dir w_bin count                    (bath histogram, 0.02 bins)
#   SUM per-class counts + residual stats
/^t=/ {
    # ... xUDD=(a,b,c) at line end (tri diag rows)
    if (match($0, /xUDD=\(([0-9.+-]+),([0-9.+-]+),([0-9.+-]+)\)/, mm)) {
        xu1 = mm[1]+0; xu2 = mm[2]+0; xd = mm[3]+0; have = 1
    }
    next
}
/^# QATOM/ {
    # # QATOM t=%.2f dir=%s w=%g e=%g i=%d Em=%.4f
    split($3, a, "="); t = a[2]+0
    split($4, a, "="); dir = a[2]
    split($5, a, "="); w = a[2]+0
    split($7, a, "="); ci = a[2]+0
    if (t < tmin) next
    if (i0 > 0 && ci >= i0 && ci <= i0+2 && have) {
        x = (ci == i0) ? xu1 : (ci == i0+1) ? xu2 : xd
        base = (dir == "DF") ? 2.9 : 1.65
        wp = base / (1 + 1.2*x)
        printf "EV %.2f %d %s %.6f %.6f %+.6f\n", t, ci, dir, w, wp, w-wp
        nev[dir]++; sr[dir] += (w-wp); sr2[dir] += (w-wp)*(w-wp)
    } else {
        b = int(w/0.02)
        hb[dir "_" b]++
        nb_[dir]++
    }
    next
}
END {
    for (k in hb) {
        split(k, a, "_")
        printf "HB %s %.3f %d\n", a[1], (a[2]+0.5)*0.02, hb[k]
    }
    for (d in nev) {
        m = sr[d]/nev[d]; v = sr2[d]/nev[d] - m*m; if (v<0) v=0
        printf "SUM obj %s n=%d resid_mean=%+.6f resid_rms=%.6f\n", \
               d, nev[d], m, sqrt(v)
    }
    for (d in nb_) printf "SUM bath %s n=%d\n", d, nb_[d]
}
