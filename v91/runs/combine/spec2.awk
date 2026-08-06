# COMBINE spectrometer — tri2-aware spec.awk (two xUDD tuples per diag
# row, cells i0..i0+3*ntri-1; UUD roles: kk=0,1 U, kk=2 D per triangle).
# Usage: gawk -v i0=2700 -v ntri=2 -v tmin=200 -f spec2.awk run.log
# Output:
#   EV  t k ci dir w_meas w_pred resid     (k = triangle index)
#   HB  dir w_bin count                    (bath histogram, 0.02 bins)
#   SUM obj T<k> <dir> n resid_mean resid_rms
#   SUM pooled <dir> n resid_mean resid_rms
#   SUM bath <dir> n
BEGIN { if (!ntri) ntri = 1 }
/^t=/ {
    s = $0; nt = 0
    while (match(s, /xUDD=\(([0-9.+-]+),([0-9.+-]+),([0-9.+-]+)\)/, mm)) {
        if (nt < ntri) { for (v = 0; v < 3; v++) xt[nt, v] = mm[v+1] + 0 }
        nt++
        s = substr(s, RSTART + RLENGTH)
    }
    if (nt >= ntri) have = 1
    next
}
/^# QATOM/ {
    split($3, a, "="); t = a[2] + 0
    split($4, a, "="); dir = a[2]
    split($5, a, "="); w = a[2] + 0
    split($7, a, "="); ci = a[2] + 0
    if (t < tmin) next
    if (i0 > 0 && ci >= i0 && ci < i0 + 3 * ntri && have) {
        k = int((ci - i0) / 3); kk = (ci - i0) % 3
        x = xt[k, kk]
        base = (dir == "DF") ? 2.9 : 1.65
        wp = base / (1 + 1.2 * x)
        printf "EV %.2f %d %d %s %.6f %.6f %+.6f\n", t, k, ci, dir, w, wp, w - wp
        key = k "_" dir
        nev[key]++; sr[key] += (w - wp); sr2[key] += (w - wp) * (w - wp)
        np[dir]++; pr[dir] += (w - wp); pr2[dir] += (w - wp) * (w - wp)
    } else {
        b = int(w / 0.02)
        hb[dir "_" b]++
        nb_[dir]++
    }
    next
}
END {
    for (k in hb) {
        split(k, a, "_")
        printf "HB %s %.3f %d\n", a[1], (a[2] + 0.5) * 0.02, hb[k]
    }
    for (key in nev) {
        split(key, a, "_")
        m = sr[key] / nev[key]; v = sr2[key] / nev[key] - m * m; if (v < 0) v = 0
        printf "SUM obj T%s %s n=%d resid_mean=%+.6f resid_rms=%.6f\n", \
               a[1], a[2], nev[key], m, sqrt(v)
    }
    for (d in np) {
        m = pr[d] / np[d]; v = pr2[d] / np[d] - m * m; if (v < 0) v = 0
        printf "SUM pooled %s n=%d resid_mean=%+.6f resid_rms=%.6f\n", d, np[d], m, sqrt(v)
    }
    for (d in nb_) printf "SUM bath %s n=%d\n", d, nb_[d]
}
