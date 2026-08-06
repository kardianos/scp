# COMBINE door analyzer — tri2-aware door.awk (6 object voices in two
# triangles). Adds ILAG: inter-triangle door lag (FD on triangle k vs
# most recent DF on the OTHER triangle) — the event-grain meter for
# cross-object exchange through the medium.
# Usage: gawk -v i0=2700 -v ntri=2 -v tmin=1000 -f door2.awk run.log
# Output: VOICE k v ... | LAG k v ... | XLAG (same-tri cross-voice) |
#         ILAG (cross-tri) | BATH
BEGIN { cap = 2.5; tmax = 0; if (!ntri) ntri = 1 }
/^t=/ {
    tt = $0; sub(/^t= */, "", tt); t = tt + 0
    if (t > tmax) tmax = t
    s = $0; nt = 0
    while (match(s, /xUDD=\(([0-9.+-]+),([0-9.+-]+),([0-9.+-]+)\)/, mm)) {
        if (nt < ntri) { for (v = 0; v < 3; v++) xv[nt, v] = mm[v+1] + 0 }
        nt++
        s = substr(s, RSTART + RLENGTH)
    }
    if (nt >= ntri && t >= tmin) {
        for (k = 0; k < ntri; k++) for (v = 0; v < 3; v++) xs[k, v] += xv[k, v]
        nx++
    }
    next
}
/INIT NC=/ { if (match($0, /INIT NC=([0-9]+)/, mm)) nc = mm[1] + 0; next }
/^# QATOM/ {
    split($3, a, "="); t = a[2] + 0
    split($4, a, "="); dir = a[2]
    split($6, a, "="); e = a[2] + 0
    split($7, a, "="); ci = a[2] + 0
    if (t < tmin) next
    if (ci >= i0 && ci < i0 + 3 * ntri) {
        k = int((ci - i0) / 3); v = (ci - i0) % 3
        if (dir == "DF") {
            ndf[k, v]++; edf[k, v] += e; lastDF[k, v] = t
            if (!((k, "any") in lastDFk) || t > lastDFk[k, "any"]) lastDFk[k, "any"] = t
        } else {
            nfd[k, v]++; efd[k, v] += e
            if ((k, v) in lastDF) lag[k, v, ++nlag[k, v]] = t - lastDF[k, v]
            best = -1
            for (u = 0; u < 3; u++) if (u != v && ((k, u) in lastDF) && lastDF[k, u] > best) best = lastDF[k, u]
            if (best >= 0) xlag[k, ++nxlag[k]] = t - best
            ko = 1 - k
            if (ntri == 2 && ((ko, "any") in lastDFk)) ilag[k, ++nilag[k]] = t - lastDFk[ko, "any"]
        }
    } else nbath++
    next
}
function q(arr, n, f,   j) { j = int(f * n); if (j < 1) j = 1; return arr[j] }
END {
    W = tmax - tmin
    if (W <= 0 || nx == 0) { print "ERROR: empty window"; exit 1 }
    for (k = 0; k < ntri; k++) for (v = 0; v < 3; v++) {
        xa = xs[k, v] / nx; emfl = xa * cap
        gross = (edf[k, v] + efd[k, v]) / W; net = (efd[k, v] - edf[k, v]) / W
        tau = gross > 0 ? emfl / gross : -1
        printf "VOICE T%d v%d nDF=%d nFD=%d rateDF=%.4f rateFD=%.4f eDFr=%.5f eFDr=%.5f gross=%.5f net=%+.5f x_avg=%.4f Emfl=%.4f tau_turn=%.1f\n", \
            k, v, ndf[k, v], nfd[k, v], ndf[k, v]/W, nfd[k, v]/W, edf[k, v]/W, efd[k, v]/W, gross, net, xa, emfl, tau
    }
    for (k = 0; k < ntri; k++) for (v = 0; v < 3; v++) {
        n = nlag[k, v]
        if (n > 0) {
            delete s2; m = 0
            for (j = 1; j <= n; j++) { s2[j] = lag[k, v, j]; m += lag[k, v, j] }
            asort(s2)
            printf "LAG  T%d v%d n=%d median=%.2f p25=%.2f p75=%.2f mean=%.2f\n", \
                k, v, n, q(s2, n, 0.5), q(s2, n, 0.25), q(s2, n, 0.75), m / n
        }
    }
    for (k = 0; k < ntri; k++) {
        n = nxlag[k]
        if (n > 0) {
            delete s2
            for (j = 1; j <= n; j++) s2[j] = xlag[k, j]
            asort(s2)
            printf "XLAG T%d n=%d median=%.2f\n", k, n, q(s2, n, 0.5)
        }
        n = nilag[k]
        if (n > 0) {
            delete s2
            for (j = 1; j <= n; j++) s2[j] = ilag[k, j]
            asort(s2)
            printf "ILAG T%d n=%d median=%.2f p25=%.2f p75=%.2f\n", k, n, q(s2, n, 0.5), q(s2, n, 0.25), q(s2, n, 0.75)
        }
    }
    if (nc > 3 * ntri) printf "BATH nev=%d rate_per_cell=%.5f\n", nbath, nbath / (W * (nc - 3 * ntri))
}
