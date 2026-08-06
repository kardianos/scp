# D-A1 — the identity lane's design envelope, measured at the door
# from a standing full-rate QATOM raw (ASTRO §6, analysis-only).
# Usage: gawk -v i0=2700 -v tmin=1000 -f door.awk run.log.raw
#   i0   = first object cell id (tri: NC-3)
#   tmin = start of settled window (end = last diag t seen)
# Output:
#   VOICE v nDF nFD rateDF rateFD eDFr eFDr gross net x_avg Emfl tau_turn
#     rates per t.u.; eDFr/eFDr = energy out/in per t.u.; gross = eDFr+eFDr;
#     net = eFDr-eDFr; Emfl = x_avg*cap (held Em+flload, cap=2.5);
#     tau_turn = Emfl/gross = time to cycle the held energy through the door
#     (the lifetime a parcel gid must survive to be load-bearing there)
#   LAG  v n median p25 p75 mean  (same-voice DF->FD lag: each FD vs that
#     voice's most recent DF — the door re-eating its own shine, per-event)
#   XLAG v n median               (control: FD vs most recent DF on any
#     OTHER object voice)
#   BATH nev rate_per_cell        (bath door rate baseline, events/t.u./cell)
BEGIN { cap = 2.5; tmax = 0 }
/^t=/ {
    tt = $0; sub(/^t= */, "", tt); t = tt + 0
    if (t > tmax) tmax = t
    if (match($0, /xUDD=\(([0-9.+-]+),([0-9.+-]+),([0-9.+-]+)\)/, mm)) {
        for (v = 0; v < 3; v++) xv[v] = mm[v+1] + 0
        if (t >= tmin) { for (v = 0; v < 3; v++) xs[v] += xv[v]; nx++ }
    }
    if (match($0, /INIT NC=([0-9]+)/, mm)) nc = mm[1] + 0
    next
}
/INIT NC=/ { if (match($0, /INIT NC=([0-9]+)/, mm)) nc = mm[1] + 0; next }
/^# QATOM/ {
    split($3, a, "="); t = a[2] + 0
    split($4, a, "="); dir = a[2]
    split($6, a, "="); e = a[2] + 0
    split($7, a, "="); ci = a[2] + 0
    if (t < tmin) next
    if (ci >= i0 && ci <= i0 + 2) {
        v = ci - i0
        if (dir == "DF") {
            ndf[v]++; edf[v] += e; lastDF[v] = t
        } else {
            nfd[v]++; efd[v] += e
            if (v in lastDF) { lag[v, ++nlag[v]] = t - lastDF[v] }
            best = -1
            for (u = 0; u < 3; u++) if (u != v && (u in lastDF) && lastDF[u] > best) best = lastDF[u]
            if (best >= 0) xlag[v, ++nxlag[v]] = t - best
        }
    } else nbath++
    next
}
function q(arr, n, f,   k) { k = int(f * n); if (k < 1) k = 1; return arr[k] }
END {
    W = tmax - tmin
    if (W <= 0 || nx == 0) { print "ERROR: empty window"; exit 1 }
    for (v = 0; v < 3; v++) {
        xa = xs[v] / nx; emfl = xa * cap
        gross = (edf[v] + efd[v]) / W; net = (efd[v] - edf[v]) / W
        tau = gross > 0 ? emfl / gross : -1
        printf "VOICE %d nDF=%d nFD=%d rateDF=%.4f rateFD=%.4f eDFr=%.5f eFDr=%.5f gross=%.5f net=%+.5f x_avg=%.4f Emfl=%.4f tau_turn=%.1f\n", \
            v, ndf[v], nfd[v], ndf[v]/W, nfd[v]/W, edf[v]/W, efd[v]/W, gross, net, xa, emfl, tau
    }
    for (v = 0; v < 3; v++) {
        n = nlag[v]
        if (n > 0) {
            delete s; m = 0
            for (k = 1; k <= n; k++) { s[k] = lag[v, k]; m += lag[v, k] }
            asort(s)
            printf "LAG  %d n=%d median=%.2f p25=%.2f p75=%.2f mean=%.2f\n", \
                v, n, q(s, n, 0.5), q(s, n, 0.25), q(s, n, 0.75), m / n
        }
        n = nxlag[v]
        if (n > 0) {
            delete s
            for (k = 1; k <= n; k++) s[k] = xlag[v, k]
            asort(s)
            printf "XLAG %d n=%d median=%.2f\n", v, n, q(s, n, 0.5)
        }
    }
    if (nc > 3) printf "BATH nev=%d rate_per_cell=%.5f\n", nbath, nbath / (W * (nc - 3))
}
