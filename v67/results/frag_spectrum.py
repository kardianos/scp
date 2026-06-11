#!/usr/bin/env python3
"""Fragment charge spectrum analysis from sfa_qball_track TSVs (v67 cf runs)."""
import sys, statistics

BINS = [0, 10, 30, 87, 170, 300, 600, 1200, 1e18]
BLAB = ["<10", "10-30", "30-87", "87-170", "170-300", "300-600", "600-1200", ">1200"]

def load(path):
    frames = {}
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split('\t')
            fi = int(p[0]); t = float(p[1])
            frames.setdefault((fi, t), []).append(
                dict(nvox=int(p[4]), mass=float(p[5]), Q=float(p[6]),
                     rho2_peak=float(p[10]), rms=float(p[11])))
    return frames

def spectrum(clusters, Qbox, qmin=87.0):
    # separate percolating web (rms_size comparable to box half-width)
    web = [c for c in clusters if c['rms'] > 10.0]
    frags = [c for c in clusters if c['rms'] <= 10.0]
    absQ = sorted(abs(c['Q']) for c in frags)
    out = {}
    out['n_frag'] = len(frags)
    out['n_web'] = len(web)
    out['web_Q'] = sum(c['Q'] for c in web)
    out['sumQ'] = sum(c['Q'] for c in frags)
    out['sumAbsQ'] = sum(absQ)
    out['n_neg'] = sum(1 for c in frags if c['Q'] < 0)
    if absQ:
        out['medQ'] = statistics.median(absQ)
        out['maxQ'] = absQ[-1]
        out['minQ'] = absQ[0]
        out['meanQ'] = out['sumAbsQ']/len(absQ)
        out['qwmean'] = sum(q*q for q in absQ)/out['sumAbsQ']  # charge-weighted mean
        above = [q for q in absQ if q >= qmin]
        out['n_above'] = len(above)
        out['Q_above'] = sum(above)
        out['frac_above_clusterQ'] = sum(above)/out['sumAbsQ'] if out['sumAbsQ'] else 0.0
        out['frac_above_boxQ'] = sum(above)/Qbox
        hist = [0]*(len(BINS)-1)
        for q in absQ:
            for i in range(len(BINS)-1):
                if BINS[i] <= q < BINS[i+1]:
                    hist[i] += 1; break
        out['hist'] = hist
    out['cluster_Q_frac_of_box'] = (out['sumQ']+out['web_Q'])/Qbox
    return out

def report(run, path, Qbox, want_times):
    frames = load(path)
    keys = sorted(frames)
    print(f"\n=== {run}  ({path})  Qbox={Qbox:.0f} ===")
    print("t\tN_clusters\tN_frag(rms<=10)\tN_web\tsum|Q|_frag\tQ_web")
    for (fi, t) in keys:
        s = spectrum(frames[(fi, t)], Qbox)
        print(f"{t:.0f}\t{len(frames[(fi,t)])}\t{s['n_frag']}\t{s['n_web']}\t{s['sumAbsQ']:.1f}\t{s['web_Q']:.1f}")
    for tw in want_times:
        key = min(keys, key=lambda k: abs(k[1]-tw))
        s = spectrum(frames[key], Qbox)
        print(f"\n--- {run} frame t={key[1]:.0f} spectrum ---")
        if s['n_frag'] == 0:
            print("no discrete fragments"); continue
        print(f"N_frag={s['n_frag']}  (web clusters: {s['n_web']}, Q_web={s['web_Q']:.0f})")
        print(f"|Q|: min={s['minQ']:.2f} median={s['medQ']:.2f} mean={s['meanQ']:.2f} "
              f"max={s['maxQ']:.1f} charge-weighted-mean={s['qwmean']:.1f}")
        print(f"signed: n_neg={s['n_neg']}/{s['n_frag']}  sumQ={s['sumQ']:.1f} sum|Q|={s['sumAbsQ']:.1f}")
        print(f">=Qmin87: n={s['n_above']}  Q={s['Q_above']:.1f}  "
              f"frac_of_fragQ={s['frac_above_clusterQ']:.3f}  frac_of_boxQ={s['frac_above_boxQ']:.4f}")
        print("hist |Q|: " + "  ".join(f"{l}:{n}" for l, n in zip(BLAB, s['hist'])))
        print(f"total tracked Q (frag+web) / box = {s['cluster_Q_frac_of_box']:.3f}")

if __name__ == "__main__":
    R = "/home/d/code/scp/v67/results"
    for thr in ["0.10", "0.15", "0.20"]:
        print(f"\n################ threshold {thr} ################")
        report("cf2_a40_e0",   f"{R}/frag_cf2_thr{thr}.tsv", 75168.8, [100, 200])
        report("cf1_a40_e05",  f"{R}/frag_cf1_thr{thr}.tsv", 75168.8, [200])
        report("cf4_a585_e05", f"{R}/frag_cf4_thr{thr}.tsv", 170656.0, [300])
