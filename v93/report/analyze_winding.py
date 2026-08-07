#!/usr/bin/env python3
"""QUENCH-3 offline winding meters (QUENCH.md §7.2).
Input: `fcsdump -mode cells run.fcs` text on stdin or as file arg.
Per frame: W_A (field-phase winding), W_th2 (live-matter clock winding),
m-spectrum R_m = |<e^{i(th2 - m*phi)}>|, radial live-cell profile.
Windings are annulus loop-sums of wrapped phase increments about the
Em/Ee-weighted centroid (advects with the packet/cloud) and the fixed
beam origin (14, 32)."""
import sys, math
import numpy as np

SRC = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (14.0, 32.0)
KX  = (float(sys.argv[4]) if len(sys.argv) > 4 else 2.0)               # seed carrier: phases read modulo the -kx*dx ramp
EM_LIVE = 0.05          # par_hi * cap: the episode-mint criterion
EE_LIT  = 1e-3

def wrap_pi(a):
    return (a + math.pi) % (2*math.pi) - math.pi

def winding(xs, ys, th, cx, cy, r1, r2):
    """loop winding number of phase field th on annulus [r1,r2] about (cx,cy)"""
    dx, dy = xs - cx, ys - cy
    r = np.hypot(dx, dy)
    m = (r >= r1) & (r < r2)
    if m.sum() < 8:
        return None
    phi = np.arctan2(dy[m], dx[m])
    order = np.argsort(phi)
    t = th[m][order]
    d = np.diff(np.concatenate([t, t[:1]]))
    d = (d + np.pi) % (2*np.pi) - np.pi
    return d.sum() / (2*np.pi)

def m_spectrum(xs, ys, th, cx, cy, r1=2.0, r2=16.0, ms=range(-4,5)):
    dx, dy = xs - cx, ys - cy
    r = np.hypot(dx, dy)
    sel = (r >= r1) & (r < r2)
    if sel.sum() < 12:
        return {}
    phi = np.arctan2(dy[sel], dx[sel])
    t = th[sel]
    return {m: abs(np.exp(1j*(t - m*phi)).mean()) for m in ms}

def main():
    f = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    header = f.readline()
    frames = {}
    for line in f:
        p = line.split()
        if len(p) < 14: continue
        frames.setdefault(float(p[0]), []).append(
            (float(p[2]), float(p[3]), float(p[7]), float(p[8]),
             float(p[11]), float(p[12]), float(p[13])))
    print("t nlive | field: cx cy W_A(c) W_A(src) | matter: cx cy W_th2(c) W_th2(src) | R2 R0 Rm_peak")
    for t in sorted(frames):
        a = np.array(frames[t])
        x, y, em, ee, fa1, fa2, th2 = a[:,0], a[:,1], a[:,2], a[:,3], a[:,4], a[:,5], a[:,6]
        # field side
        lit = ee > EE_LIT
        line = f"{t:7.1f} "
        if lit.sum() > 20:
            w = ee[lit]
            fcx, fcy = np.average(x[lit], weights=w), np.average(y[lit], weights=w)
            argA = np.arctan2(fa2, fa1) + KX * (x - SRC[0])   # carrier demod
            wa  = [winding(x[lit], y[lit], argA[lit], fcx, fcy, r1, r2)
                   for r1, r2 in [(1.5,4),(4,7),(7,11)]]
            was = [winding(x[lit], y[lit], argA[lit], *SRC, r1, r2)
                   for r1, r2 in [(1.5,4),(4,7),(7,11)]]
            wa  = [v for v in wa if v is not None]
            was = [v for v in was if v is not None]
            fW  = np.median(wa) if wa else float('nan')
            fWs = np.median(was) if was else float('nan')
            specA = m_spectrum(x[lit], y[lit], argA[lit], fcx, fcy, r1=1.0, r2=12.0)
            RA2 = specA.get(2, float('nan'))
            pkA = max(specA, key=specA.get) if specA else 0
        else:
            fcx = fcy = fW = fWs = RA2 = float('nan'); pkA = 0
        live = em >= EM_LIVE
        nl = int(live.sum())
        line += (f"{nl:5d} | F: {fcx:5.1f} {fcy:5.1f} {fW:+7.3f} {fWs:+7.3f} "
                 f"RA2={RA2:.3f} pkA=m{pkA:+d} | ")
        if nl > 20:
            w = em[live]
            mcx, mcy = np.average(x[live], weights=w), np.average(y[live], weights=w)
            wt  = [winding(x[live], y[live], th2[live], mcx, mcy, r1, r2)
                   for r1, r2 in [(1.5,4),(4,7),(7,11)]]
            wts = [winding(x[live], y[live], th2[live], *SRC, r1, r2)
                   for r1, r2 in [(1.5,4),(4,7),(7,11)]]
            wt  = [v for v in wt if v is not None]
            wts = [v for v in wts if v is not None]
            mW  = np.median(wt) if wt else float('nan')
            mWs = np.median(wts) if wts else float('nan')
            spec = m_spectrum(x[live], y[live], th2[live], mcx, mcy)
            th2d = th2 + KX * (x - SRC[0])                    # carrier demod
            specd = m_spectrum(x[live], y[live], th2d[live], mcx, mcy)
            r2v, r0v = spec.get(2, float('nan')), spec.get(0, float('nan'))
            r2d = specd.get(2, float('nan'))
            pk = max(spec, key=spec.get) if spec else 0
            pkd = max(specd, key=specd.get) if specd else 0
            line += (f"M: {mcx:5.1f} {mcy:5.1f} {mW:+7.3f} {mWs:+7.3f} | "
                     f"R2={r2v:.3f} R0={r0v:.3f} pk=m{pk:+d} | "
                     f"R2d={r2d:.3f} pkd=m{pkd:+d}")
        else:
            line += "M: (no live matter)"
        print(line)
    # final-frame radial profile of live matter about source
    t = max(frames)
    a = np.array(frames[t])
    x, y, em = a[:,0], a[:,1], a[:,2]
    live = em >= EM_LIVE
    dx, dy = x[live]-SRC[0], y[live]-SRC[1]
    r = np.hypot(dx, dy)
    hist, edges = np.histogram(r, bins=np.arange(0, 26, 2))
    print("# n(r) live cells about src @ t=%.0f:" % t,
          " ".join(f"[{edges[i]:.0f}-{edges[i+1]:.0f}]:{hist[i]}" for i in range(len(hist))))

if __name__ == "__main__":
    main()
