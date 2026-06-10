#!/usr/bin/env python3
"""Characterize the eta-mediated (curl/theta) parametric instability of the
uniform condensate for amplitudes that are phi-only stable (kappa A^6 > 1/2),
plus eta-dependence at A=0.40. Uses mi_band.full_growth."""
import numpy as np
from mi_band import full_growth, band_max, background, KHATS

def band(A, eta_, khat, kmax=2.5, nk=501, thresh=1e-6):
    ks = np.linspace(1e-4, kmax, nk)
    gs = np.array([full_growth(A,k,khat,eta_) for k in ks])
    pos = gs > thresh
    if not pos.any(): return None
    i = int(np.argmax(gs))
    klo = ks[pos][0]; khi = ks[pos][-1]
    return klo, khi, ks[i], gs[i]

if __name__ == "__main__":
    kh = KHATS["1m10"]   # worst direction for A=0.585
    print("eta scaling of max growth (direction (1,-1,0)/sqrt2):")
    print(f"{'A':>6} {'eta':>5} {'k_lo':>7} {'k_hi':>7} {'k*':>7} {'sigma*':>9} {'local exp':>9}")
    for A in (0.50, 0.585):
        prev = None
        for e_ in (0.125, 0.25, 0.5, 0.75, 1.0):
            r = band(A, e_, kh)
            if r is None:
                print(f"{A:6.3f} {e_:5.3f}   stable (no growth > 1e-6)"); prev=None; continue
            klo,khi,kst,sst = r
            ex = ""
            if prev is not None:
                ex = f"{np.log(sst/prev[1])/np.log(e_/prev[0]):9.2f}"
            print(f"{A:6.3f} {e_:5.3f} {klo:7.4f} {khi:7.4f} {kst:7.4f} {sst:9.5f} {ex:>9}")
            prev = (e_, sst)
        print()
    print("A=0.40: eta dependence of the dominant (phi-sector) mode, worst direction:")
    for e_ in (0.0, 0.25, 0.5, 0.75):
        best = max(((n,)+band(0.40,e_,k) for n,k in KHATS.items() if band(0.40,e_,k)),
                   key=lambda t: t[4])
        print(f"  eta={e_:4.2f}: dir={best[0]:4s} band=({best[1]:.3f},{best[2]:.3f}) k*={best[3]:.4f} sigma*={best[4]:.5f}")
