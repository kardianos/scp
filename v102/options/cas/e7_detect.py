#!/usr/bin/env python3
"""E7 — the door (../PREREG.md §3 E7).

(a) D0 arithmetic on a lock cell (the kernel's own rule, freecell.c:2548,
    2739): econd = 0 where the work function sees matter (pattern_only:
    space deficit e_s0 - Es >= wf_floor), wf_far = 99 elsewhere; per step
    d1 = f_conv (Ee - econd).  Steps (dt = 0.02) to halve a lock cell's
    field, in units of the lock's 1/J.
(c) Detection under D2 (S3 only): a photon on a 1-D lead (t_out = 0.9) ends
    at site o; the link (o, screen cell) carries S3 with an oscillator at
    Omega.  Screen cell holds m = 1 or 2 quanta.  The oscillator's space is
    F0 (no flow), F1 (1-D space lead, t_s = s_k = 0.06) or R-c (t_s = 1.0,
    a broad continuum).  Exact absorption probability (1 - |S|^2 summed)
    against the photon's kinetic energy across the band, polarization
    averaged.  g/Delta = 0.1 and 0.2.
"""
import sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
TOUT = 0.9
J01 = 0.009126                       # the n = 6 lock's J at g/Delta = 0.1 (E1, T0)
print("# E7(a) D0 on a lock cell: the door's arithmetic")
Ee = G.EF[2]; dt = 0.02
for (label, econd) in (("work function open (space deficit >= wf_floor): econd = 0", 0.0),
                       ("work function closed (no deficit): econd = wf_far = 99", 99.0)):
    if Ee <= econd:
        print(f"  {label}: Ee = {Ee:.4f} <= econd: the door does not fire. The lock is exempt without a label.")
        continue
    k = 0; e = Ee
    while e > Ee/2:
        e -= 0.25*(e - econd); k += 1
    print(f"  {label}: halves the cell's field in {k} steps = {k*dt:.3f} time units = {k*dt*J01:.2e}/J "
          f"(the lock exchanges once per 1/J = {1/J01:.0f} time units)")

def sigma(w, t):
    w = np.asarray(w, dtype=float); inside = np.abs(w) < 2*t
    out = np.empty(w.shape, complex)
    out[inside] = (w[inside] - 1j*np.sqrt(4*t*t - w[inside]**2))/2
    wo = w[~inside]; out[~inside] = (wo - np.sign(wo)*np.sqrt(wo**2 - 4*t*t))/2
    return out

print("\n# E7(c) detection under D2: absorption of a lone photon at a screen cell")
for gr_ in (0.1, 0.2):
    g = gr_*G.DELTA
    for m in (1, 2):
        # P1: screen states (n1 = m..0) x photon pol s;  P2: depleted screen (m-1 quanta) with c_b = 1
        scr = [(a, m - a) for a in range(m, -1, -1)]
        dep = [(a, m - 1 - a) for a in range(m - 1, -1, -1)]
        e_scr = 0.0                                           # all P1 states degenerate (single cell)
        e_dep = G.OMEGA - G.EF[1] - (G.EF[m] - G.EF[m - 1])   # P2 relative to P1 (photon E_f(1) gone, screen loses one)
        V = np.zeros((len(dep), 2*len(scr)))
        for k, (a, b) in enumerate(scr):
            for s in (0, 1):
                for al in (0, 1):
                    M = G.WEDGE[al, s]
                    if M == 0: continue
                    n_al = (a, b)[al]
                    if n_al == 0: continue
                    nc = [a, b]; nc[al] -= 1
                    V[dep.index(tuple(nc)), 2*k + s] += g*M*np.sqrt(n_al)
        for F, ts in (("F0", 0.0), ("F1", 0.06), ("R-c", 1.0)):
            eps = np.linspace(-2*TOUT + 1e-3, 2*TOUT - 1e-3, 3601)
            absb = np.zeros(len(eps))
            for q, E in enumerate(eps):
                S1 = sigma(E - e_scr, TOUT)                    # all channels share the kinetic energy
                S2 = sigma(E - e_dep, ts) if ts > 0 else 0.0
                D2 = 1.0/(E - e_dep - S2 + 1e-13j)
                A = (E - e_scr - S1)*np.eye(2*len(scr)) - D2*(V.T@V)
                Gam = -2*S1.imag
                Gf = np.linalg.inv(A)
                S = np.eye(2*len(scr)) - 1j*Gam*Gf
                absb[q] = np.mean(1 - (np.abs(S)**2).sum(0))
            eres = e_dep
            peak = absb.max(); epk = eps[absb.argmax()]
            width = (absb > peak/2).sum()*(eps[1] - eps[0]) if peak > 1e-12 else 0.0
            print(f"  g/Delta={gr_} m={m} {F:4s}: eps_res={eres:+.4f}  band-averaged absorption={absb.mean():.4f}  "
                  f"peak={peak:.4f} at eps={epk:+.4f}  FWHM={width:.4f}  absorption at the band bottom={absb[0]:.2e}", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
