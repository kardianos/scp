#!/usr/bin/env python3
"""WP-B numerics: winding quantization table, polariton branch steps,
and second-order Q-ball frequency shift in a helical (Beltrami) theta background.

Theory inputs (Maxima-verified, theta_winding.mac 20/20 PASS):
  - traveling circular theta wave: tu = B(cos psi, sin psi, 0), psi=K z - K t (m_th=0)
    energy density e_bath = K^2 B^2;  curl tu = -K tu (Beltrami, lambda_B=-K)
  - induced phi condensate: c = eta K B / m^2 (on shell)
  - second-order mean potential shift density on a ball with profile f(r):
        d<Vt> = c^2 g(f),  g(f) = (mu/2) f^4 (1 - kap f^6) / (1 + kap f^6)^3
  - ball shift: dE = c^2 * I_g,  I_g = int 4 pi r^2 g(f) dr
    domega(at fixed Q) = d(dE)/dQ  ~ [dE(w2)-dE(w1)]/[Q(w2)-Q(w1)]
"""
import numpy as np

m2, mu, kap, eta = 2.25, -41.345, 50.0, 0.5
L = 25.0                      # box side 2L = 50
dK = np.pi / L                # winding step

def Vtp(s):  return (mu/2.0)/(1.0+kap*s)**2

def rhs(r, f, fp, w2):
    return (m2 - w2)*f + 2.0*Vtp(f**6)*f**5

def shoot(f0, w2, rmax=30.0, dr=2e-3):
    """RK4 integrate; return (status, r_arr, f_arr). status: +1 overshoot (f<0),
    -1 undershoot (turnaround with f>0), 0 reached rmax."""
    n = int(rmax/dr)
    r = np.empty(n+1); F = np.empty(n+1)
    f, fp = f0, 0.0
    r[0], F[0] = 0.0, f0
    for i in range(n):
        ri = i*dr
        def acc(rr, ff, gg):
            fric = 2.0/rr*gg if rr > 1e-12 else 0.0
            return rhs(rr, ff, gg, w2) - fric
        k1f, k1g = fp, acc(ri, f, fp)
        k2f, k2g = fp+0.5*dr*k1g, acc(ri+0.5*dr, f+0.5*dr*k1f, fp+0.5*dr*k1g)
        k3f, k3g = fp+0.5*dr*k2g, acc(ri+0.5*dr, f+0.5*dr*k2f, fp+0.5*dr*k2g)
        k4f, k4g = fp+dr*k3g,     acc(ri+dr,     f+dr*k3f,     fp+dr*k3g)
        f  += dr*(k1f+2*k2f+2*k3f+k4f)/6.0
        fp += dr*(k1g+2*k2g+2*k3g+k4g)/6.0
        r[i+1], F[i+1] = ri+dr, f
        if f < 0:                 return +1, r[:i+2], F[:i+2]  # overshoot
        if f > 1.02*f0:           return +2, r[:i+2], F[:i+2]  # runaway (f0 past barrier)
        if fp > 0 and f < 0.9*f0: return -1, r[:i+2], F[:i+2]  # undershoot
    return 0, r, F

def profile(w):
    w2 = w*w
    lo, hi = 0.05, 1.2
    # bracket: find lo undershoot, hi overshoot
    for _ in range(80):
        mid = 0.5*(lo+hi)
        s, r, F = shoot(mid, w2)
        if s in (+1, +2): hi = mid
        else:             lo = mid
    s, r, F = shoot(lo, w2)
    # truncate at the minimum |f| before divergence/turnaround
    i = np.argmin(np.abs(F[len(F)//2:])) + len(F)//2
    return r[:i+1], F[:i+1]

def ball_integrals(w):
    r, f = profile(w)
    dr = r[1]-r[0]
    fp = np.gradient(f, dr)
    s = f**6
    Vt = (mu/2.0)*s/(1.0+kap*s)
    g  = (mu/2.0)*f**4*(1.0-kap*s)/(1.0+kap*s)**3
    vol = 4.0*np.pi*r**2
    Q  = np.trapz(vol*3.0*w*f**2, r)
    E  = np.trapz(vol*(1.5*(w*w*f**2 + fp**2 + m2*f**2) + Vt), r)
    Ig = np.trapz(vol*g, r)
    return f[0], Q, E, Ig

print("=== Winding quantization (L=25, box side 50) ===")
print(f"dK = pi/L = {dK:.6f}  (Beltrami eigenvalue lambda_B = -K_n, K_n = n*dK)")
print("bath band |k| in [1.1,1.7]  ->  n in [%.2f, %.2f]  ->  n = 9..13" % (1.1/dK, 1.7/dK))

print("\n=== Polariton branches at K_n (m_th=0, eta=0.5; helicity-degenerate, W4) ===")
def branches(k):
    b = 2*k*k+m2
    d = np.sqrt(m2*m2 + 4*eta*eta*k*k)
    return np.sqrt((b-d)/2), np.sqrt((b+d)/2)
print(" n    K_n     w_-      w_+     (dw_- per n)")
prev = None
for n in range(8, 15):
    k = n*dK
    wm, wp = branches(k)
    step = '' if prev is None else f"{wm-prev:.5f}"
    prev = wm
    print(f"{n:2d}  {k:.4f}  {wm:.5f}  {wp:.5f}   {step}")
vph = np.sqrt(1-eta*eta/m2)
print(f"small-k lower-branch v_phase = sqrt(1-eta^2/m^2) = {vph:.6f}  (v50: 0.9428)")

print("\n=== Induced condensate amplitude c = eta K B / m^2 ===")
e_bath = 0.2
print("fixed e_bath=0.2:  c = eta*sqrt(e_bath)/m^2 = %.6f  (n-INDEPENDENT)" %
      (eta*np.sqrt(e_bath)/m2))
B14 = np.sqrt(e_bath)/1.4
print(f"fixed B (B={B14:.4f}, i.e. e_bath=0.2 at K=1.4):  c_n = {eta*B14*dK/m2:.6f} * n  (linear in n)")
for n in (9, 11, 13):
    print(f"   n={n:2d}: K={n*dK:.4f}  c_n={eta*B14*n*dK/m2:.6f}  e_bath(n)={(n*dK*B14)**2:.4f}")

print("\n=== Q-ball profiles and second-order shift dE = c^2 * I_g ===")
c2 = (eta*np.sqrt(e_bath)/m2)**2     # fixed e_bath normalization
rows = []
for w in (1.39, 1.41, 1.44):
    f0, Q, E, Ig = ball_integrals(w)
    rows.append((w, f0, Q, E, Ig))
    print(f"w={w:.2f}: f(0)={f0:.4f}  Q={Q:7.1f}  E={E:7.1f}  I_g={Ig:8.2f}"
          f"  dE=c^2*I_g={c2*Ig:+.4f}  dE/E={c2*Ig/E:+.2e}")
print("(v66 cross-check at w=1.39: f(0)=0.64, Q=482, E=692)")

w1, _, Q1, E1, Ig1 = rows[0][0], *rows[0][1:]
w2_, _, Q2, E2, Ig2 = rows[1][0], *rows[1][1:]
domega = c2*(Ig2-Ig1)/(Q2-Q1)
print(f"\nfixed-Q frequency shift d(omega) = d(dE)/dQ ~ c^2*(I_g(1.41)-I_g(1.39))/(Q(1.41)-Q(1.39))")
print(f"  = {domega:+.2e}   (relative: {domega/1.40:+.2e})")
print(f"crude bound |dE|/Q = {abs(c2*Ig1)/Q1:.2e}")

print("\n=== n-stepping of the ball shift ===")
print("fixed e_bath: dE independent of n (c = eta sqrt(e_bath)/m^2, K cancels) [thm]")
print("fixed B:      dE(n) = (eta B dK/m^2)^2 * n^2 * I_g  -> QUADRATIC steps in n")
for n in (9, 11, 13):
    cn2 = (eta*B14*n*dK/m2)**2
    print(f"   n={n:2d}: dE = {cn2*Ig1:+.4f}")

print("\n=== AB-analog scattering phase (NOT topological) ===")
# component-2 mass shift from W7 Hessian, time/space-averaged: dm2 = mu*c^4/8
dm2 = mu*c2**2/8.0
kpkt = 1.4
ell = 50.0
print(f"d(m^2)_Phi2 = mu*c^4/8 = {dm2:+.3e}; path phase over l=50 at k=1.4:"
      f" dphi = -dm2*l/(2k) = {-dm2*ell/(2*kpkt):+.3e} rad (continuous in B, NOT n-quantized)")
