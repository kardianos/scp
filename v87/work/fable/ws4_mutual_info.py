#!/usr/bin/env python3
"""ws4 — FABLE seat, v87 B0.

The p-family scan: S_max(p) and the mutual-information cost I(p).

E_p(theta) = -1 + ((p+1)/pi) B((p+2)/2,1/2) Int_0^theta sin^p phi dphi   (ws2)

S_max(p): maximise S over coplanar settings (alpha, alpha', beta, beta').
I(p): the model's rho(lambda|a,b) = c_p |lambda.b|^p depends only on the
B-side setting axis; with settings uniform over the four CHSH pairs the
B-setting is uniform over {b, b'} at relative angle phi_opt(p), so
  I(lambda : a,b) = I(lambda : b)
    = (1/2) Int [rho_b log2(rho_b/m) + rho_b' log2(rho_b'/m)] dOmega,
  m = (rho_b + rho_b')/2.

Purpose (the Tsirelson question): if S_max(p) passes SMOOTHLY through 2sqrt2
at p = 1 with nothing structural distinguishing p = 1, then the Tsirelson
bound is IMPORTED by the choice p = 1, not derived from the relaxation.
"""
import numpy as np
from scipy import integrate
from scipy.optimize import minimize
from scipy.special import beta as Bfun

def E_p(pval, theta):
    theta = abs(theta) % (2*np.pi)
    if theta > np.pi:
        theta = 2*np.pi - theta
    I = integrate.quad(lambda ph: np.sin(ph)**pval, 0, theta,
                       epsabs=1e-12, limit=200)[0]
    return -1 + (pval+1)/np.pi * Bfun((pval+2)/2, 0.5) * I

def S_of(pval, ang):
    al, alp, be, bep = ang
    return (E_p(pval, al-be) + E_p(pval, al-bep)
            + E_p(pval, alp-be) - E_p(pval, alp-bep))

def smax(pval, ntries=40, seed=7):
    rng = np.random.default_rng(seed)
    best = (-1e9, None)
    for _ in range(ntries):
        x0 = rng.uniform(0, 2*np.pi, 4)
        r = minimize(lambda x: -S_of(pval, x), x0, method='Nelder-Mead',
                     options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 4000})
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
    return best

def mi(pval, phi):
    """I(lambda:b) in bits for b-settings {z-axis, axis at angle phi}, uniform."""
    cp = (pval + 1) / (4*np.pi)
    def f(vv, uu):
        cb = np.cos(uu)
        cbp = np.cos(uu)*np.cos(phi) + np.sin(uu)*np.cos(vv)*np.sin(phi)
        r1 = cp * abs(cb)**pval
        r2 = cp * abs(cbp)**pval
        m = 0.5*(r1 + r2)
        out = 0.0
        if r1 > 0:
            out += 0.5 * r1 * np.log2(r1/m)
        if r2 > 0:
            out += 0.5 * r2 * np.log2(r2/m)
        return out * np.sin(uu)
    val, err = integrate.dblquad(f, 0, np.pi, 0, 2*np.pi,
                                 epsabs=1e-10, epsrel=1e-10)
    return val, err

print(" p      S_max(p)    b^b' angle(deg)   I(lambda:settings) bits")
rows = []
for pval in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0):
    Sm, ang = smax(pval)
    al, alp, be, bep = ang
    phi = abs(be - bep) % (2*np.pi)
    if phi > np.pi:
        phi = 2*np.pi - phi
    Ival, Ierr = mi(pval, phi) if pval > 0 else (0.0, 0.0)
    rows.append((pval, Sm, np.degrees(phi), Ival))
    print(f"{pval:5.2f}   {Sm:.8f}   {np.degrees(phi):8.3f}       {Ival:.6f} (+-{Ierr:.1e})")

np.save('/home/d/code/scp/v87/work/fable/ws4_rows.npy', np.array(rows))

# sanity anchors
S0 = smax(0.0)[0]; S1 = smax(1.0)[0]
print(f"\nanchors: S_max(0) = {S0:.9f} (must be 2), "
      f"S_max(1) = {S1:.9f} (must be 2sqrt2 = {2*np.sqrt(2):.9f})")
assert abs(S0 - 2) < 1e-6
assert abs(S1 - 2*np.sqrt(2)) < 1e-6
print("ws4 OK.")
