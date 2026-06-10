#!/usr/bin/env python3
"""
v66 Q-ball numerical checks: existence window, radial profile (shooting),
Q(omega), E(omega), stability indicators.

Radial ODE (symmetric ansatz Phi_a = f(r) e^{i w t}, Theta = 0, eta decoupled):
    f'' + (2/r) f' = (m2 - w^2) f + 2 Vt'(f^6) f^5,   Vt'(s) = (mu/2)/(1+kappa s)^2
Effective potential (mechanical analogy: particle in -U_w with friction 2/r):
    U_w(f) = (1/2)(m2-w^2) f^2 + (1/3) Vt(f^6)
Existence: U_w(f) < 0 for some f  =>  w^2 > wmin2 = m2 + (mu/9)(2/kappa)^(2/3)

Observables:
    Q(w) = 3 w * int 4 pi r^2 f^2 dr
    E(w) = int 4 pi r^2 [ (3/2)(w^2 f^2 + f'^2 + m2 f^2) + Vt(f^6) ] dr
Checks: E/Q < m (evaporation bound), dE/dw = w dQ/dw, ODE residual.
"""
import numpy as np

m2, mu, kap = 2.25, -41.345, 50.0
m = np.sqrt(m2)

def Vt(s):  return 0.5*mu*s/(1.0+kap*s)
def Vtp(s): return 0.5*mu/(1.0+kap*s)**2
def Uw(f, w2):  return 0.5*(m2-w2)*f**2 + Vt(f**6)/3.0
def Uwp(f, w2): return (m2-w2)*f + 2.0*Vtp(f**6)*f**5

# ---- existence window (closed form, cross-check) ----
wmin2 = m2 + (mu/9.0)*(2.0/kap)**(2.0/3.0)
fstar = (2.0/kap)**(1.0/6.0)
# brute-force minimum of m2 + (mu/3) f^4/(1+kap f^6)
fg = np.linspace(1e-4, 3.0, 400000)
brute = m2 + (mu/3.0)*fg**4/(1.0+kap*fg**6)
i = np.argmin(brute)
print(f"closed form: wmin2={wmin2:.6f} wmin={np.sqrt(wmin2):.6f} f*={fstar:.6f}")
print(f"brute force: wmin2={brute[i]:.6f}            f*={fg[i]:.6f}")
assert abs(brute[i]-wmin2) < 1e-6 and abs(fg[i]-fstar) < 1e-4
print(f"existence window: omega in ({np.sqrt(wmin2):.4f}, {m:.4f})\n")

# ---- shooting solver ----
def integrate(f0, w2, R=200.0, dr=2e-3, full=False):
    """RK4 from r~0 in long double. Classify:
       'over'  : f crosses 0 (too much energy -> f0* below mid)
       'right' : f runs away above f0 (started past the dip -> f0* below mid)
       'under' : f' > 0 while 0 < f < f0 (friction won -> f0* above mid)."""
    ld = np.longdouble
    f, fp, r, dr = ld(f0), ld(0.0), ld(1e-9), ld(dr)
    if full: rs, fs, fps = [r], [f], [fp]
    def deriv(r, f, fp):
        return fp, Uwp(f, w2) - (ld(2.0)/r)*fp
    n = int(R/float(dr))
    tag = 'end'
    for _ in range(n):
        k1 = deriv(r, f, fp)
        k2 = deriv(r+dr/2, f+dr/2*k1[0], fp+dr/2*k1[1])
        k3 = deriv(r+dr/2, f+dr/2*k2[0], fp+dr/2*k2[1])
        k4 = deriv(r+dr,   f+dr*k3[0],   fp+dr*k3[1])
        f  += dr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        fp += dr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        r  += dr
        if full: rs.append(r); fs.append(f); fps.append(fp)
        if f < 0:            tag = 'over';  break
        if f > 1.05*f0:      tag = 'right'; break
        if fp > 0 and f < 0.99*f0:
            tag = 'under'; break
    if full:
        return tag, np.array(rs, float), np.array(fs, float), np.array(fps, float)
    return tag, float(r)

def find_f0(w2):
    """Bisect f0 in (f1, f_dip): f1 = zero crossing of U_w into U<0,
       f_dip = bottom of the dip (U'=0). Long-double bisection (thin wall)."""
    ld = np.longdouble
    fg = np.linspace(1e-4, 2.0, 200000)
    U = Uw(fg, w2)
    neg = np.where(U < 0)[0]
    assert len(neg) > 0, "no negative U region: omega outside window"
    f1 = fg[neg[0]]
    # precise dip: bisect U'=0 between f1 and the upper end of the U<0 region
    a, b = ld(f1), ld(fg[neg[-1]])
    for _ in range(200):
        mfd = (a+b)/2
        if Uwp(mfd, ld(w2)) < 0: a = mfd
        else: b = mfd
    fdip = a            # just below the true dip (U' < 0 side)
    lo, hi = ld(f1)*(1+ld(1e-9)), fdip
    if integrate(float(lo), w2)[0] != 'under':
        raise RuntimeError("lo bracket not under")
    for _ in range(110):
        mid = (lo+hi)/2
        tag, _ = integrate(mid, w2)
        if tag == 'under': lo = mid
        else: hi = mid          # 'over'/'right'/'end' -> reduce f0
    return float(lo)

def observables(w):
    w2 = w*w
    f0 = find_f0(w2)
    tag, rs, fs, fps = integrate(f0, w2, full=True)
    # truncate at breakdown point, attach exponential tail f ~ C e^{-k r}/r
    rc, fc, fpc = rs[-1], max(fs[-1], 0.0), fps[-1]
    k = np.sqrt(m2 - w2)
    # numerical part
    r2 = rs*rs
    I2  = np.trapz(4*np.pi*r2*fs**2, rs)                       # int f^2 dV
    Egr = np.trapz(4*np.pi*r2*1.5*fps**2, rs)
    Epot= np.trapz(4*np.pi*r2*(1.5*m2*fs**2 + Vt(fs**6)), rs)
    Ekin= 1.5*w2*I2
    # tail (f = fc*(rc/r) e^{-k(r-rc)}): int f^2 dV = 4 pi fc^2 rc^2/(2k)
    tail2 = 4*np.pi*fc**2*rc**2/(2*k) if fc > 0 else 0.0
    I2 += tail2; Ekin += 1.5*w2*tail2; Epot += 1.5*m2*tail2
    Q = 3*w*I2
    E = Ekin + Egr + Epot
    # radii (first downward crossing of f0/2)
    below = np.where(fs <= 0.5*f0)[0]
    rhalf = rs[below[0]] if len(below) else np.nan
    rQ = np.sqrt(np.trapz(4*np.pi*r2*r2*fs**2, rs)/I2)         # charge rms radius
    # central energy density
    e0 = 1.5*(w2+m2)*f0**2 + Vt(f0**6)
    # ODE residual check (midpoints, away from r=0)
    sel = (rs > 0.5) & (fs > 1e-6*f0)
    fpp = np.gradient(fps[sel], rs[sel])
    res = fpp + 2/rs[sel]*fps[sel] - Uwp(fs[sel], w2)
    resmax = np.max(np.abs(res))/max(np.max(np.abs(Uwp(fs[sel], w2))), 1e-30)
    return dict(w=w, f0=f0, Q=Q, E=E, EoQ=E/Q, rhalf=rhalf, rQ=rQ, e0=e0,
                s0=f0**6, tag=tag, resmax=resmax, tailfrac=tail2/I2)

ws = [1.315, 1.33, 1.36, 1.39, 1.42, 1.44, 1.45, 1.47, 1.49]
rows = []
print(f"{'omega':>6} {'f0':>8} {'s0=f0^6':>9} {'e(0)':>8} {'Q':>10} {'E':>10} "
      f"{'E/Q':>7} {'E-mQ':>9} {'r_half':>7} {'r_Q':>6} {'res':>8} {'tail':>8}")
for w in ws:
    o = observables(w)
    rows.append(o)
    print(f"{o['w']:6.3f} {o['f0']:8.5f} {o['s0']:9.5f} {o['e0']:8.4f} "
          f"{o['Q']:10.3f} {o['E']:10.3f} {o['EoQ']:7.4f} {o['E']-m*o['Q']:9.3f} "
          f"{o['rhalf']:7.3f} {o['rQ']:6.3f} {o['resmax']:8.1e} {o['tailfrac']:8.1e}")

# ---- stability indicators ----
print("\ndQ/domega (finite diff) and dE/dw - w*dQ/dw consistency:")
for a, b in zip(rows[:-1], rows[1:]):
    dQ = (b['Q']-a['Q'])/(b['w']-a['w'])
    dE = (b['E']-a['E'])/(b['w']-a['w'])
    wm = 0.5*(a['w']+b['w'])
    print(f"  w~{wm:5.3f}: dQ/dw={dQ:10.2f}  dE/dw={dE:10.2f}  w*dQ/dw={wm*dQ:10.2f}"
          f"  rel.err={(dE-wm*dQ)/abs(dE):+8.1e}  {'STABLE(dQ/dw<0)' if dQ<0 else 'unstable'}")

print("\nevaporation bound E < m*Q per row:",
      ["%.3f:%s" % (o['w'], 'OK' if o['E'] < m*o['Q'] else 'VIOLATED') for o in rows])
# locate the E = m*Q crossing
for a, b in zip(rows[:-1], rows[1:]):
    ga, gb = a['E']-m*a['Q'], b['E']-m*b['Q']
    if ga < 0 <= gb:
        wc = a['w'] + (b['w']-a['w'])*(-ga)/(gb-ga)
        print(f"E = m*Q crossing (evaporation threshold): omega ~ {wc:.3f}")
qmin = min(rows, key=lambda o: o['Q'])
print(f"minimum charge on computed grid: Q_min ~ {qmin['Q']:.1f} at omega = {qmin['w']}")
print("NOTE: omega=1.315 row is thin-wall-limited (bisection at long-double "
      "precision floor; see resmax/tailfrac) — treat as approximate.")
