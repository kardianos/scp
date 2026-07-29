#!/usr/bin/env python3
"""PLAST pair model — coupled (theta drift, load, d flow) for a 2-cell pair.

Reduced kernel dynamics (laws_V2g constants), Euler at the kernel's dt=0.02:
  pitch      w_i = w2/(1+q x_i)
  gates      psi_f = wrap(th1 - w1 d/C - th2), psi_b = wrap(th2 - w2 d/C - th1)
             G(psi) = ((1+cos psi)/2)^p_gate
  exchange   F_ij = k2 * geo(d) * res(det) * G(psi_f) * head_j * sqrt(m_i m^_j)
  roughness  rho = rough_k * 2|det|G_r/(det^2+G_r^2) of each delivery -> radiated
  entrain    th_j pulled toward retarded source tail at kappa_lock * mix
  PLASTICITY dd/dt = -kappa_p * Phi * dV/dd,
             V = 1 - G(psi_f) G(psi_b),  dV/dd = (1/C)[w1 G'(f)G(b) + w2 G(f)G'(b)]
             Phi = k2 * geo(d) * res(det) * sqrt(m_i m^_j)   (gate-free urge;
             exactly the S2 reactive prefactor 'base*Sm' of cellfab.c pass 2)

Vacuum safety is structural: m_i = 0 and m_j = 0  =>  Phi = 0 exactly.

Outputs: (A) one-sided vs two-sided annealing across the rung-offset axis,
(B) anneal time + stability vs kappa_p, (D) the comma paid in energy vs
geometry (shed comparison).
"""
import numpy as np

# laws_V2g
W2, Q, CAP, C = 2.9, 1.2, 2.5, 1.0
PG, GM = 8.0, 0.10
KDEP, KDEPM = 1.2, 2.0
KLOCK, LOCKF = 1.0, 0.005 * CAP
ROUGHK, GR = 0.35, 0.5
MOBF = 0.004 * CAP
R0 = 0.85
DT = 0.02

def wrap(a): return (a + np.pi) % (2*np.pi) - np.pi
def gate(p): return (0.5*(1+np.cos(p)))**PG
def dgate(p):
    h = 0.5*(1+np.cos(p))
    return -0.5*PG*h**(PG-1)*np.sin(p)

def geo(d):
    """lens overlap conductance (equal radii R0), kernel pass-2 form"""
    ri = rj = R0
    if d >= ri+rj: return 0.0
    a2 = (4*d*d*ri*ri - (d*d)**2) / (4*d*d)   # t = d^2 for equal radii
    A = np.pi*a2 if a2 > 0 else 0.0
    Aref, dref = np.pi*R0*R0, 2*R0
    return (A/Aref)*(dref/d)

def om(x): return W2/(1.0 + Q*x)

def run(d0, x1, x2, kp, T=400.0, th_lock=True, seed=1):
    rng = np.random.default_rng(seed)
    th1, th2 = rng.uniform(0, 2*np.pi), 0.0
    d = d0
    m1, m2 = x1*CAP, x2*CAP
    E0 = m1 + m2
    if th_lock:
        th2 = (th1 - om(x1)*d/C) % (2*np.pi)   # forward-gate seed (kernel practice)
    nst = int(T/DT)
    t_lock = np.nan
    for s in range(nst):
        w1, w2_ = om(m1/CAP), om(m2/CAP)
        det = w1 - w2_
        res = (GM*GM)/(GM*GM + det*det)
        pf = wrap(th1 - w1*d/C - th2)
        pb = wrap(th2 - w2_*d/C - th1)
        gf, gb = gate(pf), gate(pb)
        g0 = geo(d)
        m1h, m2h = max(m1, MOBF), max(m2, MOBF)
        me1, me2 = np.sqrt(max(m1,0)*m2h), np.sqrt(max(m2,0)*m1h)
        h1 = max(0.0, 1.0 - m1/CAP); h2 = max(0.0, 1.0 - m2/CAP)
        base = KDEP*KDEPM*g0*res
        F12 = base*gf*h2*me1          # per unit time
        F21 = base*gb*h1*me2
        rho = ROUGHK * 2*abs(det)*GR/(det*det + GR*GR)
        dm1 = (-F12 + (1-rho)*F21)*DT
        dm2 = (-F21 + (1-rho)*F12)*DT
        # entrainment (deposit+delivery folded): receiver pulled to retarded tail
        if F12 > 0:
            mix = F12*DT/(F12*DT + m2 + LOCKF)
            th2 += KLOCK*mix*wrap(th1 - w1*d/C - th2)
        if F21 > 0:
            mix = F21*DT/(F21*DT + m1 + LOCKF)
            th1 += KLOCK*mix*wrap(th2 - w2_*d/C - th1)
        # plasticity
        if kp > 0:
            Phi = base*np.sqrt(me1*me2)
            dV = (w1*dgate(pf)*gb + w2_*gate(pf)*dgate(pb))/C
            d += -kp*Phi*dV*DT
            if d < 0.5: d = 0.5
        m1 += dm1; m2 += dm2
        th1 = (th1 + w1*DT) % (2*np.pi)
        th2 = (th2 + w2_*DT) % (2*np.pi)
        if np.isnan(t_lock):
            dl = wrap((w1+w2_)*d/C)
            if abs(dl) < 0.10 and gf*gb > 0.9:
                t_lock = s*DT
    w1, w2_ = om(m1/CAP), om(m2/CAP)
    delta = wrap((w1+w2_)*d/C)
    pf = wrap(th1 - w1*d/C - th2); pb = wrap(th2 - w2_*d/C - th1)
    return dict(t_lock=t_lock, delta=delta, gg=gate(pf)*gate(pb),
                ret=(m1+m2)/E0, d=d, x=(m1+m2)/(2*CAP))

# reference: x=0.32 pair, rung at d* = pi C/om
XREF = 0.32
DSTAR = np.pi*C/om(XREF)
print(f"# pair model: x_ref={XREF} om={om(XREF):.4f} d*={DSTAR:.4f} "
      f"(foam dbar=1.505), pinch at d=2r={2*R0}")

print("\n# ==== A: the rung-offset axis, frozen (kp=0) vs plastic ====")
print("# d0/d*    delta0   | kp=0:  gg    ret    delta_f | kp=0.5: gg    ret   delta_f  t_lock")
for f in (0.82, 0.88, 0.94, 1.00, 1.06, 1.12, 1.18):
    d0 = f*DSTAR
    delta0 = wrap(2*om(XREF)*d0/C)
    a = run(d0, XREF+0.02, XREF-0.02, 0.0)
    b = run(d0, XREF+0.02, XREF-0.02, 0.5)
    print(f"#  {f:.2f}   {delta0:+.3f}  |  {a['gg']:.3f}  {a['ret']:.3f}  {a['delta']:+.3f}"
          f"   |  {b['gg']:.3f}  {b['ret']:.3f}  {b['delta']:+.3f}   "
          f"{'%.0f' % b['t_lock'] if not np.isnan(b['t_lock']) else '-':>5}")

print("\n# ==== B: anneal time + stability vs kappa_p (d0 = 1.12 d*, the frozen side) ====")
print("# kappa_p   t_lock   gg_final  ret   d_final/d*")
for kp in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
    b = run(1.12*DSTAR, XREF+0.02, XREF-0.02, kp, T=600.0)
    tl = '%.0f' % b['t_lock'] if not np.isnan(b['t_lock']) else '  -'
    print(f"#  {kp:5.2f}   {tl:>6}    {b['gg']:.3f}   {b['ret']:.3f}   {b['d']/DSTAR:.4f}")

print("\n# ==== D: the comma paid — energy vs geometry (overloaded side, delta0<0) ====")
# overloaded pair: x0 above the tuning-curve load for its d (E8 protocol)
d0 = DSTAR*1.0
x_over = XREF + 0.18       # overloaded: om smaller, delta = 2 om d - 2pi < 0
delta0 = wrap(2*om(x_over)*d0/C)
a = run(d0, x_over+0.02, x_over-0.02, 0.0, T=400)
b = run(d0, x_over+0.02, x_over-0.02, 0.5, T=400)
print(f"# delta0={delta0:+.3f} (overloaded by {x_over-XREF:.2f})")
print(f"#   frozen (kp=0):  shed={1-a['ret']:.3f}  gg={a['gg']:.3f}  delta_f={a['delta']:+.3f}  x_f={a['x']:.3f}")
print(f"#   plastic(kp=.5): shed={1-b['ret']:.3f}  gg={b['gg']:.3f}  delta_f={b['delta']:+.3f}  x_f={b['x']:.3f}  d_f/d*={b['d']/DSTAR:.4f}")

print("\n# ==== C: basin — random phases, near the anti-phase saddle ====")
nl, ok = 0, 0
for s in range(24):
    b = run(1.10*DSTAR, XREF+0.02, XREF-0.02, 0.5, th_lock=False, seed=100+s)
    nl += 1
    if b['gg'] > 0.9: ok += 1
print(f"# random-phase starts at d0=1.10 d*: locked {ok}/{nl} by T=400 (kernel e6: acquisition ~10 t.u.)")
