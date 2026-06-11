#!/usr/bin/env python3
"""v68 IN-CORE DRESSED-STATE THEORY — Q-ball immersed in an isotropic theta bath.

Mean-field model: the bath (energy density e_bath) drives a Phi condensate with
per-component variance n = n0_per_e * e_bath (FRAME band average, vacuum detuning),
locally ENHANCED inside the ball because the in-core Hessian shift W(r) = g(f^2)
lowers the detuning:  n(r) = n0 * Enh(f^2),  Enh(X) = [M2t/(M2t+g(X))]^2.

The dressed radial problem replaces the vacuum potential force g(f^2) f by the
Gaussian-averaged (Wick) force around the BALL background:

    f'' + (2/r) f' = (m^2 - omega^2) f + G_eff(f^2; n) f

Models:
  M1 (first order, exact O(n)):  G_eff = g(X) + n(X) * c(X),
      c(X) = 2 mu X (2 t^2 - 6 t + 1)/(1+t)^4,  t = kappa X^3      [thm, checked here]
      (FRAME A.5 used only the partial slope g'(X) = 2 mu X (1-2t)/(1+t)^3 — that
       misses the a-component cross term  -2 mu t X/(1+t)^3  and the second-order
       <dw^2> Hessian-trace terms  [6 t^2 + 2t(2t-4)] mu X/(1+t)^4.)
  M2 (resummed): full 6D Gauss-Hermite average of the coherent-projected force
      <2 Vt'(w0 w1 w2) w1 w2 (f + da_re)>,  w_b = (f db_a)^2..., circular Gaussian
      per complex component with <|d|^2> = n(X).  Exact in n at mean-field level.

Outputs:
  - CHECK A: M2 slope at n->0 equals c(X)            [PASS/FAIL]
  - CHECK B: vacuum branch reproduces v66 table      [PASS/FAIL]
  - CHECK C: vacuum window edge omega_min = 1.308700 [PASS/FAIL]
  - dressed branch omega(Q=300; e) for e in {0, .05, .2, .8}; matched-Q shifts
  - window edge omega_min(e); n needed to put the edge below 1.294
  - fitted n_eff per rung to reproduce measured shifts (+0.012, -0.025, -0.120)
  - crossover scan delta-omega(e) for the vacuum-drive n(e) = n0*e
"""
import numpy as np, time, sys

t_start = time.time()
m2, mu, kap, eta = 2.25, -41.345, 50.0, 0.5
m = np.sqrt(m2)
n0_per_e = 0.01462          # FRAME band-avg driven variance per component per e_bath

# band-averaged vacuum detuning M2t = m2 + K^2 - W^2(K), photon branch, K in [1.1,1.7]
def W2_photon(K): return (2*K**2 + m2 - np.sqrt(m2**2 + 4*eta**2*K**2))/2
M2t = float(np.mean([m2 + K**2 - W2_photon(K) for K in np.linspace(1.1, 1.7, 13)]))

def Vt(s):  return 0.5*mu*s/(1.0+kap*s)
def Vtp(s): return 0.5*mu/(1.0+kap*s)**2

def g(X):
    t = kap*np.asarray(X, float)**3
    return mu*np.asarray(X, float)**2/(1+t)**2

def gp(X):  # FRAME partial slope g'(X)
    t = kap*X**3
    return 2*mu*X*(1-2*t)/(1+t)**3

def c_exact(X):  # exact O(n) Gaussian slope [derived here, verified CHECK A]
    X = np.asarray(X, float)
    t = kap*X**3
    return 2*mu*X*(2*t**2 - 6*t + 1)/(1+t)**4

def Enh(X):  # in-core drive enhancement [estimate]: detuning lowered by Hessian g(X)
    return (M2t/(M2t + g(X)))**2

# ---------------------------------------------------------------- M2 (resummed)
def G_M2_table(Xg, n0, use_enh=True, H=8):
    """Gaussian-averaged coherent force / f on grid Xg. Returns G_eff(Xg)."""
    F = np.sqrt(Xg)
    nloc = n0*(Enh(Xg) if use_enh else np.ones_like(Xg))
    sig = np.sqrt(np.maximum(nloc, 0.0)/2.0)       # std per real quadrature
    x, wq = np.polynomial.hermite.hermgauss(H)
    d = np.sqrt(2.0)*x                              # delta = d_i * sig
    # precompute w_b for the two spectator components for each (i,j) node pair
    Wb = {}
    for i in range(H):
        for j in range(H):
            Wb[(i, j)] = (F + d[i]*sig)**2 + (d[j]*sig)**2
    acc = np.zeros_like(F)
    for i0 in range(H):
        a0 = d[i0]*sig
        for j0 in range(H):
            w0 = (F + a0)**2 + (d[j0]*sig)**2
            w00 = wq[i0]*wq[j0]
            for i1 in range(H):
                for j1 in range(H):
                    w1 = Wb[(i1, j1)]
                    w01 = w00*wq[i1]*wq[j1]
                    for i2 in range(H):
                        for j2 in range(H):
                            w2 = Wb[(i2, j2)]
                            s = w0*w1*w2
                            acc += (w01*wq[i2]*wq[j2]) * (2*Vtp(s)*w1*w2*(F + a0))
    acc /= np.pi**3
    Geff = np.empty_like(F)
    Geff[1:] = acc[1:]/F[1:]
    Geff[0] = 2*Geff[1] - Geff[2]                   # linear extrapolation to X=0
    return Geff

# ---------------------------------------------------------------- tables
NX, XMAX = 2400, 1.2
Xg = np.linspace(0.0, XMAX, NX)
dXg = Xg[1] - Xg[0]

def cumint(y, x):
    out = np.empty_like(y)
    out[0] = 0.0
    out[1:] = np.cumsum(0.5*(y[1:]+y[:-1])*np.diff(x))
    return out

class Tables:
    """G_eff(X) lookup + Up(X)=0.5 int G dY for brackets / energy / window edge."""
    def __init__(self, Geff):
        self.G = np.asarray(Geff, float)
        self.Up = cumint(self.G, Xg)*0.5
        with np.errstate(divide='ignore', invalid='ignore'):
            q = 2*self.Up/Xg
        q[0] = 0.0
        self.wmin2 = m2 + q.min()
        self.wmin = np.sqrt(max(self.wmin2, 0.0))
    def Gi(self, X):                                 # scalar/array linear interp
        X = np.minimum(np.maximum(X, 0.0), XMAX - 1e-12)
        i = (X/dXg).astype(int) if hasattr(X, '__len__') else int(X/dXg)
        fr = X/dXg - i
        return self.G[i]*(1-fr) + self.G[np.minimum(i+1, NX-1)]*fr
    def Upi(self, X):
        X = np.minimum(np.maximum(X, 0.0), XMAX - 1e-12)
        i = (X/dXg).astype(int) if hasattr(X, '__len__') else int(X/dXg)
        fr = X/dXg - i
        return self.Up[i]*(1-fr) + self.Up[np.minimum(i+1, NX-1)]*fr
    def U(self, f, w2):                              # effective potential U_w(f)
        X = f*f
        return 0.5*(m2 - w2)*X + self.Upi(X)

def make_tables(model, e=0.0, n_override=None, use_enh=True, H=8):
    n0 = n_override if n_override is not None else n0_per_e*e
    if n0 == 0.0:
        return Tables(g(Xg))
    if model == 'M1':
        nl = n0*(Enh(Xg) if use_enh else 1.0)
        return Tables(g(Xg) + nl*c_exact(Xg))
    return Tables(G_M2_table(Xg, n0, use_enh=use_enh, H=H))

# ---------------------------------------------------------------- shooting
DR = 4e-3
def rk4_batch(tb, f0s, w2, R=120.0):
    """classify each candidate: 1=under, 2=over/right/end."""
    B = len(f0s)
    f = f0s.copy(); fp = np.zeros(B)
    tags = np.zeros(B, int)
    r = 1e-9
    n = int(R/DR)
    act = np.ones(B, bool)
    for _ in range(n):
        def deriv(rr, ff, ffp):
            return ffp, (m2 - w2)*ff + tb.Gi(ff*ff)*ff - (2.0/rr)*ffp
        k1 = deriv(r, f, fp)
        k2 = deriv(r+DR/2, f+DR/2*k1[0], fp+DR/2*k1[1])
        k3 = deriv(r+DR/2, f+DR/2*k2[0], fp+DR/2*k2[1])
        k4 = deriv(r+DR,   f+DR*k3[0],   fp+DR*k3[1])
        df  = DR/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        dfp = DR/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        f  = np.where(act, f + df, f)
        fp = np.where(act, fp + dfp, fp)
        r += DR
        over  = act & ((f < 0) | (f > 1.05*f0s))
        under = act & (fp > 0) & (f > 0) & (f < 0.99*f0s)
        tags[over] = 2; tags[under] = 1
        act &= (tags == 0)
        if not act.any():
            break
    tags[tags == 0] = 2                              # 'end' -> reduce f0
    return tags

def find_f0(tb, w2, B=20, rounds=6):
    Uv = tb.U(np.sqrt(Xg[1:]), w2)
    neg = np.where(Uv < 0)[0]
    if len(neg) == 0:
        return None
    f1 = np.sqrt(Xg[1:][neg[0]])
    fdip = np.sqrt(Xg[1:][neg[0] + np.argmin(Uv[neg[0]:neg[-1]+1])])
    lo, hi = f1*(1+1e-9), fdip
    for _ in range(rounds):
        cand = np.linspace(lo, hi, B)
        tags = rk4_batch(tb, cand, w2)
        if tags[0] != 1:
            return None                              # bracket failure
        k = np.argmax(tags != 1) if (tags != 1).any() else B-1
        lo, hi = cand[k-1], cand[k]
    return lo

def profile_obs(tb, w):
    w2 = w*w
    f0 = find_f0(tb, w2)
    if f0 is None:
        return None
    f, fp, r = f0, 0.0, 1e-9
    rs, fs, fps = [r], [f], [fp]
    nmax = int(120.0/DR)
    for _ in range(nmax):
        def deriv(rr, ff, ffp):
            return ffp, (m2 - w2)*ff + tb.Gi(ff*ff)*ff - (2.0/rr)*ffp
        k1 = deriv(r, f, fp)
        k2 = deriv(r+DR/2, f+DR/2*k1[0], fp+DR/2*k1[1])
        k3 = deriv(r+DR/2, f+DR/2*k2[0], fp+DR/2*k2[1])
        k4 = deriv(r+DR,   f+DR*k3[0],   fp+DR*k3[1])
        f  += DR/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        fp += DR/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        r += DR
        rs.append(r); fs.append(f); fps.append(fp)
        if f < 0 or f > 1.05*f0 or (fp > 0 and f < 0.99*f0):
            break
    rs = np.array(rs); fs = np.clip(np.array(fs), 0, None); fps = np.array(fps)
    r2 = rs*rs
    I2 = np.trapz(4*np.pi*r2*fs**2, rs)
    Egr = np.trapz(4*np.pi*r2*1.5*fps**2, rs)
    Epot = np.trapz(4*np.pi*r2*(1.5*m2*fs**2 + 3*tb.Upi(fs**2)), rs)
    Ekin = 1.5*w2*I2
    keff2 = m2 + tb.G[0] - w2
    rc, fc = rs[-1], fs[-1]
    tail2 = 4*np.pi*fc**2*rc**2/(2*np.sqrt(max(keff2, 1e-6))) if fc > 0 else 0.0
    I2 += tail2; Ekin += 1.5*w2*tail2; Epot += 1.5*m2*tail2
    Q = 3*w*I2
    E = Ekin + Egr + Epot
    return dict(w=w, f0=f0, Q=Q, E=E, rs=rs, fs=fs, I2=I2, tailfrac=tail2/I2)

def omega_at_Q(tb, Qtgt=300.0, w_guess=1.410):
    """invert Q(w) = Qtgt on the dQ/dw<0 branch (secant on log Q)."""
    wlo = tb.wmin + 0.015
    w1 = max(w_guess, wlo + 1e-3)
    o1 = profile_obs(tb, w1)
    if o1 is None:
        w1 = wlo + 0.05; o1 = profile_obs(tb, w1)
    if o1 is None:
        return None, None
    w2_ = w1 + (0.008 if o1['Q'] > Qtgt else -0.008)
    o2 = profile_obs(tb, w2_)
    if o2 is None:
        return None, None
    for _ in range(24):
        l1, l2 = np.log(o1['Q']/Qtgt), np.log(o2['Q']/Qtgt)
        if abs(l2) < 2e-4:
            return w2_, o2
        if l2 == l1:
            break
        wn = w2_ - l2*(w2_ - w1)/(l2 - l1)
        wn = min(max(wn, wlo), 1.499)
        w1, o1 = w2_, o2
        w2_, o2 = wn, profile_obs(tb, wn)
        if o2 is None:
            wn = 0.5*(w1 + wlo); w2_, o2 = wn, profile_obs(tb, wn)
    return w2_, o2

# ================================================================ CHECK A
print("=== CHECK A: exact O(n) Gaussian slope c(X) vs resummed M2 numeric ===")
print("c(X) = 2 mu X (2t^2-6t+1)/(1+t)^4, t=kappa X^3 ; blue band where 2t^2-6t+1<0:")
t1, t2 = (3-np.sqrt(7))/2, (3+np.sqrt(7))/2
X1, X2 = (t1/kap)**(1/3.), (t2/kap)**(1/3.)
print(f"  t in ({t1:.4f},{t2:.4f})  =>  X=f^2 in ({X1:.4f},{X2:.4f}), f in ({np.sqrt(X1):.4f},{np.sqrt(X2):.4f})")
ok = True
nt = 1e-5
Xs = np.array([0.05, 0.15, 0.25, 0.342, 0.410, 0.55])
Gt = G_M2_table(np.concatenate([[0.0], Xs]), nt, use_enh=False, H=8)[1:]
slope = (Gt - g(Xs))/nt
for X, sl in zip(Xs, slope):
    ce = c_exact(X)
    rel = abs(sl - ce)/max(abs(ce), 1e-3)
    ok &= rel < 2e-2
    print(f"  X={X:.3f}: M2 slope {sl:+9.4f}  c_exact {ce:+9.4f}  g'(FRAME) {gp(X):+9.4f}  rel {rel:.1e}")
print("CHECK A:", "PASS" if ok else "FAIL")
print(f"  at ball core X=0.410: c={c_exact(0.410):+.4f} (RED) vs FRAME g'={gp(0.410):+.4f} (blue) — "
      "the FRAME blue-core term is cancelled by the Wick cross terms")

# ================================================================ CHECK B/C
print("\n=== CHECK B/C: vacuum branch and window edge ===")
tb0 = make_tables('vac')
print(f"omega_min(0) = {tb0.wmin:.6f}  (closed form 1.308700)",
      "PASS" if abs(tb0.wmin - 1.308700) < 2e-4 else "FAIL")
okB = True
for wref, Qref in ((1.39, 482.2), (1.42, 209.6)):
    o = profile_obs(tb0, wref)
    rel = abs(o['Q'] - Qref)/Qref
    okB &= rel < 0.03
    print(f"  w={wref}: Q={o['Q']:.1f} (v66 {Qref}) E={o['E']:.1f} f0={o['f0']:.4f} rel {rel:.1e}")
print("CHECK B:", "PASS" if okB else "FAIL")
w0_300, o0_300 = omega_at_Q(tb0, 300.0)
print(f"vacuum branch omega(Q=300) = {w0_300:.4f}  [measured eta=0.5 control: 1.4142; "
      f"offset {1.4142-w0_300:+.4f} = eta-dressing systematic, differenced out below]")

# ================================================================ dressed branches
print("\n=== Dressed branches (M2 resummed, in-core enhancement) ===")
meas = {0.05: +0.012, 0.20: -0.025, 0.80: -0.120}
rows = []
for e in (0.05, 0.20, 0.80):
    tb = make_tables('M2', e=e, H=8)
    we, oe = omega_at_Q(tb, 300.0, w_guess=w0_300)
    dw = we - w0_300
    rows.append((e, n0_per_e*e, tb.wmin, we, dw))
    print(f"  e={e:4.2f} n={n0_per_e*e:.5f}: omega_min={tb.wmin:.5f} ({tb.wmin-tb0.wmin:+.5f})  "
          f"omega(Q=300)={we:.5f}  d_omega={dw:+.5f}   [measured {meas[e]:+.3f}]")

# M1 comparison at one rung (linearized vs resummed)
tbM1 = make_tables('M1', e=0.80)
wM1, _ = omega_at_Q(tbM1, 300.0, w_guess=w0_300)
print(f"  (M1 first-order check at e=0.8: d_omega={wM1-w0_300:+.5f}, "
      f"omega_min={tbM1.wmin:.5f} — O(n^2) size on shift: {abs(wM1-rows[-1][3]):.1e})")

# ================================================================ linear response
print("\n=== Matched-Q linear response (M1 analytic, small n) vs energy envelope ===")
for n in (5e-4, 1e-3, 2e-3):
    tb = make_tables('M1', n_override=n, use_enh=True)
    w, _ = omega_at_Q(tb, 300.0, w_guess=w0_300)
    print(f"  M1 n={n:g}: d_omega = {w-w0_300:+.6f}   slope d(dw)/dn = {(w-w0_300)/n:+.3f}")
# energy-envelope route (NON-variational averaged EOM => expected to disagree):
Yg = np.linspace(0, 1.2, 24001); dYg = Yg[1]-Yg[0]
CEn = cumint(c_exact(Yg)*Enh(Yg), Yg)
env = []
for w in (1.401, 1.403, 1.405, 1.407, 1.409):
    o = profile_obs(tb0, w)
    ii = np.minimum((o['fs']**2/dYg).astype(int), len(Yg)-2)
    dE = 0.5*np.trapz(4*np.pi*o['rs']**2*CEn[ii], o['rs'])
    env.append((o['Q'], dE))
env = np.array(env)
print(f"  energy-envelope d(dE)/dQ near Q=300 = {np.polyfit(env[:,0],env[:,1],1)[0]:+.4f} per n")
print("  (mismatch vs full solve = the Wick-averaged EOM is NOT variational: the")
print("   coherent-projection cross term in c(X) has no counterpart in <E>)")

# ================================================================ window edge vs n
print("\n=== Window edge omega_min(n) — can ANY Gaussian n reach the measured 1.294? ===")
edge_ok = False
for n in (0.001, 0.003, 0.01, 0.03, 0.06, 0.10, 0.15, 0.20):
    tb = make_tables('M2', n_override=n, H=6)
    closed = tb.wmin2 >= m2 - 1e-4
    print(f"  n={n:5.3f}: omega_min = {tb.wmin:.5f}  ({tb.wmin-tb0.wmin:+.5f})"
          + ("   [window CLOSED — well smeared away]" if closed else ""))
    edge_ok |= tb.wmin < 1.294
print("RESULT: " + ("some n reaches 1.294" if edge_ok else
      "NO Gaussian n reaches 1.294 — the edge moves UP monotonically until the well "
      "is smeared away. The neutral-Gaussian mean-field CANNOT produce sub-window states."))

# ================================================================ max model shift vs n
print("\n=== Matched-Q shift d_omega(n) over the full Gaussian range (M2+enh, H=6) ===")
def dw_of_n(n, w_guess):
    tb = make_tables('M2', n_override=n, use_enh=True, H=6)
    w, _ = omega_at_Q(tb, 300.0, w_guess=w_guess)
    return (None, None, tb.wmin) if w is None else (w - w0_300, w, tb.wmin)
wg = w0_300
for n in (3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 4e-2, 8e-2):
    d, w, wm = dw_of_n(n, wg)
    if d is None:
        print(f"  n={n:6.4f}: branch lost (no Q=300 solution; omega_min={wm:.4f})")
        continue
    wg = w
    print(f"  n={n:6.4f}: d_omega = {d:+.6f}   omega_min = {wm:.5f}")
print("targets: +0.012 (e=0.05), -0.025 (e=0.2), -0.120 (e=0.8)")

# ================================================================ crossover scan
print("\n=== Crossover scan: d_omega(e) with vacuum-drive n=n0*e (M2+enh) ===")
es = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2]
prev = None
wg = w0_300
for e in es:
    tb = make_tables('M2', e=e, H=6)
    w, _ = omega_at_Q(tb, 300.0, w_guess=wg)
    if w is None:
        print(f"  e={e:4.2f}: branch lost"); continue
    wg = w
    d = w - w0_300
    print(f"  e={e:4.2f}: d_omega = {d:+.6f}   omega_min = {tb.wmin:.5f}")
    if prev is not None and prev[1]*d < 0:
        ec = prev[0] + (e - prev[0])*abs(prev[1])/(abs(prev[1]) + abs(d))
        print(f"    -> sign change: crossover e* ~ {ec:.3f}")
    prev = (e, d)

# ================================================================ coherent alternatives
print("\n=== Coherent (non-Gaussian) alternatives — required scales [estimate] ===")
print("(a) +/-omega anti-charged admixture: apparent rotation rate of")
print("    f+ e^{+iwt} + f- e^{-iwt} is w_app = w (f+^2 - f-^2)/(f+^2 + f-^2);")
print("    contamination needed for each measured red shift (w=1.414):")
for e, dwm in ((0.2, -0.025), (0.8, -0.120)):
    x = -dwm/(2*1.4142 + dwm)         # f-^2/f+^2 from w_app = w(1-x)/(1+x)
    print(f"      e={e}: f-^2/f+^2 = {x*100:.2f} %")
print("(b) driven (pumped) state below the window: existence deficit at omega=1.294:")
print(f"    omega_min^2 - omega^2 = {1.308700**2 - 1.294**2:+.4f}  — a coherent co-rotating")
print("    drive of this size in U_eff (pump power, not a free eigenstate) is required;")
print("    free dressed branches cannot reach it (see edge scan above).")
print("(c) Mathieu n=1 sum-beat pumping (THETA_DYNAMICS WP-C): net growth at the rungs")
print("    +0.017/+0.074/+0.21 per t.u. — rank-orders with the measured |red| shifts.")

print(f"\n[total runtime {time.time()-t_start:.0f} s]")

print(f"\n[total runtime {time.time()-t_start:.0f} s]")
