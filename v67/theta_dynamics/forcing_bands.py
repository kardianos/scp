#!/usr/bin/env python3
"""forcing_bands.py — WP-C: resonance-tongue tables, in-band overlap, rung
thresholds, Floquet validation of the Mathieu analysis (mathieu_forcing.mac,
17/17 PASS).

Model: ball mode q(t):  q'' + 2*g*q' + [w0^2 + eps*cos(W*t)] q = 0.
Tongues at W = 2*w0/n, n=1,2,3,...  (n=1 width eps/w0, growth eps/(4 w0);
n=2 width eps^2/(4 w0^3), growth O(eps^2)).
"""
import numpy as np

# ---------------- parameters ----------------
m2, mu, kap, eta = 2.25, -41.345, 50.0, 0.5
omega   = 1.39          # ball internal U(1) frequency (v66 recommended seed)
omega_att = 1.44        # v66 dressed attractor (E/Q ~ 1.497 parked edge)
f0      = 0.6405        # core amplitude at omega=1.39 (v66 THEORY table)
band    = (1.1, 1.7)    # running v67 bath band, massless theta: W = |k|
rungs   = [0.05, 0.2, 0.8]   # e_bath ladder of the running experiments

# measured drain (v66 FINDINGS: dE/dt=-2.4, E=692 at eta=0.5) -> amplitude damping
gamma = 0.5 * 2.4 / 692.0    # [estimate] gamma = 1.73e-3 per t.u.

# ---------------- condensate-Hessian machinery (v67 sec 2.1) ----------------
def Msym(A):
    return 4*mu*A**4*(1 - 2*kap*A**6)/(1 + kap*A**6)**3

def dMsym(A, h=1e-6):
    return (Msym(A+h) - Msym(A-h))/(2*h)

M0   = Msym(f0)                  # curvature of symmetric amplitude channel at core
Mp   = dMsym(f0)                 # dM/dA at core
Om0  = np.sqrt(M0 + 4*omega**2)  # k=0 corotating amplitude-mode frequency

print("=== core-mode parameters (omega=%.2f ball) ===" % omega)
print("f0=%.4f  M_sym(f0)=%+.4f  M'(f0)=%+.4f  Om0=sqrt(M+4w^2)=%.4f  gamma=%.2e"
      % (f0, M0, Mp, Om0, gamma))

# ---------------- epsilon(B, eta) chain [estimate] ----------------
# bath amplitude per complex component: e_bath = 3 W^2 B^2  (Maxima C4e)
# curl force (corotating):  F0 = eta * k * B = eta * W * B = eta*sqrt(e_bath/3)
# driven amplitude response at beat Ob:  |X| = F0/|Om0^2 - Ob^2|   (C4a, g<<1)
# parametric coefficient:   eps = |M'(f0)| * |X|                   (C4c chain)
def F0(e_bath):
    return eta*np.sqrt(e_bath/3.0)          # W-independent (massless)

def eps_chain(e_bath, Ob):
    X = F0(e_bath)/abs(Om0**2 - Ob**2)
    return abs(Mp)*X

# ---------------- Floquet integrator ----------------
def floquet_growth(w0, eps, W, g=0.0, nsteps=4000):
    """Largest Floquet exponent of q'' + 2g q' + (w0^2+eps cos(Wt)) q = 0."""
    T = 2*np.pi/W
    h = T/nsteps
    def deriv(t, y):
        return np.array([y[1], -(w0**2 + eps*np.cos(W*t))*y[0] - 2*g*y[1]])
    M = np.zeros((2, 2))
    for col, y0 in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        y = np.array(y0); t = 0.0
        for _ in range(nsteps):
            k1 = deriv(t, y); k2 = deriv(t+h/2, y+h/2*k1)
            k3 = deriv(t+h/2, y+h/2*k2); k4 = deriv(t+h, y+h*k3)
            y = y + h/6*(k1+2*k2+2*k3+k4); t += h
        M[:, col] = y
    ev = np.linalg.eigvals(M)
    return np.log(max(abs(ev)))/T

# ---------------- Floquet validation vs Maxima formulas ----------------
print("\n=== Floquet validation (w0=1.39) ===")
npass = nfail = 0
def chk(name, got, want, tol):
    global npass, nfail
    ok = abs(got - want) <= tol
    print("%s %-58s got=%.6g want=%.6g" % ("PASS" if ok else "FAIL", name, got, want))
    npass, nfail = npass + ok, nfail + (not ok)

w0 = omega
e_t = 0.10
chk("n=1 center growth = eps/(4 w0)",
    floquet_growth(w0, e_t, 2*w0), e_t/(4*w0), 2e-4)
chk("n=1 edge W=2w0+eps/(2w0): growth ~ 0",
    floquet_growth(w0, e_t, 2*w0 + e_t/(2*w0)), 0.0, 2e-3)
chk("just outside (1.2x half-width): growth = 0",
    floquet_growth(w0, e_t, 2*w0 + 1.2*e_t/(2*w0)), 0.0, 1e-9)
# n=2 tongue: scan W around w0, compare width to eps^2/(4 w0^3)
e_t2 = 0.45
Ws = np.linspace(w0 - 0.06, w0 + 0.03, 1801)
gr = np.array([floquet_growth(w0, e_t2, W, nsteps=2000) for W in Ws])
inside = Ws[gr > 1e-6]
width_num = inside.max() - inside.min() if len(inside) else 0.0
chk("n=2 tongue width = eps^2/(4 w0^3) (eps=0.45)",
    width_num, e_t2**2/(4*w0**3), 0.25*e_t2**2/(4*w0**3))
mu2_max = gr.max()
print("     n=2 max growth (eps=%.2f): %.4e /t.u.  [~eps^2/(16 w0^3)=%.4e]"
      % (e_t2, mu2_max, e_t2**2/(16*w0**3)))
chk("damped n=1 threshold eps_c = 4 g w0 (g=0.01)",
    floquet_growth(w0, 4*0.01*w0, 2*w0, g=0.01), 0.0, 5e-4)

# ---------------- tongue tables ----------------
print("\n=== Table 1: naive lab-frame tongues W_n = 2*w0/n  (band %.2f-%.2f) ===" % band)
for w0c, label in [(1.39, "omega_0 = 1.39 (internal U(1) rotation)"),
                   (1.44, "omega_0 = 1.44 (v66 dressed attractor)"),
                   (0.94, "omega_0 = 0.94 (v66 lattice breathing)")]:
    print(label)
    for n in range(1, 5):
        Wn = 2*w0c/n
        tag = "IN BAND" if band[0] <= Wn <= band[1] else "out"
        print("  n=%d  W_n = %.4f   %s" % (n, Wn, tag))

print("\n=== Table 2: beat-channel resonances (parametric drive = s-modulation at")
print("    Ob = |W-omega| or W+omega; tongue n: Ob = 2*w0/n; omega=1.39) ===")
for w0c, label in [(1.39, "w0=1.39"), (0.94, "w0=0.94")]:
    hits = []
    for n in range(1, 5):
        Ob = 2*w0c/n
        for Wlab, ch in [(omega + Ob, "|W-om| diff"), (omega - Ob, "|W-om| diff"),
                         (Ob - omega, "W+om sum")]:
            if Wlab > 0:
                tag = "IN BAND" if band[0] <= Wlab <= band[1] else "out"
                hits.append((n, Wlab, ch, tag))
    print(label)
    for n, Wlab, ch, tag in sorted(hits, key=lambda x: x[1]):
        print("  n=%d  W_lab = %.4f  (%s)   %s" % (n, Wlab, ch, tag))

# ---------------- rung thresholds for the running bath runs ----------------
print("\n=== Table 3: rung ladder e_bath = 0.05/0.2/0.8, in-band drive W=omega=1.39 ===")
print("eps chains [estimate]: (a) sum-beat Ob=W+omega=2.78 -> n=1 tongue on a")
print("    corotating 1.39-mode (existence [open]); (c) naive n=2 at W=1.39.")
eps_c1 = 4*gamma*omega                    # n=1 threshold
print("thresholds: n=1 eps_c = 4*g*w0 = %.3e ;  n=2 eps_c (Floquet growth=gamma):" % eps_c1)
# numeric n=2 threshold: find eps where max n=2 growth = gamma
from bisect import bisect
def n2_max_growth(eps):
    Ws2 = np.linspace(w0 - 1.2*5*eps**2/(24*w0**3) - 0.005,
                      w0 + 1.2*eps**2/(24*w0**3) + 0.005, 121)
    return max(floquet_growth(w0, eps, W, nsteps=2000) for W in Ws2)
lo, hi = 0.05, 1.2
for _ in range(30):
    mid = 0.5*(lo+hi)
    if n2_max_growth(mid) > gamma: hi = mid
    else: lo = mid
eps_c2 = 0.5*(lo+hi)
print("    eps_c2 = %.4f" % eps_c2)
print("%-8s %-9s %-9s %-22s %-22s" % ("e_bath", "F0", "eps(a)", "n=1 net growth /t.u.",
                                      "eps(c) -> n=2 verdict"))
for e in rungs:
    Ob_a = omega + omega                    # 2.78
    ea = eps_chain(e, Ob_a)
    g1 = ea/(4*omega) - gamma
    # naive channel (c): same chain magnitude, drive at W=1.39 acting as n=2
    ec = ea
    g2 = n2_max_growth(ec) - gamma if ec > 0.05 else -gamma
    v2 = "ABOVE thr (net %.1e)" % g2 if ec > eps_c2 else "below thr (eps<%.2f)" % eps_c2
    print("%-8.2f %-9.4f %-9.4f %-22s %s" % (e, F0(e), ea,
          "%.3e %s" % (g1, "GROW" if g1 > 0 else "damped"), "eps=%.3f %s" % (ec, v2)))

# ---------------- direct (non-parametric) forcing ----------------
print("\n=== Direct forcing / interference (item 4) ===")
print("in-band corotating beats |W-omega|: %.3f .. %.3f (soft modes only)"
      % (0.0, max(abs(band[0]-omega), abs(band[1]-omega))))
print("two-bath-mode beats |W1-W2|: 0 .. %.3f" % (band[1]-band[0]))
for w0c in (0.94, 1.39):
    reach = band[1]-band[0]
    print("  beat hits mode %.2f? %s" % (w0c, "YES" if w0c <= reach else "NO (max beat %.2f)" % reach))
print("Lorentzian response |X| = F0/sqrt((w0^2-Ob^2)^2+(2 g Ob)^2), Ob=|W-omega|:")
for W in (1.1, 1.25, 1.39, 1.55, 1.7):
    Ob = abs(W-omega)
    X = F0(0.2)/np.sqrt((0.94**2-Ob**2)**2 + (2*gamma*Ob)**2)
    print("  W=%.2f  Ob=%.2f  |X|(vs 0.94-mode, e=0.2) = %.4f" % (W, Ob, X))

# ---------------- broadband dilution [estimate] ----------------
print("\n=== Broadband dilution (running baths spread power over dW=%.1f) ===" % (band[1]-band[0]))
print("only in-tongue power pumps: eps_eff = eps_mono*sqrt(min(1, dW_tongue/dW_band))")
dWb = band[1] - band[0]
print("channel (a) n=1 (tongue width eps/w0), fixed point:")
for e in rungs:
    em = eps_chain(e, 2*omega)
    ee = em
    for _ in range(200):
        ee = em*np.sqrt(min(1.0, (ee/omega)/dWb))
    g1 = ee/(4*omega) - gamma
    print("  e_bath=%.2f  eps_mono=%.3f -> eps_eff=%.3f  net growth %.3e %s"
          % (e, em, ee, g1, "GROW" if g1 > 0 else "damped"))
print("channel (c) n=2 (tongue width eps^2/(4 w0^3)): fixed point eps_eff ="
      "\n  eps_mono*eps_eff/sqrt(4 w0^3 dW) -> collapses to 0 for eps_mono <"
      " sqrt(4 w0^3 dW) = %.3f" % np.sqrt(4*omega**3*dWb))
print("  => broadband n=2 pumping NEGLIGIBLE at all rungs; needs monochromatic drive.")

print("\nFloquet validation: %d PASS, %d FAIL" % (npass, nfail))
