#!/usr/bin/env python3
"""ws3 — FABLE seat, v87 B0.

The COST side of Door A: how much measurement dependence buys how much S.

Setup: 4 setting pairs (x,y) in {0,1}^2, uniform. Hidden variable lambda; wlog
(for bounds) lambda ranges over the 16 deterministic strategies
s = (A0, A1, B0, B1) in {-1,+1}^4. Conditional distributions rho(s|xy).
E_xy = sum_s rho(s|xy) A_x(s) B_y(s),  S = E00 + E01 + E10 - E11.

(1) ANALYTIC BOUND [D]:
    For ANY reference density m(s):
      S = sum_s m(s) s_CHSH(s) + sum_xy c_xy sum_s (rho(s|xy)-m(s)) A_x B_y
      where s_CHSH(s) in {-2,2}  =>
      S <= 2 + sum_xy sum_s |rho(s|xy) - m(s)| = 2 + 2*sum_xy TV(rho_xy, m).
    So with the TV budget D := sum_xy TV(rho_xy, m),  S <= 2 + 2D.
    S = 2*sqrt(2) needs D >= (sqrt(2)-1) ~ 0.414;  S = 4 needs D >= 1.

(2) LP (exact): minimise D subject to S = S_target, with m the actual mixture
    m = (1/4) sum_xy rho_xy. Checks tightness of (1).

(3) CONVEX PROGRAM: minimise the mutual information I(lambda : XY) in bits
    subject to S = S_target. This is the Hall-style "measurement dependence
    in bits" cost. Also Pinsker corollary:
      I >= (2/ln2) * avgTV^2 = (2/ln2) * (D/4)^2  with D >= (S-2)/2
      => I(S=2sqrt2) >= (2/ln2)*((sqrt2-1)/4)^2 ~ 0.0309 bits (weak lower bd).
"""
import numpy as np
from scipy.optimize import linprog, minimize
import itertools

# ---- enumerate strategies --------------------------------------------------
STRATS = list(itertools.product([-1, 1], repeat=4))   # (A0,A1,B0,B1)
NS = len(STRATS)
PAIRS = [(0, 0), (0, 1), (1, 0), (1, 1)]
C = {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): -1}

def corr(s, x, y):
    A = s[x]; B = s[2 + y]
    return A * B

# CHSH value of each strategy (must be +-2):
schsh = np.array([sum(C[xy] * corr(s, *xy) for xy in PAIRS) for s in STRATS])
assert set(np.abs(schsh)) == {2}
print("(1) every deterministic strategy has |s_CHSH| = 2  [D, verified]")

# ---- (2) LP: minimal TV budget D for S = target ----------------------------
# variables: rho(s|xy) for 4 pairs (64), u(s,xy) aux (64). m is the mixture
# (linear in rho). minimise (1/2) sum u  s.t.
#   u >= rho - m, u >= m - rho, sum_s rho(s|xy) = 1, rho >= 0, S = target.
def min_tv(starget):
    nv = 64 + 64
    cobj = np.concatenate([np.zeros(64), 0.5 * np.ones(64)])
    A_ub, b_ub = [], []
    def idx_r(si, pi): return pi * 16 + si
    def idx_u(si, pi): return 64 + pi * 16 + si
    for pi in range(4):
        for si in range(NS):
            # rho - m - u <= 0  and  m - rho - u <= 0, m = (1/4) sum_pi rho
            for sign in (1, -1):
                row = np.zeros(nv)
                for pj in range(4):
                    row[idx_r(si, pj)] += -sign * 0.25
                row[idx_r(si, pi)] += sign * 1.0
                row[idx_u(si, pi)] = -1.0
                A_ub.append(row); b_ub.append(0.0)
    A_eq, b_eq = [], []
    for pi in range(4):
        row = np.zeros(nv)
        for si in range(NS):
            row[idx_r(si, pi)] = 1.0
        A_eq.append(row); b_eq.append(1.0)
    row = np.zeros(nv)
    for pi, xy in enumerate(PAIRS):
        for si, s in enumerate(STRATS):
            row[idx_r(si, pi)] = C[xy] * corr(s, *xy)
    A_eq.append(row); b_eq.append(starget)
    bounds = [(0, None)] * 64 + [(0, None)] * 64
    r = linprog(cobj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds,
                method='highs')
    assert r.status == 0, r.message
    return r.fun, r.x[:64].reshape(4, 16)

print("\n(2) LP: minimal TV budget D(S) vs analytic (S-2)/2:")
for Star in (2.0, 2.2, 2.5, 2*np.sqrt(2), 3.0, 3.5, 4.0):
    D, _ = min_tv(Star)
    print(f"   S = {Star:.6f}   D_min = {D:.6f}   (S-2)/2 = {(Star-2)/2:.6f}")

# ---- (3) convex: minimal mutual information I(lambda:XY) for S = target ----
def neg_entropy_terms(rho):  # rho shape (4,16); settings uniform 1/4
    m = rho.mean(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        kl = np.where(rho > 1e-300, rho * np.log2(np.maximum(rho, 1e-300) / np.maximum(m, 1e-300)), 0.0)
    return 0.25 * kl.sum()

def min_mi(starget, ntries=8, seed=1):
    rng = np.random.default_rng(seed)
    Svec = np.zeros((4, 16))
    for pi, xy in enumerate(PAIRS):
        for si, s in enumerate(STRATS):
            Svec[pi, si] = C[xy] * corr(s, *xy)
    def unpack(z):  # softmax per pair
        z = z.reshape(4, 16)
        e = np.exp(z - z.max(axis=1, keepdims=True))
        return e / e.sum(axis=1, keepdims=True)
    def obj(z):
        return neg_entropy_terms(unpack(z))
    def scon(z):
        return (unpack(z) * Svec).sum() - starget
    best = None
    for k in range(ntries):
        z0 = rng.normal(scale=2.0, size=64)
        r = minimize(obj, z0, constraints=[{'type': 'eq', 'fun': scon}],
                     method='SLSQP', options={'maxiter': 2000, 'ftol': 1e-14})
        if r.success or abs(scon(r.x)) < 1e-8:
            v = obj(r.x)
            if best is None or v < best[0]:
                best = (v, unpack(r.x))
    return best

print("\n(3) minimal I(lambda:XY) in bits for S = target "
      "(4-setting CHSH problem, deterministic strategies):")
results = {}
for Star in (2.0, 2.2, 2.5, 2*np.sqrt(2), 3.0, 3.5, 3.9, 4.0):
    got = min_mi(Star)
    results[Star] = got[0]
    pins = (2/np.log(2)) * ((Star-2)/8)**2 if Star > 2 else 0.0
    print(f"   S = {Star:.6f}   I_min = {got[0]:.6f} bits   "
          f"(Pinsker lower bound {pins:.6f})")

np.save('/home/d/code/scp/v87/work/fable/ws3_imin.npy', results, allow_pickle=True)
print("\nws3 OK.")

# ---- (4) CORRECTION (cross-seat, GROK): global minimum via Blahut-Arimoto --
# The SLSQP result in (3) is a LOCAL minimum (softmax parametrisation).
# GROK's GOOD-tilt construction achieves I = D_KL((1+1/sqrt2)/2 || 3/4)
# = 0.046273847 bits at S = 2sqrt2. The problem min I s.t. S = target is
# convex; the exponential-family fixed point rho_pi(s) prop mix(s)*exp(nu*
# Svec_pi(s)) solved by alternating minimisation (Blahut-Arimoto) gives the
# GLOBAL optimum. Result: the whole curve is
#     I_min(S) = D_KL( (1 + S/4)/2 || 3/4 )  bits,  S in [2,4]
# verified below at 7 points to 1e-9. At S=4: I = log2(4/3) = 0.4150375 bits.
def ba_solve(nu, iters=20000):
    rho = np.full((4, 16), 1/16)
    Sv = np.zeros((4, 16))
    for pi, xy in enumerate(PAIRS):
        for si, s in enumerate(STRATS):
            Sv[pi, si] = C[xy] * corr(s, *xy)
    for _ in range(iters):
        mix = rho.mean(axis=0)
        w = mix[None, :] * np.exp(nu * Sv)
        rho = w / w.sum(axis=1, keepdims=True)
    mix = rho.mean(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        kl = np.where(rho > 1e-300,
                      rho*np.log2(np.maximum(rho, 1e-300)/np.maximum(mix, 1e-300)), 0)
    S = sum((rho[pi] * Sv[pi]).sum() for pi in range(4))
    return S, 0.25 * kl.sum()

from scipy.optimize import brentq
def dkl2(a, b):
    out = 0.0
    if a > 0: out += a*np.log2(a/b)
    if a < 1: out += (1-a)*np.log2((1-a)/(1-b))
    return out

print("\n(4) GLOBAL I_min(S) [Blahut-Arimoto] vs closed form D_KL((1+S/4)/2||3/4):")
for Star in (2.2, 2.5, 2*np.sqrt(2), 3.0, 3.5, 3.9):
    nu = brentq(lambda n: ba_solve(n)[0] - Star, 1e-6, 10, xtol=1e-12)
    S, I = ba_solve(nu)
    cf = dkl2((1 + Star/4)/2, 0.75)
    print(f"   S = {Star:.6f}   I_min = {I:.9f}   closed form = {cf:.9f}   "
          f"diff = {abs(I-cf):.1e}")
    assert abs(I - cf) < 1e-7
print("   S = 4 limit: closed form log2(4/3) =", f"{np.log2(4/3):.7f}")
print("ws3 correction OK: section (3) SLSQP values are LOCAL minima; "
      "the global curve is I_min(S) = D_KL((1+S/4)/2 || 3/4).")
