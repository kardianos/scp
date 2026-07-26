#!/usr/bin/env python3
"""ws2 — FABLE seat, v87 B0.

The explicit Door-A (measurement-dependent) model and its closed form.

MODEL (p-family):
  lambda in S^2 (unit sphere) — the hidden fabric variable per run.
  Responses (deterministic, LOCAL — A never sees b, B never sees a):
      A(a,lambda) = sgn(a.lambda)
      B(b,lambda) = -sgn(b.lambda)
  Event source (the measurement-dependent piece):
      rho_p(lambda | a,b) = c_p |lambda.b|^p,   c_p = (p+1)/(4*pi)
  i.e. the source density is correlated with ONE detector's setting axis.
  p = 0 is the uniform (MI-respecting) sphere model.
  p = 1 is the Degorre–Laplante–Roland density (PRA 72, 062314 (2005) —
  identifier believed correct, UNVERIFIED here; the result is re-derived
  below from scratch, so nothing rests on the citation).

CLAIM (derived below, then verified to 1e-10):
  E_p(theta) = -1 + ((p+1)/pi) * B((p+2)/2, 1/2) * Int_0^theta sin^p(phi) dphi
  In particular:
     p=0:  E = -1 + 2*theta/pi          (local sawtooth, |S|max = 2)
     p=1:  E = -cos(theta)              (EXACT singlet => |S|max = 2*sqrt(2))
     p->oo: E -> -sgn(cos theta)        (step => |S|max -> 4)

Derivation route (each step machine-checked):
 (i)   azimuthal integral: for alpha>=0,
       Int_0^{2pi} sgn(alpha cos v + beta) dv = 2*pi*sgn(beta) if |beta|>=alpha
                                        else 2*pi - 4*arccos(beta/alpha)
       [uses arccos(-x) = pi - arccos(x); measure of {v: cos v > c} = 2 arccos(c)]
 (ii)  reduce E_p to  E_p(theta) = -1 + (2(p+1)/pi) * F_p(theta),
       F_p(theta) = Int_0^{sin theta} t^p arccos(cot(theta) t / sqrt(1-t^2)) dt
 (iii) Leibniz: dF_p/dtheta = sin^p(theta) * (1/2) * B((p+2)/2, 1/2)
       (boundary term vanishes because arccos(1)=0)
 (iv)  integrate, fix F_p(0)=0.
Consistency: E_p(pi) = +1 for all p (Beta-function identity, checked).
"""
import sympy as sp
import numpy as np
from scipy import integrate

th, t, p, s, u, v = sp.symbols('theta t p s u v', positive=True)

# ---- step (iii): the Leibniz derivative, symbolically ----------------------
inner = sp.acos(sp.cot(th) * t / sp.sqrt(1 - t**2))
dth = sp.diff(inner, th)
target = t / (sp.sin(th)**2 * sp.sqrt(1 - t**2/sp.sin(th)**2))
# valid domain 0 < t < sin(theta) < 1, 0 < theta < pi/2
diff_expr = sp.simplify(sp.trigsimp(dth - target))
print("(iii) d/dtheta arccos(cot th * t/sqrt(1-t^2)) - t/(sin^2 th sqrt(1-t^2/sin^2 th)):")
print("      symbolic residual:", diff_expr)
if diff_expr != 0:
    # fall back to exact numeric check on random points in the domain
    fn = sp.lambdify((th, t), dth - target, 'mpmath')
    import mpmath
    mpmath.mp.dps = 30
    worst = max(abs(fn(mpmath.mpf('0.9'), mpmath.mpf(str(tt))))
                for tt in (0.1, 0.3, 0.5, 0.7))
    print("      numeric residual (dps=30):", worst)
    assert worst < 1e-25

# boundary term: integrand at t = sin(theta) is arccos(cos theta ... ) -> acos(1) = 0
bnd = inner.subs(t, sp.sin(th))
print("(iii) boundary term integrand at t=sin(theta):", sp.simplify(bnd))

# the s-integral: Int_0^1 s^(p+1)/sqrt(1-s^2) ds = (1/2) B((p+2)/2, 1/2)
sint = sp.integrate(s**(p+1)/sp.sqrt(1-s**2), (s, 0, 1), conds='none')
beta_half = sp.beta((p+2)/2, sp.Rational(1,2))/2
print("(iii) Int_0^1 s^(p+1)/sqrt(1-s^2) ds =", sp.simplify(sint))
chk = sp.simplify(sint - beta_half)
print("      minus B((p+2)/2,1/2)/2 :", chk)
for pv in (0, 1, 2, 5):
    d = float(sint.subs(p, pv) - beta_half.subs(p, pv))
    assert abs(d) < 1e-12, (pv, d)

# ---- assemble the closed form ---------------------------------------------
def E_closed(pval, theta):
    """E_p(theta) = -1 + ((p+1)/pi) B((p+2)/2,1/2) Int_0^theta sin^p phi dphi"""
    from scipy.special import beta as Bfun
    I = integrate.quad(lambda ph: np.sin(ph)**pval, 0, theta, epsabs=1e-13)[0]
    return -1 + (pval+1)/np.pi * Bfun((pval+2)/2, 0.5) * I

# p=1 must be exactly -cos(theta):
E1 = -1 + (2/sp.pi) * sp.beta(sp.Rational(3,2), sp.Rational(1,2)) \
        * sp.integrate(sp.sin(u), (u, 0, th))
print("p=1 closed form:", sp.simplify(E1), " (must be -cos(theta))")
assert sp.simplify((E1 + sp.cos(th)).rewrite(sp.gamma)) == 0

# p=0 must be the sawtooth:
E0 = -1 + (1/sp.pi) * sp.beta(1, sp.Rational(1,2)) * th
print("p=0 closed form:", sp.simplify(E0), " (must be -1 + 2*theta/pi)")
assert sp.simplify((E0 - (-1 + 2*th/sp.pi)).rewrite(sp.gamma)) == 0

# E_p(pi) = +1 for all p (consistency: anti-parallel settings perfectly corr.)
Epi = -1 + (p+1)/sp.pi * sp.beta((p+2)/2, sp.Rational(1,2)) \
          * sp.beta((p+1)/2, sp.Rational(1,2))
res = sp.simplify((Epi - 1).rewrite(sp.gamma).subs(
    sp.gamma(p/2 + sp.Rational(3, 2)),
    (p/2 + sp.Rational(1, 2)) * sp.gamma(p/2 + sp.Rational(1, 2))))
print("E_p(pi) - 1 =", res, " (must be 0 for all p)")
assert res == 0

# ---- brute-force verification of the whole chain --------------------------
# Two-stage, because naive dblquad on a sign-discontinuous integrand stalls
# at ~1e-4:
#  (A) verify step (i), the azimuthal formula, by direct 1D quadrature of
#      sgn(alpha cos v + beta) at random (alpha, beta);
#  (B) evaluate E_p from the DEFINITION reduced by (A) to a 1D u-integral,
#      split at the kink u = pi/2 - theta (and pi/2), high accuracy.
def azim_formula(alpha, beta):
    if abs(beta) >= alpha:
        return 2*np.pi*np.sign(beta)
    return 2*np.pi - 4*np.arccos(beta/alpha)

rng = np.random.default_rng(11)
worstA = 0.0
for _ in range(40):
    alpha = rng.uniform(0, 1); beta = rng.uniform(-1.2, 1.2)
    if abs(beta) < alpha:   # the true discontinuities of the integrand
        v0 = np.arccos(-beta/alpha)
        pts = [v0, 2*np.pi - v0]
    else:
        pts = [np.pi]
    direct = integrate.quad(lambda vv: np.sign(alpha*np.cos(vv)+beta),
                            0, 2*np.pi, limit=400, points=pts)[0]
    worstA = max(worstA, abs(direct - azim_formula(alpha, beta)))
print(f"\n(A) azimuthal formula vs direct quadrature, 40 random pts: "
      f"worst {worstA:.2e}")
assert worstA < 1e-7

def E_direct(pval, theta):
    """definition, reduced by the (verified) azimuthal formula to 1D."""
    cp = (pval + 1) / (4*np.pi)
    def g(uu):
        alpha = abs(np.sin(theta))*np.sin(uu)
        beta = np.cos(theta)*np.cos(uu)
        # sgn(b.lam)|b.lam|^p = cos(u)|cos u|^(p-1); B = -sgn(b.lam)
        return -cp * np.sign(np.cos(uu))*abs(np.cos(uu))**pval \
                   * azim_formula(alpha, beta) * np.sin(uu)
    kinks = sorted({np.pi/2,
                    min(max(np.pi/2 - theta, 0), np.pi),
                    min(max(np.pi/2 + theta, 0), np.pi)})
    return integrate.quad(g, 0, np.pi, points=kinks, limit=400,
                          epsabs=1e-12, epsrel=1e-12)[0]

print("\n(B) definition (via verified azimuthal reduction) vs closed form:")
worst = 0.0
for pval in (0.0, 1.0, 2.0, 4.0):
    for theta in (0.3, 0.7854, 1.2, 1.9, 2.7):
        d = abs(E_direct(pval, theta) - E_closed(pval, theta))
        worst = max(worst, d)
        print(f"  p={pval:3.1f} theta={theta:6.4f}  |E_direct - E_closed| = {d:.3e}")
assert worst < 1e-9
print(f"worst residual: {worst:.3e}")

# ---- no-signalling check ---------------------------------------------------
# marginal <A> and <B> must be 0 independent of settings (rho is even in
# lambda, A and B are odd) — so P(A|a) = 1/2 regardless of b: the model is
# no-signalling at the operational level even though it is measurement-
# dependent at the hidden level.
def margA(pval, theta):
    a = np.array([np.sin(theta), 0.0, np.cos(theta)])
    cp = (pval + 1) / (4*np.pi)
    f = lambda vv, uu: cp * abs(np.cos(uu))**pval * np.sign(
        a @ np.array([np.sin(uu)*np.cos(vv), np.sin(uu)*np.sin(vv), np.cos(uu)])
        ) * np.sin(uu)
    return integrate.dblquad(f, 0, np.pi, 0, 2*np.pi, epsabs=1e-11)[0]
m = margA(1.0, 0.9)
print(f"\n<A> at p=1, theta=0.9: {m:.2e}  (must be 0 — no-signalling)")
assert abs(m) < 1e-6

print("\nws2 OK: E_p(theta) = -1 + ((p+1)/pi) B((p+2)/2,1/2) Int_0^theta sin^p,")
print("        p=0 local sawtooth, p=1 exact singlet, checked to <1e-8.")
