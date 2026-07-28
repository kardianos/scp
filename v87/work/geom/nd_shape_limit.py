#!/usr/bin/env python3
"""v87 B0-G addendum: is the CHSH ceiling an INSCRIBED-SHAPE limit, and does it
change with the dimension N of the body?

Model family (continuous motion on an N-dimensional body):
    lambda uniform on the unit sphere S^{N-1} in R^N,
    responses  A = sgn(a.lambda),  B = -sgn(b.lambda)          (local, deterministic)
    sampling   rho_p(lambda|b) proportional to |lambda.b|^p     (p=0 -> no bias)

Questions turned into arithmetic:
  Q1  Does E(theta) depend on N?
  Q2  What is E(theta) for p = 0, 1, 2, ... and what waveform is it?
  Q3  What is max|S| for each, and is it 4 x (a shape factor)?
  Q4  Is the ceiling an inscribed-polygon factor? (chained Bell, n settings)
  Q5  Is the sqrt(2) the l1/l2 ratio of a unit vector -- and in what space?
"""
import numpy as np
import sympy as sp

TAU = 2 * np.pi
rng = np.random.default_rng(20260726)


# ---------------------------------------------------------------- Q1 (symbolic)
def q1_symbolic():
    print("=" * 72)
    print("Q1  Does the correlation depend on the dimension N of the body?")
    print("=" * 72)
    print("""
For lambda uniform on S^{N-1}, write the first two coordinates in polar form
(lambda_1, lambda_2) = (r cos phi, r sin phi). The marginal density of the first
two coordinates is

    dP  proportional to  (1 - r^2)^{(N-4)/2}  r dr dphi        (N >= 3)

Both responses and the weight depend ONLY on phi:

    sgn(a.lambda) sgn(b.lambda)  depends only on phi (a, b span the 1-2 plane)
    |lambda.b|^p = r^p |cos phi|^p

so the correlation is

              Int r^p (1-r^2)^{(N-4)/2} r dr  x  Int_phi |cos phi|^p sgn.sgn dphi
    E(theta) = --------------------------------------------------------------
              Int r^p (1-r^2)^{(N-4)/2} r dr  x  Int_phi |cos phi|^p dphi

THE RADIAL INTEGRAL IS IDENTICAL IN NUMERATOR AND DENOMINATOR AND CANCELS.
Therefore E(theta) is INDEPENDENT OF N, exactly, for every p.
""")
    # verify the cancellation symbolically for a few N
    r, p = sp.symbols('r p', positive=True)
    for N in (3, 4, 5, 10):
        integ = sp.integrate(r**p * (1 - r**2)**sp.Rational(N - 4, 2) * r, (r, 0, 1))
        print("   N=%2d : radial factor = %s   (identical top and bottom -> cancels)"
              % (N, sp.simplify(integ)))
    print()


def q1_numeric():
    print("   Monte-Carlo check (p=1, theta=60 deg; exact answer -cos60 = -0.5):")
    theta = np.pi / 3
    for N in (2, 3, 4, 5, 10, 25):
        n = 4_000_000
        lam = rng.normal(size=(n, N))
        lam /= np.linalg.norm(lam, axis=1, keepdims=True)
        b = np.zeros(N); b[0] = 1.0
        a = np.zeros(N); a[0] = np.cos(theta); a[1] = np.sin(theta)
        w = np.abs(lam @ b) ** 1          # p = 1
        val = -np.sign(lam @ a) * np.sign(lam @ b)
        E = np.average(val, weights=w)
        print("      N=%2d : E = %+.5f   (err %+.5f)" % (N, E, E + np.cos(theta)))
    print()


# ---------------------------------------------------------------- Q2/Q3
def E_of_theta(theta, p, ngrid=200001):
    """E(theta) = -Int |cos phi|^p sgn(cos(phi-theta)) sgn(cos phi) dphi / Int |cos phi|^p"""
    phi = np.linspace(0, TAU, ngrid)
    w = np.abs(np.cos(phi)) ** p
    f = np.sign(np.cos(phi - theta)) * np.sign(np.cos(phi))
    return -np.trapz(w * f, phi) / np.trapz(w, phi)


def chsh_max(p, coarse=121, refine=True):
    """max |S| over all four coplanar settings.

    BUG FIXED: the settings must range over the FULL circle. A first version
    searched b, b' in [0, pi] only, which EXCLUDES the CHSH optimum (which needs
    b' at -pi/4) and reported 2.50 instead of 2sqrt2 for the exact cosine. Any
    restricted box is a search artifact, not a bound.
    a is fixed to 0 WLOG: E depends only on angle differences."""
    Ecache = {}

    def E(t):
        k = round(t % TAU, 9)
        if k not in Ecache:
            Ecache[k] = E_of_theta(t, p)
        return Ecache[k]

    def scan(centers, half, npts):
        best, arg = 0.0, None
        for ap in np.linspace(centers[0] - half, centers[0] + half, npts):
            for b in np.linspace(centers[1] - half, centers[1] + half, npts):
                for bp in np.linspace(centers[2] - half, centers[2] + half, npts):
                    S = E(-b) + E(-bp) + E(ap - b) - E(ap - bp)
                    if abs(S) > best:
                        best, arg = abs(S), (0.0, ap, b, bp)
        return best, arg

    best, arg = scan((np.pi, np.pi, np.pi), np.pi, coarse)
    if refine:
        for half, npts in ((0.10, 21), (0.01, 21)):
            b2, a2 = scan(arg[1:], half, npts)
            if b2 > best:
                best, arg = b2, a2
    return best, arg


def q23():
    print("=" * 72)
    print("Q2/Q3  What waveform is E(theta), and what ceiling does it give?")
    print("=" * 72)
    ths = np.linspace(0, np.pi, 7)
    print("  p      E(0)    E(30)   E(60)   E(90)  E(120)  E(150)  E(180)   waveform")
    for p in (0, 0.5, 1, 2, 4, 8, 64):
        vals = [E_of_theta(t, p) for t in ths]
        # identify shape by comparing with the three canonical candidates
        tri = [-1 + 2 * t / np.pi for t in ths]
        cos = [-np.cos(t) for t in ths]
        sq = [-np.sign(np.cos(t)) if abs(np.cos(t)) > 1e-12 else 0.0 for t in ths]
        d = {"triangle": np.max(np.abs(np.array(vals) - tri)),
             "cosine": np.max(np.abs(np.array(vals) - cos)),
             "square": np.max(np.abs(np.array(vals) - sq))}
        tag = min(d, key=d.get)
        lbl = "%s (max dev %.1e)" % (tag, d[tag]) if d[tag] < 2e-3 else "between"
        print("  %-5g " % p + " ".join("%+.4f" % v for v in vals) + "   " + lbl)
    print()
    print("  CHSH ceiling for each p (FULL-circle search, refined):")
    print("   p        max|S|     max|S|/4    note")
    for p in (0, 0.5, 1, 2, 4, 8, 64):
        s, arg = chsh_max(p)
        note = ""
        if abs(s - 2) < 3e-3: note = "= 2         classical  (triangle wave)"
        if abs(s - 2 * np.sqrt(2)) < 3e-3: note = "= 2sqrt2    TSIRELSON  (cosine wave)"
        if abs(s - 4) < 3e-2: note = "-> 4        PR max     (square wave)"
        print("   %-5g   %8.5f   %8.5f    %s" % (p, s, s / 4, note))
    print()
    print("   exact targets:  2 = 2.00000   2sqrt2 = %.5f   4 = 4.00000"
          % (2 * np.sqrt(2)))
    print()


# ---------------------------------------------------------------- Q4
def q4():
    print("=" * 72)
    print("Q4  Is the ceiling an INSCRIBED-POLYGON factor?  (chained Bell, n settings)")
    print("=" * 72)
    print("""
The chained Bell inequality with n settings per side has
    local bound      L(n) = 2n - 2
    quantum bound    Q(n) = 2n cos(pi/(2n))
n = 2 reproduces CHSH: L = 2, Q = 4 cos(pi/4) = 2 sqrt2.

cos(pi/m) is exactly the INSCRIBED-POLYGON EFFICIENCY: the projection you keep
when you must use one of m equally spaced directions instead of the exact one.
So Q(n) is '2n x polygon efficiency', and the circle is the m -> infinity limit.
""")
    print("   n     L=2n-2      Q=2n cos(pi/2n)    Q/L      cos(pi/2n)")
    for n in (2, 3, 4, 5, 8, 16, 64, 1024):
        L = 2 * n - 2
        Q = 2 * n * np.cos(np.pi / (2 * n))
        print("   %-5d %8.4f    %14.6f  %8.5f   %10.7f"
              % (n, L, Q, Q / L, np.cos(np.pi / (2 * n))))
    print("""
   -> the RATIO Q/L is maximal at n = 2 (= sqrt2) and DECREASES with n.
      More settings do not buy more relative non-classicality; CHSH is the
      extremal case. The circle limit is classical-looking.
""")


# ---------------------------------------------------------------- Q5
def q5():
    print("=" * 72)
    print("Q5  Where does the sqrt(2) live -- body dimension, or setting count?")
    print("=" * 72)
    print("""
CHSH:  S = a.(b+b') + a'.(b-b')  <=  ||b+b'|| + ||b-b'||.
Put ||b+b'|| = 2 cos B, ||b-b'|| = 2 sin B:
    S <= 2(cos B + sin B) = 2 sqrt2 sin(B + pi/4)  <=  2 sqrt2.

cos B + sin B is the l1 norm of a UNIT l2 vector in the plane. The general fact:

    max { ||x||_1 : ||x||_2 = 1, x in R^d }  =  sqrt(d)

attained at the fully diagonal direction. d here is the number of INDEPENDENT
TERMS being combined (2), NOT the dimension N of the body lambda lives on.
Geometrically this is the cross-polytope / sphere ratio: the l1 ball inscribed
in the l2 ball touches it only on the axes, and the worst case is the diagonal.
""")
    d = sp.symbols('d', positive=True, integer=True)
    print("   max l1 of a unit l2 vector in R^d = sqrt(d):")
    for dd in (1, 2, 3, 4, 8):
        x = np.ones(dd) / np.sqrt(dd)
        print("      d=%d :  l1 = %.6f   sqrt(d) = %.6f" % (dd, x.sum(), np.sqrt(dd)))
    print("""
   CHSH combines d = 2 such terms, so the factor is sqrt2 and S <= 2 sqrt2.
   A hypothetical inequality combining d genuinely independent terms of this
   shape would carry sqrt(d) -- but that changes the INEQUALITY, not the
   physics, and the chained family (Q4) shows the ratio does not improve.
""")


if __name__ == "__main__":
    q1_symbolic()
    q1_numeric()
    q23()
    q4()
    q5()
