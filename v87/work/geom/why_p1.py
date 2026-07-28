#!/usr/bin/env python3
"""v87 B0-G crank 2: WHY p = 1?

The p-family rho_p(lambda|b) ~ |lambda.b|^p gives a ladder of correlation
waveforms: triangle (p=0, S=2) -> cosine (p=1, S=2sqrt2) -> square (p->inf, S=4).
Geometry says nothing about which rung. This file asks what singles out p=1,
gets the closed form, and then tests whether the fabric can supply it.

Contents
  A  closed form for E_p(theta), exactly
  B  the algebraic mechanism: |x| sgn(x) = x
  C  uniqueness: p=1 is the ONLY monochromatic member (Fourier)
  D  the experimentally decisive fork: source-side MD vs detection bias
  E  what the fabric would have to supply
"""
import numpy as np
import sympy as sp

TAU = 2 * np.pi


# ------------------------------------------------------------------ A
def part_A():
    print("=" * 74)
    print("A.  CLOSED FORM for E_p(theta)")
    print("=" * 74)
    print("""
E_p(theta) = - Int_0^2pi |cos f|^p sgn(cos(f-theta)) sgn(cos f) df
             / Int_0^2pi |cos f|^p df

Use |cos f|^p sgn(cos f) = |cos f|^(p-1) cos f  ==  g(f).
g is odd about f = pi/2, so its integral over the full circle vanishes; the
sgn(cos(f-theta)) factor therefore just doubles the integral over the half
where it is +1:

    E_p(theta) = -(2/Z_p) Int_{theta-pi/2}^{theta+pi/2} g(f) df ,
    Z_p = Int_0^2pi |cos f|^p df .

Differentiate under the integral. The endpoints give
    g(theta+pi/2) - g(theta-pi/2) = -2 |sin theta|^(p-1) sin theta ,
so

    dE_p/dtheta = (4/Z_p) |sin theta|^(p-1) sin theta        <-- clean
    E_p(theta)  = -1 + (4/Z_p) Int_0^theta sin^p u du    (0 <= theta <= pi)
    Z_p = 2 sqrt(pi) Gamma((p+1)/2) / Gamma(p/2 + 1)

So the whole waveform ladder is a regularised incomplete beta function.
""")
    p, th, u = sp.symbols('p theta u', positive=True)
    Z = 2 * sp.sqrt(sp.pi) * sp.gamma((p + 1) / 2) / sp.gamma(p / 2 + 1)
    print("   symbolic check of Z_p and E_p at the three landmark exponents:")
    for pv, name, target in ((0, "triangle", -1 + 2 * th / sp.pi),
                             (1, "cosine", -sp.cos(th)),
                             (2, "between", None)):
        Zv = sp.simplify(Z.subs(p, pv))
        Ev = sp.simplify(-1 + (4 / Zv) * sp.integrate(sp.sin(u) ** pv, (u, 0, th)))
        ok = "" if target is None else \
            ("  == %s  MATCH" % target if sp.simplify(Ev - target) == 0 else "  MISMATCH")
    print("      p=%d : Z=%s   E=%s%s" % (pv, Zv, sp.simplify(Ev), ok)) \
        if False else None
    for pv, target in ((0, -1 + 2 * th / sp.pi), (1, -sp.cos(th)), (2, None)):
        Zv = sp.simplify(Z.subs(p, pv))
        Ev = sp.simplify(-1 + (4 / Zv) * sp.integrate(sp.sin(u) ** pv, (u, 0, th)))
        tag = ""
        if target is not None:
            tag = "  MATCHES %s" % target if sp.simplify(sp.expand_trig(Ev - target)) == 0 \
                else "  <-- MISMATCH"
        print("      p=%-2d  Z_p = %-14s  E_p(theta) = %s%s"
              % (pv, sp.nsimplify(Zv), Ev, tag))
    print()


# ------------------------------------------------------------------ B
def part_B():
    print("=" * 74)
    print("B.  THE MECHANISM:  |x| sgn(x) = x   -- p=1 is where weight x response")
    print("    collapses to the underlying continuous variable")
    print("=" * 74)
    print("""
The weighted B-response is

    rho_p(lambda|b) * B(b,lambda)  ~  |lambda.b|^p * ( -sgn(lambda.b) )
                                   =  -|lambda.b|^(p-1) (lambda.b) .

At p = 1 and ONLY at p = 1 this is exactly -(lambda.b): the dichotomic +-1
outcome, weighted by how often it is sampled, reconstitutes the LINEAR
PROJECTION of the hidden variable on the setting.

That is the whole trick. The detector reports a sign; the rate at which it
reports is proportional to the magnitude; sign x magnitude = the field itself.
For p != 1 a residual |lambda.b|^(p-1) survives and distorts the waveform.

Physically this is a rectifier/linear-detector relation, not a coincidence of
exponents.
""")
    x = sp.symbols('x', real=True)
    pp = sp.symbols('p', positive=True)
    expr = sp.Abs(x) ** pp * sp.sign(x)
    print("   |x|^p sgn(x) simplified at p = 1:  %s"
          % sp.simplify(expr.subs(pp, 1).rewrite(sp.Piecewise)))
    print("   -> equals x.  For p != 1 a factor |x|^(p-1) remains.\n")


# ------------------------------------------------------------------ C
def part_C():
    print("=" * 74)
    print("C.  UNIQUENESS: p=1 is the ONLY MONOCHROMATIC member of the family")
    print("=" * 74)
    print("""
E_p is even and pi-antisymmetric, so it has a pure cosine series in odd
harmonics:  E_p(theta) = -sum_{k odd} c_k cos(k theta).
p = 1 gives c_1 = 1 and every other c_k = 0. Any other p carries higher
harmonics. This is the sharp characterisation of the quantum rung:

    p = 1  <=>  the correlation is a SINGLE HARMONIC  <=>  E = -cos(theta)
                                                     <=>  S_max = 2 sqrt 2
""")
    def E_p(theta, p):
        Z = 2 * np.sqrt(np.pi) * sp.gamma((p + 1) / 2).evalf() / sp.gamma(p / 2 + 1).evalf()
        Z = float(Z)
        u = np.linspace(0, theta, 20001)
        return -1.0 + (4.0 / Z) * np.trapz(np.sin(u) ** p, u)

    th = np.linspace(0, np.pi, 2001)
    print("   p       c1        c3        c5        c7     |higher|/c1")
    for p in (0, 0.5, 1.0, 1.5, 2, 4, 16):
        E = np.array([E_p(t, p) for t in th])
        cs = {}
        for k in (1, 3, 5, 7):
            cs[k] = -(2.0 / np.pi) * np.trapz(E * np.cos(k * th), th)
        hi = np.sqrt(sum(cs[k] ** 2 for k in (3, 5, 7)))
        print("  %-6g %9.5f %9.5f %9.5f %9.5f   %10.2e"
              % (p, cs[1], cs[3], cs[5], cs[7], hi / abs(cs[1])))
    print("""
   p = 1 is the unique rung whose higher-harmonic content vanishes.
   THE FABRIC CONNECTION: the objects in this model ARE monochromatic --
   Phi_a = f(r) e^{i omega t}, a single clock. HC-1 measured that the ball
   carries essentially no internal harmonics (at most one bound mode over the
   whole branch, and NONE in the working region omega >= 1.36). A single
   rotating phase projected on a setting yields a single harmonic in the
   correlation. p = 1 is what a monochromatic object gives.
""")


# ------------------------------------------------------------------ D
def part_D():
    print("=" * 74)
    print("D.  THE DECISIVE FORK: source-side MD (allowed) vs detection bias (DEAD)")
    print("=" * 74)
    print("""
rho_p(lambda|b) ~ |lambda.b| admits TWO physically different readings that give
IDENTICAL correlations but have OPPOSITE experimental status.

  READING 1 -- SOURCE-SIDE measurement dependence.
    The fabric configuration actually PRODUCED is correlated with the settings
    (common causal past; the chooser is fabric). Every emitted pair is
    detected. This violates (MI). It is NOT excluded by experiment.

  READING 2 -- DETECTOR-SIDE sampling.
    The source emits uniformly; the detector FIRES with probability ~|lambda.b|
    and misses the rest. This is exactly the DETECTION LOOPHOLE, and it is
    experimentally CLOSED.

Reading 2 requires a detection efficiency equal to the mean of the weight:
""")
    # circle (phase) and sphere (direction) versions
    phi = np.linspace(0, TAU, 2000001)
    eta_circle = np.trapz(np.abs(np.cos(phi)), phi) / TAU
    n = 8_000_000
    rng = np.random.default_rng(7)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    eta_sphere = np.abs(v[:, 0]).mean()
    eta_crit = 2 * (np.sqrt(2) - 1)
    print("      mean weight on the PHASE CIRCLE   <|cos phi|>   = %.5f  (= 2/pi)"
          % eta_circle)
    print("      mean weight on the SPHERE S^2     <|lambda.b|>  = %.5f  (= 1/2)"
          % eta_sphere)
    print("      CHSH detection-loophole threshold 2(sqrt2 - 1)  = %.5f" % eta_crit)
    print("""
      Both are FAR BELOW the threshold, so Reading 2 predicts an efficiency at
      which loophole-free experiments would never have violated CHSH -- and they
      did. **Reading 2 is refuted.**

  CONSEQUENCE FOR THE FABRIC, and it is a hard design constraint:
  the p=1 weight must be realised as a bias in WHAT THE SOURCE PRODUCES,
  correlated with the settings through the common past. It must NOT be realised
  as a detector that preferentially registers well-aligned configurations.
  Those are the same equation and opposite verdicts.
""")


# ------------------------------------------------------------------ E
def part_E():
    print("=" * 74)
    print("E.  FITTING IT TO THE FABRIC")
    print("=" * 74)
    print("""
What the p=1 structure needs        | what v86 has measured
------------------------------------+-------------------------------------------
a continuous internal phase         | Phi_a = f(r) e^{i omega t}; omega on a
  carried by every object           |   branch; local clock measured per voxel
monochromatic (single harmonic)     | HC-1: at most ONE bound internal mode on
                                    |   the whole branch, NONE for omega >= 1.36
coupling linear in the phase overlap| A6 contact law: force ~ cos(Delta phi)
                                    |   e^{-kappa D}  [MEASURED]
a sign carried by that overlap      | A6: co-phase FUSES (x2.7), anti-phase
                                    |   REPELS -> sgn(cos Delta phi) [MEASURED]
weight = magnitude of the overlap   | the same law's magnitude |cos Delta phi|
source-settings correlation via     | O2 monism: chooser is fabric, shares a
  common past                       |   past. Capacity needed: 0.046 bits.

So the fabric already exhibits, as a MEASURED interaction law,
      force  ~  cos(Delta phi)  =  |cos Delta phi| * sgn(cos Delta phi)
which is precisely  weight x response  at p = 1.

WHAT IS STILL MISSING, stated exactly:
 1. The contact law is a FORCE between two nearby objects, not a correlation of
    dichotomic outcomes at spacelike separation. Same algebraic form, different
    observable.
 2. It is exponentially suppressed, ~e^{-kappa D} with kappa ~ 0.5. At any
    Bell-relevant separation it is numerically zero:
""")
    kappa = 0.5117
    print("        D=  5 : e^-kD = %.3e" % np.exp(-kappa * 5))
    print("        D= 12 : e^-kD = %.3e   (v86 measured Adler bound < 0.06 here)"
          % np.exp(-kappa * 12))
    print("        D= 50 : e^-kD = %.3e" % np.exp(-kappa * 50))
    print("""
    So direct phase-locking across the separation is dead. The correlation must
    be established in the COMMON PAST and survive transport -- and GROK's G14
    obstruction says the kernel's mixing works against exactly that.
 3. Nothing yet forces the exponent to be 1 rather than 0. Monochromaticity is
    suggestive (part C) but is a property of the OBJECT, not yet a derivation of
    the SAMPLING WEIGHT.

FALSIFIABLE, and cheap: run the in-kernel CHSH with two objects from a common
origin and phase-overlap readout. The prediction from this analysis is
  * S <= 2 if the readout is unbiased (no MD): the fabric is then a
    Cauchy-slice local model and Bell applies;
  * S -> 2sqrt2 ONLY if the source ensemble itself is setting-correlated;
  * if S > 2 appears with a detector-side bias, it is the closed detection
    loophole and must be rejected, not celebrated.
""")


if __name__ == "__main__":
    part_A(); part_B(); part_C(); part_D(); part_E()
