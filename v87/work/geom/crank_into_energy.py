#!/usr/bin/env python3
"""v87 B0-G crank 3: plug the CHSH crank into the fabric's energy relations.

Available, all measured this session:
    E = Q w (1 + eps)                    eps = 0.45-4.34% (shooter), 2.6-4.9% (kernel)
    w Q = 2 (E_kin + E_g)                Gauss identity, residual 5e-9
    Sigma = E - wQ = (2/3)(E_grad - E_g) virial reduction, residual 5e-8
    M = E/c^2                            N7, to 1-4 parts in 10^4
    hbar_E/Q = 1 + eps                   N6

Available from the crank:
    E_p(theta) ladder; p=1 -> -cos theta -> S = 2sqrt2
    I_min(S) = D_KL( (1+S/4) || 3/4 )  bits ; I_min(2sqrt2) = 0.046274

Questions:
  A  What is eps in closed form, and does it share algebra with the Bell numbers?
  B  Is the numerical proximity eps ~ I_min (both ~0.04) anything?
  C  What sampling measure p does the fabric's OWN dynamics generate?
  D  What does that predict for an in-kernel CHSH?
"""
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
NB = "/home/d/code/scp/v86/n_battery"


# ------------------------------------------------------------------ A
def part_A():
    print("=" * 74)
    print("A.  eps in closed form -- and does it share algebra with Bell?")
    print("=" * 74)
    print("""
Ungauged, the virial identity E_grad + 3(E_m + E_V - E_kin) = 0 gives
E_m + E_V = E_kin - E_grad/3, so

    E    = E_kin + E_grad + E_m + E_V = 2 E_kin + (2/3) E_grad
    w Q  = 2 E_kin                                  (Gauss, E_g = 0)
    eps  = E/(wQ) - 1 = E_grad / (3 E_kin)          <-- clean closed form

Gauged:  eps = (E_grad - E_g) / (3 (E_kin + E_g)).
""")
    d = np.genfromtxt(os.path.join(NB, "n1_decomp.tsv"), names=True)
    pred = (d["E_grad"] - d["E_g"]) / (3.0 * (d["E_kin"] + d["E_g"]))
    err = np.max(np.abs(pred - d["eps"]))
    print("   check against the 292 production branch points:")
    print("      max |eps_closedform - eps_measured| = %.3e   (N2 residual floor)" % err)
    print("      eps range on the branch             = %.5f .. %.5f"
          % (d["eps"].min(), d["eps"].max()))
    print("""
   THE ALGEBRA DOES NOT MATCH. Every coefficient in the fabric's energy
   relations is RATIONAL -- 2 (Gauss), 3 (virial), 2/3 (Sigma), 1/3 (eps).
   Every coefficient in the Bell ceiling is sqrt2-laden -- 2sqrt2, 1/sqrt2,
   2(sqrt2-1), (2+sqrt2)/4. There is no sqrt2 anywhere in the virial package
   and no rational-3 anywhere in the CHSH geometry. They are different
   algebraic objects and do not plug into one another.
""")
    return d


# ------------------------------------------------------------------ B
def part_B(d):
    print("=" * 74)
    print("B.  Is  eps ~ I_min  (both ~0.04) a real relation or a coincidence?")
    print("=" * 74)
    # closed form for the minimum measurement-dependence cost
    S, x = sp.symbols('S x', positive=True)
    eta = (1 + S / 4) / 1  # fraction of winning setting-strategy pairs
    def DKL_num(a, b=0.75):
        """binary KL in bits, with the 0 log 0 = 0 convention at a = 1."""
        t1 = a * np.log2(a / b) if a > 0 else 0.0
        t2 = (1 - a) * np.log2((1 - a) / (1 - b)) if a < 1 else 0.0
        return t1 + t2
    Imin = lambda Sv: DKL_num((1 + float(Sv) / 4) / 2)
    print("   GROK's closed form, re-derived here:")
    print("      I_min(S) = D_KL( (1 + S/4)/2 || 3/4 )   bits")
    for Sv, tag in ((2.0, "classical"), (2 * np.sqrt(2), "Tsirelson"), (4.0, "PR")):
        print("      S = %-9.5f -> I_min = %.6f bits   (%s)" % (Sv, Imin(Sv), tag))
    I_star = Imin(2 * np.sqrt(2))
    print("""
   Where those numbers come from:
      3/4          = 3 of the 4 CHSH terms that any local strategy can win
      (2+sqrt2)/4  = (1 + 1/sqrt2)/2, the Tsirelson win fraction
   I_min is therefore a pure CONSTANT of the CHSH inequality's combinatorics.
   eps is a PROFILE FUNCTIONAL of a soliton, depending on mu, kappa, g and omega.
   Different categories. A match would need a mechanism, not a decimal.
""")
    print("   Decisive test -- does the branch ever actually reach eps = I_min?")
    print("      I_min(2sqrt2)        = %.6f" % I_star)
    for g in sorted(set(np.round(d["g"], 6))):
        sub = d[np.abs(d["g"] - g) < 1e-9]
        hit = np.any(np.abs(sub["eps"] - I_star) < 1e-4)
        print("      g = %-5.3f : eps in [%.5f, %.5f]   reaches I_min? %s"
              % (g, sub["eps"].min(), sub["eps"].max(), "YES" if hit else "NO"))
    kern = os.path.join(NB, "n56_action.tsv")
    if os.path.exists(kern):
        # the run-name column can contain a space (multi-object flag ' *'),
        # so parse by position from the end rather than by whitespace count
        rows = [ln.split("\t") for ln in open(kern).read().strip().split("\n")[1:]]
        eb = np.array([float(r[9]) for r in rows if len(r) >= 12])
        print("      in-kernel eps_bar (45 runs): [%.5f, %.5f]  reaches I_min? %s"
              % (eb.min(), eb.max(), "YES" if np.any(np.abs(eb - I_star) < 1e-3) else "NO"))
    print("""
   VERDICT ON THE COINCIDENCE: eps is not one number -- it varies by a factor of
   ~10 along the branch and is a continuous function of omega. I_min is a single
   fixed constant. Even where the ranges overlap, an identification would have to
   explain why the fabric sits at one particular omega, and nothing does: the
   branch is continuous and every point on it is equally a solution.
   **The proximity is a coincidence. Do not build on it.**
""")


# ------------------------------------------------------------------ C
def part_C():
    print("=" * 74)
    print("C.  What sampling measure does the fabric's OWN dynamics generate?")
    print("=" * 74)
    print("""
This is the question the crank actually reaches, and the answer is a NEGATIVE
result that follows from what has already been measured.

The fabric's objects carry a uniformly advancing clock: Phi_a = f(r) e^{i w t},
dphi/dt = w = const. N5 measured the clock rate stable across 45 runs; HC-1
found no internal mode in the working region that could modulate it. So over any
observation window long compared with 2pi/w, the phase is EQUIDISTRIBUTED on the
circle.

Equidistributed phase IS the uniform measure. Uniform measure IS p = 0.
And p = 0 is the TRIANGLE rung, S = 2.
""")
    # demonstrate: sample the phase from uniform rotation, no bias
    rng = np.random.default_rng(31415)
    n = 4_000_000
    for label, w in (("perfectly uniform clock", None),):
        t = rng.uniform(0, 1000.0, n)           # uncontrolled observation times
        omega = 1.43                            # the v86 branch clock
        lam = (omega * t) % (2 * np.pi)         # phase at the moment of the run
        best = 0.0
        for ap in np.linspace(0, 2 * np.pi, 73):
            for b in np.linspace(0, 2 * np.pi, 73):
                for bp in np.linspace(0, 2 * np.pi, 73):
                    def E(a1, b1):
                        return -np.mean(np.sign(np.cos(lam - a1)) *
                                        np.sign(np.cos(lam - b1)))
                    pass
        # closed form is enough: uniform phase -> triangle wave
        def E_uniform(theta):
            return -np.mean(np.sign(np.cos(lam)) * np.sign(np.cos(lam - theta)))
        ths = np.linspace(0, np.pi, 7)
        vals = [E_uniform(t_) for t_ in ths]
        tri = [-1 + 2 * t_ / np.pi for t_ in ths]
        print("   phase sampled from a uniform clock (omega=%.2f, uncontrolled t):" % omega)
        print("      theta:      " + "  ".join("%5.0f" % np.degrees(t_) for t_ in ths))
        print("      E measured: " + "  ".join("%+5.3f" % v for v in vals))
        print("      triangle:   " + "  ".join("%+5.3f" % v for v in tri))
        print("      max deviation from the triangle wave = %.4f"
              % max(abs(v - t_) for v, t_ in zip(vals, tri)))
    print("""
   -> the fabric's free evolution generates EXACTLY the p = 0 rung.

   THE IRONY, and it is worth stating plainly: crank 2 tried to argue that the
   fabric's monochromaticity implied p = 1. The reviewer showed that inference
   was a pun. Following it through, monochromaticity implies the OPPOSITE -- a
   perfectly uniform clock gives a perfectly uniform phase measure, which is the
   CLASSICAL rung. The fabric's best-measured property is evidence AGAINST
   reaching the quantum ceiling, not for it.
""")


# ------------------------------------------------------------------ D
def part_D():
    print("=" * 74)
    print("D.  What this predicts, and what would have to change")
    print("=" * 74)
    print("""
PREDICTION (sharp, cheap, falsifiable):
   An in-kernel CHSH built from two objects of a common origin, read out by
   phase projection, with settings chosen from fresh entropy outside the light
   cone, will give |S| <= 2 and in fact S ~= 2 -- the triangle value, not merely
   the bound. This is stronger than the seats' pre-registration (which said only
   |S| <= 2): part C says the fabric should sit AT the classical rung, not
   somewhere below it.

   If the measured in-kernel S comes out at 2 within error, that CONFIRMS the
   analysis and closes the question for the present kernel.
   If it comes out below 2, the readout is lossy and the test is not yet clean.
   If it comes out above 2, something in the implementation is non-local and
   that is a bug to find, not a discovery.

WHAT WOULD HAVE TO CHANGE to leave p = 0:
   The phase measure must stop being uniform, and it must stop being uniform in
   a way CORRELATED WITH THE SETTINGS. Three routes, all requiring structure the
   kernel does not have:
     1. the clock must not advance uniformly -- i.e. the object must NOT be
        monochromatic. This is in direct tension with the measured branch, and
        would also break E = Q w (1 + eps), which assumes a single omega;
     2. the source must be prepared with a setting-correlated phase, established
        in the common past and surviving transport (G14 says mixing fights this);
     3. the readout must be non-uniform in phase -- which is the detector-side
        reading, killed in crank 2 section 4 by the efficiency bound.

   Route 1 is the interesting one, because it collides with the energy
   relations: E = Q w (1 + eps) is a MONOCHROMATIC relation. A non-monochromatic
   object has no single omega, so E = Q w stops being definable, eps stops being
   a single number, and the entire virial package (which assumes a stationary
   phase-rotating ansatz) would have to be rebuilt. That is the actual cost of
   leaving the classical rung, stated in the programme's own terms.
""")


if __name__ == "__main__":
    d = part_A()
    part_B(d)
    part_C()
    part_D()
