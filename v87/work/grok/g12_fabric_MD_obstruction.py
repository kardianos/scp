#!/usr/bin/env python3
"""G12: When does a local hyperbolic fabric yield measurement independence?

Setup (abstract, no kernel run):
  Source region S at time t0 prepares λ = functional of field in S.
  Chooser A in region A, chooser B in region B, spacelike separated from
  each other at measurement time t1, each with a causal past Past(A), Past(B).
  MI says: ρ(λ | a, b) = ρ(λ), i.e. settings independent of source variable.

In a local theory with finite signal speed c_sig:
  - If Past(A) ∩ Past(B) ∩ S is empty of variables that affect both settings
    and λ in a correlated way, MI can hold.
  - O2 says choosers are fabric: their settings a,b are functionals of the
    field in Past(A), Past(B). The source λ is a functional of the field in S.
  - If there is a common past C ⊂ Past(A) ∩ Past(B) ∩ causal_past(S) that
    influences both λ and (a,b), then generically I(Λ : A,B) > 0.

Quantitative cartoon (bits):
  Suppose common past carries a classical bit X that is copied into λ and
  into a soft bias of the choosers: P(a=0|X) = (1+ε_X)/2, etc.
  Then I(Λ:A,B) can be made equal to any value in (0, H(X)] by choosing ε.

  For CHSH we need only ~0.046 bits — a *tiny* classical correlation in the
  common past suffices *if structured correctly*.

Obstruction to derivation:
  1. Existence of some correlation is easy (almost inevitable under O2).
  2. Existence of the *right* correlation (exactly the Hall tilt for every
     setting pair, and not the PR tilt) is hard — requires a dynamical
     reason that the joint (λ, a, b) lands on the quantum manifold in
     correlator space, not the PR point or a local point.
  3. Local PDE evolution from generic Cauchy data does not know about
     "CHSH angles" or "Tsirelson". Without a selection principle tied to
     instrument geometry, the predicted E(a,b) is model-dependent junk.

Conclusion [D from locality + O2, C on genericity]:
  O2 makes strict MI *unnatural* (common past exists).
  Locality makes *arbitrary* MD *unnatural* (only common-past-mediated).
  The intersection — MD of Hall type and not PR type — is a narrow target
  that has not been derived from the Cosserat PDE.

No simulation is run (brief: do not run kernel). This is a structural argument.
"""
import math

def cartoon_I_vs_bias(eps: float) -> float:
    """Common past bit X fair; setting a biased: P(a=0|X=0)=(1+eps)/2, etc.
    Simplified: I(X:a) for binary symmetric channel-like bias.
    I = 1 - h2((1+eps)/2) roughly when X determines soft a.
    """
    def h2(p):
        if p <= 0 or p >= 1:
            return 0.0
        return -p*math.log2(p) - (1-p)*math.log2(1-p)
    # X fair, A|X: P(A=X)=(1+eps)/2
    # I(X:A) = H(A) - H(A|X) = 1 - h2((1+eps)/2)  if marginal A fair
    return 1.0 - h2((1 + eps) / 2)

if __name__ == "__main__":
    print("=== G12 fabric MD obstruction (structural) ===")
    print("I needed for CHSH Tsirelson (g5/g10) = 0.046274 bits")
    for eps in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]:
        I = cartoon_I_vs_bias(eps)
        print(f"  soft bias eps={eps:.2f} → I(X:setting)≈{I:.4f} bits")
    # invert: eps such that I≈0.046
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if cartoon_I_vs_bias(mid) < 0.046274:
            lo = mid
        else:
            hi = mid
    print(f"  eps for I≈0.046 bits: {lo:.4f}")
    print("  → a few-percent soft correlation with a common-past bit is enough *in capacity*.")
    print("  The unsolved problem is *structure* (quantum not PR), not *capacity*.")
    print("G12 OK: O2 motivates MD existence; locality restricts channel; quantum structure not derived.")
