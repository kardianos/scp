#!/usr/bin/env python3
"""G8: Door B sketch — nonlocal element with preferred foliation.

Under Door B we keep (R) and (MI), drop (L): outcomes may depend on the
distant setting, A = A(a, b, λ) or through a nonlocal mediating field.

de Broglie double solution / Bohm: guidance equation is nonlocal in
configuration space; preferred foliation required for multi-time consistency.

SCP contact (from THEORY_v86 A2, not re-measured here):
  lattice preferred frame at 1–5% group-velocity level  [M, v86]
  Under Door B this is *required structure*, not a pure defect.

This worksheet does NOT claim the current SCP kernel violates Bell —
the kernel's matter+gauge evolution is *local* (link-covariant derivatives,
local potential). Door B requires *additional* nonlocal structure.

Explicit toy model (Bell 1964 style nonlocal det. model):
  λ ~ Unif[0, 2π)
  A(a, b, λ) = sgn( cos(λ − a) )           # still local in form
  B(a, b, λ) = −sgn( cos(λ − a) ) * f     # wait — true nonlocal:

Standard explicit nonlocal model achieving singlet:
  A(a, λ) = sgn(a · λ)          # local
  B(b, a, λ) = −sgn( b · λ_eff(a) )  with nonlocal dependence
or simply the PR-style:
  A(a,λ)=λ, B(b,a,λ) = λ * g(a,b) with g the desired product.

Clean nonlocal model for E = −cos(a−b):
  λ ∈ [0, 2π)
  A(a, λ) = +1 if λ ∈ [0, π) else −1     # ignore a? need a-dependence

Better known: use λ = (λ1, λ2) and
  A = sgn(sin(λ1 − a)), B = sgn(sin(λ1 − a)) * (−cos(a−b))/|...|  — messy.

Simplest exact nonlocal deterministic model for quantum CHSH:
  λ ~ Unif{±1}
  A(a, λ) = λ                    # independent of a (trivial local)
  B(b, a, λ) = λ * t(a,b)        # NONLOCAL: depends on a
  where t(a,b) ∈ {±1} is the target sign of E, and we mix for magnitude.

For continuous E = −cos(φ):
  Need stochastic or richer λ. With λ~Unif[−1,1]:
  A = sgn(λ), B = sgn( −cos(φ) − λ ) or use the Gisin nonlocal box...

Actually the cleanest:
  Keep local responses in appearance but let λ-distribution be prepared
  with knowledge of both settings AFTER they are chosen — that's Door A.
  Door B: settings free (MI), λ prior fixed, but B = B(b, a, λ).

Model NB1 (exact quantum correlators, nonlocal):
  λ ~ Unif[0,1]
  A(a, λ) = +1 if λ < 0.5 else −1
  B(b, a, λ) = A(a,λ) * sgn_target if |E|=1; for |E|=c use threshold:

  Let u ~ Unif[0,1] independent (part of λ).
  A = +1
  B = +1 if u < (1+E_target)/2 else −1
  Then E = E_target, but A is constant — marginals wrong (always A=+1).

Proper marginals ±1 fair:
  λ = (α, u) with α∈{±1} fair, u∈[0,1]
  A(a,λ) = α
  B(b,a,λ) = α if u < (1+E(a,b))/2 else −α
  ⇒ A B = +1 with prob (1+E)/2, so ⟨AB⟩=E.
  B depends on a: NONLOCAL. λ independent of (a,b): MI holds.
  Deterministic given λ: R holds.
  L fails.

CHSH: plug E = quantum targets → |S|=2√2.
No upper bound from this construction alone → can set E=PR targets → S=4.
So Door B unrestricted also overshoots; need a nonlocality principle
(Tsirelson from operator algebras / no-signalling + something).

No-signalling check for NB1:
  P(B=+1|a,b) = P(α B' = +1) ... B = α with prob (1+E)/2, else −α
  P(B=+1) = 1/2 always if α fair. Good.
  P(A=+1)=1/2. Good.
  So NS holds while L (parameter independence of the *response function*) fails
  — wait: parameter independence says P(A|a,b,λ)=P(A|a,λ). Here A doesn't
  depend on b, but B depends on a. So outcome-independence / parameter
  independence for Bob fails.

Cost of Door B: preferred foliation for the nonlocal influence (who depends
on whose setting). Signature: frame-dependent causal order of the nonlocal
link; in SCP, the lattice rest frame already preferred at 1–5% [M,v86].
"""
import math
import numpy as np

C = 1/math.sqrt(2)
# Quantum CHSH targets for S = +2√2:
# E00=E01=E10=+C, E11=−C
TARGETS = {
    (0, 0): +C,
    (0, 1): +C,
    (1, 0): +C,
    (1, 1): -C,
}

def mc_nonlocal(N=200_000, seed=7):
    rng = np.random.default_rng(seed)
    # λ = (α, u); settings random
    Es = {s: [] for s in TARGETS}
    for _ in range(N):
        a = int(rng.integers(0, 2))
        b = int(rng.integers(0, 2))
        alpha = 1 if rng.random() < 0.5 else -1
        u = rng.random()
        E_t = TARGETS[(a, b)]
        A = alpha
        B = alpha if u < (1 + E_t) / 2 else -alpha
        Es[(a, b)].append(A * B)
    Em = {s: float(np.mean(v)) for s, v in Es.items()}
    S = Em[(0,0)] + Em[(0,1)] + Em[(1,0)] - Em[(1,1)]
    return Em, S

def mc_pr_nonlocal(N=100_000, seed=8):
    """Same structure with PR targets → S=4 (failure mode)."""
    PR = {(0,0): +1, (0,1): +1, (1,0): +1, (1,1): -1}
    rng = np.random.default_rng(seed)
    Es = {s: [] for s in PR}
    for _ in range(N):
        a = int(rng.integers(0, 2)); b = int(rng.integers(0, 2))
        alpha = 1 if rng.random() < 0.5 else -1
        u = rng.random()
        E_t = PR[(a, b)]
        A = alpha
        B = alpha if u < (1 + E_t) / 2 else -alpha
        Es[(a, b)].append(A * B)
    Em = {s: float(np.mean(v)) for s, v in Es.items()}
    S = Em[(0,0)] + Em[(0,1)] + Em[(1,0)] - Em[(1,1)]
    return Em, S

if __name__ == "__main__":
    print("=== G8 Door B nonlocal toy model ===")
    Em, S = mc_nonlocal()
    print("E:", {k: round(v, 5) for k, v in Em.items()})
    print("S =", round(S, 5), "  target", round(2*math.sqrt(2), 5))
    assert abs(S - 2*math.sqrt(2)) < 0.03
    Em4, S4 = mc_pr_nonlocal()
    print("PR via same Door-B structure: S =", round(S4, 5))
    assert abs(S4 - 4) < 0.03
    print("G8 OK: Door B reaches 2√2 and can overshoot to 4 without extra principle [D/M]")
    print("Cost: nonlocal response B=B(b,a,λ); preferred foliation for 'a before b' influence.")
    print("SCP kernel as-is: local field dynamics — Door B structure is NOT present [M/context].")
    print("SCP preferred frame 1–5% g.v.: necessary for relativistic Door B, not sufficient.")
