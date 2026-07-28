#!/usr/bin/env python3
"""G4: Door A — local deterministic model with partial measurement dependence
achieving the Tsirelson CHSH value S = 2√2, with mutual information quantified.

Construction (explicit, self-contained; not a verbatim copy of Hall's paper):

Hidden variable λ = (A0, A1, B0, B1) ∈ {±1}^4  — 16 deterministic local strategies.
Response functions (R)+(L):
  A(a, λ) = A_a   for a ∈ {0,1}
  B(b, λ) = B_b   for b ∈ {0,1}

Without MD, any ρ(λ) gives |S|≤2 (G1).
With MD: for each setting pair s=(a,b) we allow a different distribution ρ_s(λ).

Target quantum CHSH correlators (Tsirelson-optimal signs):
  E00 = E01 = E10 = +c,  E11 = −c,  with c = 1/√2
  ⇒ S = 4c = 2√2.

For each s, we realise E_s = t_s by a two-point mixture on the set of λ
with product A_a B_b = +1 vs −1, using a *fixed reference* local deterministic
strategy list so that responses stay local.

Minimal-I construction used here (analytical, not globally optimised):
  For each setting pair (a,b) with target t ∈ {+c, −c}:
    put mass (1+|t|)/2 on the single strategy that realises product = sign pattern
    maximizing agreement, and (1−|t|)/2 on one that realises the opposite product,
    choosing those two strategies from a *shared* 4-element menu so the marginal
    support of λ is small (reduces H(Λ)).

We then compute:
  I(Λ : S) = H(Λ) − Σ_s P(s) H(Λ|s)     bits,  P(s)=1/4
and compare to the brief's Hall singlet figure ≈0.0663 bits (full singlet is
a stronger requirement than four-setting CHSH alone).

Also: Monte Carlo recovery of S.
"""
from __future__ import annotations
import itertools
import math
import numpy as np

# ---- strategy encoding ----
# λ index 0..15: bits (A0,A1,B0,B1) with +1 -> bit 1, -1 -> bit 0
def decode(idx: int):
    bits = [(idx >> k) & 1 for k in range(4)]
    vals = [1 if b else -1 for b in bits]
    return dict(A0=vals[0], A1=vals[1], B0=vals[2], B1=vals[3])

def product(idx: int, a: int, b: int) -> int:
    d = decode(idx)
    A = d['A0'] if a == 0 else d['A1']
    B = d['B0'] if b == 0 else d['B1']
    return A * B

ALL = list(range(16))

# ---- targets ----
C = 1.0 / math.sqrt(2.0)
TARGETS = {
    (0, 0): +C,
    (0, 1): +C,
    (1, 0): +C,
    (1, 1): -C,
}

def entropy(ps: np.ndarray) -> float:
    p = ps[ps > 0]
    return float(-np.sum(p * np.log2(p)))

def two_point_dist(a: int, b: int, t: float) -> np.ndarray:
    """Distribution on 16 strategies with <A_a B_b> = t, minimal binary support.

    Pick one λ+ with product +1 and one λ− with product −1.
    Mass: P(+) = (1+t)/2, P(−) = (1−t)/2.
    Choice of which λ±: fix a canonical pair depending only on (a,b) so that
    the four conditionals share as much support as possible.
    """
    plus = [i for i in ALL if product(i, a, b) == +1]
    minus = [i for i in ALL if product(i, a, b) == -1]
    # Canonical: smallest index in each class (deterministic choice)
    ip = min(plus)
    im = min(minus)
    rho = np.zeros(16)
    rho[ip] = (1.0 + t) / 2.0
    rho[im] = (1.0 - t) / 2.0
    return rho

def correlator(rho: np.ndarray, a: int, b: int) -> float:
    return float(sum(rho[i] * product(i, a, b) for i in ALL))

def mutual_information(rhos: dict) -> dict:
    """I(Λ;S) with S uniform on 4 setting pairs."""
    settings = list(TARGETS.keys())
    # conditional entropies
    H_cond = 0.0
    marg = np.zeros(16)
    for s in settings:
        rho = rhos[s]
        H_cond += 0.25 * entropy(rho)
        marg += 0.25 * rho
    H_marg = entropy(marg)
    I = H_marg - H_cond
    # also average KL form check
    I_kl = 0.0
    for s in settings:
        rho = rhos[s]
        for i in range(16):
            if rho[i] > 0 and marg[i] > 0:
                I_kl += 0.25 * rho[i] * math.log2(rho[i] / marg[i])
    return dict(I_bits=I, I_kl=I_kl, H_Lambda=H_marg, H_Lambda_given_S=H_cond, marg=marg)

def main():
    rhos = {s: two_point_dist(s[0], s[1], t) for s, t in TARGETS.items()}

    print("=== G4 Door A: two-point MD model for quantum CHSH ===")
    E = {}
    for s, t in TARGETS.items():
        e = correlator(rhos[s], *s)
        E[s] = e
        print(f"  E{s} = {e:.10f}  (target {t:.10f}, err {abs(e-t):.2e})")

    S = E[(0, 0)] + E[(0, 1)] + E[(1, 0)] - E[(1, 1)]
    print(f"  S = {S:.10f}  (Tsirelson 2*sqrt(2) = {2*math.sqrt(2):.10f})")
    assert abs(S - 2 * math.sqrt(2)) < 1e-12

    mi = mutual_information(rhos)
    print(f"  H(Λ)           = {mi['H_Lambda']:.6f} bits")
    print(f"  H(Λ|S)         = {mi['H_Lambda_given_S']:.6f} bits")
    print(f"  I(Λ:S)         = {mi['I_bits']:.6f} bits")
    print(f"  I via KL check = {mi['I_kl']:.6f} bits")
    print(f"  support marg   = {int(np.sum(mi['marg'] > 0))} of 16 strategies")
    print()
    print("  NOTE: this two-point-per-setting construction is FAR from MI-optimal.")
    print("  It achieves S=2√2 with local det. responses, but I(Λ:S) is O(1) bits.")
    print("  Hall's ~0.066 bit figure is for an optimised / continuous model of the")
    print("  *full singlet*, not this crude finite menu. See g5 for a tighter CHSH MI.")

    # ---- Monte Carlo ----
    rng = np.random.default_rng(42)
    N = 200_000
    # sample setting then λ then outcomes
    counts = {s: [] for s in TARGETS}
    settings_list = list(TARGETS.keys())
    for _ in range(N):
        s = settings_list[int(rng.integers(0, 4))]
        rho = rhos[s]
        idx = int(rng.choice(16, p=rho))
        counts[s].append(product(idx, *s))
    print(f"=== Monte Carlo N={N} ===")
    Emc = {}
    for s in settings_list:
        Emc[s] = float(np.mean(counts[s]))
        se = float(np.std(counts[s], ddof=1) / math.sqrt(len(counts[s])))
        print(f"  E{s}_MC = {Emc[s]:+.5f} ± {se:.5f}  (analytic {E[s]:+.5f})")
    Smc = Emc[(0, 0)] + Emc[(0, 1)] + Emc[(1, 0)] - Emc[(1, 1)]
    # rough SE on S
    print(f"  S_MC = {Smc:.5f}  (analytic {S:.5f})")
    assert abs(Smc - S) < 0.05, Smc
    print("G4 OK: symbolic S=2√2 recovered; MC matches within sampling error [M]")
    return dict(S=S, I_bits=mi['I_bits'], Emc=Emc, Smc=Smc)

if __name__ == "__main__":
    main()
