#!/usr/bin/env python3
"""G5: tighten mutual information I(Λ:S) for Door-A CHSH at S=2√2.

Two constructions:

(A) Shared-support mixture (analytical):
    Use only 4 deterministic strategies that form a "CHSH cycle":
      λ0: (A0,A1,B0,B1) = (+,+,+,+)   products: all +
      λ1: (+,+,-,+)                     00:- 01:+ 10:- 11:+  ... let's build carefully
    Better: use the 4 strategies that are the deterministic vertices extremal
    for the CHSH expression, and mix with a setting-dependent tilt.

(B) Numerical minimisation of I over conditionals on all 16 strategies
    subject to E_s = target_s (SLSQP / projected).

Also compute the *no-signalling* / PR distance and note that I→ larger as S→4.

Reference anchors (identifiers only; numbers recomputed here):
  Hall PRL 105, 250404 (2010); PRA 84, 022102 (2011)
  Barrett & Gisin PRL 106, 100406 (2011)
  Brief's ~0.0663 bits is for the *full singlet* (continuum of settings),
  not four-setting CHSH; CHSH-only MI can be smaller.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import minimize

# ---- encode strategies ----
def signs_from_idx(i):
    return [1 if (i >> k) & 1 else -1 for k in range(4)]  # A0,A1,B0,B1

def prod(i, a, b):
    A0, A1, B0, B1 = signs_from_idx(i)
    return (A0 if a == 0 else A1) * (B0 if b == 0 else B1)

C = 1 / math.sqrt(2)
TARGETS = [(0, 0, +C), (0, 1, +C), (1, 0, +C), (1, 1, -C)]
SETTINGS = [(0, 0), (0, 1), (1, 0), (1, 1)]


def entropy(p):
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log2(p)))


def I_from_rhos(rhos):
    # rhos: list of 4 arrays length 16
    marg = sum(rhos) / 4.0
    Hc = sum(entropy(r) for r in rhos) / 4.0
    return entropy(marg) - Hc, entropy(marg), Hc, marg


# ========== Construction A: analytical 2-strategy-per-pair with shared pool ==========
# Pool of strategies that individually give CHSH contribution ±2
# s(λ) = A0B0 + A0B1 + A1B0 - A1B1 ∈ {−2,+2}

def chsh_value(i):
    A0, A1, B0, B1 = signs_from_idx(i)
    return A0*B0 + A0*B1 + A1*B0 - A1*B1

GOOD = [i for i in range(16) if chsh_value(i) == 2]
BAD  = [i for i in range(16) if chsh_value(i) == -2]
print("Strategies with s(λ)=+2:", GOOD, "count", len(GOOD))
print("Strategies with s(λ)=-2:", BAD, "count", len(BAD))

# For correlators separately: construction using one global "good" mixture
# plus setting-dependent correction.
# Simpler analytical model achieving exact quantum correlators with small support:
#
# For each setting pair, mix the uniform distribution over the 8 strategies with
# product = sign(target) and the 8 with opposite, with weights (1±c)/2.
# That gives correct E but large I. Shrink by using single representatives.

def construction_A():
    """Per setting: mass (1+t)/2 on argmin-index with product +1 under sign(t),
    (1-t)/2 on opposite — same as g4. Report I."""
    rhos = []
    Es = []
    for a, b, t in TARGETS:
        r = np.zeros(16)
        plus = min(i for i in range(16) if prod(i, a, b) == 1)
        minus = min(i for i in range(16) if prod(i, a, b) == -1)
        # want <prod> = t: P(+)=(1+t)/2
        r[plus] = (1 + t) / 2
        r[minus] = (1 - t) / 2
        rhos.append(r)
        Es.append(sum(r[i] * prod(i, a, b) for i in range(16)))
    I, H, Hc, marg = I_from_rhos(rhos)
    S = Es[0] + Es[1] + Es[2] - Es[3]
    return dict(name="A_two_point", S=S, E=Es, I=I, H=H, Hc=Hc, support=int(np.sum(marg > 0)))


# ========== Construction B: single shared λ with setting-dependent bias on a 4-cycle ==========
def construction_B():
    """Use exactly four strategies that realise the four extremal CHSH patterns.

    λ++ : A0=A1=B0=B1=+1                     → all products +
    λ+- : need products for E01 focus, etc.

    More systematic: take four strategies L_s, one per setting pair s, such that
    prod(L_s, a, b) = sign(target_s) for that pair (and arbitrary elsewhere).
    Then for setting s: ρ = (1+|t|)/2 on L_s_agree + (1−|t|)/2 on a global
    disagree strategy D_s.

    To minimise I, pick the L_s to coincide as much as possible.
    Observation: no single λ gives all four targets' signs simultaneously
    (that would be s(λ)=4, impossible). Max is s(λ)=2, i.e. three agreements
    and one disagreement with the (+,+,+,−) pattern — wait:
    target signs of products: (+,+,+,−). That is EXACTLY s(λ)= A0B0+A0B1+A1B0−A1B1,
    and max is 2, not 4. So no deterministic λ realises all four quantum *signs
    with magnitude 1*. Quantum has magnitude 1/√2 < 1, so mixtures work.

    Construction: let U = uniform on GOOD (s= +2). Under U, what are the E's?
    """
    # average correlators under uniform on GOOD
    E_u = []
    for a, b, t in TARGETS:
        vals = [prod(i, a, b) for i in GOOD]
        E_u.append(sum(vals) / len(GOOD))
    print("  Uniform on GOOD strategies: E =", E_u, "S =", E_u[0]+E_u[1]+E_u[2]-E_u[3])

    # Under uniform on all 16, E=0.
    # Under uniform on GOOD: each correlator that appears with + in CHSH should
    # be positively biased.
    return E_u


# ========== Construction C: numerical MI minimisation ==========
def construction_C(maxiter=400):
    """Minimise I(Λ:S) over 4 distributions on 16 atoms (softmax parameterisation)."""
    # params: 4 * 16 raw logits
    rng = np.random.default_rng(0)
    x0 = rng.normal(scale=0.5, size=4 * 16)

    def unpack(x):
        rhos = []
        for k in range(4):
            logits = x[k * 16:(k + 1) * 16]
            logits = logits - np.max(logits)
            e = np.exp(logits)
            rhos.append(e / e.sum())
        return rhos

    def correlators(rhos):
        out = []
        for k, (a, b, _) in enumerate(TARGETS):
            out.append(sum(rhos[k][i] * prod(i, a, b) for i in range(16)))
        return out

    def objective(x):
        rhos = unpack(x)
        I, _, _, _ = I_from_rhos(rhos)
        return I

    def constraints():
        cons = []
        for k, (a, b, t) in enumerate(TARGETS):
            def make(k=k, a=a, b=b, t=t):
                def c(x):
                    rhos = unpack(x)
                    E = sum(rhos[k][i] * prod(i, a, b) for i in range(16))
                    return E - t
                return c
            cons.append({'type': 'eq', 'fun': make()})
        return cons

    # multi-start
    best = None
    for seed in range(12):
        rng = np.random.default_rng(seed)
        x0 = rng.normal(scale=1.0, size=4 * 16)
        # warm start toward two-point solution
        if seed == 0:
            x0 = np.full(4 * 16, -5.0)
            for k, (a, b, t) in enumerate(TARGETS):
                plus = min(i for i in range(16) if prod(i, a, b) == 1)
                minus = min(i for i in range(16) if prod(i, a, b) == -1)
                x0[k * 16 + plus] = math.log((1 + t) / 2 + 1e-6) + 5
                x0[k * 16 + minus] = math.log((1 - t) / 2 + 1e-6) + 5
        res = minimize(objective, x0, method='SLSQP', constraints=constraints(),
                       options={'maxiter': maxiter, 'ftol': 1e-12, 'disp': False})
        rhos = unpack(res.x)
        Es = correlators(rhos)
        I, H, Hc, marg = I_from_rhos(rhos)
        ok = all(abs(Es[k] - TARGETS[k][2]) < 1e-5 for k in range(4))
        S = Es[0] + Es[1] + Es[2] - Es[3]
        rec = dict(I=I, H=H, Hc=Hc, E=Es, S=S, ok=ok, support=int(np.sum(marg > 1e-6)), success=res.success)
        if ok and (best is None or I < best['I']):
            best = rec
        print(f"  seed {seed}: success={res.success} ok={ok} I={I:.6f} S={S:.6f} supp={rec['support']}")
    return best


# ========== Construction D: closed-form small-I model for CHSH ==========
def construction_D():
    """Analytical model with I = h((1+c)/2) − something... 

    Use ONE shared bit λ∈{±1} plus measurement dependence via rejection that
    is equivalent to a setting-dependent reweight.

    Model (finite, exact):
      Let the deterministic responses be the 'optimal classical' ones for CHSH:
        A0=+1, A1=+1, B0=+1, B1=−1   always  (one fixed strategy λ*)
      Then products: E00=+1, E01=−1, E10=+1, E11=+1  → S = 1−1+1−1 = 0. Bad.

      Instead use mixture of two strategies with s=+2:
        λp: products matching (+1,+1,+1,−1) as well as possible
        Find λ with (p00,p01,p10,p11) = (+,+,+,−): 
        A0B0=+1, A0B1=+1, A1B0=+1, A1B1=−1.
        From A0B0=A0B1 ⇒ B0=B1, but A1B0=−A1B1 ⇒ B0=−B1. Contradiction.
        So NO strategy realises the quantum sign pattern at magnitude 1.
        Max agreements: 3 of 4. Those are exactly the GOOD strategies (s=2).

    Take uniform mixture over GOOD. Compute E's (construction_B).
    Then tilt each conditional slightly from that common distribution toward
    the needed correlator — small KL cost.
    """
    nG = len(GOOD)
    base = np.zeros(16)
    for i in GOOD:
        base[i] = 1.0 / nG
    # E under base
    E0 = [sum(base[i] * prod(i, a, b) for i in range(16)) for a, b, _ in TARGETS]
    S0 = E0[0] + E0[1] + E0[2] - E0[3]
    print(f"  base on GOOD: E={E0}, S={S0}")

    # For each setting, reweight by factor (1 + α t_s * product) normalised (linear tilt)
    # Choose α so that E hits C.
    # Under reweight w_i ∝ base_i * (1 + α * t_hat * prod_i), with t_hat = sign target
    def tilted(a, b, t, alpha):
        w = np.zeros(16)
        for i in GOOD:
            w[i] = (1.0 / nG) * (1.0 + alpha * t * prod(i, a, b))
        # may go negative if alpha large; clip domain
        if np.any(w < -1e-15):
            return None
        w = np.maximum(w, 0)
        w /= w.sum()
        E = sum(w[i] * prod(i, a, b) for i in range(16))
        return w, E

    # For a single correlator under GOOD uniform: compute mean m and max
    # E(alpha) = (m + alpha t <prod^2>) / (1 + alpha t m) roughly
    # prod^2=1 always, so E = (m + alpha t)/(1 + alpha t m)
    # Set equal to t (target value, |t|=c): solve for alpha.
    rhos = []
    Es = []
    alphas = []
    for a, b, t in TARGETS:
        m = sum(base[i] * prod(i, a, b) for i in range(16))
        # E = (m + alpha*t*1)/(1 + alpha*t*m) = t
        # m + alpha t = t + alpha t^2 m
        # m - t = alpha t^2 m - alpha t = alpha t (t m - 1)
        # alpha = (m - t) / (t (t*m - 1)) = (t - m) / (t (1 - t m))
        alpha = (t - m) / (t * (1.0 - t * m))
        alphas.append(alpha)
        w, E = tilted(a, b, t, alpha)
        assert w is not None and abs(E - t) < 1e-9, (E, t, alpha, m)
        rhos.append(w)
        Es.append(E)
    I, H, Hc, marg = I_from_rhos(rhos)
    S = Es[0] + Es[1] + Es[2] - Es[3]
    return dict(name="D_tilt_GOOD", S=S, E=Es, I=I, H=H, Hc=Hc,
                alphas=alphas, support=int(np.sum(marg > 1e-12)))


def construction_E_binary_entropy_bound():
    """Information-theoretic lower sketch for this finite setting scenario.

    For a single pair of binary outcomes with correlation c, the minimum
    mutual information between λ and the setting pair, over local det models,
    is not forced to be large — but for the full CHSH *sign pattern* simultaneously
    under one family of conditionals, Hall-type bounds apply.

    Binary entropy h2(p) = -p log2 p - (1-p)log2(1-p).
    Note h2((1+1/√2)/2) = h2(0.853553) ≈ 0.600876 bits.
    Brief's 0.0663 ≈ 1/15 is much smaller — it is I(Λ:settings), not H(outcome).
    """
    def h2(p):
        if p <= 0 or p >= 1:
            return 0.0
        return -p*math.log2(p) - (1-p)*math.log2(1-p)
    p = (1 + C) / 2
    return dict(h2_of_corr=h2(p), p=p, brief_singlet_I=1/15)


if __name__ == "__main__":
    print("=== Construction A (two-point) ===")
    print(construction_A())
    print("=== Construction B (uniform GOOD) ===")
    construction_B()
    print("=== Construction D (tilt on GOOD) ===")
    d = construction_D()
    print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in d.items() if k != 'alphas'})
    print("  alphas:", d['alphas'])
    print("=== Construction E (entropy notes) ===")
    print(construction_E_binary_entropy_bound())
    print("=== Construction C (numerical min I) ===")
    best = construction_C()
    print("BEST:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in best.items() if k != 'E'},
          "E", best['E'])
