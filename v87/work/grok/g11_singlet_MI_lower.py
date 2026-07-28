#!/usr/bin/env python3
"""G11: Toward Hall's full-singlet mutual information (continuum of settings).

For the *full* singlet E(a,b)=−a·b for all directions, the MD cost is larger
than the 4-setting CHSH value I_*=0.046 bits.

Construction (circle / planar spins, continuum of angles):
  Responses: A(α,λ)=sgn cos(λ−α), B(β,λ)=−sgn cos(λ−β)
  For each pair (α,β) with θ=|α−β|, tilt level sets so ⟨AB⟩=−cos θ
  (same as g7).

Settings measure: take α,β independent uniform on [0,π) (planar projective),
or on a finite ε-net and extrapolate.

I(Λ : Α,Β) = H(Λ) − H(Λ|Α,Β) with continuous Λ on the circle.
Equivalently average D_KL(ρ_{αβ} || ρ_marginal).

We compute:
  (1) mean D_KL(ρ_θ || uniform) over θ-distribution induced by random angles
  (2) discretised I with N_ang setting bins
  (3) comparison to CHSH-only 0.046 and brief's ~0.066

Note: this is still a *model* MI, not a model-independent lower bound.
Hall's optimisation may use a different response family and get a smaller I
for the same correlators. Our number is an *achievable upper bound on the
minimum* (we achieve singlet planar correlators at this I), hence
  I_min ≤ I_our.
A matching lower bound needs a different argument (Hall PRA 2011 style).
"""
from __future__ import annotations
import math
import numpy as np

def dkl_levelset(theta: float) -> float:
    """D_KL(ρ_θ || Unif) in bits for tilt to E=-cos θ; θ in (0,π)."""
    if theta < 1e-9 or abs(theta - math.pi) < 1e-9:
        return 0.0
    Eq = -math.cos(theta)
    pp = (1.0 + Eq) / 2.0
    pm = 1.0 - pp
    mup = theta / math.pi
    mum = 1.0 - mup
    d = 0.0
    if pp > 0 and mup > 0:
        d += pp * math.log2(pp / mup)
    if pm > 0 and mum > 0:
        d += pm * math.log2(pm / mum)
    return d

def mean_dkl_random_angles(n_samples=200_000, seed=0):
    """α,β ~ Unif[0,π); θ = min(|α-β|, π-|α-β|) in [0, π/2]? 
    For E=-cos(α-β) with signed angle difference reduced to [0,π]."""
    rng = np.random.default_rng(seed)
    a = rng.uniform(0, math.pi, n_samples)  # projective planar
    b = rng.uniform(0, math.pi, n_samples)
    th = np.abs(a - b)
    th = np.minimum(th, math.pi - th)  # [0, π/2]
    # For hemisphere model, use full [0,π] separation on circle of period 2π
    # with directions defined mod π for spin-1/2... keep th in [0,π/2] and
    # note cos is even about 0; E=-cos θ with θ in [0,π/2] covers [−1,0].
    # For full coverage use θ in [0,π] without folding:
    th_full = np.abs(a - b)  # [0,π)
    dkls = np.array([dkl_levelset(float(t)) if t > 1e-12 else 0.0 for t in th_full])
    return float(np.mean(dkls)), float(np.std(dkls))

def chsh_only_dkl():
    # four CHSH angles all have θ=π/4 or 3π/4; both give same |cos|=1/√2
    return dkl_levelset(math.pi / 4)

def discretised_I(n_ang=12, n_grid=2000):
    """Finite net of angles; full I(Λ:S) with S=(α,β)."""
    angles = np.linspace(0, math.pi, n_ang, endpoint=False)
    grid = np.linspace(0, 2 * math.pi, n_grid, endpoint=False)
    rhos = []
    for a in angles:
        for b in angles:
            th = abs(a - b)
            if th > math.pi:
                th = 2 * math.pi - th
            if th > math.pi:  # noqa — th in [0,π]
                pass
            # reduce to [0,π]
            th = min(th, 2 * math.pi - th) if th > math.pi else th
            # actually a,b in [0,π), th in [0,π)
            A = np.sign(np.cos(grid - a))
            B = -np.sign(np.cos(grid - b))
            sig = A * B
            sig[sig == 0] = 1
            Eq = -math.cos(th) if th <= math.pi else -math.cos(2*math.pi - th)
            # use th as constructed
            th_use = abs(a - b)
            Eq = -math.cos(th_use)
            pp = (1 + Eq) / 2
            pm = 1 - pp
            mup = np.mean(sig > 0)
            mum = 1 - mup
            dens = np.where(sig > 0, pp / (mup + 1e-30), pm / (mum + 1e-30))
            p = dens / n_grid
            p = p / p.sum()
            rhos.append(p)
    marg = sum(rhos) / len(rhos)
    def H(p):
        p = p[p > 0]
        return float(-np.sum(p * np.log2(p)))
    Hc = sum(H(r) for r in rhos) / len(rhos)
    return H(marg) - Hc, len(rhos)

if __name__ == "__main__":
    print("=== G11 singlet / continuum MI (achievable with circle tilt) ===")
    print("CHSH-only per-setting D_KL(π/4) =", chsh_only_dkl())
    mean_d, std_d = mean_dkl_random_angles()
    print(f"Mean D_KL(ρ||unif) over random planar angles = {mean_d:.6f} ± (std {std_d:.4f})")
    print("  (This averages single-pair KL to *uniform*, not true I to marginal.)")
    for n in [4, 6, 8, 12]:
        I, nset = discretised_I(n_ang=n)
        print(f"  discretised I with {n}×{n}={nset} setting pairs: I≈{I:.6f} bits")
    print()
    print("Brief Hall singlet figure ~0.0663 bits (full singlet, optimised model).")
    print("Our CHSH-only min I=0.0463 bits; continuum mean D_KL differs.")
    print("G11: achievable MD cost for planar −cos correlators is O(0.05) bits;")
    print("  model-independent lower bound not re-derived here (see Hall PRA 2011).")
