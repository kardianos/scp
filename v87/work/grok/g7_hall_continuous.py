#!/usr/bin/env python3
"""G7: Continuous Door-A model (Hall-style sphere) for singlet correlators.

Model (local deterministic + measurement dependence):
  λ ∈ S² (unit sphere).
  A(a, λ) = sgn(a · λ)
  B(b, λ) = −sgn(b · λ)

Under uniform ρ (full MI): E(a,b) = −1 + (2/π) θ, θ = ∠(a,b) ∈ [0,π].
  CHSH max = 2.

To reproduce the quantum singlet E_qm = −a·b, we allow ρ(λ|a,b).

Exact construction achieving E = −a·b (closed form):
  Use the 'rejection / distance' reweighting of Degorre–Laplante–Roland
  type, specialised to a conditional density on the sphere.

  One explicit density that works (verified numerically below):
    Consider the 1-D reduction along the plane of (a,b). Let φ be the polar
    angle of λ projected; the classical hemisphere model is 1-D.
    We seek ρ(φ|θ) on [0, 2π) such that
      ∫ ρ(φ|θ) A(0,φ) B(θ,φ) dφ = −cos θ
    with A = sgn cos φ, B = −sgn cos(φ−θ).

  Closed-form solution via Fourier / rearrangement:
    Let the product function σ_θ(φ) = A B ∈ {±1}.
    We need ⟨σ_θ⟩_ρ = −cos θ.
    Minimum relative entropy to the uniform prior subject to this one
    correlator constraint is an exponential tilt:
      ρ(φ|θ) ∝ exp( κ(θ) σ_θ(φ) )
    with κ chosen so the correlator matches. Since σ=±1,
      ρ = e^{κσ}/(2 cosh κ) on the two level sets, and
      ⟨σ⟩ = tanh κ  ⇒  κ = artanh(−cos θ)   wait target is −cos θ.
      Actually ⟨σ⟩ = tanh(κ) if ρ ∝ e^{κ σ}.
      Set tanh(κ) = −cos θ,  κ = artanh(−cos θ).

  Mutual information for a single pair (θ fixed settings):
    I = D_KL(ρ || uniform) = log 2 + ⟨log ρ⟩ ... for two-level:
    With p = (1 + ⟨σ⟩)/2 = (1 − cos θ)/2 on σ=+1? 
    ⟨σ⟩ = p_{+}(+1) + p_{-}(−1) = 2p_{+} − 1 = −cos θ
    ⇒ p_{+} = (1 − cos θ)/2
    D_KL to uniform-on-circle: the density is constant on {σ=+1} and on {σ=−1},
    with total masses p_{+}, p_{-}. Uniform has mass = measure fraction.
    Measure of {σ=+1}: depends on θ. For the hemisphere product,
    fraction where σ = +1 equals (1 + E_cl)/2 with E_cl = −1+2θ/π,
    so μ_{+} = θ/π, μ_{-} = 1 − θ/π.

  Careful: exponential tilt on values of σ reweights the two level sets
  keeping *conditional* uniformity inside each, so
    ρ = (p_{+}/μ_{+}) on {σ=+1}, (p_{-}/μ_{-}) on {σ=−1}.
  Then D_KL(ρ||u) = p_{+} log(p_{+}/μ_{+}) + p_{-} log(p_{-}/μ_{-}).

  For the CHSH experiment we have four setting pairs with four angles.
  Settings chosen uniformly → average the KL's (if the prior is the mixture
  marginal, the true I is not exactly the average KL to uniform, but average
  KL to the mixture marginal). We compute both.

This is an *explicit* Door-A model with closed-form κ and E, and a computed I.
It is [P] postulated structure (the tilt κ(θ)), not [D] from fabric dynamics.
"""
from __future__ import annotations
import math
import numpy as np

def E_classical(theta):
    return -1.0 + 2.0 * theta / math.pi

def mu_plus(theta):
    """Fraction of λ with σ = A*B = +1 under hemisphere model.
    E_cl = 2 μ_{+} − 1 ⇒ μ_{+} = (1 + E_cl)/2 = θ/π.
    """
    return theta / math.pi

def rho_masses_for_quantum(theta):
    """Masses p± on {σ=±1} so that ⟨σ⟩ = E_qm = −cos θ."""
    Eq = -math.cos(theta)
    # ⟨σ⟩ = p+ − p− = 2p+ − 1 = Eq ⇒ p+ = (1+Eq)/2
    pp = (1.0 + Eq) / 2.0
    pm = 1.0 - pp
    return pp, pm, Eq

def dkl_to_uniform_levelsets(theta):
    """D_KL(ρ_θ || uniform) in bits for the level-set tilt model."""
    pp, pm, Eq = rho_masses_for_quantum(theta)
    mup = mu_plus(theta)
    mum = 1.0 - mup
    # guard endpoints
    if theta < 1e-12 or abs(theta - math.pi) < 1e-12:
        return 0.0, Eq
    dkl = 0.0
    if pp > 0 and mup > 0:
        dkl += pp * math.log2(pp / mup)
    if pm > 0 and mum > 0:
        dkl += pm * math.log2(pm / mum)
    return dkl, Eq

def chsh_angles():
    # a=0, a'=π/2, b=π/4, b'=-π/4
    # separations:
    thetas = {
        (0, 0): abs(0 - math.pi/4),
        (0, 1): abs(0 - (-math.pi/4)),
        (1, 0): abs(math.pi/2 - math.pi/4),
        (1, 1): abs(math.pi/2 - (-math.pi/4)),
    }
    # reduce to [0,π]
    for k, t in list(thetas.items()):
        t = t % (2 * math.pi)
        if t > math.pi:
            t = 2 * math.pi - t
        thetas[k] = t
    return thetas

def mutual_info_chsh():
    """Approximate I by average KL to uniform (upper bound if prior=uniform;
    exact I uses KL to marginal — we also compute that by discretisation)."""
    thetas = chsh_angles()
    kls = {}
    Es = {}
    for s, th in thetas.items():
        dkl, Eq = dkl_to_uniform_levelsets(th)
        kls[s] = dkl
        Es[s] = Eq
    # Quantum CHSH with E = −cos θ: need sign pattern for S
    # E00 = −cos(π/4) = −1/√2
    # E01 = −cos(π/4) = −1/√2
    # E10 = −cos(π/4) = −1/√2
    # E11 = −cos(3π/4) = −(−1/√2) = +1/√2
    # S = E00+E01+E10 − E11 = −3/√2 − 1/√2 = −4/√2 = −2√2, |S|=2√2
    S = Es[(0, 0)] + Es[(0, 1)] + Es[(1, 0)] - Es[(1, 1)]
    I_avg_uniform = sum(kls.values()) / 4.0
    return dict(thetas=thetas, Es=Es, kls=kls, S=S, I_avg_KL_uniform=I_avg_uniform)

def monte_carlo_tilt(N=300_000, seed=1):
    """MC: for each setting pair, sample φ from the tilted level-set density,
    evaluate A,B, estimate E and S."""
    rng = np.random.default_rng(seed)
    thetas = chsh_angles()
    # Precompute for each θ: sample from {σ=+1} or {σ=-1} with probs p±,
    # then uniform in that level set.
    # Level sets for A=sgn cos φ, B=-sgn cos(φ-θ), σ=A*B:
    # Use rejection from fine grid for simplicity and exactness of measure.
    grid = np.linspace(0, 2*math.pi, 10000, endpoint=False)

    def sample_phi(theta, n):
        A = np.sign(np.cos(grid))
        B = -np.sign(np.cos(grid - theta))
        sig = A * B
        # fix sign(0)
        sig[sig == 0] = 1
        plus_idx = np.where(sig > 0)[0]
        minus_idx = np.where(sig < 0)[0]
        pp, pm, _ = rho_masses_for_quantum(theta)
        n_plus = rng.binomial(n, pp)
        n_minus = n - n_plus
        i_plus = rng.choice(plus_idx, size=n_plus, replace=True)
        i_minus = rng.choice(minus_idx, size=n_minus, replace=True)
        phis = grid[np.concatenate([i_plus, i_minus])]
        rng.shuffle(phis)
        return phis

    results = {}
    for s, th in thetas.items():
        phis = sample_phi(th, N // 4)
        A = np.sign(np.cos(phis))
        B = -np.sign(np.cos(phis - th))
        AB = A * B
        AB[AB == 0] = 1
        results[s] = float(np.mean(AB))
    S = results[(0, 0)] + results[(0, 1)] + results[(1, 0)] - results[(1, 1)]
    return results, S

def discrete_I_exact(n_grid=4000):
    """Exact I(Λ:S) for discretised circle under the four tilted densities."""
    thetas = chsh_angles()
    grid = np.linspace(0, 2 * math.pi, n_grid, endpoint=False)
    dphi = 2 * math.pi / n_grid
    rhos = []
    for s, th in thetas.items():
        A = np.sign(np.cos(grid))
        B = -np.sign(np.cos(grid - th))
        sig = A * B
        sig[sig == 0] = 1
        pp, pm, _ = rho_masses_for_quantum(th)
        mup = np.mean(sig > 0)
        mum = 1 - mup
        dens = np.where(sig > 0, pp / (mup + 1e-30), pm / (mum + 1e-30))
        # dens is probability mass per grid point times n_grid? 
        # probability of bin i: dens[i] * (1/n_grid) if dens is pdf * 2π? 
        # Level-set construction: P(bin i) = (p_level / μ_level) * (1/n_grid)
        # since uniform inside level set. dens as defined is p/μ, so
        p = dens / n_grid
        p = p / p.sum()
        rhos.append(p)
    marg = sum(rhos) / 4.0
    def H(p):
        p = p[p > 0]
        return float(-np.sum(p * np.log2(p)))
    Hc = sum(H(r) for r in rhos) / 4.0
    I = H(marg) - Hc
    return I, H(marg), Hc

if __name__ == "__main__":
    print("=== G7 continuous Hall-style tilt model ===")
    info = mutual_info_chsh()
    print("angles (rad):", {k: round(v, 4) for k, v in info['thetas'].items()})
    print("E_qm:", {k: round(v, 6) for k, v in info['Es'].items()})
    print("S =", info['S'], "  |S| =", abs(info['S']), "  2√2 =", 2*math.sqrt(2))
    print("per-setting D_KL(ρ||uniform) bits:", {k: round(v, 6) for k, v in info['kls'].items()})
    print("avg KL to uniform:", round(info['I_avg_KL_uniform'], 6), "bits")
    I, H, Hc = discrete_I_exact()
    print(f"I(Λ:S) discretised = {I:.6f} bits  (H={H:.4f}, H|S={Hc:.4f})")
    print()
    Emc, Smc = monte_carlo_tilt()
    print("MC E:", {k: round(v, 5) for k, v in Emc.items()})
    print("MC S:", round(Smc, 5), "  analytic", round(info['S'], 5))
    assert abs(abs(info['S']) - 2*math.sqrt(2)) < 1e-9
    assert abs(Smc - info['S']) < 0.03
    print("G7 OK: continuous MD model hits |S|=2√2; I computed [D/M]")
    # Compare to classical at same angles
    Scl = 0
    for s, th in info['thetas'].items():
        e = E_classical(th)
        if s == (1, 1):
            Scl -= e
        else:
            Scl += e
    print("Classical (MI) S at same angles:", Scl)
