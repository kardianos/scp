# V4 Extension — Self-Confinement Models

## Conceptual Shift

Q-balls confine waves via a potential well (mass gap). The wave sits IN a pre-existing well.

The alternative: the wave's own nonlinear self-interaction creates localization. No external well needed. The interaction IS the confinement. The wave IS the well.

In 3+1D, Derrick's theorem prevents static self-confined solutions for scalar fields. But OSCILLATING configurations can evade Derrick. The question: do any concrete, massless, self-interacting field theories support long-lived oscillating localized solutions?

Key distinction: massive fields give exponentially localized, eternal solutions. Massless fields give power-law tails and quasi-stable solutions (finite but potentially very long lifetime). For the purpose of this investigation, a quasi-breather lasting ≫ 100 oscillation periods counts as POSITIVE.

---

## Model A: SU(2) Yang-Mills Quasi-Breather

**Lagrangian**: L = -1/(4g²) F^a_μν F^{μν a}, with F = dA + g[A,A]

**Key property**: Completely massless. The self-coupling [A,A] provides the ONLY nonlinearity. No potential, no mass gap in the Lagrangian. Scale-free classically.

**Hedgehog ansatz**: A_i^a = ε_{ija} x_a/r² (1 - w(r,t))

After angular integration, the dynamics reduces to a single function w(r,t):
- PDE: ẅ = w'' - w(w²-1)/r²
- Energy: E = ∫ [ẇ² + w'² + (w²-1)²/(2r²)] dr
- Note: NO 2w'/r friction term (r² already absorbed by angular integration)

**Effective potential**: V(w) = (w²-1)²/(4r²) has degenerate minima at w = ±1 (both vacuum) and an unstable maximum at w = 0. A "bubble" of w = -1 embedded in w = +1 vacuum will oscillate.

**BCs**: w(0) = 1 (gauge regularity at origin), w(∞) = 1 (vacuum).
**IC**: w(r,0) = 1 - A·exp(-(r-R₀)²/σ²), creating a Gaussian shell perturbation.

**Expected**: Quasi-breather for A ≈ 2 (w dips to -1). The bubble oscillates between w = ±1 at r ≈ R₀. Energy slowly radiates outward. Lifetime depends on amplitude and width.

**Implementation**: Leapfrog PDE solver, 1+1D.

---

## Model B: Born-Infeld Perturbations on BIon

**Lagrangian**: L_BI = -b²[√(1 - E²/b²) - 1]

**Key property**: Massless EM field with saturating nonlinearity. The BIon (static monopole) has E_r(r) = Q/(4πr²)/√(1+(Q/(4πbr²))²), which saturates at E = b near r = 0.

**Self-trapping mechanism**: On the BIon background, the phase velocity for transverse perturbations is v_⊥ = c·(1-E₀²/b²)^{1/2}. Near the BIon core (E₀ → b), v_⊥ → 0. The wave slows to a halt — a self-consistent potential well.

**Approach**: Linearize BI Maxwell equations on BIon background. Seek quasi-normal modes with angular momentum ℓ ≥ 1. If a narrow resonance exists (Im(ω) ≪ Re(ω)), the BIon supports a long-lived oscillating halo.

**Caveat**: §11.4 showed m²_eff = 0 on UNIFORM backgrounds. But the BIon background is non-uniform. The position-dependent v(r) might create a potential well.

**Implementation**: Sturm-Liouville eigenvalue problem, 1D.

---

## Model C: Pure DBI Massless Oscillon

**Lagrangian**: L = -b²[√(1 + (∂σ)²/b²) - 1], NO potential

**Key property**: Massless scalar, saturating kinetic term. At low amplitude → free wave. At high amplitude → DBI saturation slows propagation.

**PDE**: ∂_t[σ̇/Γ] - (1/r²)∂_r[r²σ'/Γ] = 0, where Γ = √(1 + (σ'² - σ̇²)/b²)

**Expected**: NULL. Without a potential, there's no restoring force. The field can take any value, and the DBI kinetic term only affects dynamics (speed), not statics (equilibrium). A Gaussian pulse will disperse, albeit slower than linear theory due to DBI saturation.

**Implementation**: 1+1D PDE (spherically symmetric).

---

## Model D: Boson Star (Gravity as Self-Confinement)

**Lagrangian**: L = |∂φ|² - μ²|φ|² + R/(16πG) (Einstein-Hilbert + massive complex scalar)

**Key property**: Gravity IS the self-interaction that confines. The scalar field's own energy density curves spacetime, which traps the field. No self-interaction potential needed beyond the mass term.

**This is a KNOWN positive result**. Boson stars are proven to exist and be stable (Kaup 1968, Ruffini & Bonazzola 1969). They directly realize "the interaction constrains it to itself" where the interaction = gravity.

**Equations**: Coupled ODEs (TOV-like):
- dM/dr = 4πr²ρ
- dΦ/dr = (M + 4πr³p)/(r²(1-2M/r))
- φ'' + [2/r + Φ' + M'_correction]φ' + [ω²e^{-2Φ} - μ²]φ/(1-2M/r) = 0

**Shooting**: Fix μ=1, vary φ(0) and ω until φ→0 at r→∞.

**Implementation**: Radial ODE shooting (like Q-ball solver).

---

## Model E: Soler Model (Nonlinear Dirac, Spin-½)

**Lagrangian**: L = ψ̄(iγ^μ∂_μ - m)ψ + λ(ψ̄ψ)²

**Key property**: Self-interacting spinor field. The (ψ̄ψ)² term creates a self-consistent potential well. Naturally spin-½. Has a mass parameter m but confinement is from self-interaction.

**Ansatz**: ψ = [f(r)χ_κ, ig(r)σ·r̂ χ_κ] e^{-iωt}

**EOM**: iγ∂ψ = (m - 2λ(ψ̄ψ))ψ, i.e. m_eff = m - 2λ·S where S = ψ̄ψ = f² - g²

**Coupled ODEs** (with λ_code = 2λ_Lagrangian, κ = -1 ground state):
  f' = (m + ω - λ(f²-g²)) g
  g' + 2g/r = (m - ω - λ(f²-g²)) f

Note: the coupling is MINUS λS (reduces effective mass), not plus.
f = upper/large component (f(0) = f₀, shooting parameter), g = lower/small (g(0) = 0).

**Result**: POSITIVE. Solutions exist for ω ∈ [ω_min, m) with ω_min ≈ 0.13m at λ=1.
Q(ω) has minimum at ω ≈ 0.95m. Soliton most compact at ω ≈ 0.7m.

**Implementation**: RK4 + bisection on f(0), with tail clipping.

---

## Execution Order

1. Model A (YM quasi-breather) — most novel, uncertain
2. Model B (BI on BIon) — concrete, testable
3. Model C (DBI massless) — likely null, quick check
4. Model D (Boson star) — known positive, for completeness
5. Model E (Soler) — known positive, for spin-½

Stop at first positive result (quasi-stable oscillating localized solution).

---

## Summary of Results

| Model | Confinement? | Mechanism | Gravity? | Spin |
|-------|:---:|---|:---:|:---:|
| A: YM breather | **NULL** | Gauge [A,A] too weak — 0% core retention at t=30 | — | 1 |
| B: BI on BIon | **NULL** | Refractive well exists but m²_eff=0 → leaks to ∞ | — | 1 |
| C: DBI massless | **NULL** | No potential → no restoring force → disperses | — | 0 |
| D: Boson star | **YES** | Gravity IS the self-interaction that traps | built-in | 0 |
| E: Soler | **YES** | Scalar self-coupling (ψ̄ψ)² reduces m_eff → binds | no | ½ |

### Key findings

1. **Massless fields cannot self-confine without gravity or a mass gap.** Models A–C all fail. The Yang-Mills gauge nonlinearity (Model A) disperses completely. Born-Infeld creates a graded-index well (Model B) but lacks the mass gap to prevent tunneling to infinity. Pure DBI (Model C) has no restoring force at all.

2. **Self-confinement requires either gravity (Model D) or a mass term + self-interaction (Model E).** In both cases, something prevents the field from dispersing: curved spacetime in Model D, the fermion mass gap in Model E. The DBI Q-ball (Phase 3, not listed here) also confines via mass gap + refractive index.

3. **Gravity remains the unsolved problem.** Models D and E both confine, but:
   - Model D assumes GR (Poisson equation for Φ) — it doesn't explain gravity, it uses it.
   - Model E has no long-range force at all — all interactions are Yukawa-screened by the mass m.
   - The confinement-range dilemma (v2/v3/Phase 5) persists: the mass that traps the wave also screens the force.

4. **Model E (Soler) is the most interesting for particle physics.** It gives spin-½ without assuming it — the Dirac equation's spinor structure is intrinsic. The soliton family is parameterized by ω/m, with a minimum at ω_min ≈ 0.13m. The particle number Q has a characteristic U-shape with minimum Q ≈ 47 near ω = 0.95m.

### What works, what to forward

**For the "particles = self-trapped waves" program:**
- DBI Q-ball (Phase 3) + Soler model (Model E) together capture the essential ingredients: nonlinear self-trapping, mass gap, spin-½ structure.
- A Dirac field with DBI-type kinetic nonlinearity would combine both: spin-½ + refractive self-trapping. This is the natural next model.

**For gravity:**
- Model D (boson star) proves gravity CAN provide self-confinement. But it assumes GR rather than deriving it.
- No model in the v4 program generates 1/r gravity from the soliton itself. The structural obstruction (confinement-range dilemma) blocks all attempts.
- The only remaining avenue: a massless, decoupled mediator (dilaton) with fine-tuned coupling α ≈ 10⁻¹⁸ — which IS the hierarchy problem.
