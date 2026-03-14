# V22 Plan: Massless Scalar Gravity Coupled to Three-Body Oscillon

## Physics

Massless scalar field Φ sourced by the oscillon's FULL energy density ρ:

    □Φ = -α ρ         (Φ equation)
    ρ = Σ_a [½v_a² + ½|∇φ_a|² + ½m_eff²φ_a²] + V(P)   (FULL energy density)
    m_eff² = m₀² - β Φ    (mass modulation by gravity)

Backreaction on each φ_a:

    ∂²φ_a/∂t² = ∇²φ_a - m_eff²·φ_a - ∂V/∂φ_a

Sign conventions:
- α > 0, ρ > 0 → static Φ = -(α/4π)∫ρ/|r-r'|d³r' < 0 near source
- β > 0, Φ < 0 → -βΦ > 0 → m_eff > m₀ near oscillon (deeper gap = stabilizing)
- Far away: Φ → -(αE_osc)/(4πr) → 1/r Newtonian tail

**Key insight (monopole suppression)**: Using the FULL energy density (including
kinetic ½v²) as source means ∫ρd³x = E_total = const (by energy conservation).
The monopole moment is constant → no monopole radiation. Previously, excluding
kinetic energy caused Q to oscillate 1→130 (100× variation), launching spurious
monopole waves. With full ρ, Q varies only ~10% (tracking boundary losses).

## Implementation Phases

### Phase 1: Static Φ from equilibrium oscillon — DONE ✓

Poisson solver (`src/poisson_phi.c`) confirms perfect 1/r tail.
Q = 85.0, Φ(0) = -2.82 at α=1. r·Φ converges to const at r>15.

### Phase 2: Dynamic Φ + backreaction — DONE ✓

`src/triad3d_grav.c` implements full coupled system (Verlet, OpenMP, Poisson init).

**Stability scan** (N=100, t=300, full ρ source):

| α, β | E(300) | fc(300) | ω | Q variation | Alive? |
|------|--------|---------|---|-------------|--------|
| 0.01, 0.01 | 86.6 | 0.989 | 0.948 | 1.1× | YES |
| 0.1, 0.1 | 89.1 | 0.988 | 0.954 | 1.1× | YES |
| 0.2, 0.2 | 90.4 | 0.986 | 0.966 | — | YES |
| 0.3, 0.3 | 96.8 | 0.970 | 0.984 | — | YES (marginal) |
| 0.5, 0.5 | — | — | — | — | DEAD (ω > m) |

Stability boundary at α≈0.4 (ω crosses mass gap). Physical α ~ 10⁻³⁸ is safe.

### Phase 3: Self-consistent bootstrap (m₀ → 0) — DONE ✗ (NEGATIVE)

Scanned m₀ from 1.0 down to 0.4 with α=0.1, β=0.1 (and β=0.5 for m₀=0.6).

**Result**: Bootstrap fails. The vacuum mass gap is m₀, not m_eff(core).
Reducing m₀ closes the gap faster than backreaction can compensate.
At m₀=0.8: ω=0.821 > m=0.800 → above gap → radiates → dies.
At m₀<0.62: false vacuum → runaway instability.
Stronger β destabilizes the oscillon before providing enough mass.

**Conclusion**: m₀ must be a fundamental parameter. The mass gap cannot
emerge purely from gravitational self-interaction in this model.

### Phase 4: Adjoint optimization

Extend adjoint3d.c to include Φ dynamics:
- Extra tangent arrays for Φ (3 arrays × NPAR = 12 arrays at N=100 = 96 MB)
- Loss = -fc + w₁·(deviation from 1/r) + w₂·(monopole amplitude)
- Optimize: α, β, m₀, μ, κ, A, σ (7 params, still feasible with forward-mode)

## Critical Issues

### 1. Source definition — RESOLVED ✓
ρ_source = full energy density of φ_a fields (kinetic + gradient + mass + potential).
Does NOT include Φ's own energy ½(∂Φ)² (avoids self-sourcing runaway).
At α ≤ 0.1, self-sourcing from kinetic energy is negligible.

### 2. Monopole radiation — RESOLVED ✓
Using full ρ (with kinetic energy), ∫ρd³x = E_total = const by conservation.
Monopole oscillation suppressed by 100× (from 130× variation to 1.1×).
Residual ~10% variation tracks slow energy loss through absorbing boundary.

### 3. Sign of mass correction
Φ < 0 near source. For m_eff² = m₀² - βΦ with β > 0:
- Core: m_eff > m₀ (deeper gap, stabilizing) ✓
- Far field: m_eff → m₀ (unchanged) ✓
For m_eff² = γΦ with m₀=0 and γ > 0:
- Core: m_eff² = γΦ < 0 (TACHYONIC!) ✗
- Must use γ < 0: m_eff² = |γ|·|Φ| > 0 in core. But far field: m_eff → 0. ✓
- This means exterior is gapless — oscillon can radiate.

### 4. G coefficient
Static: Φ ~ -(αE)/(4πr). Gravitational acceleration: g = -∇Φ = -αE/(4πr²).
Compare with Newton: g = -GM/r² = -GE/(c²r²).
Match: α/(4π) = G/c² → α = 4πG/c² ≈ 9.3e-27 m³/(kg·s²) × ...
In code units: α_code ≈ (G_N × m_code²) / (ℏc) ≈ 10⁻³⁸.
This is tiny — backreaction is perturbative. Good for stability, bad for observability.

### 5. Initialization transient
If Φ starts at 0 instead of Poisson solution, the system rings violently.
Always initialize Φ = Φ_static(ρ_initial).

## Next Steps

1. **Phase 3**: Self-consistent bootstrap (m₀ → 0). Highest risk, highest reward.
   Need to check whether m_eff(vacuum) = 0 is fatal (gapless exterior radiates).

2. **Phase 4**: Adjoint optimization with Φ. Extend adjoint3d.c with Φ tangent vectors.
   Optimize α, β jointly with μ, κ, A, σ.

3. **Two-oscillon interaction**: Place two oscillons at separation D, measure force.
   Verify attractive 1/r² force law from Φ mediation. Requires larger box (L~50+).
