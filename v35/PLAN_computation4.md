# Computation 4: θ Self-Interaction — Can θ Form Solitons?

## Goal

Add a nonlinear self-coupling to the θ equation so θ can form its
own localized structures (θ-braids / θ-solitons). If they exist,
their mass IS m_eff. Combined with ℏ_eff = E_braid/ω ≈ 23,000,
this determines the Bohr radius with zero free parameters.

## The Modification

Current θ equation (linear, no self-interaction):

    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a + η×curl(φ)_a

Modified (with triple-product self-coupling, same form as φ):

    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a + η×curl(φ)_a - ∂V_θ/∂θ_a

    V_θ(Q) = (μ_θ/2) Q² / (1 + κ_θ Q²)
    Q = θ₀ θ₁ θ₂  (angle triple product)

This mirrors V(P) for the position fields. The θ field now has
its own binding potential — if the parameters are right, θ can
form helical solitons just like φ does.

## Key question: what μ_θ and κ_θ?

For φ: μ=-41.345, κ=50. These were found by CMA-ES search (V28).
For θ: we need to scan. The θ field amplitudes are ~100× smaller
than φ (θ_rms ≈ 0.025 vs φ_rms ≈ 0.8), so μ_θ may need to be
much larger to compensate.

Scaling argument: V_θ should be comparable to the θ kinetic energy.
θ kinetic ~ ½ θ_rms² ω² ≈ 0.5 × 0.025² × 0.22² ≈ 1.5e-5.
V_θ ~ (μ_θ/2) Q² where Q ~ θ_rms³ ≈ 1.6e-5.
For V_θ ~ θ kinetic: μ_θ ≈ 2 × 1.5e-5 / (1.6e-5)² ≈ 120,000.

So μ_θ needs to be ~100,000× larger than μ_φ to have comparable
effect at the θ amplitude scale. This is large but not unreasonable.

## Experiments

### Phase A: θ self-interaction parameter scan

Modify v33_cosserat to add V_θ(Q). Scan:
- μ_θ: {0, -1000, -10000, -100000, -500000}
- κ_θ: {50, 500, 5000, 50000}
- Keep η=0.5, m_θ²=0, m²=2.25

For each: single φ-braid, N=80, L=25, T=50.
Measure: does θ develop localized structure? θ_rms? θ triple product Q?

### Phase B: θ-braid initialization

If Phase A shows θ self-organization, try initializing a θ-braid:
- Same helical init as φ but in the θ fields
- Lower amplitude (A_θ = 0.01-0.1)
- WITH a φ-braid present at some distance

Does the θ-braid survive? Does it orbit the φ-braid?

### Phase C: θ-braid mass measurement

If θ-braids exist:
- Measure E_θ-braid (energy bound in the θ soliton)
- m_eff = E_θ-braid / c²
- Compute Bohr ratio: a₀ = ℏ_eff² / (m_eff × α_eff)
- Does it match ~53,000?
