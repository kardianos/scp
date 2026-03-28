# Derived Parameters for SCP Field Theory — Revision 06

## Direction: Analytical First, Simulate Second

The V45 null result and the -05 energy decomposition revealed a precise
diagnosis: V(P) is attractive (ΔE_pot = -146 at D=0) but field distortion
cost (ΔE_kin + ΔE_mass = +1155) overwhelms it 8:1. Before running more
simulations, we should determine analytically WHETHER binding is possible
under the current equation, and if so, under what conditions.

This revision proposes using Maxima (computer algebra) to answer three
specific mathematical questions that directly determine the experimental
path.

---

## 1. Derrick's Theorem for the Full Cosserat Lagrangian

### The question
Does the η curl coupling between φ and θ have the right scaling properties
to overcome the gradient energy penalty that prevents binding in 3D?

### The calculation
The full Lagrangian density is:

    L = ½(∂_t φ)² - ½(∇φ)² - ½m²φ² - V(P)
      + ½(∂_t θ)² - ½(∇θ)² - ½m_θ²θ²
      + η φ·(∇×θ)

where P = φ₀φ₁φ₂ and V(P) = (μ/2)P²/(1+κP²).

Under spatial rescaling x → λx (Derrick scaling), each term transforms as:

| Term | Scaling | Type |
|------|---------|------|
| ½(∇φ)² | λ^(d-2) = λ | Gradient (E₂) |
| ½m²φ² | λ^d = λ³ | Mass |
| V(P) | λ³ | Potential |
| ½(∇θ)² | λ | Gradient |
| ½m_θ²θ² | λ³ | Mass |
| **η φ·(∇×θ)** | **λ^?** | **Cross-coupling** |

The critical question is how η φ·(∇×θ) scales. The curl ∇×θ has one
spatial derivative, so:

    φ·(∇×θ) ~ φ × (1/λ) × θ × λ³ = λ² × (φ·θ content)

Wait — more carefully: under x → λx with φ(x) → φ(λx):
- φ(λx) has amplitude unchanged but compressed spatially
- ∇×θ → (1/λ)∇×θ (derivative brings down 1/λ)
- Integration d³x → λ³ d³x

So: ∫ φ·(∇×θ) d³x → λ³ × (1/λ) × ∫ φ·(∇×θ) d³x = **λ² × E_coupling**

This is DIFFERENT from both E₂ (scales as λ) and E_mass (scales as λ³).
It scales as λ², which is intermediate.

### Derrick's equilibrium condition

At equilibrium (dE/dλ = 0 at λ=1):

    1·E_grad + 2·E_coupling + 3·(E_mass + E_pot) = 0

Or equivalently:

    E_grad + 2·E_coupling + 3·E_mass + 3·E_pot = 0

For a bound state of TWO solitons at separation D, the interaction energy
ΔE(D) = E(D) - E(∞) must have a minimum at some finite D. The η coupling
term scales as λ², which means it can provide a NEGATIVE contribution to
the virial that partially cancels the positive E_grad term.

**This is the key result**: the curl coupling has λ² scaling, which sits
between gradient (λ) and potential (λ³). It is the ONLY term in the
Lagrangian with this intermediate scaling. In the Skyrme model, the
analogous role is played by the E₄ Skyrme term (also λ² in 3D), which
is what enables Skyrmion binding.

### Maxima calculation needed
Verify this scaling argument symbolically. Then compute: for what value of
η does the virial condition admit a bound two-soliton solution? This gives
η_crit analytically.

---

## 2. Time-Averaged Effective Potential

### The question
What is the effective inter-baryon potential ⟨V_eff(D)⟩ after averaging
over the fast carrier oscillation?

### The calculation
Each braid has carrier wave φ_a = A_a(r) × cos(kz + δ_a + Δ_carrier).
The triple product is:

    P = φ₀φ₁φ₂ = A₀A₁A₂ × cos(kz+δ₀)cos(kz+δ₁)cos(kz+δ₂)

Using the product-to-sum identity for three cosines:

    cos(α)cos(β)cos(γ) = ¼[cos(α+β+γ) + cos(α+β-γ) + cos(α-β+γ) + cos(-α+β+γ)]

Time-averaging over one carrier cycle kills all oscillating terms, leaving:

    ⟨P⟩ = A₀A₁A₂ × C(δ₀,δ₁,δ₂)

where C depends on the specific phase offsets. For δ = {0, 3.0005, 4.4325}:

    C = ¼[cos(δ₀+δ₁+δ₂) + cos(δ₀+δ₁-δ₂) + cos(δ₀-δ₁+δ₂) + cos(-δ₀+δ₁+δ₂)]

This determines whether the effective potential is attractive or repulsive
on average. If C > 0, ⟨V(P)⟩ provides a net attractive well.

For two baryons at separation D, the effective interaction is:

    ⟨V_int(D)⟩ = μ ∫ [⟨P_total²⟩ - ⟨P₁²⟩ - ⟨P₂²⟩] / [1 + κ⟨P_total²⟩] dV

The cross-terms in P_total = P₁ + P₂ depend on the RELATIVE carrier
phases between the two baryons. If they are phase-locked (same Δ_carrier),
the cross-terms are constructive → deeper potential well. If random phases,
the cross-terms average to zero → no interaction.

### Maxima calculation needed
Compute C(δ₀,δ₁,δ₂) for the empirical δ values. Compute the cross-term
structure for two baryons as a function of their relative carrier phase.
This determines whether the gen_deuterium seeds (which set carrier phases
deterministically) are in a constructive or destructive configuration.

---

## 3. Adler Synchronization Threshold

### The question
How strong must the inter-baryon breathing coupling be for phase-locking
to occur?

### The calculation
Model each baryon as a nonlinear oscillator with natural frequency
ω₀ = 2π/150t (breathing period). Two oscillators coupled with strength ε:

    d²x₁/dt² + ω₀²x₁ = ε × f(x₂)
    d²x₂/dt² + ω₀²x₂ = ε × f(x₁)

The Adler equation for the phase difference Δψ = ψ₁ - ψ₂ is:

    dΔψ/dt = Δω + ε_eff × sin(Δψ)

where Δω is the natural frequency difference (zero for identical baryons)
and ε_eff is the effective coupling, which depends on the field overlap
at separation D:

    ε_eff(D) ~ η² × ∫ |θ₁(x)| × |θ₂(x-D)| d³x

For identical oscillators (Δω=0), locking is ALWAYS achieved eventually
— the question is how long it takes. The locking time is:

    T_lock ~ 1 / ε_eff(D)

If ε_eff is small (weak coupling, large D), locking takes many cycles.
The V45 runs at T=500 ≈ 3 breathing cycles may simply be too short.

### Maxima calculation needed
Compute ε_eff(D) from the measured θ radial profile (V34 braid_hires data:
θ_rms(r) peaks at r≈1.25, decays to 0.037 at r=25). Estimate T_lock for
D=40. If T_lock > 500, that explains V45's null result without requiring
new physics.

---

## 4. Carried Forward from -05

### η=0 baseline (Priority 1)
Still the most important single experiment. Removes θ entirely. If φ-only
interaction is STILL repulsive (expected from Derrick), confirms that
binding MUST come from the θ sector.

### η sweep (Priority 2)
The Derrick analysis above predicts that η_crit exists (the curl coupling
has the right scaling). The sweep finds it numerically. If Maxima gives
an analytical η_crit, we can validate it with one or two targeted runs
instead of a full sweep.

### T=5000 long run at D=40 (Priority 3)
Tests breathing synchronization. The Adler analysis above predicts
T_lock — if T_lock ≈ 1000-2000, a T=5000 run should show locking.
If T_lock > 10,000, we need stronger coupling (higher η).

### Massive θ at D=40 (Priority 4)
Still a valid test. A massive mediator creates Yukawa binding regardless
of breathing synchronization. This is the "brute force" path — add the
physics if the existing equation can't produce it.

### Controlled quench (Priority 6)
Worth trying after the analytical work clarifies what parameters to use.
Running a quench with parameters that Maxima shows can't bind is wasted
GPU time.

---

## 5. The Analytical Program

Before any v46 simulation, run these Maxima calculations:

### Calculation A: Derrick virial for Cosserat
Input: Full Lagrangian with V(P) = (μ/2)P²/(1+κP²), curl coupling η.
Output: Virial identity. Confirm λ² scaling of coupling term.
Derive: η_crit as function of soliton profile integrals.

### Calculation B: Phase-averaged effective potential
Input: Triple product with carrier phases δ = {0, 3.0005, 4.4325}.
Output: ⟨P⟩ and ⟨P²⟩ averaged over one carrier cycle.
Derive: Whether cross-terms between two baryons are constructive or
destructive for gen_deuterium carrier phase choices.

### Calculation C: Adler coupling estimate
Input: θ_rms(r) profile from V34, breathing period 150t.
Output: ε_eff(D=40), T_lock estimate.
Derive: Whether T=500 is sufficient for synchronization.

### Calculation D: Massive θ Yukawa binding
Input: Modified Lagrangian with m_θ > 0.
Output: Yukawa potential depth as function of m_θ and η.
Derive: m_θ_crit where binding appears.

If calculations A-D show binding is possible within the current equation
(just at different η or longer T), the experimental path is clear.
If they show it's mathematically impossible, we need new physics
(massive θ, higher-order terms, or the determinant potential from -03 §F).

---

## 6. Summary

The shift from -05: **do the math before burning GPU hours.**

The λ² scaling of the curl coupling is the most important theoretical
finding in this series. It means the η coupling has the same role as the
Skyrme E₄ term — the ONLY term that enables soliton binding in 3D. The
question is whether η=0.5 is above or below the binding threshold. Maxima
can answer this in hours. If the answer is "η is too low," we know exactly
what to change. If the answer is "η is fine but T is too short," we know
to run longer. Either way, we stop guessing.
