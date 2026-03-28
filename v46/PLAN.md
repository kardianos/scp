# V46 — Analytical Foundations: Derrick Scaling, Emergent Gauss's Law, Binding Threshold

## Goal

Before further GPU simulation, establish the analytical foundations:
1. Derrick's virial identity for the full Cosserat Lagrangian
2. The λ² scaling of the η curl coupling (binding enabler)
3. Monopole radiation from composite baryons (emergent Gauss's law)
4. Adler synchronization threshold for breathing-mode locking
5. η_crit for nuclear binding (if it exists)

Tools: Maxima (computer algebra), existing V44/V45 numerical data.

## Analysis Program

### A. Derrick Virial for Cosserat (Maxima)

Compute the energy functional scaling under x → λx for:
    L = ½(∂φ)² - ½(∇φ)² - ½m²φ² - V(φ₀φ₁φ₂)
      + ½(∂θ)² - ½(∇θ)² + η φ·(∇×θ)

Confirm each term's λ-exponent. Derive the virial identity:
    a·E_grad + b·E_coupling + c·(E_mass + E_pot) = 0

If b=2 (as predicted in -06), the curl coupling has Skyrme-like scaling
and CAN enable binding. Derive η_crit from the virial balance.

### B. Phase-Averaged Effective Potential (Maxima)

Compute ⟨P²⟩ averaged over one carrier cycle for:
    P = A₀cos(kz+δ₀) × A₁cos(kz+δ₁) × A₂cos(kz+δ₂)
    δ = {0, 3.0005, 4.4325}

Then compute the cross-term structure for two baryons:
    P_total = P₁(x) + P₂(x-D)
    ⟨P_total²⟩ - ⟨P₁²⟩ - ⟨P₂²⟩ = 2⟨P₁ P₂⟩

Determine whether this cross-term is constructive or destructive
for the gen_deuterium carrier phase choices.

### C. Monopole Radiation Derivation (Maxima + numerical)

Compute the monopole moment of three orthogonal magnetic dipoles:
- Dipole 1 along x: B₁ ∝ (2cos²θ - sin²θ)/r³ pattern
- Dipole 2 along y: B₂ rotated 90°
- Dipole 3 along z: B₃ rotated 90°

Show |B₁ + B₂ + B₃|² has nonzero l=0 (monopole) component.
Compute Q_eff and compare with V44 OQ1 measurement (C_0 = 0.000489 at r=20).

### D. Adler Synchronization (Maxima + V34 data)

Using θ_rms(r) profile from V34 braid_hires:
    θ_rms peaks at r≈1.25, decays to 0.037 at r=25

Compute the overlap integral for two baryons at D=40:
    ε_eff(D) ~ η² ∫ θ₁(r) θ₂(|r-D|) d³r

Estimate T_lock = 2π / ε_eff. If T_lock >> 500 (V45 run time),
the null result is explained by insufficient time, not missing physics.

### E. Targeted Simulations (after analytical results)

Based on analytical findings:
- If η_crit > 0.5: run η sweep to find binding
- If T_lock > 500: run D=40 at T=5000
- If cross-terms destructive: redesign seeds with correct carrier phases
- If monopole Q_eff matches OQ1: emergent Gauss's law confirmed

## Deliverables

- v46/derrick_virial.mac — Maxima worksheet for Derrick scaling
- v46/phase_average.mac — Maxima worksheet for ⟨P²⟩ cross-terms
- v46/monopole_radiation.mac — Maxima worksheet for composite monopole
- v46/adler_coupling.mac — Maxima worksheet for synchronization
- v46/RESULTS.md — analytical results and implications
