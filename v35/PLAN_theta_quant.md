# V35: θ and ℏ Computation Plan

Three parallel paths to determine the two free parameters (ℏ_sim, m_eff)
from the field dynamics. If BOTH can be derived (not fit), the Bohr
radius becomes a zero-free-parameter prediction.

---

## Computation 1: θ Frequency Spectrum (cheapest — analysis only)

**Goal**: Find natural resonant frequencies of the θ field around the braid.
If the spectrum has discrete peaks, those are the orbital energy levels.

**Method**:
- Read 264 frames from V34 sfa_hires.sfa (N=80, L=25, T=50, dt≈0.19)
- At each radial shell (r=5, 10, 15, 20), extract θ_φ(t) time series
- FFT → power spectrum P(ω)
- Also: θ_r(t), θ_z(t) spectra for comparison

**What we learn**:
- Discrete peaks → natural modes → energy levels WITHOUT needing ℏ
- Peak frequencies → ω_n → energy spacing ΔE = ω₁ - ω₀
- The ratio ΔE/E₀ compared to hydrogen (should be 3/4 for n=1→2)
- If continuous spectrum: no natural quantization, ℏ must be imposed

**Compute time**: Minutes (pure Python FFT on existing data)

## Computation 2: θ Dispersion Relation (medium — small simulations)

**Goal**: Determine the effective mass of θ perturbations near the braid.

**Method**:
- Initialize small-amplitude θ plane waves at various k values
  in the Cosserat simulation (with a settled braid present)
- Measure the propagation frequency ω for each k
- Plot ω(k) → the dispersion relation
- In free space: ω² = k²c² + m_θ²c⁴ = k² (since m_θ=0, c=1)
- Near the braid: ω² = k² + V_eff(r) → effective mass from the well

**What we learn**:
- m_eff from the band bottom: m_eff = ℏ ω_min / c²
- Whether the braid creates a mass gap for θ perturbations
- The group velocity of θ waves at different k

**Two approaches**:
A. Analytical: linearize the Cosserat equation around the settled braid,
   compute the eigenfrequencies of the linearized θ operator.
   This is a 3D eigenvalue problem but can be reduced to 1D radial.

B. Numerical: inject θ waves, measure their behavior.
   Simpler but requires multiple simulation runs.

**Compute time**: 1-2 hours (several short sims or one eigenvalue computation)

## Computation 3: Field Quantization ε (medium — code mod + simulations)

**Goal**: Add the minimum field quantum and determine ℏ_sim.

**Method**:
- Modify the Verlet step: after updating φ and θ, round each value
  to the nearest multiple of ε
- ℏ_sim = ε × dx (minimum action = field quantum × grid spacing)
- Run at ε = 0.0001, 0.001, 0.01, 0.1
- Measure: braid survival, θ_rms, energy conservation

**What we learn**:
- At what ε does the braid break? (stability threshold)
- Does quantization create discrete energy levels in θ?
- The natural ℏ_sim = ε_max × dx where ε_max is the largest stable ε
- Combined with m_eff from Computation 2: Bohr radius = ℏ²/(m_eff × α)

**Compute time**: 1-2 hours (4 simulations at different ε)

## Computation 4: θ Self-Interaction (longer — code mod + search)

**Goal**: Can θ form its own solitons? If yes, their mass IS m_eff.

**Method**:
- Add V_θ = (μ_θ/2)(θ₀θ₁θ₂)²/(1+κ_θ(θ₀θ₁θ₂)²) to the θ equation
- Scan μ_θ and κ_θ for parameter regimes where θ-braids form
- If θ-braids exist: measure their mass, size, stability
- Place a θ-braid near a φ-braid: does it orbit?

**What we learn**:
- Whether the electron is a θ-soliton
- m_eff from the θ-braid mass directly
- Whether θ-braids can orbit in the φ-braid's potential well

**Compute time**: Days (parameter search required)

---

## Execution Order

1. **Computation 1** (frequency spectrum) — start NOW, free
2. **Computation 3** (field quantization) — start in parallel, small code mod
3. **Computation 2** (dispersion) — after 1 reveals the frequency structure
4. **Computation 4** (θ self-interaction) — after 2 determines if m_eff needs solitons

If Computation 1 shows discrete peaks: ℏ may not be needed explicitly —
the quantization is emergent from the braid's geometry (like a drum's
resonant modes). This would be the best outcome.
