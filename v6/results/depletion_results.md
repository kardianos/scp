# V6 Depletion Mechanism — Numerical Results

## Summary

The density conservation + depletion mechanism was tested numerically with
THREE models for the density equation.

**What works:**
- Self-consistent density profiles converge (B=1 hedgehog with ρ-dependent E₂)
- Universality formula ΔQ = E₂/(2αρ₀) verified to 0.3% accuracy
- Core depletion confirmed: ρ(0)/ρ₀ = 0.985 at α=100 (soliton depletes density)
- Poisson integral gives 1/r (mathematical identity, from localized deficit)
- **MASSLESS WAVE MODEL gives δρ ~ 1/r (LONG-RANGE!)** — confirmed numerically

**What fails with the algebraic/fluid models:**
- Algebraic density δρ(r) = -|ω|²/(4α) is LOCAL → decays as 1/r⁶
- Fluid wave (continuity + Euler) converges to algebraic → also 1/r⁶
- The Gauss obstruction: source ∇²|ω|² has zero monopole → no 1/r

**What the massless wave model fixes:**
- Equation: □δρ = -½|ω|² (massless, causal at speed c)
- Source ½|ω|² has NONZERO monopole Q_eff = E₂/ρ₀ = 51.3
- Solution: δρ(r,t) = -(Q_eff/4πr)·θ(ct - r) → **1/r behind wavefront**
- Power-law fit: exponent = -0.989 (essentially -1.0)
- Bypasses Gauss obstruction because source is |ω|², NOT ∇²|ω|²

**Remaining problems:**
- G_eff/G_Newton ≈ 3.6 × 10⁴⁰ — coupling too strong
- Conservation: ∫δρ grows as t² (compensated by wavefront spike at r=ct)
- Scalar, not tensor (no gravitational waves)

## Code

`src/depletion.c` — self-contained solver (~800 lines)
- Tests 1-5: Static depletion (algebraic equilibrium, Poisson, universality)
- Test 6 (`-wave`): Causal wave equation with Models A (massless) and B (fluid)

Build: `cd src && make depletion`
Run: `../bin/depletion [-alpha A] [-scan] [-wave]`

## Key Numbers — Algebraic Model

| Quantity | Value | Meaning |
|----------|-------|---------|
| a (slope) | 1.421 | Nearly unchanged from σ-model (1.420) |
| E₂ | 51.30 | Twist energy (ρ-weighted) |
| E₄ | 51.64 | Skyrme energy (ρ-independent) |
| ΔQ | 0.257 | Total density deficit |
| ρ(0)/ρ₀ | 0.985 | Core depletion (1.5% at α=100) |
| δρ exponent | -5.95 | Decays as 1/r⁶ (not 1/r) |
| Φ exponent | -1.00 | Poisson integral is 1/r (mathematical) |
| E_overlap exponent | -5.93 | Physical interaction is 1/r⁶ |

## Key Numbers — Massless Wave Model (ds = c·dt)

| Quantity | Value | Meaning |
|----------|-------|---------|
| Q_eff = E₂/ρ₀ | 51.3 | Monopole of source ½|ω|² |
| δρ exponent | -0.989 | 1/r (LONG-RANGE!) |
| δρ/prediction ratio | 1.003 | Matches -Q_eff/(4πr) at r>5 |
| Total deficit | grows as t² | ½c²Q_eff t² verified to 0.4% |
| G_eff | 1/(16πρ₀²) ≈ 0.02 | From Q_eff²/(4πM²) |
| G_eff/G_Newton | 3.6×10⁴⁰ | Still too strong (nuclear scale) |

## Three Models Compared

| Model | Equation | Source monopole | δρ decay | Conservation |
|-------|----------|----------------|----------|-------------|
| Algebraic | δρ = -\|ω\|²/(4α) | N/A (local) | 1/r⁶ | N/A |
| Fluid wave | ∂²δρ/∂t² = c_s²∇²δρ + ½∇²\|ω\|² | 0 (Gauss) | 1/r⁶ | ∫δρ = 0 ✓ |
| **Massless wave** | **□δρ = -½\|ω\|²** | **E₂/ρ₀ ≠ 0** | **1/r** | **∫δρ grows** |
