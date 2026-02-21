# Gravity No-Go: Monopole Coupling Analysis

## 1. Goal

Provide a concise, definitive structural argument for why the Skyrme model (L₂+L₄) cannot produce 1/R gravity. The argument: for 1/R potential between two solitons of mass M, you need a MASSLESS field that couples to the MONOPOLE (L=0) of T₀₀. Systematically test all propagating channels.

## 2. Foundation References

**FOUNDATION.md §18 — Gravity:**
> "The field's energy-momentum density could modulate the effective propagation speed for null-rotors."

**Caveat #1 — Gravity has no quantitative support:**
> All tested paths give nuclear-scale or null results.

**Depends on**:
- `tasks/done/gravity-null-rotor-metric.md` (BLV metric: nuclear Yukawa, not 1/r)
- `tasks/done/depletion-gravity.md` (inter-soliton: dipole 1/R³, not 1/R)

This task synthesizes both prior results into a single structural argument.

## 3. Technical Requirements

### The no-go argument (analytical):

For 1/R gravitational potential: need mediator φ with:
(a) Massless propagator: G(r) = 1/(4πr)
(b) Monopole source: ∫ φ·T₀₀ d³x ≠ 0

Test all Skyrme model channels:

| Channel | Massless? | Monopole? | Status |
|---------|-----------|-----------|--------|
| Pion (π_a) | Yes (m_π=0) | No (isovector dipole) | → 1/R³ |
| Breathing (σ) | No (m_σ→∞) | Yes | → Yukawa |
| Weight (J,P) | N/A (ε₀²=0) | N/A | Non-dynamical |
| Spin-2 tensor | Does not exist | N/A | SU(2) has j≤1 |

### Numerical verification:

Compute from the B=1 hedgehog profile:
1. Energy monopole M = E₂ + E₄ (confirm = total soliton mass)
2. Pion dipole moment D = (4π/3)ρ₀ ∫ sin(f) r² dr
3. Pion monopole coupling = 0 (isovector × spherical = 0, verify numerically)
4. Breathing source integral g_σ = 2E₂/ρ₀ (monopole coupling strength)
5. Breathing mass m_σ = √(8λ)ρ₀ for various λ
6. Hypothetical Planck mass if σ were massless: M_P_eff = √(4π)ρ₀
7. Ratio to real Planck mass → factor ~10⁵³ too strong

### Code to write:
- `src/monopole_coupling.c` — single file, 1D profile integrals only (no 3D grid)

### Input data:
- `data/profile_sigma_e1.dat` — massless B=1 hedgehog

## 4. Success Criteria

- Pion monopole coupling verified = 0 (to numerical precision)
- Breathing source integral g_σ = M/ρ₀ ± 1% (using virial E₂=E₄=M/2)
- Hypothetical Planck mass ratio computed, confirming ~10⁵³ mismatch
- Clean summary table of all channels

## 5. Null Result Criteria

**Expected null result**: No massless monopole mediator exists in the Skyrme model.

This provides the structural explanation for all prior negative gravity results:
- BLV metric (nuclear Yukawa) → because breathing mode is massive
- Depletion (1/R³ dipole) → because pion is dipole-coupled, not monopole
- WZW / constraint / spectral (all null) → because no additional massless fields

## 6. Implementation Plan

### Step 1: Compute energy integrals from profile
### Step 2: Compute pion multipoles (dipole moment, verify monopole = 0)
### Step 3: Compute breathing mode coupling and hypothetical gravity
### Step 4: Print summary table
### Step 5: Documentation — move to done/, update RESULTS.md

## 7. Reusable Infrastructure

- Profile loading from `nullrotor_metric.c` / `depletion_gravity.c`
- Known: c₄ = 2ρ₀²/e², profile cols 1-3 only
- **Bug found**: For 1D radial E₄ integral, use (1/e²) NOT c₄. The c₄ = 2ρ₀²/e² includes a factor for the 3D commutator normalization.

## 8. Results — NEGATIVE (No-Go Confirmed)

**Status**: COMPLETE. Gravity is structurally absent from L₂+L₄.

### Core Argument

For 1/R gravity between solitons, need a mediator satisfying BOTH:
- (a) MASSLESS: propagator G(r) = 1/(4πr)
- (b) MONOPOLE-COUPLED: source ∫ φ·T₀₀ d³x ≠ 0

| Channel | Massless? | Monopole? | Gravity? | Interaction |
|---------|-----------|-----------|----------|-------------|
| Pion (π) | YES | NO (dipole) | NO | 459/R³ (dipole-dipole) |
| Breathing (σ) | NO (frozen) | YES | NO | Yukawa (range→0) |
| Weight (J,P) | N/A | N/A | NO | Non-dynamical |
| Spin-2 | N/A | N/A | NO | Doesn't exist |

No channel satisfies both (a) and (b). Gravity is absent.

### Quantitative Results

| Quantity | Value | Verification |
|----------|-------|-------------|
| M (soliton mass) | 103.14 code = 938 MeV | E/E_FB = 1.232 ✓ |
| E₂/E₄ (virial) | 0.9992 | Should be 1.000 ✓ |
| A (tail coefficient) | 4.276 | Consistent with Part A ✓ |
| g_σ (breathing source) | 103.10 = 2E₂/ρ₀ | = M/ρ₀ at virial ✓ |
| g_σ/M | 0.9996 | Confirms g_σ = M/ρ₀ ✓ |

### Pion Channel: Massless but Dipole

- Pion field: π_a = ρ₀ sin(f) x_a/r (isovector, L=1)
- Monopole coupling: ∫ π_a dΩ = 0 (odd angular integral, exact)
- Dipole coefficient: C = 8πρ₀²A² = 459.5
- Interaction: E_int = C·g(α)/R³ (orientation-dependent, ⟨SU(2)⟩=0)

### Breathing Channel: Monopole but Massive

- Source integral: g_σ = 2E₂/ρ₀ = M/ρ₀ (couples to total energy)
- E₄ is ρ-independent (right-current A_μ = q̃∂q/|q|² cancels ρ)
- Mass: m_σ = √(8λ)ρ₀ → ∞ in σ-model limit
- At minimum viable λ=8000: m_σ=253 code, range=0.002 fm
- **Hypothetical (if massless)**:
  - M_P_eff = 3.55 code = 32.3 MeV
  - M_P_real/M_P_eff = 3.8 × 10²⁰
  - Force ratio: (M_P/M_P_eff)² = 1.4 × 10⁴¹ (times too strong)
  - AND: making σ massless (λ→0) destroys the soliton (Derrick)

### Explanatory Power

This single argument explains ALL prior negative gravity results:
- BLV metric (Yukawa, not 1/r) → σ is massive
- Depletion (1/R³, not 1/R) → pion is dipole-coupled
- WZW/constraint/spectral (null) → no additional massless fields
- L₆ sextic (nuclear-scale) → same mechanism, no new massless channel
- All v2 gravity avenues (null) → g_top is free, no mechanism to fix it

### Conclusion

Gravity from L₂+L₄ is structurally impossible. Any viable gravity mechanism requires physics beyond the Skyrme Lagrangian: either new massless fields, nonperturbative effects (instanton tunneling, S~32), or coupling to additional sectors not present in Cl⁺(3,0,1).
