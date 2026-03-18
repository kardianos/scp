# T10 Field Study Results

## T10B: Effective Metric Extraction — MIXED SIGNAL

### Method
Computed the linearized mass matrix M²_ab(x) = m²δ_ab + ∂²V/∂φ_a∂φ_b
on the equilibrated braid background. Eigenvalues give effective mass for
each mode; refractive index n = sqrt(k²+M²)/sqrt(k²+m²_vac).

### Key Results
- **Center (r<2.5)**: ALL eigenvalues ABOVE vacuum → ALL modes slower → ATTRACTIVE
  - M² = (2.74, 2.78, 3.88) vs vacuum 2.25
  - n_eff = (1.013, 1.014, 1.044) → 1-4% speed reduction
- **Shell (r=3-5)**: Minimum eigenvalue DROPS BELOW vacuum → one mode FASTER → REPULSIVE
  - M²_min ≈ 1.3-1.5 at r≈4 → n_min ≈ 0.97 for that mode
  - Other two modes still above vacuum
- **Far field (r>8)**: Returns to vacuum M²=2.25 everywhere
- **Anisotropy**: 41% at center → genuine tensor structure (3 different eigenvalues)
- **Time-averaged** (20 samples): Same qualitative structure, slightly smoothed

### Interpretation
NOT like gravity. Gravity requires ALL modes to slow down uniformly (scalar
potential Φ/c² < 0 everywhere). Here we have:
- Birefringence: different polarizations propagate at different speeds
- Sign change: attractive at center, one mode repulsive at r≈4
- Short range: Yukawa decay (range ~6), not 1/r
- The triple-product Hessian has INHERENT anisotropy (3 distinct fields)

## T10A: Wave-Wave Scattering — NEGATIVE

### Method
Collide wave packets head-on through the triple-product coupling.
Test phase dependence, amplitude dependence, and cross-field scattering.

### Key Results
- **Same-component collision** (P=0): Pure superposition, no scattering. ✓
- **Three-field collision**: |P| ≈ 10⁻¹⁷ even at collision → effectively ZERO
  - L/R asymmetry = -1/3 exactly → just momentum ratio (2 left, 1 right)
  - NO phase dependence: asymmetry varies by <0.05% across 5 phase offsets
  - NO amplitude dependence up to amp=1.0 (then very slight at amp=2.0)
- **Perpendicular collision**: 50/50 split → no preferred scattering angle

### Interpretation
The triple-product coupling is TOO WEAK for wave-wave scattering.
Plane waves don't overlap long enough to accumulate significant interaction.
The coupling needs the PERSISTENT, STRUCTURED overlap of a braid.
→ EM at the field level is NOT present in this model.

## T10C: Dispersion Relation — NO TACHYONIC INSTABILITY

### Method
Computed Hessian eigenvalues analytically for uniform condensate φ_a = v.

### Key Results
- At v=0: M² = (2.25, 2.25, 2.25) → standard massive KG
- At v=0.3: M² = (0.82, 2.56, 2.56) → large mass reduction for one mode
  - This is the SYMMETRIC mode (all fields in phase)
  - The two ANTISYMMETRIC modes have INCREASED mass
- At v≥0.5: κ regularization kicks in, M² returns toward vacuum
- NO TACHYONIC INSTABILITY at any amplitude
  - Minimum M²_min = 0.82 at v≈0.3, never goes negative
  - κ=50 prevents the instability that μ=-41.3 alone would cause

### Mode Structure
- Symmetric mode (1,1,1)/√3: softest, mass drops with v
- Two antisymmetric modes: mass increases with v
- The symmetric mode is the one that forms braids (all fields in phase, rotating)

### Interpretation
Mass does NOT emerge from band structure in this model. The m=1.5 mass
term is irreducible — it's a Lagrangian parameter, not an emergent quantity.
However, the effective mass CAN be modified by ~60% inside a condensate,
which explains why the braid has different internal dynamics.

## T10E: Impulse Response (Green's Function) — CHANNEL MIXING DETECTED

### Method (v2 — differential)
Run two copies of braid simulation (with/without impulse), measure difference.
This cleanly isolates the impulse response from braid dynamics.

### Key Results
- **At braid center**: 3.7% channel mixing (impulse in φ₀ excites φ₁,φ₂)
- **At r=5**: 1.9% mixing (weaker, as expected)
- **At r=10 (vacuum)**: 0.04% mixing (numerical noise only)
- **Signal speed**: ≈c (lattice artifacts at high k give apparent v>c)
- **Isotropy**: vz/vx ≈ 1.0 at center (cylindrical symmetry preserved)

### Interpretation
The braid acts as a MODE CONVERTER: an excitation in one field partially
converts to the other fields as it passes through the braid. This is the
microscopic mechanism behind the "torsion" diagnostic — it's not EM, but it
IS a field-mediated interaction. The conversion rate (3.7%) is small but
nonzero, and localized to the braid core.

## Overall Assessment

| Property | Status | Evidence |
|----------|--------|----------|
| Gravity (1/r potential) | ✗ NOT PRESENT | Yukawa range, birefringent, sign-changing |
| EM (phase-dependent force) | ✗ NOT PRESENT | No wave-wave scattering |
| Mass emergence | ✗ NOT PRESENT | Mass is a Lagrangian parameter |
| Stable soliton | ✓ PRESENT | fc=0.96, survives perturbation |
| Mode coupling | ✓ PRESENT | 3.7% channel mixing in braid |
| Tensor structure | ✓ PRESENT | 41% anisotropy in effective metric |

## What Would Need to Change

For gravity: Need a potential where the Hessian is positive semi-definite
EVERYWHERE (all eigenvalues > vacuum). The triple-product V(P) inherently
creates sign-changing, anisotropic corrections. A potential of the form
V(ρ) = V(φ₀²+φ₁²+φ₂²) would give isotropic, same-sign corrections, but
wouldn't support the triple-product topology.

For EM: Need a gauge symmetry (U(1) or similar). The three real scalar
fields have no gauge invariance. Could upgrade to complex fields or
add a gauge field.

For mass emergence: Need a mechanism where m² can be driven to zero in some
limit, then regenerated by the field configuration. The κ-regularized
triple product doesn't provide this.

## T10F: Massless Mediator Test — NEGATIVE (mass essential for confinement)

### Method
Test mixed-mass configurations to see if a massless mode can mediate 1/r forces.

### Key Results
| Config | m₀ | m₁,m₂ | fc | Far-field power | Notes |
|--------|-----|--------|------|-----------------|-------|
| standard | 1.5 | 1.5 | 0.88 | ~exponential | control |
| m0_zero | 0.0 | 1.5 | 0.39 | r^(-2.5) | braid partially lost |
| m0_small | 0.1 | 1.5 | 0.25 | r^(-3.4) | worse than massless |
| m0_half | 0.75 | 1.5 | 0.53 | r^(-1.9) | intermediate |
| all_small | 0.3 | 0.3 | 0.13 | nearly flat | field fills box |

### Interpretation
Making one field massless DESTROYS the braid (fc: 0.88 → 0.39). The massless
field freely radiates away. The far field decays as r^(-2.5), not 1/r.
**Mass is essential for confinement. This is the fundamental tension: mass confines
the soliton but prevents 1/r gravity.**

## T10G: Complex Field Upgrade — POSITIVE (charge-carrying soliton!)

### Method
Upgrade from 3 real to 3 complex scalars: ψ_a = (φ_a + iχ_a)/√2.
Potential V = (μ/2)|ψ₀ψ₁ψ₂|²/(1+κ|ψ₀ψ₁ψ₂|²).
This has U(1)² global symmetry with two conserved charges.

### Key Results
- **Braid survives**: fc = 0.91 (vs 0.88 for real case)
- **Carries charge**: Q₀ ≈ 3, Q₁ ≈ 3, Q₂ ≈ 14
- **Charges localized**: core charges are larger than total (edge radiation carries opposite charge)
- **Charges not individually conserved**: only Q₀+Q₁+Q₂ constraint from triple product

### Interpretation
The complex upgrade is viable! The braid naturally becomes a charge-carrying soliton.
This is the first step toward EM: if we promote the global U(1) to a gauge symmetry
(add A_μ field), we get a charged soliton coupled to EM.

But: this still doesn't give us gravity. And adding a gauge field is "putting EM in by hand."

## T10D: Kibble-Zurek Defect Formation — NEGATIVE (no spontaneous braids)

### Method
Random initial field with controlled amplitude, cooled via viscous damping.
36 configurations tested: 4 damping rates x 3 amplitudes x 3 random seeds.
Grid: N=96, L=20, T_cool=500.

### Key Results
- ALL configurations: fc → 0, n_peaks → 0. No localized structures form.
- Energy GROWS despite damping (mu<0 drives triple-product instability)
- Eventually saturates at uniform high-energy state (E ≈ 3.6M at gamma=0.05)
- No Kibble-Zurek defects at any cooling rate or amplitude

### Interpretation
Braids do NOT form spontaneously from random initial conditions.
The triple-product coupling drives the field toward a uniform high-energy state,
not toward localized structures. Braids require specific initialization (helical
geometry with matched phases). This means braids are NOT natural topological
defects of this field theory — they are metastable artifacts of careful preparation.

## T10H: Two Charged Complex Braids — NEGATIVE (charges don't matter)

### Method
Two complex braids at D=15 in L=30 box. Three configs: same charge, opposite charge,
no charge (real only). Track separation over T=200.

### Key Results
| Config | ΔD | Direction |
|--------|-----|-----------|
| Same charge | -5.71 | attraction |
| Opposite charge | -6.41 | attraction |
| No charge | -5.21 | attraction |

All three attract similarly (within ~20%). The U(1) charges play NO role in the force.
The attraction is purely from Yukawa field overlap (massive fields).

### Interpretation
Global U(1) symmetry does NOT produce EM-like forces. A GAUGE field (A_μ) is needed
to mediate long-range charge-dependent forces. The complex upgrade gives the braid
charges but no mechanism to transmit those charges at long range.

## Next Steps
1. T9 (substrate): What micro-structure produces V(P)?
2. Model modifications:
   - Complex fields + GAUGE field for EM (A_μ needed, not just global U(1))
   - Need fundamentally new mechanism for gravity (massless mediator)
   - Consider whether V6/V7 density wave approach can be combined with braid model
3. Consider entirely different framework: gauge-Higgs system that PRODUCES braids
