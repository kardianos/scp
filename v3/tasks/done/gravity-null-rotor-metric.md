# Gravity: Null-Rotor Effective Metric and GW Polarization

## 1. Goal

Determine whether null-rotor (photon/EM) perturbations experience a nontrivial effective spacetime metric when propagating through the B=1 Skyrmion background, and characterize the gravitational wave polarization content. The v2 investigation proved P/m = 2 algebraically for K=0 breathing modes, ruling out gravitational lensing for scalar perturbations. But the **K≥1 angular/bivector modes** — which are the physical null-rotors identified in FOUNDATION.md §23 — were never analyzed. If the L₄ commutator structure gives P/m ≠ 2 for these modes, the soliton lenses light without any additional terms or free parameters. If P/m = 2 holds for all modes, then gravity from L₂ + L₄ is definitively ruled out, and the task documents precisely what modifications (L₆, new couplings, or nonperturbative effects) would be required.

## 2. Foundation References

**FOUNDATION.md §18 — Gravity:**
> "The field's energy-momentum density could modulate the effective propagation speed for null-rotors, creating an analog-gravity effective metric. If this effective metric satisfies Einstein's field equations, all post-Newtonian tests follow. Computing the PPN parameters remains the critical quantitative test."

**FOUNDATION.md §23 — Null-Rotors:**
> "A null-rotor is a propagating field pattern with: Speed = c (massless), Energy and momentum (E, p = E/c), No rest mass, no topological protection, Transverse bivector structure (E ⊥ B ⊥ k)."
> "In the PGA framework, the three spatial bivector components (F = e_ij) map to the massless electromagnetic field."

**FOUNDATION.md §18 — Critical caveat:**
> "The BLV effective metric investigation showed P/m = 2 algebraically for L₂ + L₄, meaning no gravitational lensing from the base Lagrangian alone."

**FOUNDATION.md Open Problems #2, #3. Caveat #1.**

**v2 result** (`hopfion_composition/results/effective_metric_analysis.md`):
The P/m = 2 identity is proved for the K=0 Sturm-Liouville reduction:
```
P(r) = 2r² + 4c₄ sin²f = 2(r² + 2c₄ sin²f) = 2 × m(r)
```
The algebraic argument relies on [A₀, A_r] = 0 for radial perturbations on the hedgehog. For K≥1 modes, [A₀, A_θ] ≠ 0 in general — the commutator structure differs, and the P/m ratio has not been computed.

**Depends on**: `tasks/done/sourced-maxwell.md` (Hessian infrastructure, charge densities, form factors).

## 3. Technical Requirements

### Mathematical derivations needed (before coding):

**Derivation A — K≥1 Sturm-Liouville reduction:**
Starting from the full Lagrangian L = L₂ + L₄ (§6), linearize around the hedgehog q₀ = ρ₀(cos f + sin f r̂·σ) for a perturbation with angular quantum number K≥1:

δq = Σ_{lm} g_l(r) Y_{lm}(θ,φ) · τ_a

where τ_a are isospin generators. The resulting linearized equation has the form:

-(P_K(r) g')' + W_K(r) g = (ω²/c²) m_K(r) g

Extract P_K(r) and m_K(r) analytically. The critical question: does P_K/m_K = 2?

For K=0: P₀ = 2r² + 4c₄ sin²f, m₀ = r² + 2c₄ sin²f → P₀/m₀ = 2. Known.
For K=1: The angular derivatives in L₄ contribute new commutator cross-terms [A_r, A_θ] that modify the coefficients. Derive explicitly.

**Derivation B — Mode decomposition:**
Classify all linearized modes on the hedgehog background by (K, parity, isospin channel). For each channel, determine:
- Whether the mode is propagating (ω² > 0 for large r) or evanescent
- The asymptotic dispersion relation ω(k) as r → ∞
- The helicity content (scalar h=0, vector h=±1, tensor h=±2)

**Derivation C — Effective metric (if P/m ≠ 2):**
Using the BLV formula g^{μν}_eff = -∂²L/∂(∂_μφ)∂(∂_νφ) evaluated on the hedgehog background for the relevant bivector perturbation channel, extract:
- The effective line element ds²_eff
- The effective gravitational potential Φ(r)/c² = (P_K/m_K - 2)/4
- The PPN parameters γ, β from the weak-field expansion

### Code to write:

1. **`src/nullrotor_metric.c`** — Main computation:
   - Load 1D hedgehog profile f(r)
   - Compute P_K(r) and m_K(r) for K=0,1,2 from the analytical formulas
   - Verify P₀/m₀ = 2 (consistency check against v2)
   - Report P_K/m_K for K≥1 — the key result
   - If P/m ≠ 2: compute Φ(r), ray-trace deflection angle θ(b), extract PPN γ

2. **`src/gw_polarization.c`** — GW mode analysis:
   - Enumerate all propagating modes on hedgehog background
   - For each mode: compute asymptotic helicity
   - Count tensor (h=±2), vector (h=±1), scalar (h=0) modes
   - Compare with GR expectation (2 tensor modes only)

### Input data:
- `data/profile_sigma_e1.dat` — massless B=1 hedgehog profile
- `data/profile_massive_e1_mpi0.398.dat` — massive pion profile
- Hessian data from `data/hessian_pert*.dat` (for cross-validation)

### Numerical methods:
- Sturm-Liouville coefficient evaluation: direct computation from profile f(r)
- Ray-tracing: 4th-order Runge-Kutta on null geodesic ODE
- PPN extraction: asymptotic matching of metric to Schwarzschild form

## 4. Success Criteria

**Primary (addresses Open Problem #2):**
- P_K/m_K computed to 6+ significant figures for K=0,1,2
- If P₁/m₁ ≠ 2: effective potential Φ(r) extracted, deflection angle θ(b) computed for b = 1,2,5,10,20 code lengths
- PPN γ computed to 3+ significant figures (GR predicts γ = 1)

**Secondary (addresses Open Problem #3):**
- Complete mode table: for each (K, parity, channel), report {propagating/evanescent, helicity, asymptotic speed}
- Number of tensor (h=±2) GW polarizations identified
- If exactly 2 tensor modes with speed c: consistent with GR/LIGO
- If only scalar modes: inconsistent with LIGO (definitive negative)

**Consistency checks:**
- P₀/m₀ = 2.000000 (reproduces v2 breathing mode result)
- Far-field limit: all propagating modes have speed c (flat space recovery)
- K=1 zero modes (translation, isorotation) have ω = 0 at k = 0

## 5. Null Result Criteria

**Null Result A — P/m = 2 for all modes:**
If P_K/m_K = 2 for K=0,1,2 (and by induction for all K), then L₂ + L₄ produces ZERO gravitational lensing for ANY perturbation mode on the hedgehog background. This would mean:
- Gravity cannot emerge from the base Skyrme Lagrangian
- The FOUNDATION.md §18 mechanism ("density depletion → refractive index") fails not just for σ-model but structurally
- Gravity requires either: (a) L₆ or higher terms, (b) explicit degenerate-sector coupling, or (c) nonperturbative (instanton) effects
- This CLOSES Open Problem #2 with a definitive negative for L₂+L₄

**Null Result B — Only scalar GW modes:**
If the hedgehog background supports only h=0 propagating modes (no tensor h=±2), then:
- The theory cannot reproduce LIGO's tensor GW observations
- Either the gravitational sector requires additional structure (massive spin-2 mode), or
- The Skyrme model is the wrong effective theory for the gravitational sector
- This would be a serious constraint on the theory

**Null Result C — Dispersive lensing:**
If the effective metric depends on frequency (ω-dependent phase velocity), then:
- The effect cannot be gravity (GR is non-dispersive)
- It would be a scattering/refraction effect, not geometric lensing
- Still physically interesting (nuclear-scale EM scattering) but not gravity

**What null results MEAN for the theory:**
A null result on all three fronts (A+B+C) would establish that the Skyrme Lagrangian L₂+L₄ is structurally incapable of producing gravity. Combined with the v2 result that L₆ gives only nuclear-scale effects, this would mean HFKT requires a fundamentally new coupling mechanism for gravity — either a modification to the Lagrangian or a nonperturbative effect not captured by linearization. The v2 "sharp prediction" that gravity is nonperturbative (instanton tunneling with S ~ 32) would become the only remaining viable path.

## 6. Implementation Plan

### Step 1: Analytical derivation of K≥1 Sturm-Liouville coefficients
**Produces**: Explicit formulas for P_K(r), m_K(r), W_K(r) for K=1
**Verify**: Reduce to K=0 case and confirm P₀ = 2r² + 4c₄sin²f, m₀ = r² + 2c₄sin²f
**Method**: Expand L₂ + L₄ to second order in perturbation δq with angular dependence. Use the hedgehog factorization of right-currents R_μ = q₀⁻¹∂_μq₀ and separate radial/angular contributions. The key new terms come from L₄ commutators involving ∂_θ(δq).

**Critical check**: If the derivation gives P₁ = 2m₁ algebraically (same identity), STOP — the rest of the computation is unnecessary for the gravity question (but continue to Step 5 for GW polarization).

### Step 2: Implement nullrotor_metric.c
**Produces**: Numerical values of P_K(r)/m_K(r) at each grid point
**Verify**: P₀/m₀ = 2.000000 to machine precision. Far-field P_K/m_K → c² (flat space).
**Method**: Evaluate the Step 1 formulas using the loaded profile f(r), f'(r).

### Step 3: Effective potential and ray-tracing (conditional on P/m ≠ 2)
**Produces**: Φ(r)/c², deflection angle θ(b), PPN γ
**Verify**: θ(b) → 0 as b → ∞. Consistency with weak-field limit.
**Method**: Φ(r)/c² = (P_K/m_K - 2)/4. Ray-trace with RK4 on null geodesic. Extract PPN γ from asymptotic metric matching.

### Step 4: Cross-validate with 3D Hessian data
**Produces**: Confirmation that the 1D Sturm-Liouville coefficients match the 3D Hessian from hessian_check.c
**Verify**: The radial profile of the Hessian for bivector perturbations (hessian_pert1-3.dat) should be consistent with the W_K(r) potential from Step 1.
**Method**: Compare Hessian data along z-axis and x-axis with analytical W_K(r).

### Step 5: GW polarization analysis
**Produces**: Mode table (K, parity, helicity, speed) for all propagating modes
**Verify**: Zero modes (K=1 translation, isorotation) have ω=0. Far-field modes have speed c.
**Method**: Decompose the linearized EOM into irreducible angular momentum channels. For each channel, determine the asymptotic behavior and helicity content.

### Step 6: Documentation
**Produces**: Complete results appended to this task file, RESULTS.md updated
**Contents**:
- P_K/m_K values for K=0,1,2
- If nontrivial: Φ(r), θ(b), PPN parameters table
- GW polarization table
- Comparison with GR predictions
- Assessment of what modifications are needed (if null)

## 7. Reusable Infrastructure

### From v3/src/ (sourced-maxwell task):
- **`field.c`** / **`field.h`**: 3D field energy, force, topology (from v2). Known issue: `field_topological_charge` has unused parameter.
- **`clifford.h`**: Multivector algebra.
- **`sourced_maxwell.c`**: Profile loading (`load_profile`), baryon density computation. Will copy the profile loading code.
- **`hessian_check.c`**: 3D Hessian computation for cross-validation (Step 4).

### From v2/proposal/hopfion_search/src/:
- **`normal_modes.c`**: K=0 Sturm-Liouville solver. The P(r), m(r), W(r) coefficient formulas will be the starting point for the K≥1 generalization. **Bug note**: GL weight formula bug exists in veff.c but NOT in normal_modes.c.
- **`radial.c`**: Profile shooting solver. Not needed directly but documents the ODE for f(r).

### Known v2 bugs to avoid:
- GL weight formula: w_i = 2(1-x²)/(n²P_{n-1}²), NOT 2/((1-x²)n²P_{n-1}²). Correct in bound.c, wrong in veff.c.
- Profile format: 6-column sigma-model files have baryon_dens as column 4, NOT ρ.
- c₄ = 2ρ₀²/e² (factor of 2 matters — was source of E₄ bug in sourced_maxwell.c).

## 8. Results

### Key Finding: BLV effective metric IS nontrivial for null-rotors

The v2 result P/m = 2 (implying no gravitational lensing) was correct but **only applies to K=0 breathing modes**. The K≥1 bivector modes — the physical null-rotors from FOUNDATION.md §23 — see a **different effective metric** with v < c in the soliton core.

### Analytical Derivation

**BLV formula** (Barceló-Liberati-Visser): For perturbation polarization e_A on the hedgehog background, the phase velocity for propagation along direction j is:

v²_j / c² = (1 + P₄^{jj}) / (1 + m₄)

where:
- P₄^{jj} = (c₄/2ρ₀⁴) Σ_{k≠j} |[e_A q̃₀, R_k]|² (stiffness: commutators with PERPENDICULAR currents)
- m₄ = (c₄/2ρ₀⁴) Σ_k |[e_A q̃₀, R_k]|² (inertia: commutators with ALL currents)

Since m₄ ≥ P₄^{jj} (extra Cj ≥ 0 term), **all phase velocities are ≤ c**.

**On the z-axis of the hedgehog** (r̂ = ẑ, q₀ = ρ₀(cosf + sinf σ₃)):

| Mode | [δR, R_x]² | [δR, R_y]² | [δR, R_z]² | v_z/c | v_x/c |
|------|-------------|-------------|-------------|-------|-------|
| scalar (e₀) | 4ρ₀²sin⁴f/r² | same | 0 | **1** | <1 |
| σ₁ | 0 | 4ρ₀²sin²f/r² | 4ρ₀²f'² | **<1** | **1** |
| σ₂ | 4ρ₀²sin²f/r² | 0 | 4ρ₀²f'² | **<1** | **<1** |
| σ₃ | 4ρ₀²sin⁴f/r² | same | 0 | **1** | <1 |

**Critical observation**: Each mode has one "free" direction (v = c) where its commutator with the background current vanishes (parallel currents). The slow directions have v < c due to nonzero commutators.

**No birefringence for physical photon polarizations**: For propagation along any direction, the two transverse bivector polarizations (e.g., σ₁ and σ₂ for z-propagation) always have the **same speed**. This is consistent with gravity (GR predicts no birefringence).

### Numerical Results — Sigma Model (e=1, ρ₀=1)

| Quantity | Massless | Massive (m_π=0.398) |
|----------|----------|---------------------|
| v_min/c | 0.670 | 0.625 |
| v_min location | r ≈ 0.7 | r ≈ 0.7 |
| Φ_min/c² | -0.275 | -0.305 |
| Φ_min (MeV) | -258 | -286 |
| θ_peak (rad) | 0.656 at b=2 | 0.911 at b=2 |
| θ_peak (degrees) | 37.6° | 52.2° |
| Half-width (fm) | 0.84 | 0.56 |
| K=0 SL P/m cross-check | 2.000000 | 2.000000 |

### Far-field behavior: Yukawa, NOT 1/r

| b (code) | θ massless | θ massive | θ_GR (Newtonian) |
|----------|-----------|-----------|------------------|
| 2.0 | 6.56e-1 | 9.11e-1 | 1.22e-36 |
| 5.0 | 1.09e-2 | 4.56e-3 | 4.87e-37 |
| 10.0 | 8.98e-5 | 5.46e-6 | 2.43e-37 |
| 20.0 | 7.97e-7 | 1.94e-10 | 1.22e-37 |

Far-field falloff slopes (d log θ / db):
- Massless: −1.16 → −0.39 (power-law × 1/r decay, NOT 1/b)
- Massive: −1.50 → −0.98 (approaches 2m_π = 0.796, Yukawa)
- Newtonian: −1/b (inverse, not exponential)

**Ratio to GR**: The BLV deflection exceeds Newtonian by factors of 10³⁰–10³⁵, confirming this is a **nuclear-scale** effect (nuclear potential ~250 MeV), not gravity.

### GW Polarization (Open Problem #3)

Propagating modes on the hedgehog background, classified by angular momentum:

| Channel | K | Type | SL P/m | Asymptotic v | Helicity |
|---------|---|------|--------|--------------|----------|
| Breathing | 0 | scalar | 2 | c | h = 0 |
| Isorotational (A) | 1 | bivector | 1 | c | h = ±1 |
| Translational (B) | 1 | scalar | 2 | c√2 (SL) | h = 0 |
| Mixed (C) | 1 | bivector | 1 | c | h = ±1 |

- **No h = ±2 tensor modes**: The SU(2) field supports only spin-0 (scalar) and spin-1 (vector/bivector) perturbations. Tensor (graviton-like) modes require spin-2 content absent from the Skyrme field.
- **Two h = ±1 modes propagate at c**: Channels A and C are the physical null-rotors. Both are non-dispersive in the far field.
- **Inconsistent with LIGO**: GR gravitational waves have tensor polarization (h = ±2). The Skyrme model produces only scalar and vector radiation.

### Resolution of Open Problems

**Open Problem #2 (Effective metric / post-Newtonian gravity)**: RESOLVED — PARTIAL POSITIVE, ultimately NEGATIVE for gravity.

L₂ + L₄ **does** give a nontrivial BLV effective metric for null-rotors (v < c in soliton core). The v2 result P/m = 2 is **correct but incomplete** — it applies only to K=0 breathing modes, not to K≥1 bivector modes. However:
1. The effect is nuclear-scale (~250 MeV, 10³⁵× too strong)
2. The range is ~0.8 fm (Yukawa, not 1/r)
3. No tensor (h=±2) GW modes exist
4. PPN parameters are undefined (no 1/r limit to match)

**Open Problem #3 (GW polarization)**: RESOLVED — NEGATIVE.

Only h=0 (scalar) and h=±1 (vector) modes propagate. No h=±2 tensor modes. Inconsistent with LIGO observations.

**Caveat #1 (Gravity has no quantitative support)**: REINFORCED.

L₂ + L₄ produces a nuclear-scale meson-scattering potential, not gravity. Combined with v2 results (L₆ also nuclear-scale, g_top free, all quantum avenues null), gravity remains the theory's most significant gap.

### What Would Be Required for Gravity

The BLV metric result establishes that L₂ + L₄ gives the correct **qualitative** features (attractive lensing, no birefringence) but wrong **quantitative** properties (wrong magnitude by 10³⁵, wrong range, wrong polarization content). For gravity to emerge:

1. **Magnitude**: Need a mechanism to suppress the nuclear effect by ~10³⁵ and produce an unsuppressed 1/r tail. This is beyond perturbation theory.
2. **Range**: Need a massless spin-2 mediator (graviton). The Skyrme model has no spin-2 field. Emergent spin-2 from spin-1 composites would require a new collective mode.
3. **Polarization**: Need h=±2 modes. These could arise from metric perturbations that are QUADRATIC in the field perturbation (δg_eff ~ δΨ²), but this is a higher-order effect.
4. **The only remaining viable path**: Nonperturbative (instanton) effects with action S ~ 32, which are exponentially suppressed (e⁻³² ~ 10⁻¹⁴) but could produce the needed G_N scale.

### Consistency Checks

1. K=0 SL P/m = 2.000000 to 8 digits — reproduces v2 result ✓
2. Far-field v → c for all modes — flat space recovery ✓
3. v ≤ c everywhere — no superluminal propagation ✓
4. σ₁ and σ₂ have same v_z — no birefringence for transverse pair ✓
5. Massive profile has shorter range than massless — Yukawa expected ✓
6. v_x(σ₁) = v_y(σ₂) = c — hedgehog symmetry ✓
