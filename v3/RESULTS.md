# HFKT v3 — Results Registry

This document tracks which FOUNDATION.md open problems and caveats have been addressed by completed tasks. FOUNDATION.md itself is never modified — all results are recorded here.

Each entry references the relevant FOUNDATION.md section and the completed task file in `tasks/done/`.

---

## Resolved Open Problems

### Open Problem #1: Sourced Maxwell equations and charge derivation — RESOLVED

**FOUNDATION.md refs**: §9, §15, Caveat #2, "What HFKT Does Not Yet Have" (sourced Maxwell)

**Task**: `tasks/done/sourced-maxwell.md`

**Result**: POSITIVE (standard Skyrme model result, as expected)

The sourced Maxwell equation is derived via the Adkins-Nappi-Witten (1983) gauged U(1)_EM mechanism:

1. **Charge identification**: Q_EM = B/2 + I₃, where B = topological baryon number and I₃ = quantized isospin projection. Proton Q = 1.000000, Neutron Q = 0.000000.

2. **Charge is NOT classically topological**: Electric charge requires quantum collective coordinate quantization of the SU(2) isospin zero mode. The classical hedgehog has Q_class = B/2 = 1/2; integer charge comes from quantizing A(t) ∈ SU(2).

3. **Sourced Maxwell structure**: ∇·E = ρ_EM where ρ_EM = (1/2)B⁰ + (1/2)ρ_I(r). Far-field: E(r) → Q/(4πr²) (Coulomb). The isoscalar density ρ_I(r) = M(r)/Λ is normalized to ∫ρ_I 4πr² dr = 1.

4. **Linearized equation on hedgehog background**: □δF + M²(r,θ)·δF = 0 (ungauged). The Hessian M²(r,θ) is nonzero in the soliton core and vanishes in the far field (free Maxwell). Gauging introduces the source: □δF + M²·δF = J_EM.

5. **Convergence**: Hessian verified to O(ε²) convergence for central differences across 4 perturbation types (breathing, bivector monopole/dipole/quadrupole).

**Electromagnetic properties (σ-model)**:
| Quantity | This work | ANW target | Experiment |
|----------|-----------|------------|------------|
| Q_proton | 1.000000 | 1 | 1 |
| Q_neutron | 0.000000 | 0 | 0 |
| μ_p/μ_n | -1.5000 | -1.500 | -1.46 |
| r_E (σ-model) | 1.88 fm | — | 0.84 fm |
| r_E (massive) | 0.90 fm | 0.59 fm | 0.84 fm |
| F₁(0) | 1.000000 | 1 | 1 |
| Λ (MOI) | 141.55 | 141.6 | — |

**Key structural finding**: The Hessian for bivector perturbations is NONZERO and ANISOTROPIC at the soliton core. This means the soliton scatters EM-like excitations nontrivially. The source term J in the sourced equation comes from the gauging (external input), not from the ungauged Lagrangian. This matches standard Skyrme model physics.

---

## Resolved Caveats

### Caveat #2: Q is baryon number, not electric charge — RESOLVED

**Status**: Q remains baryon number. Electric charge is derived via Q_EM = B/2 + I₃ through gauged U(1)_EM ⊂ SU(2)_I. Charge quantization is quantum-mechanical (collective coordinate quantization), not classical-topological. This is the standard and expected Skyrme model result.

### Open Problem #2: Effective metric / post-Newtonian gravity — RESOLVED (NEGATIVE for gravity)

**FOUNDATION.md refs**: §18, §23, Open Problem #2, Caveat #1

**Task**: `tasks/done/gravity-null-rotor-metric.md`

**Result**: L₂ + L₄ gives a nontrivial BLV effective metric for null-rotors, but it is a nuclear-scale scattering effect, not gravity.

The v2 result P/m = 2 was correct but **only for K=0 breathing modes**. The K≥1 bivector modes (physical null-rotors from §23) have P/m < 2, giving v < c in the soliton core:

| Profile | v_min/c | Φ_min/c² | Φ_min (MeV) | θ_peak | Range (fm) |
|---------|---------|----------|-------------|--------|------------|
| Massless | 0.670 | -0.275 | -258 | 0.656 rad | 0.84 |
| Massive (m_π=0.398) | 0.625 | -0.305 | -286 | 0.911 rad | 0.56 |

**Why NOT gravity**: (1) Magnitude 10³⁵× too strong, (2) Yukawa range ~0.8 fm not 1/r, (3) No spin-2 mediator. This is the nuclear potential for meson-nucleon scattering, a standard Skyrme model result.

**No birefringence**: The two transverse null-rotor polarizations see identical speeds, consistent with GR's prediction but insufficient for gravity identification.

### Open Problem #3: Gravitational wave polarization — RESOLVED (NEGATIVE)

**FOUNDATION.md refs**: §18, Open Problem #3

**Task**: `tasks/done/gravity-null-rotor-metric.md`

**Result**: Only h=0 (scalar) and h=±1 (vector) propagating modes exist on the hedgehog background. **No h=±2 (tensor) modes** — the SU(2) field supports only spin ≤ 1 perturbations. This is inconsistent with LIGO's observation of tensor gravitational waves.

---

### Caveat #1: Gravity has no quantitative support — CLOSED (structural no-go)

**Status**: The monopole coupling analysis provides a definitive structural argument: for 1/R gravity, need a MASSLESS field with MONOPOLE coupling to T₀₀. The pion is massless but dipole-coupled (→ 1/R³). The breathing mode is monopole-coupled but massive/frozen (→ Yukawa). No field satisfies both requirements. Even hypothetically, a massless breathing mode would give a force 10⁴¹ × too strong and would destroy the soliton (λ→0 violates Derrick). This explains all prior null results (BLV, depletion, WZW, constraint, spectral, L₆).

### Open Problem #2 (continued): Depletion gravity / inter-soliton potential — RESOLVED (NEGATIVE)

**FOUNDATION.md refs**: §6, §18, Caveat #1

**Task**: `tasks/done/depletion-gravity.md`

**Result**: The σ-model constraint |q|=ρ₀ does NOT generate a 1/R gravitational potential between solitons.

The depletion hypothesis was: one soliton's field concentration forces a compensating deficit felt by a second soliton at large separation R, producing a 1/R potential absent from single-soliton analysis. Four independent tests rule this out:

| Test | Result | Implication |
|------|--------|-------------|
| A: Tail coefficient | A=4.275, T₀₀ ~ r⁻⁶ | No long-range energy overlap |
| B: Lagrange multiplier | λ ~ r⁻⁶ (not r⁻²) | No Poisson source → no 1/r |
| C: Manton formula | E_int = C·g(α)/R³ | Dipole-dipole, orientation-dependent |
| D: 3D product ansatz | n ≈ 3.9–4.5 (corrected) | Consistent with dipole, not 1/R |

The interaction is orientation-dependent (varies by 33% with isorotation angle α) and averages to zero over SU(2), ruling out universal gravity. Combined with the BLV null result (single-soliton, Yukawa range), this closes all perturbative paths to gravity from L₂+L₄.

### Gravity No-Go: Monopole Coupling — RESOLVED (STRUCTURAL NO-GO)

**FOUNDATION.md refs**: §18, Caveat #1

**Task**: `tasks/done/monopole-gravity-nogo.md`

**Result**: Definitive structural argument for why L₂+L₄ cannot produce 1/R gravity.

For 1/R gravity between solitons of mass M, need a mediator that is (a) MASSLESS and (b) MONOPOLE-COUPLED to T₀₀. Systematic enumeration of all channels:

| Channel | Massless? | Monopole? | Interaction |
|---------|-----------|-----------|-------------|
| Pion (π) | YES | NO (dipole, I=1) | 459/R³ |
| Breathing (σ) | NO (m→∞) | YES (g=M/ρ₀) | Yukawa |
| Weight (J,P) | N/A (ε₀²=0) | N/A | non-dynamical |
| Spin-2 | N/A | N/A | doesn't exist |

No channel satisfies both (a) and (b).

**Why the pion fails**: π_a = ρ₀ sin(f) x_a/r is an isovector (I=1) with angular dependence L=1. Its monopole coupling to isoscalar T₀₀ vanishes exactly: ⟨I=0|I=1|I=0⟩ = 0.

**Why the breathing mode fails**: g_σ = 2E₂/ρ₀ = M/ρ₀ (verified numerically to 0.04%). But m_σ = √(8λ)ρ₀ → ∞ in σ-model. Even hypothetically massless: M_P_eff = 32 MeV → force 10⁴¹ × too strong. And λ→0 destroys soliton (Derrick).

This single argument explains ALL prior negative results across v2 and v3.

---

## Open Items Remaining (from FOUNDATION.md)

### Still open (not yet addressed by any task):

- **Open Problem #4**: Hopfion composites
- **Open Problem #5**: Weight sector dynamics
- **Open Problem #6**: Particle spectrum
- **Open Problem #7**: Color and SU(3)
- **Open Problem #8**: Three generations
- **Open Problem #9**: Parity violation / weak force
- **Open Problem #10**: Born rule

### Still open caveats:

- **Caveat #3**: Weight sector is non-dynamical
- **Caveat #4**: Lattice topology is fragile
- **Caveat #5**: Hopfion composites are numerically unstable
- **Caveat #6**: Topological baryon number conservation vs. baryogenesis
- **Caveat #7**: Born rule is not derived

---

## Code Inventory

| File | Purpose | Status |
|------|---------|--------|
| `src/sourced_maxwell.c` | ANW EM current, charge densities, form factors | Complete |
| `src/hessian_check.c` | L₄ Hessian linearization verification | Complete |
| `src/field.c` | 3D energy/force computation (from v2) | Complete |
| `src/field.h` | Data structures (from v2) | Complete |
| `src/clifford.h` | Multivector algebra (from v2) | Complete |
| `data/profile_sigma_e1.dat` | Massless B=1 hedgehog profile | Input |
| `data/profile_massive_e1_mpi0.398.dat` | Massive pion B=1 profile | Input |
| `data/em_charge_densities.dat` | Output: ρ_p, ρ_n, B⁰, ρ_I | Output |
| `data/em_proton_field.dat` | Output: E-field, potential, Q_enc | Output |
| `data/em_form_factors.dat` | Output: F₁(q²) proton/neutron | Output |
| `data/hessian_pert[0-3].dat` | Output: Hessian radial profiles | Output |
| `src/nullrotor_metric.c` | BLV effective metric for null-rotors | Complete |
| `data/blv_metric.dat` | Output: phase velocities for 4 polarizations | Output |
| `data/blv_deflection.dat` | Output: deflection angles vs impact parameter | Output |
| `src/depletion_gravity.c` | Depletion gravity: tail, λ(r), Manton, 3D scan | Complete |
| `data/depletion_tail.dat` | Output: tail coefficient and T₀₀ falloff | Output |
| `data/depletion_lagrange.dat` | Output: Lagrange multiplier λ(r) | Output |
| `data/depletion_manton.dat` | Output: Manton analytical potential | Output |
| `data/depletion_3d_scan.dat` | Output: 3D product ansatz energy scan | Output |
| `src/monopole_coupling.c` | Gravity no-go: monopole coupling analysis | Complete |
| `data/monopole_coupling.dat` | Output: coupling integrals and masses | Output |
