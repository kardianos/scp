# SCP Future Directions and Open Questions

**Purpose**: Tracks future research directions, unresolved questions, and
exploratory paths. Items marked [CONFIRMED] have been tested under the current
V50/C4 equations. Items marked [V44-ERA] were tested under earlier equations
and may need retesting. Maintain with CONCEPT.md and EM_THEORY.md.

**Current equation set**: V50/C4 (Cosserat strain α=0.1, curl²-hardening β=0.5).
All new experiments use pre-converged proton templates, NOT analytical seeds.

---

## Critical Priority — V52 Quantitative Measurements

These are detailed in `v52/PLAN.md`. Summary:

### F1: Force Power Law — F(D) exponent [V44-ERA, needs retest]

V33 measured F ∝ 1/D^1.8 for braids (mixed gravity+EM). Newton requires n=2.
Need proton-proton measurement with C4 equations at D = 15..50.
**V52 Test 1**: 6 runs × 10 min GPU.

### F2: Depletion Profile δρ(r) [V44-ERA, needs retest]

V34 measured δρ ∝ 1/r^1.2 for a z-aligned braid. Need isotropic proton
measurement with C4. If n→2, the force is automatically inverse-square.
**V52 Test 2**: 1 run × 35 min GPU.

### F25: Gravitational Coupling Constant C_proton [PARTIALLY CONFIRMED — V51]

V51 gradient test (C4, steep gradient 0.15/0.05) shows proton drifts +6
code units toward low ρ over T=400. Drift rate ~0.015/t. Need multiple
gradient strengths to confirm F ∝ ∇ρ linearity and extract C precisely.
**V52 Test 3**: 3 runs × 35 min GPU (or use V51 data).

### F6: Charge-Dependent Force [V44-ERA, needs retest]

V34 measured same-winding +27%, opposite -57% for braids. Need UUD+UUD vs
UUD+UDD comparison with C4.
**V52 Test 4**: 2 runs × 20 min GPU. Requires UDD template generation.

### F26: Equivalence Principle [UNTESTED]

Test whether different-mass particles (braid, proton, deuterium) experience
the same gravitational acceleration in the same gradient. If a = const
regardless of mass, equivalence holds.
**V52 Test 6**: 3 runs × 35 min GPU.

### F_NEW: Emergent Gauss' Law [UNTESTED]

Test whether ∮ S·dA is radius-independent around a proton (S = Poynting
flux from θ radiation). Equivalently, measure θ_rms²(r) on shells —
if ∝ 1/r², Gauss' law holds for the radiation field. This is the
foundation for 1/r² Coulomb force.
**V52 Test 7**: Reuses Test 2 data + analysis tool.

---

## V51 Results (March 2026)

### Proton-Proton Collision [CONFIRMED — V51]

Two pre-converged UUD protons at D=25, v=±0.1c. V50/C4 equations.
- Protons condense by t=50, maintain distinct cores
- Closest approach D≈19 at t≈125
- **Bound state formation**: protons oscillate in separation (8-25 code units)
- Slow inspiral from θ radiation losses (25→17 over T=300)
- Cosserat hardening prevents merger even at closest approach
- Details: `v51/RESULTS.md`

**Open**: A higher-energy collision (v=0.3c or 0.5c) should be tested to
find the threshold for proton disruption vs bound-state formation.

### Proton Gravitational Drift [CONFIRMED — V51]

Single proton in steep gradient (A_high=0.15, A_low=0.05). V50/C4 equations.
- Proton drifts +6 code units toward low ρ over T=400
- Energy conservation: -0.032% over T=500
- Gradient persists (verified by slab measurement)
- Consistent with V43 results under different equations
- Details: `v51/RESULTS.md`

### Analytical Seed Failure [CONFIRMED — V51]

The V44 analytical seed (gen_proton_analytical, 512³) does NOT condense
into a proton under C4 equations. The Cosserat strain and hardening terms
prevent the six raw braids from merging. All production runs must use
pre-converged templates from V43 (evolved without C4 terms).

**Open**: Can protons form spontaneously under C4? If not, what formation
pathway works? Lower α/β during condensation, then raise to production
values? This is a fundamental question about whether the C4 equations
support particle formation or only particle stability.

---

## V50 Results (March 2026)

### Polariton EM Wave [CONFIRMED — V50]

The Cosserat curl coupling produces a hybridized θ-δφ wave (polariton).
Dispersion relation proven in Lean (zero sorrys):
- v_phase = c × √(1 − η²/m²) = 0.9428c (predicted), 0.9377c (measured)
- Mixing ratio: δφ/δθ = −ηk/m² (~18% φ admixture at λ=8)
- Propagation mechanism: θ → curl(θ) → δφ → curl(δφ) → θ (E↔B cycle)
- True eigenmode is a plane wave (infinite transverse extent)
- Background A_bg must be 0 for clean propagation tests
- Details: `v50/EM_WAVE_RESULTS.md`, `lean/Polariton.lean`

### Two-Proton Bound State [CONFIRMED — V50]

V50/C4 produces stable two-body bound state with distinct cores (V44
allowed merger into a blob). Shell structure confirmed: core → hardened
shell → exterior. Details: `v50/RESULTS.md`

### Single-Pass Force Expansion [VERIFIED — V50]

The two-pass Cosserat+hardening force computation can be algebraically
expanded to a single GPU pass, eliminating 6×N³ intermediate arrays.
Verified in Maxima (10⁻¹⁸ precision) and Lean (zero sorrys).
Details: `v50/em_wave/single_pass.mac`, `lean/SinglePass.lean`
**Status**: Not yet implemented in the CUDA kernel.

---

## High Priority

### F17: Nuclear Binding Energy — ³He and ⁴He [V44-ERA]

²H confirmed (V42, V51). Next: ³He (2 UUD + 1 UDD) and ⁴He (2 UUD + 2 UDD).
Need N=512+ grids. Requires pre-converged UDD template (currently missing
for C4 equations).

### F18: Neutron Stability Under C4 [UNTESTED]

UDD behavior will change with Cosserat constraint. The V44 result (UDD
survives T=1000) needs verification with C4.

### F19: Force Equilibration [V44-ERA, needs retest]

The V42 strong/EM ratio equilibration (259:1 → 1:1) was measured with V44
equations. C4's Cosserat constraint changes θ dynamics fundamentally.

### F_NEW: High-Energy Proton Collision [UNTESTED]

V51 tested v=0.1c (kinetic/rest ~0.5%). Test v=0.3c and v=0.5c to find:
- At what energy do protons break apart permanently?
- Is there a transition from elastic bounce → bound state → disruption?
- Do the fragments resemble known particles (pions, kaons)?

### F3: Lorentz Contraction Verification [UNTESTED]

Boost a proton at v=0.1c, 0.3c, 0.5c. Verify γ contraction and time
dilation. The C4 equations are Lorentz-invariant (Lagrangian-derived),
so boosted solutions must transform correctly.

### F4: Isotropic Background [UNTESTED]

Test braids/protons in a random-phase (non-z-aligned) background.
If the theory requires the z-preferred background, there's a hidden
preferred-frame problem.

---

## Medium Priority

### F22: Breathing Mode Spectrum [V44-ERA]

FFT of per-shell ρ(t) to reveal full mode spectrum. The V50 two-proton
data shows compound breathing (4.5-unit beat pattern).

### F24: Mass Defect Measurement [V44-ERA, redesign for C4]

Differential method: same seeds, different D. The reference equation set
is now C4. Need pre-converged templates and T=500+ equilibration.

### F_NEW: Shell Thickness vs Parameters

Sweep β (hardening) at {0.1, 0.3, 0.5, 1.0, 2.0} with α=0.1 fixed.
How does shell thickness depend on β? Does it affect force law or binding?

### F_NEW: Proton Formation Pathway Under C4

The analytical seed failure (V51) raises the question: can protons form
spontaneously under C4? Test:
1. Start with V44 equations (α=β=0), form proton
2. Slowly ramp α and β to production values over T=200
3. Does the proton survive the transition? At what ramp rate?

### F21: Group Velocity Subluminal Inside Breathing Oscillator [UNTESTED]

Verify |∂φ/∂t| > c at antinodes is phase velocity, not group velocity.
Launch a δφ perturbation at the proton edge, track arrival at far side.

### F9: Analytical Effective Potential

Derive V_eff(D) from linearized perturbation theory. Would predict the
force law without simulation.

---

## Lower Priority / Speculative

### F10: Spin and Helical Handedness

Enhanced by V50: Cosserat constraint makes θ a true geometric quantity.
Spin identification more natural.

### F11: Gravitational Waves

Accelerated braid radiation — does it have spin-2 tensor structure?

### F12: Dark Matter Profiles

Does the depletion profile match NFW/Burkert halos? Quantitative comparison.

### F13: Quantization

Classical theory → QFT. Triple-product potential well-defined quantum mechanically?

### F14/F15: Cosmological Constant / Background Origin

The Cosserat strain energy E_strain = α|M|² contributes to vacuum energy.

### Speculative: Vacuum Refractive Index and Dark Matter

The polariton speed v < c gives the vacuum a refractive index n = c/v ≈ 1.06.
Near massive objects, depletion lowers A_bg → n closer to 1. The spatial
variation of n produces lensing that mimics dark matter halos. See
previous FUTURE.md for full discussion of MOND-like behavior, Cherenkov
radiation limit, chromatic micro-lensing smoking gun, and heliospheric
implications.

---

## Resolved / Abandoned

### Confirmed (no further action)
- F2: Gravity mechanism — dynamic footprint, NOT energy minimization [V33/V43/V51]
- F5: Photon speed — v = 0.9377c, matches polariton prediction [V50]
- F20: Intermediate phase group in deuterium — null result, not a bond [V43]
- Composite particles: UUD stable T=500, phase confinement works [V41]
- Nuclear binding: deuterium UUD+UDD bound T=500 [V42]
- Two-proton bound state with C4 [V50/V51]

### Abandoned
- X1: S/B two-component split (artificial)
- X2: c(ρ) speed modification (not needed)
- X3: Binding-weighted gradient coupling (net repulsive)
- X4: SPH particle-based simulation (poor conservation)
- X5: Rotating FRW expansion (no braid formation)
- X6: Asymmetric drag hypothesis (refuted by V33 drag test)
- X7: Single-particle density-dependent κ collapse (theta prevents it)
- X8: Random parameter sweep (superseded by first-principles V41)
