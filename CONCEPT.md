# SCP — A 6-Field Nonlinear Wave Simulator

This is a numerical simulation of coupled nonlinear wave equations on a 3D grid.
It produces interesting localized structures (oscillons) with some suggestive
parallels to particle physics, but no true stable particles or conserved
topological charges have been established.

---

## 1. The Equations

Six real scalar fields on a 3D grid: three "position" fields φ_a and three
"angle" fields θ_a (a = 0, 1, 2). The Cosserat equations:

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η × curl(θ)_a     (1a)

    ∂²θ_a/∂t² = ∇²θ_a + η × curl(φ)_a                          (1b)

    V(P) = (μ/2) P² / (1 + κP²)                                  (2)

    P = φ₀ φ₁ φ₂    (triple product)                              (3)

Standard parameters: m² = 2.25, μ = -41.345, κ = 50, η = 0.5

The φ fields are massive (m² = 2.25) with a saturating triple-product
potential V(P). The θ fields are massless (m_θ² = 0) and couple to φ
through the curl terms. The system is Lorentz-invariant by construction.

### Background

The simulations typically start with a uniform oscillating background:

    φ_a(x) = A_bg × cos(k·z + 2πa/3) + (localized structures)

with A_bg ≈ 0.1. This background is not required by the equations — it is
an initial condition that provides a medium for interactions. Some experiments
use A_bg = 0 (vacuum).

---

## 2. What the Simulator Produces

### Oscillons (Breathing Localized Structures)

Certain initial conditions produce localized field concentrations that persist
for hundreds of time units under periodic boundary conditions. These are
oscillons — nonlinear breather solutions that periodically concentrate and
disperse.

**Lissajous seeds** (3-axis superposition with phase offsets Δ = {0, π/3, π/3})
produce the longest-lived oscillons. The initial theta field is set to
θ = -G × curl(φ) with gain G ≈ 2.5.

**Observed behavior**:
- The localized mass oscillates by 95%+ during each breathing cycle (period ~25-50 t)
- At the breathing minimum, E_pot approaches zero — the structure has essentially
  dispersed into low-amplitude waves
- Under periodic BC, these waves wrap around and reconcentrate, giving the appearance
  of indefinite survival
- Under absorbing BC, oscillons have finite lifetimes (radiation losses eventually
  drain the structure)
- The structure often fragments into 2-3 sub-clusters that exchange mass

**The |μ|/m² < 7 condition**: Through analysis (confirmed in Lean 4), we derived
that sub-threshold breathing requires |μ|/m² < ~7. This is better understood as
a slow-dispersal condition — below this threshold, the nonlinear potential produces
breathing amplitudes that stay small enough to avoid rapid radiation loss. It does
not guarantee true stability.

### Theta Field Growth

The η × curl coupling transfers energy from φ to θ over time. Starting from
θ = 0, theta kinetic energy grows to match phi kinetic energy over T ≈ 600.
This occurs for any localized phi oscillation and is a generic consequence of
the linear coupling, not evidence of a topological mechanism.

### H_cross (Helicity Integral)

We measure H_cross = ∫ θ · curl(φ) dV as a chirality indicator. Observations:

- H_cross is consistently negative for seeds initialized with negative theta gain
- H_cross is NOT conserved — it tracks the localized mass (drops to near-zero
  when the oscillon disperses, recovers when it reconcentrates)
- Attempts to construct positive-H_cross (opposite chirality) oscillons have
  all failed — the dynamics drives any configuration toward negative H_cross
- This single-sign preference suggests H_cross reflects the equation structure,
  not a conserved topological charge

---

## 3. Gradient Response

When an oscillon is placed in a density gradient (non-uniform A_bg), it drifts.

**Observed** (V33, V43):
- In a 3-field-only simulation (no θ), localized structures drift toward low density
- The drift is proportional to the gradient strength (confirmed at multiple strengths)
- The mechanism appears to be an asymmetric Yukawa profile: the perturbation
  extends further on the low-density side

**Limitations**:
- The drifting objects are oscillons, not stable particles
- Whether this mechanism produces 1/r² force law at long range is unclear
  (measured exponent ~1.8, but contaminated by periodic BC and finite grid)
- The "gravitational" drift is small (~1-2 code units over T=200)
- No two-body orbit or sustained gravitational interaction has been demonstrated

---

## 4. Charge-Dependent Force

In the 6-field Cosserat system, the curl coupling produces winding-dependent
interactions between oscillons:

**Observed** (V34):
- Same-winding oscillons attract 27% more than the 3-field baseline
- Opposite-winding oscillons attract 57% less
- The θ field around an oscillon has a DC component that reverses with winding

This is consistent with an electromagnetic analog where winding plays the role
of charge. However:
- The interaction is wave-mediated (oscillating, period ~4t), not static 1/r
- The 0.2% DC component is small relative to the oscillation amplitude
- The oscillons involved are not stable particles

---

## 5. Composite Structures

### 3-Braid Composites

Three oscillons with carrier phase offsets Δ = {0, 2π/3, 4π/3} arranged along
orthogonal axes can be initialized as a composite. The phase structure creates
destructive interference at the triple overlap (P → 0), preventing immediate
merger. This is analogous to color confinement.

**Observed** (V41):
- UUD composites (two same-chirality + one opposite) survive T=500 under
  periodic BC with anti-phase breathing
- UDD composites are less stable (no θ synchronization, phase convergence)
- The "proton" (UUD) is more stable than the "neutron" (UDD)

**Caveat**: These composites exist under periodic BC. Their long-term fate
under absorbing BC has not been established. The individual components are
oscillons with the same dispersal-and-reconcentration dynamics.

### Deuterium Analog

A UUD + UDD pair at N=512 was simulated for T=500 (V42). The two composites
maintained separation and showed inter-baryon attraction with force equilibration
(strong:EM ratio → 1:1). The system was compacting at the end of the run.

This is an interesting dynamical result, but again under periodic BC.

---

## 6. What Has NOT Been Established

1. **True particle stability**: No simulation has produced a localized structure
   that maintains a persistent core under absorbing BC for arbitrarily long times.
   All "particles" are oscillons that breathe through near-zero mass.

2. **Conserved topological charge**: H_cross tracks mass, not topology. No
   quantized or conserved charge has been identified in the field configurations.

3. **Opposite chirality**: Only one sign of H_cross is dynamically preferred.
   A true chiral charge would be stable in both signs.

4. **Newtonian gravity**: The gradient drift is suggestive but the force law
   exponent, coupling constant, and long-range behavior are not established.
   No gravitational orbit or sustained two-body attraction has been shown.

5. **Quantitative connection to real physics**: The parallels (proton > neutron
   stability, charge-dependent force, confinement from phase cancellation) are
   qualitative. No dimensionless ratio from the simulation matches a measured
   physical constant.

---

## 7. Technical Infrastructure

### Simulation Kernel
- **CPU**: `sfa/sim/scp_sim.c` — 6-field Cosserat, config-driven, OpenMP
- **GPU**: `sfa/sim/scp_sim.cu` — CUDA port, async I/O, V100 ~50× faster than CPU
- Config: `sfa/sim/scp_config.h`, Init: `sfa/sim/scp_init.h`
- Build/run via `scp-runner` MCP server

### Output Format
- SFA format (`sfa/format/sfa.h`) — chunk-based binary with zstd compression
- Supports voxel frames (FRMD) and polynomial vector frames (FRVD)
- Viewer: `sfa/volview/main.go` (OpenGL volume renderer)

### Seed Generators
- `sfa/seed/gen_composite.c` — stamp templates into grids
- `v52/seeds/lissajous/gen_lissajous_seed.c` — 3-axis superposition
- Various others in `sfa/seed/`

### Analysis
- `sfa/analysis/sfa_particle_track.c` — flood-fill cluster detection, per-particle
  mass/centroid/H_cross/E_pot tracking across frames
- `sfa/analysis/sfa_extract.c` — field extraction and profiling

### Historical Code
- `v33/src/v33.c` — original 3-field gravity tests
- `v34/torsion_coupling/src/v33_cosserat.c` — first 6-field implementation

---

## 8. Parameter Space Explored

| Parameter | Range Tested | Notes |
|-----------|-------------|-------|
| m² | 1.0 - 9.0 | Higher m² → shorter-lived oscillons but deeper potential |
| μ | -20 to -60 | Controls nonlinear binding; |μ|/m² < 7 for slow dispersal |
| κ | 10 - 100 | Saturation parameter |
| η | 0 - 1.0 | φ-θ curl coupling strength |
| A_bg | 0 - 0.15 | Background amplitude |
| N | 64 - 512 | Grid resolution (results are resolution-independent above N=96) |
| BC | periodic, absorbing | Periodic preserves oscillons; absorbing reveals finite lifetime |

---

## 9. Version History

- **V28**: CMA-ES parameter search, bimodal seed discovery
- **V29**: Validation campaign (13 findings on oscillons and their properties)
- **V30-V32**: Failed approaches (c(ρ), FRW expansion, SPH)
- **V33**: 3-field gradient drift measurement (F ∝ ∇ρ, 1/D^1.8)
- **V34**: 6-field Cosserat, charge-dependent force, θ characterization
- **V35**: Hydrogen-like bound state (θ shell + Schrödinger wedge)
- **V36**: CUDA GPU port
- **V39**: BLV effective metric analysis (dispersive, not geometric)
- **V40**: Evolutionary composite search, stability signatures
- **V41**: Phase confinement, UUD/UDD composites
- **V42**: Deuterium analog (two-baryon binding under periodic BC)
- **V43**: Proton gradient test (gravitational drift confirmed for composites)
- **V44-V50**: Equation refinements (Cosserat strain, curl-hardening)
- **V52**: Lissajous 3-axis seed construction, chirality measurement
- **V53**: Stability condition |μ|/m² < 7 derived; N=192 grid tests confirm
  oscillons are breathing structures, not stable particles
