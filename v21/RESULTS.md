# V21 Results: Three-Body Oscillon with Mass Gap

## Summary

**First positive result in the soliton search program.** Three massive scalar fields coupled by a saturating triple-product potential form a long-lived, localized, oscillating bound state (oscillon). The oscillation frequency lies below the mass gap, preventing radiative decay. The bound state is genuinely three-body: removing any one field destroys the binding.

**Works in both 1D and 3D.** The 1D oscillon (μ=-10, κ=10) lives >10,000 time units. The 3D oscillon requires stronger coupling (μ=-20, κ=20) to compensate for spherical spreading, and has been **confirmed stable for 3000 time units** (full production run complete) with f_core > 0.98 at N=200 resolution. Energy loss rate decelerates to ~0.0004/t.u., giving extrapolated lifetime >200,000 t.u. **Gravity backreaction (β>0) tested: oscillon survives at α,β ≤ 0.1; stability boundary near α≈0.3.**

---

## Model

Lagrangian density (three scalar fields φ_a, a = 1,2,3):

    L = Σ_a [ ½(∂_t φ_a)² - ½(∂_x φ_a)² - ½m²φ_a² ] - V(P)

    V = (μ/2) P² / (1 + κ P²),   P = φ₁ φ₂ φ₃

- Mass term m²φ²/2: ensures φ=0 is the true vacuum, gives exponential far-field decay
- Triple product P: nonzero only when all three fields overlap
- Saturating potential: caps binding energy, prevents collapse
- Dispersion: ω² = k² + m², mass gap ω_min = m at k=0
- Oscillon mechanism: nonlinear frequency shift pushes ω below m → radiation forbidden

Parameters for main results: μ = -10, κ = 10, m = 1.0

**Key condition**: φ=0 must be the true vacuum (not a false vacuum). Requires:
m² > |μ| × (2/κ)^{2/3} / 9. For μ=-10, κ=10: need m > 0.62.

---

## Test Results

**Code**: `src/triad1d.c` | Grid: Nx=2000, xmax=60, dx=0.06, dt=0.031

### Test Comparison (μ=-10, κ=10, m=1.0, A=0.8, σ=3.0)

| Test | Description | E(0) | E(2000) | E retained | f_core(2000) | Peak | ω_peak |
|------|-------------|------|---------|-----------|-------------|------|--------|
| **1** | **Symmetric triad (A,A,A)** | **3.84** | **3.46** | **90%** | **0.999** | **0.72** | **0.82** |
| 2 | Asymmetric triad (A,0.8A,0.6A) | 2.93 | 1.83 | 62% | 0.996 | 0.48 | 0.90 |
| 3 | Two-field control (A,A,0) | 3.59 | 1.74* | 48% | ~0.15 | 0.19 | — |
| 4 | Single-field control (A,0,0) | 1.80 | 0.77* | 43% | ~0.15 | 0.12 | — |

(*) at t≈900, extrapolated

**The triad retains 90% of energy in the core for 2000+ time units, while controls lose 50-60%.**

### Long-Duration Test (t=10,000)

| Time range | E | f_core | dE/dt |
|-----------|---|--------|-------|
| t = 0 | 3.84 | 1.00 | — (initial shedding) |
| t = 1000 | 3.50 | 0.999 | -3.4e-4 |
| t = 5000 | 1.80 | 0.999 | -3.6e-5 |
| t = 10000 | 1.77 | 0.999 | -1.6e-5 |

After the initial shedding phase (t < 3000), the energy loss rate drops to ~10⁻⁵ per time unit.
Extrapolated lifetime at this rate: > 100,000 time units (>15,000 oscillation periods).

### Oscillation Frequency

DFT of φ₁(x=0, t) for the second half of the simulation:
- **Peak ω = 0.82** (symmetric) / **0.90** (asymmetric)
- **Mass gap m = 1.0**
- ω < m → subgap oscillon (radiation into dispersive modes forbidden)

### Mode Equalization (Test 2)

Starting from asymmetric amplitudes (0.80, 0.64, 0.48), the three fields
equalize within ~300 time units to nearly identical amplitudes (~0.58).
This is the three-body resonance sharing mechanism: the coupling P = φ₁φ₂φ₃
creates a restoring force that equalizes the modes.

---

## Why This Works (And Why v19/v20 Failed)

### v19 failure (S³ + Skyrme, no mass gap)
- Massless dispersion ω = k → all frequencies can radiate
- No trapping mechanism below a gap
- Energy escapes at speed c

### v20 failure (triple product, no mass)
- Two failure modes: collapse (weak κ) or dispersal (strong κ)
- No exponential localization (fields spread as 1/r^n)
- Derrick theorem without mass: no static solution

### v21 success: three ingredients
1. **Mass gap (m > 0)**: creates ω_min = m. Oscillon at ω < m cannot radiate.
2. **Triple-product coupling (μ < 0)**: provides attractive binding between three fields.
   Nonlinear frequency shift: ω² ≈ m² - C|μ|A₁²A₂²A₃² pushes ω below gap.
3. **Saturation (κ > 0)**: prevents collapse. Caps V at |μ|/(2κ).

The mass gap is the crucial missing ingredient. Without it, there's no frequency
below which radiation is forbidden. With it, the oscillon is protected from
decay by energy-momentum conservation (no dispersive modes at ω < m).

---

## Validation: Standard Oscillon

To confirm the code correctly captures oscillons, we tested the well-known
φ⁴-φ⁶ model (single scalar):

    V = m²φ²/2 - gφ⁴/4 + hφ⁶/6

**Code**: `src/oscillon_std.c` | Parameters: m=1, g=2, h=1, A=0.8, σ=3

Result: ω = 0.708 < m = 1.0, f_core > 0.96, energy loss 17% over 1000 time units.
This confirms the simulation framework correctly detects and sustains oscillons.

---

## Implications

### For the SCP program
This is the first demonstration that three-body resonance binding works in a
classical field theory. The key was adding a mass gap, which:
- Provides exponential localization (φ ~ e^{-mr} at large r)
- Creates a frequency band [0, m) that cannot radiate
- Enables Derrick-compatible equilibrium with the triple coupling

---

## 3D Results

### Parameter Discovery

The original 1D parameters (μ=-10, κ=10, m=1.0) fail in 3D: the oscillon dies at t≈130.
A 1D parameter scan (1280 configurations) + 3D validation identified the critical fix:
**stronger coupling with matched saturation (μ=-20, κ=20, m=1.0)** pushes ω deeper below
the mass gap, providing enough margin for 3D survival.

| Parameters | ω/m (1D) | ω/m (3D) | 3D survival |
|-----------|----------|----------|-------------|
| μ=-10, κ=10, m=1.0 (original) | 0.82 | ~1.0+ | **Dies at t≈130** |
| μ=-20, κ=20, m=1.0 (optimized) | 0.76 | 0.948 | **Alive at t=300+** |
| μ=-30, κ=20, m=1.0 | 0.74 | ~0.94 | **Alive at t=300+** |

**Key insight**: In 3D, the nonlinear frequency shift is weaker (energy spreads over 4πr²),
so the 1D gap margin (~18%) is insufficient. The optimized parameters provide a 24% gap
margin in 1D, which yields a 5% gap margin in 3D — just enough.

### 3D Production Run (μ=-20, κ=20, m=1.0, A=0.8, σ=3.0)

**Code**: `src/triad3d.c` | Grid: N=200, L=30, dx=0.30, dt=0.075

| t | E_total | f_core | E_pot | Peak amp | Notes |
|---|---------|--------|-------|----------|-------|
| 0 | 134.3 | 0.999 | -33.9 | 0.797 | Initial Gaussian |
| 30 | 131.0 | 0.893 | -19.8 | 0.666 | Shedding phase |
| 60 | 112.7 | 0.947 | -5.5 | 0.504 | Core recovering |
| 90 | 103.2 | 0.968 | -1.4 | 0.396 | Stabilizing |
| 120 | 97.9 | 0.975 | -0.5 | 0.334 | Core tightening |
| 150 | 95.3 | 0.986 | -0.2 | 0.296 | Near steady-state |
| 210 | 91.7 | 0.983 | -0.7 | 0.354 | Amplitude growing |
| 270 | 90.0 | 0.987 | -2.8 | 0.453 | Self-focusing |
| 330 | 89.0 | 0.990 | -9.6 | 0.583 | Binding deepening |
| 390 | 88.1 | 0.982 | -16.1 | 0.786 | Strong pulsation |
| 450 | 87.7 | 0.989 | -19.4 | 0.921 | Peak of breathing cycle |
| 480 | 87.3 | 0.982 | -17.1 | 0.861 | Breathing continues |
| 600 | 86.7 | 0.986 | -0.8 | 0.364 | Trough of breathing |
| 720 | 86.3 | 0.991 | -0.2 | 0.288 | Next trough |
| 840 | 86.0 | 0.990 | -17.2 | 0.842 | Second peak |
| 930 | 85.8 | 0.990 | -11.7 | 0.664 | Steady breathing |
| 1050 | 85.5 | 0.992 | -0.0 | 0.035 | Deep trough |
| **1260** | **85.3** | **0.989** | **-16.6** | **0.850** | **Confirmed >1000 t.u.** |
| 1500 | 84.9* | 0.988* | — | — | Steady breathing |
| 1920 | 84.7 | 0.985 | -17.0 | 0.882 | Steady breathing |
| 2220 | 84.6 | 0.986 | -17.2 | 0.890 | Steady breathing |
| 2700 | 84.3 | 0.990 | -0.0 | 0.009 | Deep trough |
| **3000** | **84.2** | **0.989** | **-0.0** | **0.038** | **FINAL — 3000 t.u. confirmed** |

**Final ω = 0.954 < m = 1.0 (4.6% gap margin). Production run complete.**

**Self-focusing transition**: After shedding (t<150), the core amplitude grows
from 0.30 to 0.92 over 300 time units, then enters a stable breathing mode.
The nonlinear frequency shift deepens as amplitude increases, concentrating energy.

**Energy loss rate decelerates** (NOT T⁴ thermal radiation):

| Period | dE/dt avg | Extrapolated lifetime |
|--------|-----------|----------------------|
| t=300-600 | -0.009 | ~10,000 |
| t=600-900 | -0.002 | ~40,000 |
| t=900-1260 | -0.001 | **>80,000** |
| t=1260-1920 | -0.0008 | >100,000 |
| t=1920-2220 | -0.0006 | >140,000 |
| t=2220-3000 | -0.0004 | **>200,000** |

The loss rate drops roughly exponentially, consistent with residual non-equilibrium
shedding, NOT blackbody-like T⁴ emission. No runaway instability observed.

**3D oscillation frequency**: ω = 0.942–0.948 < m = 1.0 (DFT, confirmed in multiple runs).

### 3D Control Comparison (N=100, μ=-20, κ=20, t=300)

| Metric | Triad (3 fields) | Control (1 field) |
|--------|-----------------|-------------------|
| f_core | **0.996** | 0.65 (oscillating) |
| E retained | 51% | **16%** |
| Peak amplitude | 0.596 | 0.108 |
| E_pot (coupling) | -16.2 | 0.0 |
| ω (DFT) | **0.948 < m** | 1.032 > m |
| Oscillon? | **YES** | NO |

The single-field control has ω > m (above the mass gap), confirming that the
three-body coupling is essential for pushing the frequency below the radiation threshold.

### Why Original Parameters Fail in 3D

In 3D, the field amplitude scales as A_3D ~ A_1D / r (spherical spreading).
The frequency shift Δω² ∝ A₁²A₂²A₃² ∝ A⁶ drops as r⁶.
For the oscillon to survive, need ω(A_core) < m after the shedding phase.

- Original (μ=-10): 1D ω=0.82, margin = 18%. After 3D shedding, ω rises above m → dead
- Optimized (μ=-20): 1D ω=0.76, margin = 24%. After 3D shedding, ω=0.948 < m → alive

### 3D Adjoint Optimization (Forward-Mode Tangent)

**Code**: `src/adjoint3d.c` | N=80, L=20, tfinal=150, 4 parameters (A, σ, μ, κ)

Propagates tangent vectors dφ/dθ through the Verlet integrator simultaneously with
the forward simulation. Shares Jacobian dF_a/dφ_b computation across all 4 tangent
vectors. 15 iterations at lr=0.005:

| Parameter | Start | End | dfc/dθ |
|-----------|-------|-----|--------|
| A | 0.800 | 0.804 | +0.050 |
| σ | 3.000 | 3.003 | +0.044 |
| μ | -20.0 | -20.0 | -0.009 |
| κ | 20.0 | 20.0 | -0.006 |

fc improved from 0.983 to 0.983 (negligible). **Conclusion**: (A=0.8, σ=3.0, μ=-20, κ=20)
is near a local optimum for 3D fc at t=150. The gradients confirm that slightly larger A
and σ help, while μ and κ are already well-matched.

---

## Gravity Coupling (Exploratory)

### Static Gravitational Potential (Poisson Solver)

**Code**: `src/poisson_phi.c` | Input: production run profile at t=1000

Solves Φ(r) = -α[I1(r)/r + I2(r)] using the spherical Green's function.

- **Total source**: Q = ∫ρ·4πr²dr = 85.04 (matches E_total = 85.5 from production run)
- **Φ(0)** = -2.82 (at α=1)
- **Far-field**: r·Φ converges to -6.768 for r > 15, confirming perfect 1/r tail
- Profile extends to r~10, exponential decay beyond

### Dynamic Passive Gravity Test

**Code**: `src/triad3d_grav.c` | α=0.01, β=0 (no backreaction), N=100, L=20, t=200

Three scalar fields + massless scalar Φ evolved simultaneously. Φ initialized from
Poisson solution of initial ρ (avoids transient ringing).

- Oscillon survives with fc=0.991, ω=0.942 < m (unchanged from no-gravity case)
- **Monopole radiation identified**: Q_monopole = ∫ρd³x oscillates between 1 and 130
  due to the breathing mode. This modulation radiates scalar monopole waves in Φ.
- Box too small (L=20) to see clean 1/r in dynamic test; static Poisson confirms it.

**Key issue**: Scalar gravity inherently allows monopole radiation (∂²Q/∂t² ≠ 0).
Real gravity is tensorial and forbids monopoles. This is a fundamental limitation
of the scalar proxy.

### Backreaction Test (α=0.01, β=0.01, N=100, t=300)

m_eff² = m₀² - β·Φ. Since Φ < 0 near core, m_eff > m₀ → deeper mass gap.

| t | E (gravity) | E (pure) | fc (gravity) | fc (pure) |
|---|-------------|----------|--------------|-----------|
| 0 | 134.3 | 134.3 | 0.999 | 0.999 |
| 100 | 98.1 | 97.9* | 0.984 | 0.975 |
| 200 | 89.5 | 91.7* | 0.991 | 0.983 |
| 300 | 86.6 | 89.0 | 0.989 | 0.990 |

(*pure values from N=200 run, interpolated to N=100 timing)

**Result: Oscillon survives with gravity backreaction.** Energy loss is ~3% faster
(monopole radiation in Φ carries extra energy), but fc is comparable or higher.
The gravity-induced mass enhancement (m_eff > m₀) appears to help concentrate
the core during the shedding phase.

ω = 0.948 < m = 1.0 (unchanged). Φ_c oscillates between -0.025 and +0.01 with
the breathing mode. Monopole moment Q swings from 0.5 to 124.

### Gravity Coupling Strength Scan (N=100, t=300)

| α, β | E(300) | fc(300) | ω | Φ_c range | Alive? |
|------|--------|---------|---|-----------|--------|
| 0, 0 (pure) | 89.0 | 0.990 | 0.948 | — | YES |
| 0.01, 0.01 | 86.6 | 0.989 | 0.948 | [-0.025, +0.01] | **YES** |
| 0.1, 0.1 | 87.5 | 0.985 | 0.954 | [-0.20, +0.05] | **YES** |
| 0.2, 0.2 | 90.4 | 0.986 | 0.966 | [-0.43, +0.19] | **YES** |
| 0.3, 0.3 | 96.8 | 0.970 | 0.984 | [-0.65, +0.32] | **YES** (marginal) |
| 0.5, 0.5 | — | 0.18→0.70 | — | [-1.3, +0.7] | **DEAD (t≈100)** |
| 1.0, 1.0 | — | — | — | [-12, —] | **RUNAWAY (t≈50)** |

**Stability boundary**: Between α=0.3 and α=0.5. The mechanism is clear:
increasing α,β shifts ω upward toward the mass gap m=1.0. At α=0.3, ω=0.984
(only 1.6% gap margin). At α=0.5, ω crosses above m → radiation is no longer
forbidden → oscillon dies. At α=1.0, self-sourcing runaway (positive feedback).

**Key insight**: The backreaction m_eff² = m₀² - βΦ deepens the mass gap in the
core (Φ<0), but the increased effective mass *raises* the oscillation frequency.
Too much mass enhancement pushes ω above the gap, destroying the sub-gap
trapping mechanism.

At α≤0.2: the coupled oscillon-gravity system enters a stable breathing mode.
Energy oscillates between φ_a and Φ fields. The oscillon is robust to
gravitational backreaction at moderate coupling.

### Full-ρ Source (Option A: monopole suppression)

**Fix**: Include kinetic energy ½v² in ρ_source. Since total energy is conserved,
∫ρ_full d³x = E_total = const → monopole moment is constant → no monopole radiation.

| α, β | Source | Q range | Q variation | E(300) | fc(300) |
|------|--------|---------|-------------|--------|---------|
| 0.01, 0.01 | old (no kinetic) | 1–130 | 130× | 86.6 | 0.989 |
| 0.01, 0.01 | **full ρ** | **83–91** | **1.1×** | **86.6** | **0.989** |
| 0.1, 0.1 | old (no kinetic) | ~1–130 | ~100× | 87.5 | 0.985 |
| 0.1, 0.1 | **full ρ** | **85–93** | **1.1×** | **89.1** | **0.988** |

**Monopole oscillation suppressed by ~100×.** The residual 10% Q variation tracks
slow energy loss through the absorbing boundary. The full-ρ source also retains
slightly more energy (89.1 vs 87.5 at α=0.1), confirming monopole radiation was
draining the system.

At α=0.1 with full ρ: ω=0.954, Φ_c ~ -0.25, Q ~ 85–93 (nearly constant).
The 1/r Newtonian tail is preserved (static monopole = E_total).

### Open questions for v22+
1. ~~Backreaction test~~ — **DONE**: oscillon survives at α,β ≤ 0.3; dies at 0.5.
2. **Monopole radiation**: Scalar Φ emits monopole waves from breathing mode.
   Tensor gravity would suppress this (no monopole GW). Key limitation of scalar proxy.
3. **Mass origin**: What sets m physically? Could come from vacuum structure.
4. **Lifetime**: dE/dt ≈ 0.0004 at t=3000. Is it truly asymptotic to zero?
   Higher harmonics: 2ω = 1.908 > m → second harmonic can radiate (slowly).
5. **Self-consistent bootstrap**: Can m₀→0 with mass generated entirely by
   backreaction (m_eff² = -βΦ)? Risk: vacuum is gapless → all ω can radiate.

---

## Comparison Across Versions

| Version | Mechanism | Mass gap? | Localized? | Stable? | Three-body? |
|---------|-----------|----------|-----------|---------|-------------|
| v19 | S³ + Skyrme | No | No | No | Yes (but no binding) |
| v20 | Triple product (no mass) | No | No | No (collapse/disperse) | Yes |
| **v21 (1D)** | **Triple product + mass** | **Yes** | **Yes** | **Yes (>10⁴ t.u.)** | **Yes** |
| **v21 (3D)** | **Triple product + mass** | **Yes** | **Yes** | **Yes (3000 t.u.)** | **Yes** |
| Std oscillon | φ⁴ single field | Yes | Yes | Yes | No (single field) |

---

## Files

| File | Description |
|------|-------------|
| `src/triad1d.c` | 1D three-field oscillon simulator |
| `src/triad3d.c` | 3D three-field oscillon simulator (OpenMP) |
| `src/scan1d.c` | 1D parameter scan (1280 configurations) |
| `src/adjoint1d.c` | 1D adjoint optimizer (forward-mode tangent) |
| `src/adjoint3d.c` | 3D adjoint optimizer (forward-mode tangent, OpenMP) |
| `src/oscillon_std.c` | Standard single-field oscillon (validation) |
| `src/biharm1d.c` | Biharmonic + triple product (exploratory) |
| `src/static1d.c` | Static gradient flow solver (negative: no static solution) |
| `src/oscillon1d.c` | Single-field massive + triple product (exploratory) |
| `src/triad3d_grav.c` | 3D simulator with massless scalar gravity (OpenMP) |
| `src/poisson_phi.c` | Spherical Poisson solver for oscillon potential |
| `PLAN_v22.md` | Gravity coupling roadmap (Phases 1-4) |
| `GRM.md` | Gravity mechanism discussion |
| `data/triad1d_test{N}_ts.tsv` | 1D time series for test N |
| `data/triad1d_test{N}_spectrum.tsv` | 1D DFT power spectrum |
| `data/triad3d_test{N}_ts.tsv` | 3D time series for test N |
| `data/triad3d_test{N}_profile_t{T}.tsv` | 3D radial profile at time T |
| `data/scan1d_results.txt` | 1D parameter scan results |
| `data/oscillon_std_ts.tsv` | Standard oscillon time series |
| `data/oscillon_std_spectrum.tsv` | Standard oscillon spectrum |
