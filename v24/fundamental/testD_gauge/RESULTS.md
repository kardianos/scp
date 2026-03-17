# Test D Results: Local Breathing Gauge (Ω Modulator)

## Summary

A massless scalar field Ω(x,t), sourced by the oscillon's energy density and
coupled via mass modulation m²_eff = m²(1 + g_Ω·Ω), **destabilizes the oscillon
over long timescales**. The Ω field develops an **oscillating** (not static)
profile that tracks the breathing mode, with a slowly growing DC component that
eventually pushes the effective mass out of the sub-gap regime.

**Key result**: The Ω-mediated mass modulation is DESTRUCTIVE, not stabilizing.
It acts as a parametric perturbation that erodes the oscillon over ~8000 time units
(vs >10000 for control). No static 1/r analog profile develops — the 1D massless
Green's function produces a linearly growing background that is absorbed at
boundaries, leaving an oscillating residual.

---

## Parameters

    μ = -20, κ = 20, m = 1.0, A = 0.8, σ = 3.0
    g_Ω = 0.1, g_source = 0.01
    Nx = 4000, xmax = 100, dx = 0.05, dt = 0.02
    tfinal = 10000

---

## Test 1: Single Oscillon + Ω Coupling

| t     | E_total | f_core | Ω(0)   | Ω(10)  | Ω(40)  | peak_φ |
|-------|---------|--------|--------|--------|--------|--------|
| 0     | 3.37    | 1.000  | 0.000  | 0.000  | 0.000  | 0.800  |
| 1000  | 2.77    | 0.974  | +0.73  | +0.62  | +0.20  | 0.540  |
| 2000  | 1.44    | 0.981  | +1.13  | +1.07  | +0.80  | 0.314  |
| 3000  | 1.42    | 0.981  | +1.28  | +1.21  | +0.75  | 0.305  |
| 5000  | 1.26    | 0.994  | -0.31  | -0.35  | -0.29  | 0.189  |
| 7000  | 1.34    | 0.982  | +0.77  | +0.72  | +0.50  | 0.169  |
| 8000  | 1.37    | 0.964  | +1.31  | +1.24  | +0.76  | 0.195  |
| 9000  | 1.08    | **0.087** | +0.51  | +0.50  | +0.41  | 0.078  |
| 10000 | 0.69    | **0.179** | -0.31  | -0.30  | -0.25  | 0.093  |

**Spectrum (second half DFT)**: ω = 0.888 < m = 1.0

The oscillon survives until t ≈ 8000, then rapidly delocalizes (f_core drops
from 0.96 to 0.09). By t = 10000, E = 0.69 with peak φ = 0.09 — effectively dead.

### Ω Profile: Oscillating, NOT Static

The Ω field oscillates coherently with the breathing mode period (~7.1 time units
= 2π/0.888). At any given time:
- Ω has a mild spatial gradient: Ω(0) > Ω(10) > Ω(40) (peaked at core)
- The gradient is SMALL: at t=3000, Ω ranges from 1.28 (core) to 0.75 (x=40)
- The field reverses sign: Ω swings between roughly -0.4 and +1.3
- There is NO static component building up — the absorbing boundaries drain
  outgoing Ω waves, preventing secular growth

This is **not** the 1/|x| profile one would expect from a static source.
The reason: the oscillon breathing modulates the source ρ(x,t) at frequency
2ω ≈ 1.78, which is ABOVE the Ω mass gap (m_Ω = 0), so the source radiates
freely. The dominant Ω signal is a standing wave, not a static potential.

### Frequency Shift

| Test | ω (DFT) | Δω from control |
|------|---------|-----------------|
| Test 1 (g_Ω = 0.1) | 0.888 | +0.012 (+1.4%) |
| Test 3 (g_Ω = 0.0) | 0.876 | — (reference) |

The Ω coupling shifts ω **upward** by ~1.4%. This is because Ω is on average
slightly positive (positive source), so m_eff > m on average, raising ω.
The shift is small but measurable. This reduces the sub-gap margin from 12.4%
to 11.2%, slightly weakening the oscillon protection.

---

## Test 3: Control (g_Ω = 0)

| t     | E_total | f_core | peak_φ |
|-------|---------|--------|--------|
| 0     | 3.37    | 1.000  | 0.800  |
| 1000  | 1.43    | 0.962  | 0.162  |
| 2000  | 1.51    | 0.976  | 0.507  |
| 5000  | 1.24    | 0.992  | 0.382  |
| 8000  | 1.40    | 0.974  | 0.345  |
| 10000 | 1.26    | **0.992** | 0.492  |

**Spectrum**: ω = 0.876 < m = 1.0

The control oscillon remains healthy at t = 10000 with f_core = 0.99.
The Ω field still evolves (sourced by ρ) but does NOT feed back into the
matter fields (g_Ω = 0). Energy: 1.26 vs 0.69 for test 1.

**Conclusion**: The Ω coupling causes a factor ~1.8× faster energy loss
and eventual oscillon death.

---

## Test 2: Two Oscillons (Force Test)

Separation D = 30, both at A = 0.8.

| t    | sep    | phi_L  | phi_R  | Ω_mid | E_total |
|------|--------|--------|--------|-------|---------|
| 0    | 30.0   | 0.800  | 0.800  | 0.000 | 6.73    |
| 200  | 68.3   | -0.012 | -0.012 | 2.57  | 7.03    |
| 400  | 128.5  | 0.012  | 0.012  | -0.24 | 2.81    |
| 800  | 119.2  | -0.003 | -0.003 | 2.01  | 2.74    |
| 1200 | 51.0   | -0.027 | -0.026 | 1.04  | 1.48    |
| 1600 | 84.5   | 0.033  | 0.033  | -1.67 | 0.66    |
| 2000 | 70.0   | -0.045 | -0.045 | -0.01 | 0.51    |

**Spectrum**: ω = 1.074 > m = 1.0 (ABOVE mass gap — oscillons dead)

Both oscillons effectively die within ~200 time units (phi drops from 0.8 to 0.01).
The Ω field grows to amplitude ~4 during the initial shedding phase, which
massively perturbs m_eff. The "separation" oscillates because the energy centroids
track the dispersed wave packets, not coherent solitons.

**No coherent force was observed** — the oscillons die too quickly for
a force measurement. The Ω field is dominated by radiation from the
initial transient, not by a static mediated potential.

---

## Physics Discussion

### Why Ω is Oscillating, Not Static

In the proposal's picture, a static source ρ₀(x) would produce a static
Ω(x) with □Ω = ∇²Ω = g_s·ρ₀, giving a linear-in-|x| profile (1D
Coulomb analog). However:

1. The oscillon breathes at frequency ω ≈ 0.88, so ρ(x,t) oscillates
   at frequency 2ω ≈ 1.76 (energy density = φ²-terms).
2. Since Ω is massless, all frequencies propagate freely. The oscillating
   source radiates Ω waves at 2ω.
3. The absorbing boundaries drain the outgoing radiation, preventing
   secular buildup.
4. The net result: Ω is a standing wave modulated by the breathing,
   with NO static monopole component.

A static Ω profile would require either:
- A truly static source (not an oscillon — would need a true soliton)
- A massive Ω field with m_Ω > 2ω (screens oscillating source, passes DC)

### Destabilization Mechanism

The Ω oscillation modulates m_eff at frequency 2ω. This creates parametric
resonance: the mass modulation at 2ω drives the fundamental mode at ω.
Over thousands of oscillation periods, this slowly pumps energy out of the
coherent breathing mode into dispersive radiation.

The effect is stronger for larger g_Ω: with g_Ω = 0.1 and Ω ≈ 1, the
mass modulation is δm²/m² ≈ 0.1 (10%), which is significant.

### Implications for Emergent Gravity

1. **No static gravitational potential**: The massless scalar Ω does not
   develop a 1/r (or linear-in-1D) static profile around the oscillon.
   The breathing dynamics produce radiation, not a potential well.

2. **Parametric instability**: Coupling the breathing frequency to a
   massless mediator creates a parametric resonance channel that
   destabilizes the oscillons. This is the OPPOSITE of gravitational
   self-trapping.

3. **Two-body force**: No measurable force between oscillons — the
   Ω radiation destroys them before a static potential can form.

4. **Frequency shift**: The ω shift is upward (+1.4%), reducing the
   sub-gap margin. Gravity should create a redshift (downward ω shift).

---

## Files

| File | Description |
|------|-------------|
| `src/gauge.c` | Simulation code |
| `data/gauge_test{1,2,3}_ts.tsv` | Time series |
| `data/gauge_test{1,2,3}_spectrum.tsv` | DFT power spectra |
| `data/gauge_test{1,2,3}_profile_t*.tsv` | Spatial profiles at selected times |

Compile: `gcc -O3 -Wall -o gauge src/gauge.c -lm`
