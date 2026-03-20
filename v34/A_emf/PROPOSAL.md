# Track A: Electromagnetism from the Background Mode Structure

**Motivation**: The phonon test revealed that the background carries a
massless collective mode (the amplitude phonon) that mediates gravity.
The background has THREE fields with specific phase relationships. This
means there should be MULTIPLE collective modes — not just the amplitude
mode. The other modes are candidates for electromagnetism.

---

## The Mode Structure of the Background

The background is φ_a = A_bg × cos(kz + 2πa/3), with a = 0, 1, 2.
At each point in space, the three fields form a vector in field space
(φ₀, φ₁, φ₂) that rotates as you move along z.

Perturbations of this background decompose into three types:

### Mode 1: Amplitude (all three grow/shrink together)

    δφ_a = ε(x,t) × cos(kz + 2πa/3)

This changes the LOCAL AMPLITUDE of the background without changing the
phase relationships. The triple product P changes proportionally.
This is the mode that creates depletion → the PHONON → GRAVITY.

The phonon test confirmed this mode is massless and decays as 1/r^1.2.

### Mode 2: Relative Phase (phase shifts between fields)

    δφ_a = A_bg × [cos(kz + 2πa/3 + δ_a(x,t)) - cos(kz + 2πa/3)]

where δ₀ + δ₁ + δ₂ = 0 (no net phase shift = not a spatial translation).

This changes the ORIENTATION of the field triplet in field space without
changing the total amplitude Σφ². It rotates the helical pattern without
stretching or compressing it.

**Key property**: This mode does NOT create depletion (Σφ² unchanged) →
no gravitational effect. But it DOES change P (the relative alignment
matters for the triple product). So it interacts with braids (which have
specific P structure) but not with gravity.

This is exactly what EM should be:
- Propagating wave that doesn't gravitate (photons are massless in GR)
- Interacts with charged particles (braids with specific winding)
- Has two polarizations (two independent relative phase shifts)

### Mode 3: Global Phase (spatial translation)

    δφ_a = A_bg × [cos(kz + δ + 2πa/3) - cos(kz + 2πa/3)]

This is just a shift z → z - δ/k. It's a phonon with ω = ck.
Same as Mode 1 in the long-wavelength limit.

---

## Predictions

If Mode 2 (relative phase) is the EM analog:

| Property | Standard EM | SCP Phase Mode |
|----------|------------|----------------|
| Speed | c | c (if massless) |
| Mass | 0 | 0 (if Goldstone mode) |
| Polarizations | 2 | 2 (two independent δ_a with constraint) |
| Gravitates? | Yes (E=hν) | No (Σφ² unchanged) — PROBLEM |
| Couples to | charge | braid winding (helical handedness) |
| Source | accelerating charge | braid mode transitions |

**The gravitational coupling is a problem**: In GR, photons gravitate
(light bends around stars, E=mc²). If the phase mode doesn't change
Σφ², it doesn't create depletion → doesn't gravitate. This would
violate the equivalence principle.

**Possible resolution**: The phase mode may change Σφ² at second order
(δΣφ² ∝ δ²), giving a weak gravitational coupling proportional to
the wave's energy. This would be consistent with GR where photon
gravity is proportional to energy.

---

## Experimental Plan

### Phase 1: Identify Background Modes (Analytical + Numerical)

**1a. Analytical linearization**:
Write φ_a = φ_a^(bg)(z,t) + δφ_a(x,t). Expand Eq. (1) to linear order
in δφ. The linearized equation is:

    ∂²δφ_a/∂t² = ∇²δφ_a - m²δφ_a - Σ_b M_ab(z,t) δφ_b

where M_ab = ∂²V/∂φ_a∂φ_b evaluated at the background. This is a
3×3 matrix that couples the three field perturbations.

Diagonalize M_ab → three eigenmodes with effective masses m₁, m₂, m₃.
If any m_i = 0, that's a massless mode.

**1b. Numerical dispersion measurement**:
Initialize a small-amplitude plane wave in each mode direction.
Measure the propagation speed and decay rate. Plot ω vs k.
- ω² = k² + m²  → massive (Yukawa mediator)
- ω = c|k|       → massless (photon/phonon candidate)

### Phase 2: Wave Packet Propagation

For each identified mode, initialize a localized wave packet:

    δφ_a(x,0) = ε × e_a × exp(-(x-x₀)²/2σ²) × cos(k₀·x)

where e_a is the mode eigenvector (amplitude, phase-1, or phase-2).

**Measure**:
1. Does the packet maintain shape? (non-dispersive = good)
2. Propagation speed? (c = massless)
3. Energy conservation during propagation?
4. Does it radiate or scatter off the background structure?

**Three packets to test**:
- Amplitude packet: e_a ∝ cos(kz + 2πa/3) (aligned with background)
- Phase packet 1: e_a ∝ (-sin(kz), sin(kz+2π/3)/2, sin(kz+4π/3)/2)
  (relative phase shift, δ₀=-δ₁-δ₂)
- Transverse packet: e_a = (1, 0, 0) (perpendicular to background)

### Phase 3: Braid-Wave Interaction

Fire each wave packet at a stationary braid:

**Setup**: Braid at origin, wave packet initialized at x=-20 moving
in +x direction. N=256, L=40.

**Measure**:
1. Transmission/reflection/absorption coefficients
2. Does the braid move? (momentum transfer → scattering cross-section)
3. Does the braid's internal state change? (excitation → absorption)
4. Is there a difference between amplitude and phase packets?

If amplitude packets scatter/absorb strongly but phase packets pass
through (or vice versa): the two modes interact differently with
braids → different "charges" → gravity vs EM distinction.

### Phase 4: Braid Emission Spectrum

Perturb a settled braid and analyze the radiation:

**Methods of perturbation**:
1. Velocity kick (as in drag test)
2. Amplitude compression (squeeze the Gaussian envelope)
3. Phase kick (shift one field's phase relative to others)

**Measure** the emitted radiation at r=15-30:
1. Decompose into amplitude and phase components
2. Frequency spectrum (FFT in time at fixed r)
3. Are there discrete frequencies? (→ quantized emission)
4. Angular distribution (isotropic? dipole? quadrupole?)

If the braid emits at discrete frequencies → the braid has internal
energy levels. Transitions between levels emit wave packets = photons.

---

## Implementation Notes

### Code structure
Base: v33_G.c (already has single alloc, Verlet, periodic BC)
Modifications needed:
- Wave packet initialization (Gaussian × plane wave × mode eigenvector)
- Mode decomposition diagnostic (project field onto amplitude/phase basis)
- Spectral analysis (FFT at fixed radius)

### The linearization matrix M_ab

At the background, P_bg = φ₀φ₁φ₂. The coupling matrix is:

    M_ab = ∂²V/∂φ_a∂φ_b = V''(P) × (∂P/∂φ_a)(∂P/∂φ_b) + V'(P) × ∂²P/∂φ_a∂φ_b

    ∂P/∂φ_0 = φ₁φ₂,  ∂P/∂φ_1 = φ₀φ₂,  ∂P/∂φ_2 = φ₀φ₁

    ∂²P/∂φ_0∂φ_1 = φ₂,  ∂²P/∂φ_0∂φ_2 = φ₁,  ∂²P/∂φ_1∂φ_2 = φ₀
    ∂²P/∂φ_a² = 0

So M is:
    M_00 = V'' × (φ₁φ₂)²
    M_01 = V'' × (φ₁φ₂)(φ₀φ₂) + V' × φ₂
    M_11 = V'' × (φ₀φ₂)²
    ... (symmetric, with off-diagonal V' terms)

This matrix oscillates with z (through the background). The eigenmodes
are Bloch waves — periodic in z with a band structure. The band gaps
and dispersion relations determine which modes are massive vs massless.

### Resource estimate
- Phase 1a: analytical, no compute (pen and paper + Python verification)
- Phase 1b: N=256, L=40, T=50 per mode × 3 modes = ~30 min total
- Phase 2: N=256, L=40, T=100 per packet × 3 packets = ~60 min total
- Phase 3: N=256, L=40, T=200 per test × 3 packets = ~120 min total
- Phase 4: N=256, L=40, T=500 per perturbation × 3 methods = ~300 min

Total: ~8.5 hours if run sequentially. Parallelizable to ~3 hours.

---

## What Success Looks Like

**Minimum viable result**: Identify at least TWO distinct mode types
in the background, one that couples to Σφ² (gravity) and one that
doesn't (EM candidate). Show they propagate at different speeds or
have different dispersion relations.

**Strong result**: Show that the phase mode is massless, propagates
at c, has two polarizations, and interacts differently with braids
than the amplitude mode. This would be EM emerging from the same
equation that gives gravity, with no new fields or couplings.

**Grand slam**: Show that the braid emits phase-mode radiation at
discrete frequencies when perturbed, and that these frequencies
correspond to internal mode transitions. This would be photon
emission from atomic energy levels — quantum mechanics emerging
from classical field dynamics.
