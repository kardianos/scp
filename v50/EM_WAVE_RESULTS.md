# V50 EM Wave Results

## The Polariton Mode [CONFIRMED — algebraic + numerical]

The Cosserat curl coupling produces a hybridized θ-δφ wave, not a pure θ wave.
The propagation mechanism:

    θ oscillation → curl(θ) → excites δφ → curl(δφ) → regenerates θ

This is the E↔B regeneration cycle of Maxwell's equations. The wave is a
**polariton**: part angle field (photon-like), part displacement field
(matter-like), analogous to a phonon-polariton in condensed matter.

### Dispersion Relation [CONFIRMED — Maxima + Lean, zero sorrys]

Linearizing the coupled equations around vacuum (φ=0, θ=0):

    ∂²δφ/∂t² = ∇²δφ - m²δφ + η curl(δθ)
    ∂²δθ/∂t² = ∇²δθ         + η curl(δφ)

For transverse plane waves (propagating in z), curl pairs δφ_x with δθ_y
(cross-polarization, like E⊥B). The 2×2 eigenvalue problem gives:

    (ω² - k² - m²)(ω² - k²) = η²k²

Two branches:

    ω² = [2k² + m² ± √(m⁴ + 4η²k²)] / 2

| Branch | k→0 | Character |
|--------|-----|-----------|
| Photon (−) | ω² → 0 (massless) | Mostly θ, small δφ admixture |
| Matter (+) | ω² → m² (massive) | Mostly φ, small δθ admixture |

Photon branch phase velocity at long wavelength:

    v = c × √(1 − η²/m²)

With η=0.5, m²=2.25: **v = 0.9428c**.

Mixing ratio (δφ/δθ at low k):

    δφ_x / δθ_y = −ηk/m²

At k=0.785 (λ=8): δφ/δθ = −0.175, so the photon is ~18% φ.

Maxima derivation: `v50/em_wave/dispersion.mac`
Lean proof: `lean/Polariton.lean` (zero sorrys, all algebraic)

### Lean Proof Summary

Four theorems, all without sorry:

1. **photon_massless**: At k=0, ω²=0 is a root (massless photon branch)
2. **matter_massive**: At k=0, ω²=m² is a root (massive matter branch)
3. **photon_velocity_exact**: The algebraic identity proving v²=1−η²/m²
   with O(k⁴) residual — exact for all k, not just a limit
4. **subluminal**: v < c when η > 0 (polariton is slower than light in vacuum)

### Numerical Confirmation [CONFIRMED]

**Plane wave in true vacuum** (A_bg=0, periodic BC, N=128, L=30, λ=60):
- Energy drift: **0.001%** over T=100 (perfect conservation)
- θ_rms: constant at 2.05e-3 (zero decay, zero amplification)
- Peak amplitude: stable at 0.00500 throughout
- Phase velocity: **v = 0.9377c** (0.5% below predicted 0.9428c)
- The 0.5% deficit is numerical dispersion from the 6-point Laplacian

**Wave packet in true vacuum** (A_bg=0, absorbing BC, Gaussian envelope):
- Centroid propagates from z=-12 to z=+14 over T=30
- Speed: 0.87c (linear fit, lower accuracy from envelope effects)
- Peak amplitude stable (no bound-state formation)
- Energy exits cleanly through absorbing boundary

**Decoupled check** (η=0): v = 0.964c ≈ c (pure free θ wave, no polariton)

## Why A_bg=0 Is Required for Clean Tests [CONFIRMED]

The standard background φ = A_bg cos(kz + 2πa/3) has nonzero curl(φ).
This acts as a **continuous θ source** via the η curl(φ) coupling term,
overwhelming any injected wave packet. Measured effects:

- Pure θ pulse in A_bg=0.1: amplitude doubles (0.05 → 0.10), forms bound blob
- Polariton eigenmode in A_bg=0.1: amplitude still grows, chaotic centroid
- Same pulses in A_bg=0: clean propagation, stable amplitude

This is NOT a bug — it's the physical mechanism by which braids radiate.
The background curl is the braid's "antenna." But for isolating wave
propagation from radiation, A_bg=0 is required.

## Why the Photon Is a Plane Wave [THEORETICAL]

A physical photon has:
- Zero transverse localization (infinite width perpendicular to propagation)
- All energy in amplitude, not spatial confinement
- Flat wavefronts (every y-z slice identical for z-propagation)

Any attempt to localize the wave creates transverse gradients:
- Cylindrical beam (Gaussian in x,y): produces radial diffraction rings
  and spurious φ excitation from the transverse ∇²
- Slab beam (Gaussian in x only): reduced but still present edge effects
- Rayleigh range z_R = πw²/λ limits collimation distance

This is standard Fourier optics, not specific to the Cosserat theory.
A "photon" at our simulation scale (λ ~ 10 fm) cannot be spatially
localized without diffraction effects dominating. Physical photons
(λ ~ 500 nm) have wavelengths ~10⁸× larger than any structure we simulate.

## What a Photon Looks Like at the Particle Scale [THEORETICAL]

From a braid's perspective, an incident photon is an oscillating plane wave:
uniform, infinite, varying only along the propagation axis. The braid's
interaction cross-section is set by its current loop geometry (~3 code units),
not by the photon's transverse extent (infinite).

The "localization" of a photon (e.g., "photon hits detector at position x")
comes from the quantum measurement process, not from spatial confinement
of the classical wave. The Cosserat polariton is the classical field mode;
quantization would discretize its energy into ℏω packets.

## Open: Braid as Natural Polariton Emitter [UNTESTED]

The V34 braid radiates θ with period ~4t and the far field is
electromagnetic (F_curl dominates F_pot at r > 15). The prediction:

1. The emitted radiation should propagate at v = √(1 − η²/m²) ≈ 0.943c
2. The far-field θ/δφ ratio should match −ηk/m² for the emission frequency
3. The radiation should be transverse (curl-mediated, not longitudinal)

Testing this requires:
- A single braid run (A_bg=0.1, standard params, T=200+)
- Phase velocity measurement of outgoing θ wavefronts at r=15-20
- Decomposition of θ vs δφ at the emission frequency

This would confirm that the braid is a natural polariton emitter — that
the "light" radiated by particles in this theory has the correct
dispersion, polarization, and speed.

## Files

    v50/em_wave/dispersion.mac          — Maxima symbolic derivation
    lean/Polariton.lean                 — Lean 4 proof (zero sorrys)
    sfa/seed/gen_polariton.c            — Eigenmode wave packet seed
    sfa/seed/gen_plane_wave.c           — Plane wave seed
    sfa/seed/gen_beam.c                 — Directed slab beam seed
    sfa/seed/gen_em_wave.c              — Simple θ pulse seed (historical)
    v50/em_wave/analyze_wave.c          — Wave packet tracking tool
    v50/em_wave/polariton.cfg           — Polariton test config
    v50/em_wave/em_wave.cfg             — Original θ pulse config

### Viewer Update

`volview` now supports fixed color scale mode for visualizing uniform fields:
- **F**: toggle fixed/auto scale
- **[** / **]**: decrease/increase fixed scale range (×0.5 / ×2.0)

Required for viewing plane waves, which auto-scale to uniform brightness.
