# Computation 1: θ Frequency Spectrum Results

## Method

Analyzed 264 frames from V34 sfa_hires.sfa (N=80, L=25, T=50, η=0.5,
m_θ=0, dt_frame=0.19). Decomposed θ into cylindrical components
(θ_r, θ_φ, θ_z) around the braid axis. FFT of θ_φ(t) at six radii.

## Frequency Peaks

| Radius | ω₁ (dominant) | Period | Power | Character |
|--------|--------------|--------|-------|-----------|
| 3 | 0.219 | 4.56 | 6.76e-3 | Braid breathing mode |
| 5 | 0.219 | 4.56 | 5.53e-3 | Same breathing mode |
| 8 | 0.040 | 25.1 | 1.48e-3 | Low-frequency wave |
| 12 | 0.020 | 50.1 | 8.87e-4 | Even lower frequency |
| 16 | 0.020 | 50.1 | 7.62e-4 | Same |
| 20 | 0.040 | 25.1 | 6.83e-4 | Low-frequency wave |

## Two distinct features

### 1. Braid breathing mode (ω = 0.219, T ≈ 4.6)

Present at ALL radii but strongest near the core (r=3-5). This is the
braid's internal oscillation frequency — the source that drives the
θ field. It matches the V34 finding (θ oscillation period ~4 time units).

This is the only truly DISCRETE frequency in the spectrum. It comes
from the braid's own dynamics, not from the θ potential well.

### 2. Low-frequency dispersive waves

The dominant low-frequency peak SHIFTS with radius:
- r=3-5: ω ≈ 0.06 (T=17)
- r=8: ω ≈ 0.04 (T=25)
- r=12-16: ω ≈ 0.02 (T=50)
- r=20: ω ≈ 0.04 (T=25)

Frequency depends on position → DISPERSIVE (continuous), not DISCRETE.
These are outward-propagating θ waves that slow down / lengthen with
distance, consistent with massive-like dispersion in the braid's
effective potential.

## Verdict

**The θ spectrum is CONTINUOUS.** No discrete resonant frequencies
emerge from the geometry alone. The only discrete feature is the
braid's own breathing mode (ω=0.22), which is the radiation source.

**Implication**: Quantized energy levels (electron orbitals) require
imposing ℏ externally through field quantization (ε). They do NOT
emerge spontaneously from the classical field dynamics.

This directs us to Computation 3: add the field quantum ε to the
simulation and determine ℏ_sim = ε × dx.

## Files

- `theta_frequency.py` — analysis script
- `data/theta_spectra.tsv` — full power spectra at all radii
